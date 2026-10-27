"""E4b: GENIE3 / GRNBoost2 (sklearn re-implementations) on the H3 target set. Resumable.

Usage
  python e4b_grn.py --jobs lung:T1:genie3:0 lung:T1:grnboost2:0 ...   (run listed jobs in order)
  python e4b_grn.py --time ctx:tier:method:seed --limit 50              (timing estimate, writes nothing final)

Each job = (context, tier, method, seed). Targets are processed in batches of BATCH; each finished batch is
saved to partial/<ctx>_<tier>_<method>_s<seed>/b<k>.npz, so a killed job resumes at the first missing batch.
When all batches exist the job is assembled into partial/<ctx>_<tier>_<method>_s<seed>.npz
(imp: n_targets x n_regs float32, targets, regs). See e4b_common.py for the full specification.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
from joblib import Parallel, delayed

sys.path.insert(0, str(Path(__file__).resolve().parent))
import e4b_common as C  # noqa: E402

PART = C.HERE / "partial"
BATCH = 50
N_JOBS = int(__import__("os").environ.get("E4B_NJOBS", "3"))  # <= 3


def make_model(method: str, seed: int):
    if method == "genie3":
        from sklearn.ensemble import ExtraTreesRegressor
        return ExtraTreesRegressor(n_estimators=300, max_features="sqrt", random_state=seed, n_jobs=1)
    if method == "grnboost2":
        from sklearn.ensemble import GradientBoostingRegressor
        return GradientBoostingRegressor(learning_rate=0.01, n_estimators=1000, max_depth=3, subsample=0.9,
                                         max_features=0.1, n_iter_no_change=25, validation_fraction=0.1,
                                         random_state=seed)
    raise ValueError(method)


def fit_one(Xr: np.ndarray, y: np.ndarray, self_col: int, method: str, seed: int):
    """Normalised importances of all regulators for one target; the target's own column (if it is a
    regulator) is removed from the features and gets importance 0."""
    n_reg = Xr.shape[1]
    cols = np.arange(n_reg) if self_col < 0 else np.delete(np.arange(n_reg), self_col)
    out = np.zeros(n_reg, dtype=np.float32)
    n_iter = -1
    if np.var(y) <= 0:
        return out, 0, n_iter
    m = make_model(method, seed)
    m.fit(Xr[:, cols], y)
    fi = np.asarray(m.feature_importances_, dtype=np.float64)
    s = fi.sum()
    if s > 0:
        out[cols] = (fi / s).astype(np.float32)
    if method == "grnboost2":
        n_iter = int(m.n_estimators_)
    return out, 1, n_iter


def lock_alive(bdir, k) -> bool:
    lf = bdir / f"b{k:04d}.lock"
    try:
        pid = int(lf.read_text().strip() or 0)
    except (FileNotFoundError, ValueError):
        return False
    try:
        os.kill(pid, 0); return True
    except OSError:
        return False


def claim(bdir, k) -> bool:
    """Atomically claim batch k (lock file holding our pid); a lock whose pid is dead is taken over."""
    lf = bdir / f"b{k:04d}.lock"
    for _ in range(2):
        try:
            fd = os.open(lf, os.O_CREAT | os.O_EXCL | os.O_WRONLY); os.write(fd, str(os.getpid()).encode()); os.close(fd)
            return True
        except FileExistsError:
            if lock_alive(bdir, k): return False
            lf.unlink(missing_ok=True)
    return False


def run_job(ctx, tier, method, seed, limit=None, log=print):
    sets = C.load_sets(ctx)
    regs, targets = sets["regs"], sets["targets"]
    if limit:
        targets = targets[:limit]
    expr = C.load_expr(ctx, tier)
    Xr = np.ascontiguousarray(expr[:, regs])
    rpos = {int(r): k for k, r in enumerate(regs)}
    tag = f"{ctx}_{tier}_{method}_s{seed}"
    final = PART / f"{tag}.npz"
    if final.exists() and not limit:
        log(f"[{tag}] already complete"); return
    bdir = PART / tag; bdir.mkdir(parents=True, exist_ok=True)
    nb = int(np.ceil(targets.size / BATCH))
    t0 = time.time(); done_now = 0
    with Parallel(n_jobs=N_JOBS, backend="loky") as par:
        for k in range(nb):
            bf = bdir / f"b{k:04d}.npz"
            if bf.exists() and not limit:
                continue
            if not limit and not claim(bdir, k):
                continue                    # another live process is fitting this batch
            if not limit and bf.exists():
                (bdir / f"b{k:04d}.lock").unlink(missing_ok=True); continue
            tb = targets[k * BATCH:(k + 1) * BATCH]
            res = par(delayed(fit_one)(Xr, expr[:, t].astype(np.float64), rpos.get(int(t), -1), method, seed)
                      for t in tb)
            imp = np.stack([r[0] for r in res]); fitted = np.array([r[1] for r in res], np.int8)
            nit = np.array([r[2] for r in res], np.int32)
            if not limit:
                tmp = bdir / f"tmp_b{k:04d}_{os.getpid()}.npz"
                np.savez(tmp, imp=imp, targets=tb, fitted=fitted, n_iter=nit); os.replace(tmp, bf)
                (bdir / f"b{k:04d}.lock").unlink(missing_ok=True)
            done_now += tb.size
            el = time.time() - t0
            log(f"[{tag}] batch {k + 1}/{nb}  {done_now} targets in {el:.0f}s  ({el / done_now:.2f}s/target)"
                + (f"  median n_iter {int(np.median(nit))}" if method == "grnboost2" else ""))
    if limit:
        return (time.time() - t0) / max(done_now, 1)
    while True:                             # wait for batches other live processes are still fitting
        missing = [k for k in range(nb) if not (bdir / f"b{k:04d}.npz").exists()]
        if not missing: break
        if not any(lock_alive(bdir, k) for k in missing):
            log(f"[{tag}] batches {missing[:5]}... missing with no live owner; rerun to fill"); return
        time.sleep(20)
    if final.exists():
        log(f"[{tag}] assembled by another process"); return
    parts = [np.load(bdir / f"b{k:04d}.npz") for k in range(nb)]
    imp = np.concatenate([p["imp"] for p in parts]); tt = np.concatenate([p["targets"] for p in parts])
    assert np.array_equal(tt, targets)
    tmpf = PART / f"tmp_{tag}_{os.getpid()}.npz"
    np.savez(tmpf, imp=imp, targets=targets, regs=regs,
             fitted=np.concatenate([p["fitted"] for p in parts]), n_iter=np.concatenate([p["n_iter"] for p in parts]),
             n_cells=np.int64(expr.shape[0]), seconds=np.float64(time.time() - t0))
    os.replace(tmpf, final)
    log(f"[{tag}] COMPLETE {targets.size} targets x {regs.size} regulators")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", nargs="*", default=[])
    ap.add_argument("--time", nargs="*", default=[])
    ap.add_argument("--limit", type=int, default=50)
    a = ap.parse_args()
    log = lambda s: print(time.strftime("%H:%M:%S ") + s, flush=True)  # noqa: E731
    for j in a.time:
        ctx, tier, method, seed = j.split(":")
        s = run_job(ctx, tier, method, int(seed), limit=a.limit, log=log)
        n = C.load_sets(ctx)["targets"].size
        log(f"TIMING {j}: {s:.3f} s/target x {n} targets = {s * n / 60:.1f} min per seed")
    for j in a.jobs:
        ctx, tier, method, seed = j.split(":")
        run_job(ctx, tier, method, int(seed), log=log)


if __name__ == "__main__":
    main()
