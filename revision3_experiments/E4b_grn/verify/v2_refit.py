"""Verifier check 2: refit GENIE3 / GRNBoost2 with the verifier's own code (spec hyperparameters, raw expression from
cache, regulators rebuilt in v1) and compare the per-target normalised importances with E4b's partial/*.npz.

Jobs: kidney_A T1 genie3 s0 + grnboost2 s0 (all 152 targets); lung T2 grnboost2 s1 (30 random targets + the unfitted
zero-variance targets); immune T2 genie3 s1 (30 random targets). Writes verify/v2_refit.json.
"""
import json, sys
from pathlib import Path
import numpy as np
from joblib import Parallel, delayed
from sklearn.ensemble import ExtraTreesRegressor, GradientBoostingRegressor

V = Path(__file__).resolve().parent; E = V.parent; CACHE = E.parent / "cache"


def model(method, seed):
    if method == "genie3":
        return ExtraTreesRegressor(n_estimators=300, max_features="sqrt", random_state=seed, n_jobs=1)
    return GradientBoostingRegressor(learning_rate=0.01, n_estimators=1000, max_depth=3, subsample=0.9, max_features=0.1,
                                     n_iter_no_change=25, validation_fraction=0.1, random_state=seed)


def fit(Xall, regs, t, method, seed):
    y = Xall[:, t].astype(np.float64); keep = regs != t
    imp = np.zeros(regs.size)
    if y.var() == 0: return imp
    m = model(method, seed).fit(Xall[:, regs[keep]], y)
    fi = m.feature_importances_; imp[keep] = fi / fi.sum()
    return imp


res = {}
jobs = [("kidney_A", "T1", "genie3", 0, None), ("kidney_A", "T1", "grnboost2", 0, None),
        ("lung", "T2", "grnboost2", 1, 30), ("immune", "T2", "genie3", 1, 30)]
for ctx, tier, method, seed, nsub in jobs:
    z = np.load(E / "partial" / f"{ctx}_{tier}_{method}_s{seed}.npz")
    regs, targets, imp_e4b = z["regs"], z["targets"], z["imp"].astype(np.float64)
    Xall = np.asarray(np.load(CACHE / ctx / ("expr_extraction.npy" if tier == "T1" else "expr.npy")), dtype=np.float32)
    if nsub:
        rng = np.random.default_rng(123); sel = np.sort(rng.choice(targets.size, nsub, replace=False))
        sel = np.union1d(sel, np.flatnonzero(z["fitted"] == 0))
    else:
        sel = np.arange(targets.size)
    mine = np.stack(Parallel(n_jobs=3)(delayed(fit)(Xall, regs, int(targets[k]), method, seed) for k in sel))
    ref = imp_e4b[sel]
    key = f"{ctx}_{tier}_{method}_s{seed}"
    res[key] = dict(n_targets_checked=int(sel.size), max_abs_diff=float(np.abs(mine - ref).max()),
                    rows_identical_tol1e6=int((np.abs(mine - ref).max(1) < 1e-6).sum()),
                    n_unfitted_in_check=int((z["fitted"][sel] == 0).sum()),
                    unfitted_rows_all_zero=bool(np.all(ref[z["fitted"][sel] == 0] == 0)))
    print(key, res[key], flush=True)
(V / "v2_refit.json").write_text(json.dumps(res, indent=1))
