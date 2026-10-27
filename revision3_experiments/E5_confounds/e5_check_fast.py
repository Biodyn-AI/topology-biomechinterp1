"""Check clogit_fit_fast against clogit_fit on bootstrap and permutation refits (lung L0, kidney_B L0, immune L3)."""
import sys, time, json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
import e5_lib as L
S = L.S
out = []
for ctx, layer, nb in [("immune", 3, 10), ("lung", 0, 4)]:
    C = L.load_context(ctx); emb = L.load_emb(ctx)
    tf, tg, y = C["tf"], C["tg"], C["y"].astype(np.int64)
    tfs, tfc, tgs, tgc = L.tf_target_codes(tf, tg); sp = S.StratifiedPairs(tf, tg, y)
    Xa = np.column_stack([L.zscore(L.layer_scores(emb, C["obs_idx"], layer, tf, tg)), L.zscore(C["cov_raw"])])
    fa = L.clogit_fit(L.CLogitData(Xa, y, tfc))
    for b in range(nb):
        wtf, wtg = L.boot_weights(ctx, b, tfs.size, tgs.size)
        Xe, ye, ge = L.expand(Xa, y, tfc, wtf, wtg, tgc); D = L.CLogitData(Xe, ye, ge, gw_by_group=wtf)
        t0 = time.time(); f1 = L.clogit_fit(D, beta0=fa["beta"]); t1 = time.time() - t0
        t0 = time.time(); f2 = L.clogit_fit_fast(D, beta0=fa["beta"]); t2 = time.time() - t0
        out.append(dict(ctx=ctx, kind="boot", b=b, max_abs_beta_diff=float(np.max(np.abs(f1["beta"] - f2["beta"]))),
                        t_newton=t1, t_fast=t2, iters_fast=f2["iters"]))
        yp = L.perm_labels(ctx, b, sp).astype(np.int64); D = L.CLogitData(Xa, yp, tfc)
        t0 = time.time(); f1 = L.clogit_fit(D); t1 = time.time() - t0
        t0 = time.time(); f2 = L.clogit_fit_fast(D); t2 = time.time() - t0
        out.append(dict(ctx=ctx, kind="perm", b=b, max_abs_beta_diff=float(np.max(np.abs(f1["beta"] - f2["beta"]))),
                        t_newton=t1, t_fast=t2, iters_fast=f2["iters"]))
        print(out[-2], out[-1], flush=True)
(Path(__file__).parent / "check_fast.json").write_text(json.dumps(dict(
    max_abs_beta_diff=max(o["max_abs_beta_diff"] for o in out), rows=out), indent=1))
print("max diff", max(o["max_abs_beta_diff"] for o in out))
