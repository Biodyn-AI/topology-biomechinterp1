"""Independent verifier, part 3: residualised within-TF AUROC two-way bootstrap CI for selected cells, re-implemented
without e5_lib: weighted OLS of the score on covariates + explicit TF dummies (pair weight = TF mult x target mult),
weighted Mann-Whitney by explicit comparison matrices per TF. Draws regenerated from default_rng([20261026,5,ctx_id,b]).
usage: python verify_resid_boot.py <ctx> <layer> [<ctx> <layer> ...]
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
V = Path(__file__).resolve().parent; E5 = V.parent; R3 = E5.parent
CTXS = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"


def wauc(s, y, tfc, w):
    num = den = 0.0
    for t in range(tfc.max() + 1):
        m = tfc == t; p = m & (y == 1); n = m & (y == 0)
        wp, wn = w[p], w[n]
        if wp.sum() == 0 or wn.sum() == 0: continue
        sp, sn = s[p], s[n]
        cmp = (sp[:, None] > sn[None, :]) + 0.5 * (sp[:, None] == sn[None, :])
        num += float(wp @ cmp @ wn); den += float(wp.sum() * wn.sum())
    return num / den


out = {}
a = sys.argv[1:]
for ctx, layer in zip(a[::2], map(int, a[1::2])):
    d = R3 / "cache" / ctx
    obs_idx = np.load(d / "obs_idx.npy"); P = np.load(d / "h3_pairs.npz")
    tf, tg, y = P["tf"].astype(int), P["tg"].astype(int), P["y"].astype(int)
    tfu, tfc = np.unique(tf, return_inverse=True); tgu, tgc = np.unique(tg, return_inverse=True)
    cov = np.load(V / f"cov_{ctx}.npy")
    emb = np.load(EMB / f"{ctx}.npy", mmap_mode="r")
    M = np.asarray(emb[layer], np.float64)[obs_idx]; M = M - M.mean(0); M = M / np.linalg.norm(M, axis=1, keepdims=True)
    s = (M[tf] * M[tg]).sum(1)
    D = np.zeros((y.size, tfu.size)); D[np.arange(y.size), tfc] = 1
    Xf = np.column_stack([cov, D])
    bs = []
    for b in range(500):
        rng = np.random.default_rng([20261026, 5, CTXS.index(ctx), b])
        wtf = np.bincount(rng.integers(0, tfu.size, tfu.size), minlength=tfu.size)
        wtg = np.bincount(rng.integers(0, tgu.size, tgu.size), minlength=tgu.size)
        w = (wtf[tfc] * wtg[tgc]).astype(float); keep = w > 0; sw = np.sqrt(w[keep])
        beta = np.linalg.lstsq(Xf[keep] * sw[:, None], s[keep] * sw, rcond=None)[0]
        r = s - cov @ beta[:5]
        bs.append(wauc(r, y, tfc, w))
        if b % 100 == 0: print(ctx, layer, b, bs[-1], flush=True)
    bs = np.array(bs); e = np.load(E5 / "ckpt" / f"{ctx}_L{layer:02d}.npz")["boot_auc_res"]
    ra = pd.read_csv(E5 / "residual_auroc.csv").query("context == @ctx and layer == @layer").iloc[0]
    out[f"{ctx}_L{layer}"] = dict(max_abs_draw_diff=float(np.max(np.abs(bs - e))),
                                  ci_verify=np.quantile(bs, [0.025, 0.975]).tolist(), ub95_verify=float(np.quantile(bs, 0.95)),
                                  ci_e5=[float(ra.auroc_resid_lo), float(ra.auroc_resid_hi)], ub95_e5=float(ra.auroc_resid_upper95_onesided))
    print(json.dumps(out[f"{ctx}_L{layer}"]), flush=True)
    (V / f"verify_resid_boot_{'_'.join(a)}.json").write_text(json.dumps(out, indent=1))
