"""Independent check of continuity variant C (kidney_A) point estimates: own SVD of the V-centred observed matrix
(numpy, no stats3), SV bands by 0-based column slices, band cosine, residual on covariates + TF dummies, within-TF
AUROC (rankdata). Also checks that spread hashseed_0 equals the main-run variant A/B values, and that C full_512 equals
residual_auroc.csv kidney_A."""
import json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import rankdata
V = Path(__file__).resolve().parent; E5 = V.parent; R3 = E5.parent
d = R3 / "cache" / "kidney_A"; obs_idx = np.load(d / "obs_idx.npy"); P = np.load(d / "h3_pairs.npz")
tf, tg, y = P["tf"].astype(int), P["tg"].astype(int), P["y"].astype(int)
tfu, tfc = np.unique(tf, return_inverse=True); cov = np.load(V / "cov_kidney_A.npy")
D = np.zeros((y.size, tfu.size)); D[np.arange(y.size), tfc] = 1; Xf = np.column_stack([cov, D])
emb = np.load(Path.home() / "biomechinterp-local/subproject41_r3/emb/kidney_A.npy", mmap_mode="r")
def auc(s):
    num = den = 0
    for t in range(tfu.size):
        m = tfc == t; sp, sn = s[m & (y == 1)], s[m & (y == 0)]
        r = rankdata(np.r_[sp, sn]); num += r[:sp.size].sum() - sp.size * (sp.size + 1) / 2; den += sp.size * sn.size
    return num / den
bands = {"SV5_SV7": (4, 7), "SV2_SV7": (1, 7), "SV4_SV6": (3, 6), "SV1_SV6": (0, 6)}
c = pd.read_csv(E5 / "continuity_kidneyA.csv"); c = c[c.variant == "C_paper_benchmark"]
rows = []
for L in range(12):
    X = np.asarray(emb[L], np.float64)[obs_idx]; X = X - X.mean(0)
    U, S, _ = np.linalg.svd(X, full_matrices=False); US = U * S
    for nm, (a, b) in list(bands.items()) + [("full_512", (None, None))]:
        Z = X if a is None and nm == "full_512" else US[:, a:b]
        Zn = Z / np.linalg.norm(Z, axis=1, keepdims=True); s = (Zn[tf] * Zn[tg]).sum(1)
        beta = np.linalg.lstsq(Xf, s, rcond=None)[0]; r = s - cov @ beta[:5]
        e = c[(c.band == nm) & (c.layer.astype(int) == L)].iloc[0]
        rows.append(dict(layer=L, band=nm, r_rb_verify=2 * auc(r) - 1, r_rb_e5=e.r_rb, r_rb_raw_verify=2 * auc(s) - 1, r_rb_raw_e5=e.r_rb_raw))
df = pd.DataFrame(rows); df["d"] = (df.r_rb_verify - df.r_rb_e5).abs(); df["d_raw"] = (df.r_rb_raw_verify - df.r_rb_raw_e5).abs()
df.to_csv(V / "verify_continuity_C.csv", index=False)
ra = pd.read_csv(E5 / "residual_auroc.csv").query("context=='kidney_A'").sort_values("layer")
f512 = df[df.band == "full_512"].sort_values("layer")
main = pd.read_csv(E5 / "continuity_kidneyA.csv"); main = main[main.variant.str[0].isin(["A", "B"]) & (main.layer.astype(str) != "joint_mean_over_layers")]
sp0 = pd.read_csv(E5 / "continuity_hashseed_spread/hashseed_0.csv")
m = main.assign(layer=main.layer.astype(int)).merge(sp0, on=["variant", "band", "layer"], suffixes=("_main", "_hs0"))
res = dict(max_abs_diff_C_resid=float(df.d.max()), max_abs_diff_C_raw=float(df.d_raw.max()),
           max_abs_diff_C_full512_vs_residual_auroc=float(np.max(np.abs(f512.r_rb_e5.to_numpy() - ra.r_rb_resid.to_numpy()))),
           n_A_B_rows_matched=int(len(m)), max_abs_diff_main_vs_hashseed0=float((m.r_rb_main - m.r_rb_hs0).abs().max()),
           n_hashseed_files=len(list((E5 / "continuity_hashseed_spread").glob("hashseed_*.csv"))))
(V / "verify_continuity_C.json").write_text(json.dumps(res, indent=1)); print(json.dumps(res, indent=1))
