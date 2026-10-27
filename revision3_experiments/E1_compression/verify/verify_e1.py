"""Independent verification of E1_compression headline numbers (adversarial audit).

Does NOT import stats3 or e1_compression: own SVD, own ER/PR, own observed mask, own shuffle, own matching,
own bootstrap (different seed). Compares against the implementer's CSVs and writes verify/verify_results.csv
(one row per check: name, reproduced, reported, abs_diff, tolerance, pass) and verify/verify_extra.json.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
E1 = HERE.parent
R3 = E1.parent
CACHE = R3 / "cache"
LOCAL = Path.home() / "biomechinterp-local/subproject41_r3/emb"
sys.path.insert(0, str(R3))
import contexts as X  # noqa: E402  (only for the original file paths)

CORE = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
rng0 = np.random.default_rng(777)
per = pd.read_csv(E1 / "per_context_layer.csv")
gaps = pd.read_csv(E1 / "gaps.csv")
det = pd.read_csv(E1 / "detection_corr.csv")
shuf = pd.read_csv(E1 / "shuffle.csv")
mat = pd.read_csv(E1 / "matched.csv")
zr = pd.read_csv(E1 / "zero_row_artefact.csv")
checks, extra = [], {}


def chk(name, rep, ref, tol):
    checks.append(dict(check=name, reproduced=float(rep), reported=float(ref), abs_diff=float(abs(rep - ref)),
                       tol=tol, ok=bool(abs(rep - ref) <= tol)))


def sv(M):
    M = np.asarray(M, np.float64)
    M = M - M.mean(0)
    return np.linalg.svd(M, compute_uv=False)


def er(s):
    p = s / s.sum(); p = p[p > 1e-300]
    return float(np.exp(-np.sum(p * np.log(p))))


def pr(s):
    v = s ** 2
    return float(v.sum() ** 2 / np.sum(v ** 2))


def v1(s):
    return float(s[0] ** 2 / np.sum(s ** 2))


def P(c, L, col):
    return float(per[(per.context == c) & (per.layer == L)][col].iloc[0])


# ---- 0. data integrity: local copies == originals; observed mask == cache obs_idx; ncells alignment
for c in X.CONTEXTS:
    loc = np.load(LOCAL / f"{c}.npy", mmap_mode="r")
    oi = np.load(CACHE / c / "obs_idx.npy")
    nc = np.load(CACHE / c / "ncells.npy")
    nz = [np.flatnonzero(np.abs(np.asarray(loc[L])).sum(1) > 0) for L in range(12)]
    same_all_layers = all(np.array_equal(nz[0], z) for z in nz)
    cnt_p = X.CONTEXTS[c]["emb"].parent / "layer_gene_embedding_count.npy"
    cnt = np.load(cnt_p)
    cnt0 = cnt[0] if cnt.ndim > 1 else cnt
    d = dict(obs_mask_equal_all_layers=bool(same_all_layers), obs_idx_matches=bool(np.array_equal(nz[0], oi)),
             ncells_equals_count=bool(np.array_equal(cnt0[oi], nc)), count_zero_iff_unobserved=bool(
                 np.array_equal(np.flatnonzero(cnt0 > 0), oi)), V=int(oi.size), G=int(loc.shape[1]))
    if c in ("kidney_A", "immune", "lung"):
        orig = np.load(X.CONTEXTS[c]["emb"], mmap_mode="r")
        d["local_equals_original"] = bool(np.array_equal(np.asarray(orig), np.asarray(loc)))
    extra[c] = d
    print(c, d, flush=True)

# ---- 1. ER / PR / SV1 var on observed genes, core contexts, L0 and L11; all 12 layers for monotonicity
mono = {}
for c in CORE + ["kidney_ckpt_kidney", "kidney_ckpt_brain"]:
    emb = np.load(LOCAL / f"{c}.npy", mmap_mode="r")
    oi = np.flatnonzero(np.abs(np.asarray(emb[0])).sum(1) > 0)
    ers, prs = [], []
    for L in range(12):
        s = sv(np.asarray(emb[L])[oi])
        ers.append(er(s)); prs.append(pr(s))
        if L in (0, 11):
            chk(f"{c} L{L} ER", ers[-1], P(c, L, "er"), 1e-6)
            chk(f"{c} L{L} PR", prs[-1], P(c, L, "pr"), 1e-6)
            chk(f"{c} L{L} SV1var", v1(s), P(c, L, "sv1_var"), 1e-9)
            chk(f"{c} L{L} cum10", float(np.sum(s[:10] ** 2) / np.sum(s ** 2)), P(c, L, "cum10_var"), 1e-9)
        if L == 11:
            # 3. detection Spearman of SV1, SV2 (sign-free) with log n_obs, own SVD
            Xc = np.asarray(emb[L], np.float64)[oi]; Xc = Xc - Xc.mean(0)
            U, S_, _ = np.linalg.svd(Xc, full_matrices=False)
            nc = np.load(CACHE / c / "ncells.npy").astype(float)
            for k in (1, 2):
                r = abs(spearmanr(U[:, k - 1], np.log(nc))[0])
                ref = float(det[(det.context == c) & (det.layer == 11) & (det.k == k)].abs_spearman.iloc[0])
                chk(f"{c} L11 |spearman(SV{k}, log n_obs)|", r, ref, 1e-6)
            # gap ratio point estimates k=1..3
            for k in (1, 2, 3):
                ref = float(gaps[(gaps.context == c) & (gaps.layer == 11) & (gaps.k == k)].ratio.iloc[0])
                chk(f"{c} L11 gap S{k}/S{k+1}", S_[k - 1] / S_[k], ref, 1e-9)
    mono[c] = dict(er=ers, pr=prs, er_strict_dec=bool(np.all(np.diff(ers) < 0)),
                   pr_rises=[(i, i + 1, round(prs[i], 3), round(prs[i + 1], 3)) for i in range(11)
                             if prs[i + 1] > prs[i]],
                   rho_pr=float(spearmanr(np.arange(12), prs)[0]))
extra["monotonicity"] = mono

# ---- 7. zero-row artefact (full matrix incl. zero rows), own computation
for c in ("kidney_A", "immune"):
    emb = np.load(LOCAL / f"{c}.npy", mmap_mode="r")
    G = emb.shape[1]
    ind = (np.abs(np.asarray(emb[0])).sum(1) > 0).astype(float)
    for L in (0, 11):
        Xf = np.asarray(emb[L], np.float64); Xf = Xf - Xf.mean(0)
        U, s, _ = np.linalg.svd(Xf, full_matrices=False)
        z = zr[(zr.context == c) & (zr.layer == L) & (zr.matrix == "full_incl_zero_rows")].iloc[0]
        chk(f"{c} L{L} FULL ER (zero rows incl.)", er(s), z.er, 1e-6)
        chk(f"{c} L{L} FULL SV1var", v1(s), z.sv1_var, 1e-9)
        chk(f"{c} L{L} FULL |pearson(SV1, observed)|", abs(np.corrcoef(U[:, 0], ind)[0, 1]),
            abs(z.pearson_sv1_observed_indicator), 1e-6)
        chk(f"{c} n zero rows", G - ind.sum(), z.n_zero_rows, 0)

# ---- 4. feature-shuffle: own per-column permutation loop, 5 draws, different seed
for c, L in (("lung", 11), ("kidney_A", 11), ("kidney_ckpt_kidney", 11), ("lung", 0)):
    emb = np.load(LOCAL / f"{c}.npy", mmap_mode="r")
    oi = np.load(CACHE / c / "obs_idx.npy")
    Xr = np.asarray(emb[L], np.float64)[oi]
    vals = []
    for _ in range(5):
        Y = np.empty_like(Xr)
        for j in range(Xr.shape[1]):
            Y[:, j] = Xr[rng0.permutation(Xr.shape[0]), j]
        vals.append(er(sv(Y)))
    ref = shuf[(shuf.context == c) & (shuf.layer == L)].iloc[0]
    chk(f"{c} L{L} shuffled ER mean (5 own draws)", np.mean(vals), ref.er_shuf_mean, 5 * ref.er_shuf_sd + 0.5)

# ---- 6. matched n=1300, pooled deciles of log n_obs, own implementation, 40 draws, different seed
data = {c: (np.load(LOCAL / f"{c}.npy", mmap_mode="r"), np.load(CACHE / c / "obs_idx.npy"),
            np.log(np.load(CACHE / c / "ncells.npy").astype(float))) for c in CORE}
pooled = np.concatenate([data[c][2] for c in CORE])
edges = np.unique(np.quantile(pooled, np.linspace(0, 1, 11)))
nb = edges.size - 1
binof = lambda x: np.minimum(np.maximum(np.searchsorted(edges, x, side="right") - 1, 0), nb - 1)
avail = np.array([np.bincount(binof(data[c][2]), minlength=nb) for c in CORE]).min(0)
q = np.floor(avail / avail.sum() * 1300).astype(int)
for j in np.argsort(-(avail / avail.sum() * 1300 - q))[: 1300 - q.sum()]:
    q[j] += 1
extra["matched_quota_own"] = q.tolist()
extra["matched_quota_reported"] = json.loads((E1 / "results.json").read_text())["matched"]["pooled_quota"]
for c in ("kidney_A", "lung"):
    emb, oi, lg = data[c]
    b = binof(lg)
    for L in (0, 11):
        Xr = np.asarray(emb[L], np.float64)[oi]
        vals = []
        for _ in range(40):
            idx = np.concatenate([rng0.choice(np.flatnonzero(b == k), q[k], replace=False) for k in range(nb)])
            vals.append(er(sv(Xr[idx])))
        ref = mat[(mat.context == c) & (mat.layer == L) & (mat.scheme == "pooled")].iloc[0]
        chk(f"{c} L{L} matched-1300 ER mean (40 own draws)", np.mean(vals), ref.er_mean,
            4 * ref.er_sd / np.sqrt(40) + 4 * ref.er_sd / np.sqrt(200) + 0.05)

# ---- 2. gap bootstrap near-degeneracy flags, own bootstrap (B=200, different seed)
gb = []
for c, L, k in (("immune", 0, 1), ("immune", 10, 1), ("kidney_A", 11, 2), ("lung", 11, 3), ("lung", 11, 1),
                ("kidney_A", 11, 1)):
    emb = np.load(LOCAL / f"{c}.npy", mmap_mode="r")
    oi = np.load(CACHE / c / "obs_idx.npy")
    Xr = np.asarray(emb[L], np.float64)[oi]
    s = sv(Xr); r = s[k - 1] / s[k]
    rb = []
    for _ in range(200):
        sb = sv(Xr[rng0.integers(0, oi.size, oi.size)])
        rb.append(sb[k - 1] / sb[k])
    rb = np.array(rb)
    basic_lo = 2 * r - np.quantile(rb, 0.975)
    ref = gaps[(gaps.context == c) & (gaps.layer == L) & (gaps.k == k)].iloc[0]
    gb.append(dict(context=c, layer=L, k=k, ratio=float(r), own_basic_lo=float(basic_lo),
                   own_flag=bool(basic_lo <= 1), reported_basic_lo=float(ref.basic_lo),
                   reported_flag=bool(ref.near_degenerate), own_boot_mean=float(rb.mean()),
                   reported_boot_mean=float(ref.boot_mean)))
    print(gb[-1], flush=True)
extra["gap_bootstrap_recheck"] = gb

out = pd.DataFrame(checks)
out.to_csv(HERE / "verify_results.csv", index=False)
(HERE / "verify_extra.json").write_text(json.dumps(extra, indent=1, default=float))
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
print(out.to_string())
print("ALL OK:", bool(out.ok.all()), f"{int(out.ok.sum())}/{len(out)}")
