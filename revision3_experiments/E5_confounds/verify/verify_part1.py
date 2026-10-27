"""Independent verifier, part 1 (no model fits): recompute from raw inputs, WITHOUT e5_lib, the E5 covariates,
the layer scores, the raw and residualised within-TF AUROC (all 5 contexts x 12 layers), and the stratum-weight
table; re-derive the joint permutation p and bootstrap CIs from the stored draws. Writes verify_part1.json/.csv.

Independent choices: STRING degree re-parsed from the gzipped links file; GO counts from gene2go directly; Pearson r
via np.corrcoef-equivalent on centred columns; residualisation by OLS with explicit TF dummy columns (FWL route
instead of within-centring); AUROC per TF via scipy.stats.rankdata (Mann-Whitney U), pooled num/den.
"""
import gzip, json, pickle, sys, hashlib
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import rankdata

V = Path(__file__).resolve().parent; E5 = V.parent; R3 = E5.parent; REPO = R3.parents[1]
CTXS = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"
out = {}

# ---- convention hash
h = hashlib.sha256((R3 / "CONVENTION3.md").read_bytes()).hexdigest()
frozen = (R3 / "FROZEN.txt").read_text().split()[2]
out["convention_hash_match"] = bool(h == frozen); out["convention_hash"] = h

# ---- STRING genome-wide physical degree (own parse)
info = {}
with gzip.open(REPO / "subproject_61_two_axes_revision/data/string_protein_info_v12.txt.gz", "rt") as f:
    next(f)
    for line in f:
        p = line.split("\t"); info[p[0]] = p[1]
nlines = 0; pairs = set()
with gzip.open(REPO / "subproject_61_two_axes_revision/data/string_links_v12.txt.gz", "rt") as f:
    next(f)
    for line in f:
        nlines += 1
        a, b, s = line.split()
        if int(s) < 700: continue
        ga, gb = info.get(a), info.get(b)
        if not ga or not gb or ga == gb: continue
        pairs.add((ga, gb) if ga < gb else (gb, ga))
deg = {}
for a, b in pairs:
    deg[a] = deg.get(a, 0) + 1; deg[b] = deg.get(b, 0) + 1
out["string_links_data_lines"] = nlines; out["string_pairs_ge700"] = len(pairs)
g2g = pickle.load(open(REPO / "single_cell_mechinterp/data/perturb/gene2go_all.pkl", "rb"))


def auroc_within_tf(s, y, tf, w=None):
    num = den = 0.0
    for t in np.unique(tf):
        m = tf == t; yp = y[m] == 1
        sp, sn = s[m][yp], s[m][~yp]
        if w is None:
            if sp.size == 0 or sn.size == 0: continue
            r = rankdata(np.r_[sp, sn]); U = r[:sp.size].sum() - sp.size * (sp.size + 1) / 2
            num += U; den += sp.size * sn.size
        else:
            wp, wn = w[m][yp], w[m][~yp]
            if wp.sum() == 0 or wn.sum() == 0: continue
            # weighted Mann-Whitney by explicit comparison matrix
            cmp = (sp[:, None] > sn[None, :]) + 0.5 * (sp[:, None] == sn[None, :])
            num += float(wp @ cmp @ wn); den += float(wp.sum() * wn.sum())
    return num / den


def resid_dummies(score, cov, tf):
    u, code = np.unique(tf, return_inverse=True)
    D = np.zeros((tf.size, u.size)); D[np.arange(tf.size), code] = 1
    Xf = np.column_stack([cov, D])
    beta = np.linalg.lstsq(Xf, score, rcond=None)[0]
    return score - Xf @ beta + (D @ beta[cov.shape[1]:])   # remove covariates only, keep TF means (irrelevant within TF)


rows = []; covchk = {}
for ctx in CTXS:
    d = R3 / "cache" / ctx
    genes = json.loads((d / "obs_genes.json").read_text()); obs_idx = np.load(d / "obs_idx.npy")
    P = np.load(d / "h3_pairs.npz"); tf, tg, y = P["tf"].astype(int), P["tg"].astype(int), P["y"].astype(int)
    nobs = np.load(d / "ncells.npy").astype(float)
    Ex = np.load(d / "expr_extraction.npy").astype(np.float64)
    Ec = Ex - Ex.mean(0); sd = np.sqrt((Ec ** 2).mean(0))
    r = (Ec[:, tf] * Ec[:, tg]).mean(0) / (sd[tf] * sd[tg])
    cov = np.column_stack([np.log(np.abs(r) + 1e-3), np.log(Ex.mean(0))[tg], np.log(nobs)[tg],
                           np.array([deg.get(genes[i], 0) for i in tg], float),
                           np.array([len(g2g.get(genes[i], ())) for i in tg], float)])
    # compare with E5's covariates via its cached STRING degree json + its formula (reading e5 outputs only)
    sdj = json.loads((E5 / "string_physical700_degree_genomewide.json").read_text())
    degE5 = np.array([sdj["degree"].get(genes[i], 0) for i in tg], float)
    covchk[ctx] = dict(n_pairs=int(y.size), n_pos=int(y.sum()), n_tf=int(np.unique(tf).size),
                       n_targets=int(np.unique(tg).size),
                       max_abs_deg_diff_vs_e5_json=float(np.max(np.abs(degE5 - cov[:, 3]))),
                       e5_json_n_pairs=int(sdj["n_pairs"]))
    emb = np.load(EMB / f"{ctx}.npy", mmap_mode="r")
    # zero-row check
    nz = np.linalg.norm(np.asarray(emb[0], float), axis=1) > 0
    covchk[ctx]["obs_idx_equals_nonzero_rows"] = bool(np.array_equal(np.flatnonzero(nz), np.sort(obs_idx)))
    for L in range(12):
        M = np.asarray(emb[L], np.float64)[obs_idx]; M = M - M.mean(0)
        M = M / np.linalg.norm(M, axis=1, keepdims=True)
        s = (M[tf] * M[tg]).sum(1)
        res = resid_dummies(s, cov, tf)
        rows.append(dict(context=ctx, layer=L, auroc_raw_verify=auroc_within_tf(s, y, tf),
                         auroc_resid_verify=auroc_within_tf(res, y, tf),
                         corr_score_cov_max_abs=float(np.max(np.abs([np.corrcoef(s, cov[:, j])[0, 1] for j in range(5)])))))
    print(ctx, "done", flush=True)
    np.save(V / f"cov_{ctx}.npy", cov)

df = pd.DataFrame(rows)
ra = pd.read_csv(E5 / "residual_auroc.csv")
df = df.merge(ra[["context", "layer", "auroc_raw", "auroc_resid"]], on=["context", "layer"])
e4 = pd.concat([pd.read_csv(R3 / f"E4a_edges/h3_primary_{c}.csv") for c in CTXS])
e4 = e4[e4.scorer.str.contains("PRIMARY")][["context", "layer", "auroc"]].rename(columns={"auroc": "auroc_e4a"})
df = df.merge(e4, on=["context", "layer"])
df["d_raw_vs_e5"] = df.auroc_raw_verify - df.auroc_raw; df["d_resid_vs_e5"] = df.auroc_resid_verify - df.auroc_resid
df["d_raw_vs_e4a"] = df.auroc_raw_verify - df.auroc_e4a
df.to_csv(V / "verify_part1_auroc.csv", index=False)
out["auroc"] = dict(max_abs_d_raw_vs_e5=float(df.d_raw_vs_e5.abs().max()), max_abs_d_resid_vs_e5=float(df.d_resid_vs_e5.abs().max()),
                    max_abs_d_raw_vs_e4a=float(df.d_raw_vs_e4a.abs().max()))
out["covariates"] = covchk

# ---- joint permutation and CIs re-derived from stored draws
jp = {}
for ctx in CTXS:
    Z = [np.load(E5 / "ckpt" / f"{ctx}_L{l:02d}.npz") for l in range(12)]
    obs = np.mean([z["beta_adj"][0] for z in Z]); null = np.column_stack([z["perm_adj"][:, 0] for z in Z]).mean(1)
    b = int((null >= obs).sum()); bm = np.column_stack([z["boot_adj"][:, 0] for z in Z]).mean(1)
    lo, hi = np.exp(np.quantile(bm, [0.025, 0.975]))
    per = [dict(layer=l, or_lo=float(np.exp(np.quantile(z["boot_adj"][:, 0], 0.025))),
                or_hi=float(np.exp(np.quantile(z["boot_adj"][:, 0], 0.975))),
                n_nonconv=int((~z["boot_conv"]).sum() + (~z["perm_conv"]).sum()),
                any_nan=bool(np.isnan(z["boot_adj"]).any() or np.isnan(z["perm_adj"]).any())) for l, z in enumerate(Z)]
    jp[ctx] = dict(T_obs=float(obs), or_T=float(np.exp(obs)), b=b, p=(b + 1) / (null.size + 1), N=int(null.size),
                   z=float((obs - null.mean()) / null.std(ddof=1)), or_T_boot=[float(lo), float(hi)],
                   min_or_lo=min(p["or_lo"] for p in per), total_nonconv=sum(p["n_nonconv"] for p in per),
                   any_nan=any(p["any_nan"] for p in per))
out["joint_from_ckpt"] = jp

# ---- stratum weights (own)
sw = {}
for ctx in CTXS:
    P = np.load(R3 / "cache" / ctx / "h3_pairs.npz"); tf, y = P["tf"], P["y"].astype(int)
    u = np.unique(tf); k = np.array([y[tf == t].sum() for t in u], float); n = np.array([(tf == t).sum() for t in u], float)
    w = k * (n - k); wc = w / n
    sw[ctx] = dict(top5_auroc=float(np.sort(w)[::-1][:5].sum() / w.sum()), eff_auroc=float(w.sum() ** 2 / (w ** 2).sum()),
                   top5_clogit=float(np.sort(wc)[::-1][:5].sum() / wc.sum()), eff_clogit=float(wc.sum() ** 2 / (wc ** 2).sum()))
out["stratum_weights"] = sw
(V / "verify_part1.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
