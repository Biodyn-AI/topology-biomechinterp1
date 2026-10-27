"""E5 continuity row (kidney_A): the submitted version's confound residualisation, then the paper convention.

Rows of continuity_kidneyA.csv (variant column):
  A_submitted_zero_row     exact re-run of revision2_experiments/run_confounds_sv57.py (same RNG stream, same
                           pairs, same covariates, same N = 10,000 per-layer two-sided permutation): SVD of the FULL
                           4,803 x 512 matrix including the 2,831 never-observed all-zero rows, bands SV5-SV7 and
                           SV2-SV7, pooled random negatives (4 per positive, any non-target of the TF among the 187
                           curated genes), covariates = co-expression r, log1p mean expression x2, STRING
                           (functional, >= 0.7) degree sum, GO count sum; pooled rank-biserial of OLS residuals.
  B_observed_svd_identity  identical pairs, covariates, residualisation and permutation scheme; ONLY the SVD changes
                           to the paper convention (stats3.svd_observed: observed genes V, V-centred); same band
                           indices SV5-SV7 / SV2-SV7.
  B_observed_svd_mapped    as B but with the legacy bands mapped to the observed basis by principal angles, fixed
                           from geometry alone (E9_screen/legacy_band_mapping.json, legacy_full_zero:
                           SV5-7 -> SV4-6, SV2-7 -> SV1-6).
  C_paper_benchmark        paper convention end to end: H3 pairs (TRRUST positives, matched TRRUST-negative candidate
                           pairs), the E5 covariates (log|co-expression| on the extraction cells, target log mean
                           expression, target log n_obs, genome-wide STRING-physical degree, GO count; TF-level terms
                           absorbed), within-TF residualisation and within-TF stratified AUROC, r_rb = 2*AUROC - 1,
                           two-way cluster bootstrap CI (B = 500, the E5 draws), for observed-basis bands (identity and
                           mapped) and for the full-dimensional score.
Per-layer permutation p-values are reported only for A/B, under their own (pooled, non-stratified) scheme; the
round-2 Stouffer combination is replaced by a joint permutation (one permutation for all 12 layers, statistic =
mean r_rb over layers, N = 10,000).
"""
import sys, json, time, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import e5_lib as L
S = L.S
sys.path.insert(0, str(HERE.parents[1] / "revision2_experiments"))
import common as C2          # read-only: round-2 loaders and estimators

N_PERM = 10000
BANDS_ID = {"SV5_SV7": slice(4, 7), "SV2_SV7": slice(1, 7)}
BANDS_MAP = {"SV5_SV7": ("SV4_SV6", slice(3, 6)), "SV2_SV7": ("SV1_SV6", slice(0, 6))}


def log(*a): print(time.strftime("%H:%M:%S"), *a, flush=True)


def rb(x, lab):
    r = pd.Series(x).rank().to_numpy()
    return float(2 * (r[lab == 1].mean() - r[lab == 0].mean()) / lab.size)


def main(mode="main"):
    import os
    hashseed = os.environ.get("PYTHONHASHSEED", "unset")
    t0 = time.time()
    mapping = json.loads((L.R3 / "E9_screen/legacy_band_mapping.json").read_text())
    mp = {b["legacy_band"]: b["mapped_observed_band"] for b in mapping["bands"]
          if b["context"] == "kidney_A" and b["legacy_convention"] == "legacy_full_zero"}
    assert mp["SV5-7"] == "SV4-6" and mp["SV2-7"] == "SV1-6", mp
    # ------------------------------------------------------------------ round-2 setup, verbatim RNG stream
    rng = np.random.default_rng(20260819)
    genes = C2.h5ad_genes(); emb = C2.load_emb(); nz0 = C2.nonzero_mask(emb, 0)
    sub = np.array([i for i in C2.curated_genes(genes) if nz0[i]])
    names = [genes[i] for i in sub]; loc = {g: k for k, g in enumerate(names)}
    t = pd.read_csv(C2.TRRUST_PATH, sep="\t", header=None, names=["TF", "target", "mode", "pmid"])
    tr = sorted({(a, b) for a, b in zip(t.TF.astype(str), t.target.astype(str)) if a in loc and b in loc and a != b})
    obj = json.loads(C2.STRING_CACHE.read_text())
    st_all = [(p["g1"], p["g2"], p["score"]) for p in obj["pairs"]]
    gene2go = C2.load_gene2go()
    import h5py
    with h5py.File(C2.H5AD, "r") as f:
        gs = [g.decode() if isinstance(g, bytes) else str(g) for g in f["var"]["_index"][:]]
        gi = {g: i for i, g in enumerate(gs)}
        keep = [g for g in names if g in gi]; cols = [gi[g] for g in keep]
        Xh = f["X"]
        assert isinstance(Xh, h5py.Group)
        data, indices, indptr = Xh["data"][:], Xh["indices"][:], Xh["indptr"][:]
        ncell = indptr.shape[0] - 1
        sel = np.sort(rng.choice(ncell, min(3000, ncell), replace=False))
        cpos = {c: k for k, c in enumerate(cols)}; cset = set(cols)
        E = np.zeros((sel.size, len(cols)), dtype=np.float32)
        for r, ci in enumerate(sel):
            a, b = indptr[ci], indptr[ci + 1]
            for j, v in zip(indices[a:b], data[a:b]):
                if j in cset: E[r, cpos[j]] = v
    kloc = {g: k for k, g in enumerate(keep)}
    Z = (E - E.mean(0)) / (E.std(0) + 1e-9)
    def coexp(a, b):
        if a not in kloc or b not in kloc: return 0.0
        return float(np.dot(Z[:, kloc[a]], Z[:, kloc[b]]) / Z.shape[0])
    deg = {g: 0 for g in names}
    for g1, g2, s in st_all:
        if s >= 0.7:
            if g1 in deg: deg[g1] += 1
            if g2 in deg: deg[g2] += 1
    ngo = {g: len(gene2go.get(g, set())) for g in names}
    mexp = {g: float(E[:, kloc[g]].mean()) if g in kloc else 0.0 for g in names}
    neg = []; trset = set(tr)
    for a in {x for x, _ in tr}:
        for b in names:
            if b != a and (a, b) not in trset: neg.append((a, b))
    neg = [neg[i] for i in rng.choice(len(neg), min(4 * len(tr), len(neg)), replace=False)]
    allp = tr + neg
    y = np.array([1] * len(tr) + [0] * len(neg))
    Xc = np.column_stack([[coexp(a, b) for a, b in allp],
                          [np.log1p(mexp[a]) for a, b in allp], [np.log1p(mexp[b]) for a, b in allp],
                          [deg[a] + deg[b] for a, b in allp], [ngo[a] + ngo[b] for a, b in allp]]).astype(np.float64)
    X1 = np.column_stack([np.ones(len(allp)), Xc])
    log(f"[setup] {len(names)} genes, {len(tr)} positives, {len(neg)} negatives")
    ia = np.array([loc[a] for a, b in allp]); ib = np.array([loc[b] for a, b in allp])
    # observed-basis rows of the 187 genes
    d = L.R3 / "cache" / "kidney_A"
    obs_idx = np.load(d / "obs_idx.npy"); pos_obs = {int(r): k for k, r in enumerate(obs_idx)}
    assert all(int(i) in pos_obs for i in sub)
    sub_obs = np.array([pos_obs[int(i)] for i in sub])
    emb_local = L.load_emb("kidney_A")

    def band_cos(US, cols_):
        U = US[:, cols_]; Un = U / (np.linalg.norm(U, axis=1, keepdims=True) + 1e-12)
        return np.einsum("ij,ij->i", Un[ia], Un[ib])

    def resid(cos):
        beta, *_ = np.linalg.lstsq(X1, cos, rcond=None); return cos - X1 @ beta

    rows = []
    if mode == "spread":
        # only the point estimates of A and B under this interpreter's string-hash seed (set iteration order
        # decides which negatives the round-2 sampler draws)
        for tag, cols_ in BANDS_ID.items():
            for Ly in range(12):
                rows.append(dict(hashseed=hashseed, variant="A_submitted_zero_row", band=tag, layer=Ly,
                                 r_rb=rb(resid(band_cos(C2.layer_us(emb, Ly)[sub], cols_)), y)))
        for vname, bands in [("B_observed_svd_identity", {k: (k, v) for k, v in BANDS_ID.items()}),
                             ("B_observed_svd_mapped", BANDS_MAP)]:
            for legacy, (tag, cols_) in bands.items():
                for Ly in range(12):
                    US, _, _ = S.svd_observed(emb_local, Ly, obs_idx)
                    rows.append(dict(hashseed=hashseed, variant=vname, band=tag, layer=Ly,
                                     r_rb=rb(resid(band_cos(US[sub_obs], cols_)), y)))
        out = HERE / "continuity_hashseed_spread"; out.mkdir(exist_ok=True)
        pd.DataFrame(rows).to_csv(out / f"hashseed_{hashseed}.csv", index=False)
        log("spread", hashseed, "done"); return
    # A: verbatim (consumes the round-2 RNG stream in the original order -> p-values reproduce exactly)
    resA = {}
    for tag, cols_ in [("SV5_SV7", slice(4, 7)), ("SV2_SV7", slice(1, 7))]:
        resA[tag] = []
        for Ly in range(12):
            US = C2.layer_us(emb, Ly)[sub]
            r_ = resid(band_cos(US, cols_)); resA[tag].append(r_)
            obs = rb(r_, y)
            null = np.array([rb(r_, rng.permutation(y)) for _ in range(N_PERM)])
            p, _, _ = C2.emp_p(obs, null, tail="two-sided")
            rows.append(dict(variant="A_submitted_zero_row", band=tag, layer=Ly, r_rb=obs, p_two_sided=p, n_perm=N_PERM,
                             n_pos=int(y.sum()), n_neg=int((1 - y).sum())))
        log("A", tag, [round(r["r_rb"], 4) for r in rows if r["band"] == tag and r["variant"].startswith("A")])
    # B: observed SVD, identical everything else; fresh documented RNG
    resB = {}
    for vcode, vname, bands in [(1, "B_observed_svd_identity", {k: (k, v) for k, v in BANDS_ID.items()}),
                                (2, "B_observed_svd_mapped", BANDS_MAP)]:
        for legacy, (tag, cols_) in bands.items():
            rngB = np.random.default_rng([L.SEED, 7, vcode, int(cols_.start)])
            resB[(vname, tag)] = []
            for Ly in range(12):
                US, _, _ = S.svd_observed(emb_local, Ly, obs_idx)
                r_ = resid(band_cos(US[sub_obs], cols_)); resB[(vname, tag)].append(r_)
                obs = rb(r_, y)
                null = np.array([rb(r_, rngB.permutation(y)) for _ in range(N_PERM)])
                p, _, _ = C2.emp_p(obs, null, tail="two-sided")
                rows.append(dict(variant=vname, band=tag, legacy_band=legacy, layer=Ly, r_rb=obs, p_two_sided=p,
                                 n_perm=N_PERM, n_pos=int(y.sum()), n_neg=int((1 - y).sum())))
            log(vname, tag, "done")
    S.assert_resolution(N_PERM, 0.05, 12)
    # joint permutation (replaces Stouffer) for A and B
    rngJ = np.random.default_rng([L.SEED, 77])
    sets = {("A_submitted_zero_row", k): v for k, v in resA.items()}; sets.update(resB)
    ranks = {k: [pd.Series(r_).rank().to_numpy() for r_ in v] for k, v in sets.items()}
    obsJ = {k: np.mean([rb(r_, y) for r_ in v]) for k, v in sets.items()}
    nullJ = {k: np.empty(N_PERM) for k in sets}
    n = y.size
    for i in range(N_PERM):
        yp = rngJ.permutation(y).astype(bool)
        for k, rk in ranks.items():
            nullJ[k][i] = np.mean([2 * (r[yp].mean() - r[~yp].mean()) / n for r in rk])
    for (vname, tag), o in obsJ.items():
        nl = nullJ[(vname, tag)]
        p2 = (np.sum(np.abs(nl) >= abs(o)) + 1) / (N_PERM + 1)
        rows.append(dict(variant=vname, band=tag, layer="joint_mean_over_layers", r_rb=o, p_two_sided=p2, n_perm=N_PERM,
                         null_mean=nl.mean(), null_sd=nl.std(ddof=1), null_q975=np.quantile(nl, 0.975),
                         n_pos=int(y.sum()), n_neg=int((1 - y).sum())))
    log("joint permutations done")
    # C: paper benchmark (H3 pairs, E5 covariates, within-TF)
    Cc = L.load_context("kidney_A")
    tf, tg, yy = Cc["tf"], Cc["tg"], Cc["y"]
    tfs, tfc, tgs, tgc = L.tf_target_codes(tf, tg)
    sp = S.StratifiedPairs(tf, tg, yy)
    boots = [L.boot_weights("kidney_A", b, tfs.size, tgs.size) for b in range(L.B_BOOT)]
    specs = [("SV5_SV7", slice(4, 7)), ("SV2_SV7", slice(1, 7)), ("SV4_SV6", slice(3, 6)), ("SV1_SV6", slice(0, 6)),
             ("full_512", None)]
    for Ly in range(12):
        US, _, _ = S.svd_observed(emb_local, Ly, obs_idx)
        for tag, cols_ in specs:
            if cols_ is None:
                sc = L.layer_scores(emb_local, obs_idx, Ly, tf, tg)
            else:
                U = US[:, cols_]; Un = U / (np.linalg.norm(U, axis=1, keepdims=True) + 1e-12)
                sc = np.einsum("ij,ij->i", Un[tf], Un[tg])
            r_, _ = L.residualise_within_tf(sc, Cc["cov_raw"], tfc)
            a_raw = sp.auroc(sc); a_res = sp.auroc(r_)
            bs = np.empty(L.B_BOOT)
            for b, (wtf, wtg) in enumerate(boots):
                w = wtf[tfc] * wtg[tgc]
                rb_, _ = L.residualise_within_tf(sc, Cc["cov_raw"], tfc, w); bs[b] = sp.auroc(rb_, w)
            lo, hi = np.nanquantile(bs, [0.025, 0.975])
            rows.append(dict(variant="C_paper_benchmark", band=tag, layer=Ly, r_rb=2 * a_res - 1,
                             auroc_resid=a_res, auroc_raw=a_raw, r_rb_raw=2 * a_raw - 1,
                             r_rb_ci_lo=2 * lo - 1, r_rb_ci_hi=2 * hi - 1, n_boot=L.B_BOOT,
                             n_pos=int(yy.sum()), n_neg=int((1 - yy).sum()), n_tf=int(tfs.size)))
        log("C layer", Ly, "done")
    df = pd.DataFrame(rows)
    df.insert(0, "pythonhashseed", hashseed)
    df.to_csv(HERE / "continuity_kidneyA.csv", index=False)
    repro_check(df, hashseed)


def repro_check(df, hashseed):
    """Compare variant A with the submitted round-2 CSVs (numeric layers only)."""
    chk = {}
    for tag in BANDS_ID:
        old = pd.read_csv(HERE.parents[1] / f"revision2_experiments/M1_confounds_{tag}.csv")
        new = df[(df.variant == "A_submitted_zero_row") & (df.band == tag)
                 & (df.layer.astype(str) != "joint_mean_over_layers")].copy()
        new["layer"] = new.layer.astype(int); new = new.sort_values("layer")
        d = np.abs(old.r_rb.to_numpy() - new.r_rb.to_numpy().astype(float))
        chk[tag] = dict(max_abs_r_rb_diff=float(d.max()), mean_abs_r_rb_diff=float(d.mean()),
                        max_abs_p_diff=float(np.max(np.abs(old.p.to_numpy() - new.p_two_sided.to_numpy().astype(float)))),
                        n_layers_p_lt_005_submitted=int((old.p < 0.05).sum()),
                        n_layers_p_lt_005_rerun=int((new.p_two_sided < 0.05).sum()))
    chk["pythonhashseed"] = hashseed
    chk["note"] = ("round-2 negatives are drawn from a list built by iterating a Python set of TF symbols, whose order "
                   "depends on the per-process string-hash seed; the submitted values are therefore reproducible only "
                   "under the (unrecorded) hash seed of the original run. See continuity_hashseed_spread.csv.")
    (HERE / "continuity_reproduction_check.json").write_text(json.dumps(chk, indent=1))
    log("reproduction check", chk)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "check":      # re-run only the reproduction check on the saved CSV
        _d = pd.read_csv(HERE / "continuity_kidneyA.csv"); repro_check(_d, str(_d.pythonhashseed.iloc[0]))
    else:
        main(sys.argv[1] if len(sys.argv) > 1 else "main")
