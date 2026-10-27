"""Independent verification of E3_node headline numbers (part 1). Helper functions in verify_common.py (no stats3 / e3_lib)."""
from verify_common import *
out = {}; rows = []
for ctx in ("lung", "lung_external"):
    emb, obs, genes, cnt, y, det = ctx_data(ctx)
    log(ctx, "V", len(genes), "TF", int(y.sum()))
    Xs = []
    for l in range(12):
        X = np.asarray(emb[l], float)[obs]; X = X - X.mean(0); Xs.append(X)
    s_l = [loo(X, y) for X in Xs]
    a = [strat(s, y, det) for s in s_l]; u = [mw(s, y) for s in s_l]
    ref = pd.read_csv(E3 / f"h2_primary_{ctx}.csv"); ref = ref[ref.row_type == "H2_layer"]
    for l in range(12):
        rows.append(dict(context=ctx, quantity=f"H2 strat AUROC L{l}", mine=a[l], reported=float(ref.strat_auroc.iloc[l])))
        rows.append(dict(context=ctx, quantity=f"H2 unstrat AUROC L{l}", mine=u[l], reported=float(ref.auroc.iloc[l])))
    # small independent joint within-stratum permutation null (N=300) -> z check
    rng = np.random.default_rng(99); NN = 300; nul = np.empty((NN, 12))
    Ens = [X / np.linalg.norm(X, axis=1)[:, None] for X in Xs]
    for b in range(NN):
        yp = y.copy()
        for st in np.unique(det):
            idx = np.flatnonzero(det == st); yp[idx] = yp[rng.permutation(idx)]
        for l in range(12):
            v = Ens[l][yp].sum(0); s = Ens[l] @ v; k = yp.sum(); s = np.where(yp, (s - 1) / (k - 1), s / k)
            nul[b, l] = strat(s, yp, det)
    zmine = (np.array(a) - nul.mean(0)) / nul.std(0, ddof=1)
    for l in range(12):
        rows.append(dict(context=ctx, quantity=f"z L{l} (own N=300 null)", mine=float(zmine[l]), reported=float(ref.z.iloc[l])))
    out[ctx] = dict(null_mean=nul.mean(0).tolist(), null_sd=nul.std(0, ddof=1).tolist(),
                    max_null_over_layers=float(nul.max()), min_obs=float(min(a)))
    log(ctx, "primary + null done")
    # baselines: log mean expression (extraction cells)
    ee = np.load(R3 / "cache" / ctx / "expr_extraction.npy").astype(float)
    lme = np.log(ee.mean(0) + 1e-3)
    bl = pd.read_csv(E3 / "baselines.csv"); bl = bl[bl.context == ctx].set_index("baseline")
    rows.append(dict(context=ctx, quantity="log mean expr (extraction) strat AUROC", mine=strat(lme, y, det),
                     reported=float(bl.loc["log_mean_expr_extraction", "strat_auroc"])))
    rows.append(dict(context=ctx, quantity="log n_obs unstrat AUROC", mine=mw(np.log(np.maximum(cnt, 1)), y),
                     reported=float(bl.loc["log_nobs", "auroc"])))
    # static embeddings
    for nm, path, key in (("esm2", ESM, "static_esm2"), ("scgpt", SCG, "static_scgpt")):
        z = np.load(path); gn = [str(g) for g in z["gene_names"]]; first = {}
        for i, g in enumerate(gn): first.setdefault(g, i)
        iv = np.array([k for k, g in enumerate(genes) if g in first]); E = np.asarray(z["embeddings"][[first[genes[k]] for k in iv]], float)
        keep = np.linalg.norm(E, axis=1) > 0; iv = iv[keep]; E = E[keep]; E = E - E.mean(0)
        ys, ds = y[iv], det[iv]; s = loo(E, ys)
        rows.append(dict(context=ctx, quantity=f"{key} strat AUROC (n={iv.size}, dim={E.shape[1]})", mine=strat(s, ys, ds),
                         reported=float(bl.loc[key, "strat_auroc"])))
        bp = pd.read_csv(E3 / "baselines_paired.csv"); bp = bp[(bp.context == ctx) & (bp.static == nm)]
        for l in (0, 5, 11):
            sl = loo(Xs[l][iv], ys)
            rows.append(dict(context=ctx, quantity=f"{key} minus layer L{l} delta", mine=strat(s, ys, ds) - strat(sl, ys, ds),
                             reported=float(bp[bp.layer == l].delta.iloc[0])))
    log(ctx, "static done")
    # band cell from own SVD: best cos cell and SV8 L0
    bs = pd.read_csv(E3 / f"bandscan_{ctx}.csv"); bs = bs[bs.row_type == "band_cell"]
    for (l, i, j) in ((2, 3, 12), (0, 8, 8), (2, 2, 7), (5, 5, 7)):
        U, Sv, Vt = np.linalg.svd(Xs[l], full_matrices=False); US = U * Sv
        B = US[:, i - 1:j]
        val = strat(B[:, 0], y, det) if i == j else strat(loo(B, y), y, det)
        if i == j: val = 0.5 + abs(val - 0.5)
        band = f"SV{i}-{j}" if j > i else f"SV{i}"
        rows.append(dict(context=ctx, quantity=f"band {band} L{l} cos strat AUROC", mine=val,
                         reported=float(bs[(bs.layer == l) & (bs.band == band)].strat_auroc_cos.iloc[0])))
    # logistic adjustment at L5 (own Newton-IRLS)
    def logit_ll(Xd, yy):
        b = np.zeros(Xd.shape[1])
        for _ in range(50):
            p = 1 / (1 + np.exp(-Xd @ b)); W = p * (1 - p)
            H = Xd.T @ (Xd * W[:, None]); g = Xd.T @ (yy - p); step = np.linalg.solve(H, g); b += step
            if np.abs(step).max() < 1e-10: break
        p = 1 / (1 + np.exp(-Xd @ b)); ll = float(np.sum(yy * np.log(p) + (1 - yy) * np.log(1 - p)))
        se = np.sqrt(np.diag(np.linalg.inv(Xd.T @ (Xd * (p * (1 - p))[:, None]))))
        return b, se, ll
    zs = lambda v: (v - v.mean()) / v.std()
    lg = pd.read_csv(E3 / "logistic_adjustment.csv"); lg = lg[(lg.context == ctx) & (lg.expr_covariate == "extraction")]
    for l in (5, 11):
        ones = np.ones(y.size); yy = y.astype(float)
        Xr = np.column_stack([ones, zs(np.log(np.maximum(cnt, 1))), zs(lme)]); Xf = np.column_stack([ones, zs(s_l[l]), Xr[:, 1:]])
        bf, sef, llf = logit_ll(Xf, yy); _, _, llr = logit_ll(Xr, yy)
        r = lg[lg.layer == l].iloc[0]
        rows.append(dict(context=ctx, quantity=f"logit beta(score) L{l}", mine=float(bf[1]), reported=float(r.beta_score_per_sd)))
        rows.append(dict(context=ctx, quantity=f"logit LR p L{l}", mine=float(chi2.sf(2 * (llf - llr), 1)), reported=float(r.LR_p)))
    # restricted V (n_obs >= 10), L0 and L11
    rv = pd.read_csv(E3 / "restricted_V.csv"); rv = rv[(rv.context == ctx) & (rv.row_type == "H2_layer") & (rv.min_nobs == 10)]
    sub = cnt >= 10; ys = y[sub]; ds = deciles(np.log(np.maximum(cnt[sub], 1)))
    for l in (0, 11):
        X = np.asarray(emb[l], float)[obs][sub]; X = X - X.mean(0)
        rows.append(dict(context=ctx, quantity=f"restricted n_obs>=10 L{l}", mine=strat(loo(X, ys), ys, ds),
                         reported=float(rv[rv.layer == l].strat_auroc.iloc[0])))
    log(ctx, "bands/logit/restricted done")

df = pd.DataFrame(rows); df["absdiff"] = (df.mine - df.reported).abs()
df.to_csv(HERE / "verify_e3.csv", index=False)
out["max_absdiff_excluding_z"] = float(df[~df.quantity.str.startswith("z ")].absdiff.max())
(HERE / "verify_e3.json").write_text(json.dumps(out, indent=1))
pd.set_option("display.width", 200); print(df.to_string()); print("max absdiff (non-z)", out["max_absdiff_excluding_z"])
