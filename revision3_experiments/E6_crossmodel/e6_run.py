"""E6 step 2: every estimator for every embedding of one gene set.

usage: python e6_run.py SET [EMB ...]
  SET  tierA | tierB_<ctx> | legacy186 | legacy178
  EMB  embedding ids (default: all for the set). Static: model keys; contextual: ctx@L<k>;
       projections (tierA only): esm2_uce@rp1024s<k>, scgpt_static@rp102s<k>
Writes work/<SET>/<EMB>.json (rows) and .npz (bootstrap draws, max-T nulls). Skips finished embeddings.

Estimators
  full     basis-free: rows of the gene set, mean-centred over the gene set, full dimension (primary, (i))
  SV2-SV7, SV1-SV6   spectral (ii): U*S of the SVD of the model's whole vocabulary (static) or of V (scGPT
           contextual, stats3.svd_observed), restricted to the gene set; no re-centring
  bestband (iii): best of the 78 contiguous bands of SV1..SV12, Westfall-Young max-T over the 78 bands
           (stats3.max_t_adjust) with N = 1,000 joint null draws
Tasks
  node_lambert   Lambert TF vs every other gene of the set; LOO mean cosine to the TF set (stats3.tf_scores_loo)
  node_trrust    TRRUST TF vs TRRUST target-only (global TRRUST labels, as in the submitted Table 1)
  ppi            STRING >= 700 pairs within the set vs one Maslov-Sneppen rewiring (true pairs removed); pair cosine
  edge           within-TF stratified AUROC (stats3.StratifiedPairs) of signed pair cosine, TRRUST edges vs
                 prep_pairs-rule candidate non-edges
"""
import json, sys, time, zlib
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import e6lib as L
S = L.S

B_NODE, N_NODE = 1000, 2000
B_PPI, N_PPI = 200, 1000
B_EDGE_FULL, B_EDGE_BAND, N_EDGE = 2000, 1000, 2000
N_MAXT = 1000


def sd(setname, task):
    return L.SEED + zlib.crc32(f"{setname}|{task}".encode()) % 100000


def embeddings_for(setname):
    if setname == "tierA":
        return L.STATIC_KEYS + [f"esm2_uce@rp1024s{k}" for k in range(5)] + [f"scgpt_static@rp102s{k}" for k in range(5)]
    if setname.startswith("tierB_"):
        ctx = setname[6:]; return L.STATIC_KEYS + [f"{ctx}@L{l}" for l in range(12)]
    if setname == "legacy186":
        return [k for k in L.STATIC_KEYS if k != "tgpt"] + [f"kidney_A@L{l}" for l in range(12)]
    if setname == "legacy178":
        return L.STATIC_KEYS + [f"kidney_A@L{l}" for l in range(12)]
    raise ValueError(setname)


def provide(emb_id, genes):
    """-> dict(E=set rows centred over the set, Ev=set rows centred over vocabulary/V, C=12 spectral coords or None,
    d, label, kind, layer, model)"""
    if "@L" in emb_id:
        ctx, l = emb_id.split("@L"); l = int(l)
        og, oi = L.ctx_obs(ctx); emb = L.ctx_emb(ctx); pos = {g: k for k, g in enumerate(og)}
        r = np.array([pos[g] for g in genes])
        Xv = np.asarray(emb[l], dtype=np.float64)[oi]
        US, Sv, _ = S.svd_observed(emb, l, oi)
        Ev = (Xv - Xv.mean(0))[r]; E = Xv[r] - Xv[r].mean(0)
        return dict(E=E, Ev=Ev, C=US[r, :12], d=512, label=f"scGPT contextual ({ctx}) L{l}", kind="contextual residual stream",
                    layer=l, model="scgpt_contextual", ctx=ctx)
    if "@rp" in emb_id:
        key, rest = emb_id.split("@rp"); k, s = rest.split("s"); k, s = int(k), int(s)
        g, Ef = L.load_static(key); pos = {x: i for i, x in enumerate(g)}
        Xr = Ef[[pos[x] for x in genes]].astype(np.float64)
        R = np.random.default_rng(1000 + s).standard_normal((Ef.shape[1], k)) / np.sqrt(k)
        Z = Xr @ R
        return dict(E=Z - Z.mean(0), Ev=None, C=None, d=k, label=f"{L.STATIC[key][1]} random projection {Ef.shape[1]}->{k} seed {s}",
                    kind="random projection", layer=-1, model=key, ctx="")
    key = emb_id
    g, Ef = L.load_static(key); pos = {x: i for i, x in enumerate(g)}
    r = np.array([pos[x] for x in genes])
    mu = Ef.astype(np.float64).mean(0)
    Xr = Ef[r].astype(np.float64)
    gg, US, _ = L.static_spectral(key); assert gg == g
    return dict(E=Xr - Xr.mean(0), Ev=Xr - mu, C=US[r], d=Ef.shape[1], label=L.STATIC[key][1], kind=L.STATIC[key][3],
                layer=-1, model=key, ctx="")


# ------------------------------------------------------------------------------------- node
def node_block(Eest, tf_pos, neg_pos, seed, setname, task, est, rows_meta, do_check):
    s = S.tf_scores_loo(Eest, tf_pos); auc = S.auroc(s[tf_pos], s[neg_pos])
    boot, c0 = L.node_boot(Eest, tf_pos, neg_pos, B_NODE, seed)
    if do_check: L.check_node_boot(Eest, tf_pos, neg_pos, c0, boot[0])
    nul = L.fast_tf_auroc_null(Eest, tf_pos, neg_pos, N_NODE, np.random.default_rng(seed + 1))
    p, b, N = S.emp_p(auc, nul)
    plo, phi = L.ci(boot); shift = auc - float(np.mean(boot))
    lo, hi = max(plo + shift, 0.0), min(phi + shift, 1.0)
    ns = L.null_summary(nul)
    return dict(**rows_meta, task=task, estimator=est, metric="AUROC (LOO mean cosine to TF set)",
                n_pos=int(tf_pos.size), n_neg=int(neg_pos.size), value=auc, ci_lo=lo, ci_hi=hi,
                ci_method=f"gene bootstrap B={B_NODE}, leave-gene-out centroid, percentile re-centred on estimate",
                ci_lo_pct=plo, ci_hi_pct=phi, boot_mean=float(np.mean(boot)), boot_sd=float(np.std(boot, ddof=1)),
                p=p, p_b=b, p_N=N, p_method="label permutation (stats3.tf_auroc_null draws; e6lib.fast_tf_auroc_null)",
                null_mean=ns["mean"], null_sd=ns["sd"], null_q975=ns["q975"]), boot


def node_maxt(C, tf_pos, neg_pos, seed):
    """t_obs (78,), t_null (N, 78) on joint label permutations."""
    lab = np.concatenate([tf_pos, neg_pos]); k = tf_pos.size; n = lab.size
    t_obs = np.empty(len(L.BANDS)); t_null = np.empty((N_MAXT, len(L.BANDS)))
    # the permutations tf_auroc_null draws with this seed, regenerated for the 1-D (tie-heavy) bands
    rng = np.random.default_rng(seed); perms = []
    done = 0
    while done < N_MAXT:
        b = min(500, N_MAXT - done); idx = np.argsort(rng.random((n, b)), axis=0)[:k]; perms.append(idx); done += b
    perms = np.concatenate(perms, 1)                                     # k x N (positions within lab)
    for j, (a, bb) in enumerate(L.BANDS):
        Cb = C[:, a - 1:bb]
        s = S.tf_scores_loo(Cb, tf_pos); t_obs[j] = S.auroc(s[tf_pos], s[neg_pos])
        if bb > a:
            t_null[:, j] = L.null_from_perms(Cb, tf_pos, neg_pos, perms)
            if j == 1:
                ref = L.fast_tf_auroc_null(Cb, tf_pos, neg_pos, 500, np.random.default_rng(seed))
                assert np.max(np.abs(ref - t_null[:500, j])) < 1e-12, "null_from_perms mismatch"
        else:
            Cl = Cb[lab]
            for r in range(N_MAXT):
                m = np.zeros(n, bool); m[perms[:, r]] = True
                s = S.tf_scores_loo(Cl, perms[:, r]); t_null[r, j] = S.auroc(s[m], s[~m])
    # consistency: fast null == stats3.tf_auroc_null; loop path == fast path on a continuous band
    L.check_fast_null(C[:, 0:2], tf_pos, neg_pos, seed)
    Cb = C[:, 0:2][lab]
    for r in range(2):
        m = np.zeros(n, bool); m[perms[:, r]] = True
        s = S.tf_scores_loo(Cb, perms[:, r]); assert abs(S.auroc(s[m], s[~m]) - t_null[r, 1]) < 1e-9, "node perm mismatch"
    return t_obs, t_null


# ------------------------------------------------------------------------------------- ppi
def ppi_rows(est_mats, gs, n, seed, rows_meta, C, do_bands=True):
    pi, pj, ni, nj = gs["pi"], gs["pj"], gs["ni"], gs["nj"]
    out = []; obs = {}; negsorted = {}
    if pi.size < 3:
        return out, None, None
    for est, Eest in est_mats.items():
        En = S.l2norm_rows(np.asarray(Eest, float))
        cp = L.pair_cos(En, pi, pj); cn = L.pair_cos(En, ni, nj)
        auc = S.auroc(cp, cn); obs[est] = (En, auc); negsorted[est] = np.sort(cn)
        assert abs(L.auroc_sorted(cp, negsorted[est]) - auc) < 1e-9

        def stat(nodes, a, b, cp=cp, cn=cn):
            keep = np.zeros(n, bool); keep[nodes] = True
            x = cp[keep[pi] & keep[pj]]; y = cn[keep[ni] & keep[nj]]
            return S.auroc(x, y) if x.size > 2 and y.size > 2 else np.nan
        full, lo, hi, draws = S.vertex_subsample_ci(stat, n, pi, pj, B_PPI, np.random.default_rng(seed + 2))
        out.append(dict(**rows_meta, task="ppi", estimator=est, metric="AUROC pair cosine, STRING-physical>=700 vs MS-rewired non-pairs",
                        n_pos=int(pi.size), n_neg=int(ni.size), value=auc, ci_lo=lo, ci_hi=hi,
                        ci_method=f"vertex subsampling m=n/2 B={B_PPI}", boot_sd=float(np.nanstd(draws, ddof=1))))
    if C is None or not do_bands:
        return out, None, None
    # null: N_MAXT further rewirings used as pseudo-positives (true pairs removed), same for every band / estimator
    t_obs = t_null = None
    if C is not None:
        bc_neg = L.band_cos(C, ni, nj); neg_b = [np.sort(bc_neg[:, j]) for j in range(len(L.BANDS))]
        bc_pos = L.band_cos(C, pi, pj)
        t_obs = np.array([S.auroc(bc_pos[:, j], bc_neg[:, j]) for j in range(len(L.BANDS))])
        t_null = np.empty((N_MAXT, len(L.BANDS)))
    est_null = {e: np.empty(N_MAXT) for e in est_mats}
    draws = L.ms_draws_cached(setname_g[0], pi, pj, n, N_MAXT, seed + 3)
    one_d = [lo == hi for lo, hi in L.BANDS]
    for r, (a, b) in enumerate(draws):
        bc = L.band_cos(C, a, b)
        for j in range(len(L.BANDS)): t_null[r, j] = L.auroc_sorted_fast(bc[:, j], neg_b[j], one_d[j])
    if C is not None:
        for e, (a_, b_) in L.FIXED_BANDS.items():
            if e in est_null: est_null[e] = t_null[:, L.BANDS.index((a_, b_))]
    for row in out:
        if C is not None and row["estimator"] == "full": continue
        nul = est_null[row["estimator"]]; p, b, N = S.emp_p(row["value"], nul); ns = L.null_summary(nul)
        row.update(p=p, p_b=b, p_N=N, p_method=f"{N_MAXT} further MS rewirings as pseudo-positives",
                   null_mean=ns["mean"], null_sd=ns["sd"], null_q975=ns["q975"])
    return out, t_obs, t_null


# ------------------------------------------------------------------------------------- edge
def edge_rows(est_mats, gs, seed, rows_meta, C, do_check):
    if gs["e_y"].size == 0 or gs["e_y"].sum() < 3: return [], None, None
    sp = S.StratifiedPairs(gs["e_tf"], gs["e_tg"], gs["e_y"])
    out = []
    Yp = L.edge_perms_cached(setname_g[0], sp, N_EDGE, seed + 4)
    for est, Eest in est_mats.items():
        En = S.l2norm_rows(np.asarray(Eest, float)); s = L.pair_cos(En, sp.tf, sp.tg)
        auc = sp.auroc(s); fs = L.FastStrat(sp, s)
        Bn = B_EDGE_FULL if est == "full" else B_EDGE_BAND
        rngw = np.random.default_rng(seed + 5); boot = []
        for c in range(0, Bn, 250):
            W = L.two_way_weights(sp, min(250, Bn - c), rngw)
            if do_check and c == 0: L.check_faststrat(sp, s, W, Yp)
            boot.append(fs.auroc_w(W))
        boot = np.concatenate(boot); del W
        nul = np.concatenate([fs.perm_null(Yp[:, c:c + 250]) for c in range(0, N_EDGE, 250)])
        p, b, N = S.emp_p(auc, nul); ns = L.null_summary(nul); lo, hi = L.ci(boot)
        out.append(dict(**rows_meta, task="edge", estimator=est, metric="within-TF stratified AUROC, signed cosine",
                        n_pos=int(sp.y.sum()), n_neg=int((~sp.y).sum()), n_tf=int(sp.tfs.size), value=auc,
                        ci_lo=lo, ci_hi=hi, ci_method=f"two-way TF x target bootstrap B={Bn}",
                        boot_sd=float(np.std(boot, ddof=1)), upper95_one_sided=float(np.quantile(boot, 0.95)),
                        p=p, p_b=b, p_N=N, p_method="within-TF label permutation",
                        null_mean=ns["mean"], null_sd=ns["sd"], null_q975=ns["q975"]))
    t_obs = t_null = None
    if C is not None:
        bc = L.band_cos(C, sp.tf, sp.tg); Y = Yp[:, :N_MAXT]
        t_obs = np.empty(len(L.BANDS)); t_null = np.empty((N_MAXT, len(L.BANDS)))
        for j in range(len(L.BANDS)):
            fs = L.FastStrat(sp, bc[:, j]); t_obs[j] = sp.auroc(bc[:, j])
            t_null[:, j] = np.concatenate([fs.perm_null(Y[:, c:c + 250]) for c in range(0, N_MAXT, 250)])
    return out, t_obs, t_null


def maxt_row(rows_meta, task, t_obs, t_null, metric):
    z, pa, pr = S.max_t_adjust(t_obs, t_null)
    jb = int(np.argmax(t_obs)); jz = int(np.argmax(z))
    return dict(**rows_meta, task=task, estimator="bestband", metric=metric + " (best of 78 bands SV1..SV12)",
                value=float(t_obs[jb]), best_band=L.BAND_NAMES[jb], best_band_z=float(z[jb]), p=float(pa[jb]),
                p_raw=float(pr[jb]), p_N=int(t_null.shape[0]), p_method="Westfall-Young max-T over 78 bands",
                best_z_band=L.BAND_NAMES[jz], best_z=float(z[jz]), p_best_z=float(pa[jz]),
                n_bands_padj05=int((pa < 0.05).sum()))


setname_g = [None]
SKIP_MAXT_STATIC_TIERB = True


def run(setname, emb_id):
    setname_g[0] = setname
    out = L.WORK / setname; out.mkdir(exist_ok=True)
    fj = out / f"{emb_id}.json"
    if fj.exists(): return
    t0 = time.time()
    gs = L.load_geneset(setname); genes = gs["genes"]; n = len(genes)
    P = provide(emb_id, genes)
    meta = dict(set=setname, embedding=emb_id, model=P["model"], label=P["label"], kind=P["kind"], ctx=P["ctx"],
                layer=P["layer"], d=P["d"], n_genes=n)
    est = {"full": P["E"]}
    if setname.startswith("tierB_") and P["model"] != "scgpt_contextual" and SKIP_MAXT_STATIC_TIERB:
        P["C_maxt"] = None          # compute budget: estimator (iii) only for scGPT contextual in Tier B
    else:
        P["C_maxt"] = P["C"]
    if P["C"] is not None:
        for nm, (a, b) in L.FIXED_BANDS.items(): est[nm] = P["C"][:, a - 1:b]
    small = n < 6000
    rows = []; npz = {}
    tasks = {"node_lambert": (np.flatnonzero(gs["lambert"]), np.flatnonzero(~gs["lambert"])),
             "node_trrust": (np.flatnonzero(gs["trrust_tf"]), np.flatnonzero(gs["trrust_tgt"]))}
    import pickle
    def ck(task):
        return out / f"{emb_id}.part_{task}.pkl"
    def ck_load(task):
        f = ck(task)
        if f.exists():
            r_, z_ = pickle.loads(f.read_bytes()); rows.extend(r_); npz.update(z_); return True
        return False
    def ck_save(task, r0, z0):
        f = ck(task); tmp = f.with_suffix(".tmp"); tmp.write_bytes(pickle.dumps((rows[r0:], {k: npz[k] for k in npz if k not in z0}))); tmp.rename(f)
    for task, (tp, ng) in tasks.items():
        if tp.size < 3 or ng.size < 3: continue
        if ck_load(task): continue
        r0, z0 = len(rows), set(npz)
        seed = sd(setname, task)
        for e, Eest in est.items():
            row, boot = node_block(Eest, tp, ng, seed, setname, task, e, meta, do_check=(small or e != "full"))
            if e == "full":
                if P["Ev"] is not None:
                    sv = S.tf_scores_loo(P["Ev"], tp); row["value_vocab_centred"] = S.auroc(sv[tp], sv[ng])
                if task == "node_lambert":
                    sf = S.tf_scores_loo(Eest, tp); isp = np.zeros(n, bool); isp[tp] = True
                    row["value_detection_stratified"] = S.stratified_auroc(sf, isp, gs["det_dec"])
            rows.append(row); npz[f"boot_{task}_{e}"] = boot.astype(np.float32)
        if P["C_maxt"] is not None:
            to, tn = node_maxt(P["C"], tp, ng, seed + 10)
            rows.append(maxt_row(meta, task, to, tn, "AUROC (LOO mean cosine to TF set)"))
            npz[f"tobs_{task}"] = to; npz[f"tnull_{task}"] = tn.astype(np.float32)
        ck_save(task, r0, z0)
        print(f"  {emb_id} {task} done {time.time()-t0:.0f}s", flush=True)
    if not ck_load("ppi"):
        r0, z0 = len(rows), set(npz)
        r, to, tn = ppi_rows(est, gs, n, sd(setname, "ppi"), meta, P["C"], P["C_maxt"] is not None); rows += r
        if to is not None:
            rows.append(maxt_row(meta, "ppi", to, tn, "pair-cosine AUROC")); npz["tobs_ppi"] = to; npz["tnull_ppi"] = tn.astype(np.float32)
        ck_save("ppi", r0, z0)
    print(f"  {emb_id} ppi done {time.time()-t0:.0f}s", flush=True)
    r, to, tn = edge_rows(est, gs, sd(setname, "edge"), meta, P["C_maxt"], do_check=True); rows += r
    if to is not None:
        rows.append(maxt_row(meta, "edge", to, tn, "within-TF stratified AUROC")); npz["tobs_edge"] = to; npz["tnull_edge"] = tn.astype(np.float32)
    np.savez(out / f"{emb_id}.npz", **npz)
    fj.write_text(json.dumps(dict(rows=rows, seconds=time.time() - t0), default=float))
    for f in out.glob(f"{emb_id}.part_*.pkl"): f.unlink()
    print(f"{setname} {emb_id} finished in {time.time()-t0:.0f}s", flush=True)


def job_list():
    big = ["esm2_uce", "scgpt_static", "gf_v2_316m", "tgpt", "esm2_650m", "aido_100m", "gf_v2_104m",
           "scfoundation", "gf_v1_10m", "hpa_coexp", "seq_features"]
    jobs = []
    for s_ in ["legacy186", "legacy178"]: jobs += [(s_, e) for e in embeddings_for(s_)]
    jobs += [("tierA", e) for e in big]
    for s_ in ["tierB_lung", "tierB_lung_external"]: jobs += [(s_, e) for e in embeddings_for(s_)]
    jobs += [("tierA", e) for e in embeddings_for("tierA") if "@rp" in e]
    for s_ in ["tierB_kidney_B", "tierB_immune", "tierB_kidney_A"]: jobs += [(s_, e) for e in embeddings_for(s_)]
    return jobs


def worker():
    import os
    for setname, e in job_list():
        out = L.WORK / setname; out.mkdir(exist_ok=True)
        if (out / f"{e}.json").exists(): continue
        try:
            fd = os.open(out / f"{e}.lock", os.O_CREAT | os.O_EXCL | os.O_WRONLY); os.close(fd)
        except FileExistsError:
            continue
        run(setname, e)


def prep(sets):
    import os
    for setname in sets:
        gs = L.load_geneset(setname); n = len(gs["genes"])
        if gs["pi"].size >= 3:
            L.ms_draws_cached(setname, gs["pi"], gs["pj"], n, N_MAXT, sd(setname, "ppi") + 3)
        if gs["e_y"].size and gs["e_y"].sum() >= 3:
            sp = S.StratifiedPairs(gs["e_tf"], gs["e_tg"], gs["e_y"])
            L.edge_perms_cached(setname, sp, N_EDGE, sd(setname, "edge") + 4)
        print("prepped", setname, flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "prep":
        prep(sys.argv[2:])
    elif sys.argv[1] == "worker":
        worker()
    else:
        setname = sys.argv[1]; embs = sys.argv[2:] or embeddings_for(setname)
        for e in embs:
            run(setname, e)
