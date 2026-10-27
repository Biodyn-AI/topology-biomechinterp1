"""Independent verification of E6 headline numbers (verifier; does NOT import e6lib, e6_run, stats3 or prep_pairs).

Everything is recomputed from the raw files: model .npz files, Lambert TF census, TRRUST tsv, DoRothEA tsvs,
STRING v12 links + protein info, cache/<ctx>/obs_genes.json / obs_idx.npy / ncells.npy / expr.npy, and the local
scGPT residual-stream copies. AUROCs use sklearn.metrics.roc_auc_score / scipy.stats.mannwhitneyu (independent of
stats3.auroc). Leave-one-out mean cosine is re-implemented here.

usage: python verify_e6.py part1|part2|part3      (results appended to verify/verify_results.json)
"""
import gzip, json, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
E6 = HERE.parent; R3 = E6.parent; REPO = R3.parents[1]
CACHE = R3 / "cache"; EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"
NET = REPO / "single_cell_mechinterp/external/networks"
OUT = HERE / "verify_results.json"
MODELS = {
    "scgpt_static": "subproject_53_scgpt_gpl_replication/embeddings/scgpt_gene_embeddings.npz",
    "gf_v2_316m": "subproject_54_geneformer_gpl_replication/embeddings/geneformer_gene_embeddings.npz",
    "gf_v1_10m": "subproject_55_geneformer_v1_10m_gpl_replication/embeddings/geneformer_v1_10m_gene_embeddings.npz",
    "gf_v2_104m": "subproject_56_geneformer_v2_104m_gpl_replication/embeddings/geneformer_v2_104m_gene_embeddings.npz",
    "aido_100m": "subproject_57_aido_cell_100m_gpl_replication/embeddings/aido_cell_100m_gene_embeddings.npz",
    "scfoundation": "subproject_58_scfoundation_gpl_replication/embeddings/scfoundation_gene_embeddings.npz",
    "tgpt": "subproject_59_tgpt_gpl_replication/embeddings/tgpt_gene_embeddings.npz",
    "esm2_uce": "subproject_60_uce_esm2_control_gpl_replication/embeddings/uce_esm2_human_gene_embeddings.npz",
    "esm2_650m": "subproject_61_two_axes_revision/data/derived_embeddings/esm2_650m_own.npz",
    "seq_features": "subproject_61_two_axes_revision/data/derived_embeddings/seq_features.npz",
    "hpa_coexp": "subproject_61_two_axes_revision/data/derived_embeddings/hpa_coexpression.npz",
}


def save(key, val):
    d = json.loads(OUT.read_text()) if OUT.exists() else {}
    d[key] = val; OUT.write_text(json.dumps(d, indent=1, default=float)); print(key, json.dumps(val, default=float)[:600], flush=True)


def load_model(k):
    z = np.load(REPO / MODELS[k]); g = np.array([str(x) for x in z["gene_names"]]); E = np.asarray(z["embeddings"], np.float32)
    nz = np.any(E != 0, axis=1)
    return dict(zip(g[nz], np.flatnonzero(nz))), E


def lambert():
    t = pd.read_csv(REPO / "subproject_61_two_axes_revision/data/humantfs_v1.01.txt", sep="\t")
    return set(t.loc[t["Is TF?"].astype(str).str.strip() == "Yes", "HGNC symbol"].astype(str))


def trrust():
    t = pd.read_csv(NET / "trrust_human.tsv", sep="\t", header=None, names=["TF", "target", "mode", "pmid"]).astype(str)
    return t


def loo_scores(X, is_tf):
    """centre rows of X, unit-normalise, leave-one-out mean cosine to the TF rows."""
    X = np.asarray(X, np.float64); X = X - X.mean(0)
    Xn = X / np.linalg.norm(X, axis=1, keepdims=True)
    c = Xn[is_tf].sum(0); s = Xn @ c; k = is_tf.sum()
    return np.where(is_tf, (s - 1.0) / (k - 1), s / k)


def strat_auc(s, y, strata):
    num = den = 0.0
    for st in np.unique(strata):
        m = strata == st; p = y[m].sum(); q = m.sum() - p
        if p and q:
            num += roc_auc_score(y[m], s[m]) * p * q; den += p * q
    return num / den


def my_deciles(x):
    """deciles of x by rank (own convention): stratum = floor(10 * (average rank - 0.5) / n)."""
    r = rankdata(x, method="average"); return np.minimum((10 * (r - 0.5) / len(x)).astype(int), 9)


def ctx_obs(ctx):
    return json.loads((CACHE / ctx / "obs_genes.json").read_text()), np.load(CACHE / ctx / "obs_idx.npy")


def part1():
    t0 = time.time()
    voc = {}; mats = {}
    for k in MODELS:
        pos, E = load_model(k); voc[k] = pos; mats[k] = E
    A = sorted(set.intersection(*[set(v) for v in voc.values()]))
    three = set(voc["gf_v1_10m"]) & set(voc["gf_v2_104m"]) & set(voc["esm2_uce"])
    save("vocab", dict(sizes={k: len(v) for k, v in voc.items()}, tierA=len(A), three_submitted=len(three)))
    LAM = lambert(); T = trrust(); TF = set(T.TF); TGO = set(T.target) - TF
    # ---- Tier A node tasks for ESM2-15B and scGPT static
    lam = np.array([g in LAM for g in A])
    lg, li = ctx_obs("lung"); nc = np.load(CACHE / "lung/ncells.npy").astype(float)
    lnob = dict(zip(lg, np.log(np.maximum(nc, 1))))
    # lung detection deciles over lung V (own rank convention), imputed median decile for genes outside lung V
    dl = dict(zip(lg, my_deciles(np.array([lnob[g] for g in lg]))))
    med = int(np.median([dl[g] for g in A if g in dl]))
    detA = np.array([dl.get(g, med) for g in A])
    nl = np.load(CACHE / "lung/node_labels.npz"); dl_c = dict(zip(lg, nl["det_decile"].astype(int)))
    detA_c = np.array([dl_c.get(g, int(np.median([dl_c[x] for x in A if x in dl_c]))) for g in A])
    ttf = np.array([g in TF for g in A]); ttg = np.array([g in TGO for g in A]); lab = ttf | ttg
    res = {}
    for k in ("esm2_uce", "scgpt_static", "gf_v2_104m"):
        X = mats[k][[voc[k][g] for g in A]]
        s = loo_scores(X, lam)
        s_tr = loo_scores(X[lab], ttf[lab])
        res[k] = dict(lambert_auc=roc_auc_score(lam, s), lambert_detstrat_own_deciles=strat_auc(s, lam, detA),
                      lambert_detstrat_cache_deciles=strat_auc(s, lam, detA_c), trrust_auc=roc_auc_score(ttf[lab], s_tr),
                      n_lambert=int(lam.sum()), n_trrust_tf=int(ttf.sum()), n_trrust_tgt=int(ttg.sum()))
    res["n_imputed_detection"] = int(sum(g not in dl for g in A))
    save("tierA_node", res)
    # ---- Tier B lung: contextual L1 and L11, ESM2-15B
    emb = np.load(EMB / "lung.npy", mmap_mode="r")
    Aset = set(A); posl = {g: i for i, g in enumerate(lg)}
    Bg = sorted(set(lg) & Aset); r = np.array([posl[g] for g in Bg])
    lamB = np.array([g in LAM for g in Bg]); dB = np.array([dl_c[g] for g in Bg])
    out = dict(n_tierB_lung=len(Bg))
    for L in (1, 11):
        Xv = np.asarray(emb[L], np.float64)[li]
        assert np.all(np.linalg.norm(Xv, axis=1) > 0)
        s = loo_scores(Xv[r], lamB)
        out[f"lung_L{L}_detstrat"] = strat_auc(s, lamB, dB); out[f"lung_L{L}_unstrat"] = roc_auc_score(lamB, s)
        sV = loo_scores(Xv, np.array([g in LAM for g in lg]))     # V-centred, all of V (CONVENTION3 universe)
        out[f"lung_L{L}_detstrat_on_V_all_genes"] = strat_auc(sV, np.array([g in LAM for g in lg]), nl["det_decile"])
    s = loo_scores(mats["esm2_uce"][[voc["esm2_uce"][g] for g in Bg]], lamB); out["esm2_detstrat"] = strat_auc(s, lamB, dB)
    s = loo_scores(mats["scgpt_static"][[voc["scgpt_static"][g] for g in Bg]], lamB); out["scgpt_static_detstrat"] = strat_auc(s, lamB, dB)
    save("tierB_lung_node", out)
    # ---- legacy: re-derive the 186 genes and estimator (i) / submitted subset-SVD
    e = pd.read_csv(REPO / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main/cycle1_edge_dataset.tsv", sep="\t")
    named = set(e.source.astype(str)) | set(e.target.astype(str))
    kg, ki = ctx_obs("kidney_A"); kemb = np.load(EMB / "kidney_A.npy", mmap_mode="r")
    # non-zero at L0 among the 'named' genes: use the full kidney_A gene axis via obs_idx complement check
    nz0 = np.linalg.norm(np.asarray(kemb[0], np.float64), axis=1) > 0
    assert set(np.flatnonzero(nz0)) == set(ki.tolist())
    cur = sorted(g for g in named if g in set(kg))
    leg = sorted(set(cur) & three)
    ttf_l = np.array([g in TF for g in leg]); assert all((g in TF) or (g in TGO) for g in leg)
    out = dict(n_named=len(named), n_named_nonzero_kidneyA=len(cur), n_legacy=len(leg), n_tf=int(ttf_l.sum()))
    stored = json.loads((E6 / "legacy_genes.json").read_text())["legacy186"]
    out["legacy_set_identical_to_stored"] = sorted(stored) == leg
    kpos = {g: i for i, g in enumerate(kg)}
    for k in ("kidney_A@L11", "gf_v1_10m", "gf_v2_104m", "esm2_uce"):
        if k.startswith("kidney_A"):
            X = np.asarray(kemb[11], np.float64)[ki][[kpos[g] for g in leg]]
        else:
            X = mats[k][[voc[k][g] for g in leg]].astype(np.float64)
        s = loo_scores(X, ttf_l)
        Xc = X - X.mean(0); U, S_, _ = np.linalg.svd(Xc, full_matrices=False); US = (U * S_)[:, 1:7]
        s2 = loo_scores_nocentre(US, ttf_l)
        out[k] = dict(est_i=roc_auc_score(ttf_l, s), subsetSVD_SV2_SV7=roc_auc_score(ttf_l, s2))
    save("legacy", out)
    print("part1", time.time() - t0)


def loo_scores_nocentre(X, is_tf):
    Xn = X / np.linalg.norm(X, axis=1, keepdims=True); c = Xn[is_tf].sum(0); s = Xn @ c; k = is_tf.sum()
    return np.where(is_tf, (s - 1.0) / (k - 1), s / k)


def string700(genes):
    G = set(genes); ens = {}
    with gzip.open(REPO / "subproject_61_two_axes_revision/data/string_protein_info_v12.txt.gz", "rt") as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t"); ens[p[0]] = p[1]
    pairs = set()
    with gzip.open(REPO / "subproject_61_two_axes_revision/data/string_links_v12.txt.gz", "rt") as f:
        next(f)
        for line in f:
            a, b, s = line.split()
            if int(s) < 700: continue
            ga, gb = ens.get(a), ens.get(b)
            if ga and gb and ga != gb and ga in G and gb in G: pairs.add((min(ga, gb), max(ga, gb)))
    return sorted(pairs)


def rewire(pi, pj, n, n_acc, seed):
    """own double-edge swap (degree preserving, simple graph)."""
    rng = np.random.default_rng(seed); a = np.minimum(pi, pj).tolist(); b = np.maximum(pi, pj).tolist(); m = len(a)
    ks = set(zip(a, b)); acc = 0; t = 0; E1 = []
    while acc < n_acc:
        if t >= len(E1):
            E1 = rng.integers(0, m, 200000).tolist(); E2 = rng.integers(0, m, 200000).tolist(); F = (rng.random(200000) < 0.5).tolist(); t = 0
        e1, e2, fl = E1[t], E2[t], F[t]; t += 1
        if e1 == e2: continue
        u, v = a[e1], b[e1]; x, y = a[e2], b[e2]
        if fl: x, y = y, x
        if u == y or x == v: continue
        k1 = (min(u, y), max(u, y)); k2 = (min(x, v), max(x, v))
        if k1 == k2 or k1 in ks or k2 in ks: continue
        ks.discard((a[e1], b[e1])); ks.discard((a[e2], b[e2])); ks.add(k1); ks.add(k2)
        a[e1], b[e1] = k1; a[e2], b[e2] = k2; acc += 1
    return np.array(a), np.array(b)


def part2():
    t0 = time.time()
    voc, mats = {}, {}
    for k in ("esm2_uce", "scgpt_static"):
        voc[k], mats[k] = load_model(k)
    allv = [set(load_model(k)[0]) for k in MODELS]; A = sorted(set.intersection(*allv)); pos = {g: i for i, g in enumerate(A)}
    st = string700(A); pi = np.array([pos[a] for a, b in st]); pj = np.array([pos[b] for a, b in st])
    gsA = np.load(E6 / "work/geneset_tierA.npz")
    same_pairs = bool(np.array_equal(np.sort(pi * len(A) + pj), np.sort(gsA["pi"] * len(A) + gsA["pj"])))
    ra, rb = rewire(pi, pj, len(A), 10 * pi.size, 4242)
    true = set(zip(pi.tolist(), pj.tolist())); keep = np.array([(x, y) not in true for x, y in zip(ra.tolist(), rb.tolist())])
    ni, nj = ra[keep], rb[keep]
    deg0 = np.bincount(np.r_[pi, pj], minlength=len(A)); deg1 = np.bincount(np.r_[ra, rb], minlength=len(A))
    out = dict(n_pairs=int(pi.size), same_pairs_as_E6=same_pairs, degree_preserved=bool(np.array_equal(deg0, deg1)),
               own_rewire_nonpairs=int(ni.size))
    for k in ("esm2_uce", "scgpt_static"):
        X = mats[k][[voc[k][g] for g in A]].astype(np.float64); X -= X.mean(0); Xn = X / np.linalg.norm(X, axis=1, keepdims=True)
        cp = np.einsum("ij,ij->i", Xn[pi], Xn[pj])
        cn_own = np.einsum("ij,ij->i", Xn[ni], Xn[nj]); cn_e6 = np.einsum("ij,ij->i", Xn[gsA["ni"]], Xn[gsA["nj"]])
        y = np.r_[np.ones(cp.size), np.zeros(cn_own.size)]
        out[k] = dict(ppi_auc_own_rewiring=roc_auc_score(y, np.r_[cp, cn_own]),
                      ppi_auc_E6_nonpairs=roc_auc_score(np.r_[np.ones(cp.size), np.zeros(cn_e6.size)], np.r_[cp, cn_e6]))
        # own vertex subsampling (m = n/2, B = 200): SD of subsample AUROC; delete-d SE = SD * sqrt(m/(n-m))
        rng = np.random.default_rng(99); n = len(A); m = n // 2; d = []
        for r in range(200):
            kk = np.zeros(n, bool); kk[rng.choice(n, m, replace=False)] = True
            x = cp[kk[pi] & kk[pj]]; z = cn_e6[kk[gsA["ni"]] & kk[gsA["nj"]]]
            d.append(roc_auc_score(np.r_[np.ones(x.size), np.zeros(z.size)], np.r_[x, z]))
        d = np.array(d); out[k].update(subsample_sd=float(d.std(ddof=1)), deleted_d_halfwidth_approx=float(1.96 * d.std(ddof=1) * np.sqrt(m / (n - m))),
                                       prereg_halfwidth_approx=float(1.96 * d.std(ddof=1) * np.sqrt(m / n)))
        Xn32 = Xn.astype(np.float32); out[k]["_"] = None
        # ---- edges: stored E6 pairs, own within-TF AUROC + own construction checks
        tf, tg, yy = gsA["e_tf"], gsA["e_tg"], gsA["e_y"].astype(bool)
        s = np.einsum("ij,ij->i", Xn[tf], Xn[tg]); num = den = 0.0
        for t in np.unique(tf):
            mm = tf == t; p = yy[mm].sum(); q = mm.sum() - p
            if p and q: num += roc_auc_score(yy[mm], s[mm]) * p * q; den += p * q
        out[k]["edge_within_tf_auc"] = num / den
    # pair construction rules
    T = trrust(); T = T[T.TF != T.target]; trr = set(zip(T.TF, T.target))
    D1 = pd.read_csv(NET / "dorothea_human.tsv", sep="\t"); D2 = pd.read_csv(NET / "dorothea_chipseq_human.tsv", sep="\t")
    anyreg = trr | set(zip(D1.source.astype(str), D1.target.astype(str))) | set(zip(D2.source.astype(str), D2.target.astype(str)))
    G = np.array(A); P = list(zip(G[tf[yy]], G[tg[yy]])); N = list(zip(G[tf[~yy]], G[tg[~yy]]))
    trA = {p for p in trr if p[0] in pos and p[1] in pos}
    tgt_of = {}
    for a, b in trA: tgt_of.setdefault(b, set()).add(a)
    out["edge_checks"] = dict(
        n_pos=len(P), n_neg=len(N), n_tf=int(np.unique(tf).size),
        all_pos_in_trrust=all(p in trr for p in P), any_neg_in_anyreg=any(p in anyreg for p in N),
        any_neg_reverse_in_anyreg=any((b, a) in anyreg for a, b in N),
        all_neg_target_of_other_tf=all(any(x != a for x in tgt_of.get(b, ())) for a, b in N),
        min_pos_per_tf=int(min(np.bincount(np.unique(tf[yy], return_inverse=True)[1]))),
        n_trrust_within_tierA=len(trA), dup_pairs=len(set(zip(tf.tolist(), tg.tolist()))) != tf.size)
    save("tierA_ppi_edge", out)
    print("part2", time.time() - t0)


def part3():
    """selection-bias null, own implementation (2,000 matched subsets) for scGPT static and ESM2-15B."""
    t0 = time.time()
    allv = [set(load_model(k)[0]) for k in MODELS]; A = sorted(set.intersection(*allv))
    T = trrust(); TF = set(T.TF); TGO = set(T.target) - TF
    leg = json.loads((E6 / "legacy_genes.json").read_text())["legacy186"]
    pool = [g for g in A if g in TF or g in TGO]
    lg, _ = ctx_obs("lung"); ex = np.load(CACHE / "lung/expr.npy", mmap_mode="r"); me = dict(zip(lg, np.asarray(ex, np.float64).mean(0)))
    v = np.array([me[g] for g in pool if g in me]); edges = np.quantile(v, np.linspace(0, 1, 11))
    def dec(g):
        return "NA" if g not in me else int(min(max(np.searchsorted(edges, me[g], side="right") - 1, 0), 9))
    cell = {g: ("TF" if g in TF else "TG", dec(g)) for g in set(pool) | set(leg)}
    byc = {}
    for g in pool: byc.setdefault(cell[g], []).append(g)
    need = pd.Series([cell[g] for g in leg]).value_counts().to_dict()
    rng = np.random.default_rng(31337); out = dict(pool=len(pool), n_tf_pool=int(sum(g in TF for g in pool)))
    draws = []
    for r in range(2000):
        sel = []
        for c, n in need.items(): sel += list(rng.choice(byc[c], n, replace=False))
        draws.append(sel)
    for k in ("scgpt_static", "esm2_uce"):
        pos, E = load_model(k)
        def stat(gl):
            X = E[[pos[g] for g in gl]].astype(np.float64); y = np.array([g in TF for g in gl])
            return roc_auc_score(y, loo_scores(X, y))
        obs = stat(leg); nul = np.array([stat(d) for d in draws])
        out[k] = dict(legacy=obs, null_mean=float(nul.mean()), percentile=float(np.mean(nul < obs) + 0.5 * np.mean(nul == obs)))
    save("selection_bias", out)
    print("part3", time.time() - t0)


if __name__ == "__main__":
    dict(part1=part1, part2=part2, part3=part3)[sys.argv[1]]()
