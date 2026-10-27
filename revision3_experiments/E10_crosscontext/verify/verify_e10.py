"""Independent verifier for E10 (adversarial audit). Recomputes headline numbers from RAW inputs with code that
does not import e10_common or stats3 estimators (only file paths from contexts.py):
  embeddings (local copies), h5ad var index (h5py), layer_gene_embedding_count.npy (n_obs), humantfs_v1.01.txt,
  trrust_human.tsv, STRING v12 physical links + protein info. The H3 benchmark (cache h3_pairs.npz) is the shared
  pre-registered pair set; only its AUROC is recomputed here (sklearn).
Checks (V1..V9) write verify/verify_results.json; checkpoint per check in verify/ck_<name>.json.
"""
import gzip, json, sys, time
from pathlib import Path
import numpy as np, pandas as pd, h5py
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold

HERE = Path(__file__).resolve().parent
E10 = HERE.parent; R3 = E10.parent
sys.path.insert(0, str(R3))
import contexts as X  # paths only

EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"
LAM = X.REPO / "subproject_61_two_axes_revision/data/humantfs_v1.01.txt"


def log(m):
    print(time.strftime("%H:%M:%S"), m, flush=True)


def raw_ctx(ctx):
    emb = np.load(EMB / f"{ctx}.npy", mmap_mode="r")
    with h5py.File(X.CONTEXTS[ctx]["h5ad"], "r") as f:
        genes = [g.decode() if isinstance(g, bytes) else str(g) for g in f["var"]["_index"][:]]
    nz = np.linalg.norm(np.asarray(emb[0], float), axis=1) > 0
    idx = np.flatnonzero(nz)
    cnt = np.load(Path(X.CONTEXTS[ctx]["emb"]).parent / "layer_gene_embedding_count.npy")[0][idx].astype(float)
    return emb, [genes[i] for i in idx], idx, cnt


def lambert():
    d = pd.read_csv(LAM, sep="\t")
    return set(d.loc[d["Is TF?"].astype(str).str.strip() == "Yes", "HGNC symbol"].astype(str))


def trrust_pairs():
    t = pd.read_csv(X.NET / "trrust_human.tsv", sep="\t", header=None, names=["TF", "tg", "m", "p"])
    return {(a, b) for a, b in zip(t.TF.astype(str), t.tg.astype(str)) if a != b}


def mat(emb, L, idx, centre=True):
    M = np.asarray(emb[L], float)[idx]
    if centre: M = M - M.mean(0)
    return M


def nrm(M):
    return M / np.linalg.norm(M, axis=1, keepdims=True)


def loo(En, pos_idx):
    s = En @ En[pos_idx].sum(0); k = len(pos_idx); out = s / k
    out[pos_idx] = (s[pos_idx] - 1.0) / (k - 1)
    return out


def deciles(x):
    e = np.unique(np.quantile(x, np.linspace(0, 1, 11)))
    return np.clip(np.searchsorted(e, x, side="right") - 1, 0, len(e) - 2)


def strat_auc(s, y, st):
    num = den = 0.0
    for v in np.unique(st):
        m = st == v
        if y[m].any() and (~y[m]).any():
            w = y[m].sum() * (~y[m]).sum(); num += roc_auc_score(y[m], s[m]) * w; den += w
    return num / den


def within_tf_auc(tf, y, s):
    num = den = 0.0
    for t in np.unique(tf):
        m = tf == t
        if y[m].any() and (~y[m]).any():
            w = y[m].sum() * (~y[m]).sum(); num += roc_auc_score(y[m], s[m]) * w; den += w
    return num / den


def ck(name, fn):
    f = HERE / f"ck_{name}.json"
    if f.exists(): return json.loads(f.read_text())
    t0 = time.time(); r = fn(); r["_seconds"] = round(time.time() - t0, 1)
    f.write_text(json.dumps(r, indent=1, default=float)); log(f"{name}: {json.dumps(r, default=float)[:600]}")
    return r


def e10_point(ctx, stat):
    z = np.load(E10 / "partial" / f"stats_{ctx}.npz"); st = [str(s) for s in z["stats"]]
    return z["point"][st.index(stat)]


# ---------------------------------------------------------------- V1 S2 (H2 primary statistic)
def v1():
    lam = lambert(); out = {}
    for ctx in ["immune", "lung", "lung_external", "kidney_B", "kidney_A"]:
        emb, g, idx, cnt = raw_ctx(ctx)
        y = np.array([x in lam for x in g]); st = deciles(np.log(np.maximum(cnt, 1))); P = np.flatnonzero(y)
        vals = [strat_auc(loo(nrm(mat(emb, L, idx)), P), y, st) for L in range(12)]
        e = e10_point(ctx, "S2")
        out[ctx] = dict(mine=vals, e10=e.tolist(), max_abs_diff=float(np.max(np.abs(np.array(vals) - e))),
                        mine_layer_mean=float(np.mean(vals)), e10_layer_mean=float(e.mean()), n_tf=int(y.sum()), n_V=len(g))
    return out


# ---------------------------------------------------------------- V2 S3 (H3 primary statistic)
def v2():
    out = {}
    for ctx in ["immune", "lung", "kidney_B"]:
        emb, g, idx, _ = raw_ctx(ctx)
        z = np.load(R3 / "cache" / ctx / "h3_pairs.npz"); tf, tg, y = z["tf"], z["tg"], z["y"].astype(bool)
        assert json.loads((R3 / "cache" / ctx / "obs_genes.json").read_text()) == g      # index alignment
        vals = []
        for L in range(12):
            En = nrm(mat(emb, L, idx)); vals.append(within_tf_auc(tf, y, (En[tf] * En[tg]).sum(1)))
        e = e10_point(ctx, "S3")
        out[ctx] = dict(mine=vals, e10=e.tolist(), max_abs_diff=float(np.max(np.abs(np.array(vals) - e))),
                        mine_layer_mean=float(np.mean(vals)), e10_layer_mean=float(e.mean()))
    return out


# ---------------------------------------------------------------- V3 S1 with own STRING mapping + own rewiring
def string_physical700():
    ens = {}
    with gzip.open(X.STRING_INFO, "rt") as f:
        next(f)
        for line in f:
            p = line.split("\t"); ens[p[0]] = p[1]
    E = set()
    with gzip.open(X.STRING_LINKS, "rt") as f:
        next(f)
        for line in f:
            a, b, s = line.split()
            if int(s) < 700: continue
            x, y = ens.get(a), ens.get(b)
            if x and y and x != y: E.add((min(x, y), max(x, y)))
    return E


def rewire_nx(edges, n, seed):
    import networkx as nx
    G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from(map(tuple, edges))
    m = G.number_of_edges()
    nx.double_edge_swap(G, nswap=11 * m, max_tries=200 * m, seed=seed)
    R = np.array([(min(a, b), max(a, b)) for a, b in G.edges()])
    key = set(map(tuple, edges))
    return np.array([e for e in map(tuple, R) if e not in key])


def v3():
    E = string_physical700(); out = {"n_string700_genomewide": len(E)}
    for ctx, ndraw in (("immune", 20), ("kidney_ckpt_kidney", 6)):
        emb, g, idx, _ = raw_ctx(ctx); pos = {x: i for i, x in enumerate(g)}
        ed = np.array(sorted((pos[a], pos[b]) for a, b in E if a in pos and b in pos))
        ed = np.c_[np.minimum(ed[:, 0], ed[:, 1]), np.maximum(ed[:, 0], ed[:, 1])]
        cache_e = np.load(R3 / "cache" / ctx / "edges_string.npy")
        cache_set = {(min(a, b), max(a, b)) for a, b in cache_e}
        Ens = [nrm(mat(emb, L, idx)) for L in range(12)]
        draws = []
        for d in range(ndraw):
            neg = rewire_nx(ed, len(g), seed=1000 + d)
            pa = np.r_[ed[:, 0], neg[:, 0]]; pb = np.r_[ed[:, 1], neg[:, 1]]
            yy = np.r_[np.ones(len(ed)), np.zeros(len(neg))]
            draws.append(float(np.mean([roc_auc_score(yy, (En[pa] * En[pb]).sum(1)) for En in Ens])))
        out[ctx] = dict(n_edges_mine=len(ed), n_edges_cache=len(cache_set), edges_identical=set(map(tuple, ed)) == cache_set,
                        S1_layer_mean_draws=draws, draw_mean=float(np.mean(draws)), draw_sd=float(np.std(draws, ddof=1)),
                        e10_single_draw_layer_mean=float(e10_point(ctx, "S1").mean()))
    return out


# ---------------------------------------------------------------- V4 legacy bridge (submitted zero-row SVD)
def v4():
    tsv = X.REPO / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main/cycle1_edge_dataset.tsv"
    e = pd.read_csv(tsv, sep="\t"); named = set(e["source"].astype(str)) | set(e["target"].astype(str))
    t = pd.read_csv(X.NET / "trrust_human.tsv", sep="\t", header=None, names=["TF", "tg", "m", "p"])
    tfs = set(t.TF.astype(str)); tgo = set(t.tg.astype(str)) - tfs
    ea = np.load(EMB / "kidney_A.npy", mmap_mode="r"); eb = np.load(EMB / "kidney_B.npy", mmap_mode="r")
    def genes(ctx):
        with h5py.File(X.CONTEXTS[ctx]["h5ad"], "r") as f:
            return [x.decode() for x in f["var"]["_index"][:]]
    gA, gB = genes("kidney_A"), genes("kidney_B")
    nzA = np.linalg.norm(np.asarray(ea[0], float), axis=1) > 0; nzB = np.linalg.norm(np.asarray(eb[0], float), axis=1) > 0
    namesA = [x for i, x in enumerate(gA) if x in named and nzA[i]]
    subB = [i for i, x in enumerate(gB) if x in set(namesA) and nzB[i]]; namesB = [gB[i] for i in subB]
    tf = np.array([k for k, x in enumerate(namesB) if x in tfs]); tg = np.array([k for k, x in enumerate(namesB) if x in tgo])
    vals = []
    for L in range(12):
        M = np.asarray(eb[L], float); M = M - M.mean(0)
        U, S, _ = np.linalg.svd(M, full_matrices=False)
        En = nrm((U * S)[subB][:, 1:7]); s = loo(En, tf)
        vals.append(roc_auc_score(np.r_[np.ones(tf.size), np.zeros(tg.size)], np.r_[s[tf], s[tg]]))
    k = pd.read_csv(E10 / "kidneyA_vs_B.csv")
    ref = k[k.statistic == "legacy_bridge_B_curated167_submitted_zero_row_SVD_band"].assign(L=lambda d: d.layer.astype(int)).sort_values("L").A.values
    return dict(n=len(namesB), n_tf=int(tf.size), n_tg=int(tg.size), mine=vals, e10=ref.tolist(),
                max_abs_diff=float(np.max(np.abs(np.array(vals) - ref))), manuscript_L0_L3=[0.493, 0.500, 0.529, 0.551])


# ---------------------------------------------------------------- V5 kidney matched S2b_full A-B
def v5():
    trr = trrust_pairs()
    ea, gA, ia_, _ = raw_ctx("kidney_A"); eb, gB, ib_, _ = raw_ctx("kidney_B")
    M = sorted(set(gA) & set(gB)); pa = {x: i for i, x in enumerate(gA)}; pb = {x: i for i, x in enumerate(gB)}
    ia = np.array([pa[x] for x in M]); ib = np.array([pb[x] for x in M]); V = set(M)
    pv = [(a, b) for a, b in trr if a in V and b in V]; T = {a for a, _ in pv}; G = {b for _, b in pv} - T
    tf = np.array([k for k, x in enumerate(M) if x in T]); tg = np.array([k for k, x in enumerate(M) if x in G])
    yy = np.r_[np.ones(tf.size), np.zeros(tg.size)]
    A, B = [], []
    for L in range(12):
        for e_, i_, ii, acc in ((ea, ia_, ia, A), (eb, ib_, ib, B)):
            En = nrm(mat(e_, L, i_)[ii]); s = loo(En, tf); acc.append(roc_auc_score(yy, np.r_[s[tf], s[tg]]))
    k = pd.read_csv(E10 / "kidneyA_vs_B.csv"); r = k[(k.statistic == "S2b_full") & (k.layer.astype(str) == "mean")].iloc[0]
    return dict(M=len(M), n_tf=int(tf.size), n_tg=int(tg.size), A_mean=float(np.mean(A)), B_mean=float(np.mean(B)),
                delta=float(np.mean(A) - np.mean(B)), e10=dict(A=r.A, B=r.B, delta=r.delta))


# ---------------------------------------------------------------- V6 detection-only AUROC
def v6():
    lam = lambert(); trr = trrust_pairs(); out = {}
    d = pd.read_csv(E10 / "detection.csv", low_memory=False); d = d[d.part == "a"]
    for ctx in ["immune", "lung", "kidney_B", "kidney_A"]:
        _, g, _, cnt = raw_ctx(ctx); ln = np.log(np.maximum(cnt, 1)); y = np.array([x in lam for x in g])
        V = set(g); pv = [(a, b) for a, b in trr if a in V and b in V]; T = {a for a, _ in pv}; G = {b for _, b in pv} - T
        t = np.array([x in T for x in g]); o = np.array([x in G for x in g]); m = t | o
        out[ctx] = dict(lambert=roc_auc_score(y, ln), trrust=roc_auc_score(t[m], ln[m]),
                        e10_lambert=float(d[(d.context == ctx) & (d.test == "lognobs_auroc_lambertTF_vs_nonTF")].value.iloc[0]),
                        e10_trrust=float(d[(d.context == ctx) & (d.test == "lognobs_auroc_trrustTF_vs_targetonly")].value.iloc[0]))
    return out


# ---------------------------------------------------------------- V7 transfer (immune -> novel genes)
def v7():
    lam = lambert(); ei, gi, ii, _ = raw_ctx("immune"); yi = np.array([x in lam for x in gi]); imm = set(gi)
    out = {}
    tgt = {t: raw_ctx(t) for t in ["kidney_B", "lung"]}
    per = {t: [] for t in tgt}; per_in = {t: [] for t in tgt}
    for L in range(12):
        Xi = nrm(mat(ei, L, ii)); best, bC = -1, None
        for C in (0.01, 0.1, 1.0, 10.0):          # own CV for C (different fold seed from E10)
            s = np.zeros(yi.size)
            for tr, te in StratifiedKFold(5, shuffle=True, random_state=7 + L).split(Xi, yi):
                s[te] = LogisticRegression(C=C, class_weight="balanced", max_iter=3000).fit(Xi[tr], yi[tr]).decision_function(Xi[te])
            a = roc_auc_score(yi, s)
            if a > best: best, bC = a, C
        mdl = LogisticRegression(C=bC, class_weight="balanced", max_iter=3000).fit(Xi, yi)
        for t, (e, g, idx, _) in tgt.items():
            Xt = nrm(mat(e, L, idx)); y = np.array([x in lam for x in g]); nov = np.array([x not in imm for x in g])
            per[t].append(roc_auc_score(y[nov], mdl.decision_function(Xt[nov])))
            s = np.zeros(y.size)
            for tr, te in StratifiedKFold(5, shuffle=True, random_state=11 + L).split(Xt, y):
                s[te] = LogisticRegression(C=bC, class_weight="balanced", max_iter=3000).fit(Xt[tr], y[tr]).decision_function(Xt[te])
            per_in[t].append(roc_auc_score(y[nov], s[nov]))
        log(f"V7 L{L} C={bC} " + " ".join(f"{t}={per[t][-1]:.3f}/{per_in[t][-1]:.3f}" for t in tgt))
    tr_ = pd.read_csv(E10 / "transfer.csv")
    for t in tgt:
        r = tr_[(tr_.layer.astype(str) == "mean") & (tr_.test == t) & (tr_.subset == "novel") & (tr_.centred == True)].iloc[0]
        out[t] = dict(transfer_mean=float(np.mean(per[t])), incontext_mean=float(np.mean(per_in[t])),
                      n_layers_transfer_below=int(np.sum(np.array(per[t]) < np.array(per_in[t]))),
                      e10_transfer=float(r.auroc), e10_incontext=float(r.incontext_auroc))
    return out


# ---------------------------------------------------------------- V8 varcomp re-derivation from the stored bootstrap
def v8():
    out = {}
    groups = {"resample_immune": ["immune", "immune_seed43", "immune_seed44"], "hvg_axis": ["kidney_A", "kidney_B"],
              "tissue": ["immune", "kidney_B", "lung", "lung_external"], "checkpoint": ["kidney_ckpt_whole_human", "kidney_ckpt_kidney", "kidney_ckpt_brain"]}
    v = pd.read_csv(E10 / "varcomp.csv")
    for stat in ["S1", "S2", "S3"]:
        for g, cs in groups.items():
            th, bt = [], []
            for c in cs:
                z = np.load(E10 / "partial" / f"stats_{c}.npz"); st = [str(s) for s in z["stats"]]; i = st.index(stat)
                th.append(z["point"][i].mean()); bt.append(z["boot"][i].mean(0))
            th = np.array(th); bt = np.array(bt); bt = bt[:, np.all(np.isfinite(bt), 0)]
            # (a) DerSimonian-Laird with diagonal (independent-error) variances
            vi = bt.var(1, ddof=1); w = 1 / vi; mu = (w * th).sum() / w.sum(); Q = (w * (th - mu) ** 2).sum(); k = th.size
            tau2_dl = max(0.0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum()))
            # (b) the full-covariance moment estimator, re-derived by Monte Carlo identity: E[var(e)] from bootstrap columns
            dev = bt - bt.mean(1, keepdims=True)
            samp = float(np.mean(np.var(dev, axis=0, ddof=1)))
            tau2_mom = float(np.var(th, ddof=1) - samp)
            r = v[(v.stat == stat) & (v.group == g) & (v.layer.astype(str) == "mean")].iloc[0]
            out[f"{stat}|{g}"] = dict(tau2_DL_diag=tau2_dl, tau2_mom_mine=tau2_mom, tau2_e10=float(r.tau2), Q_diag=float(Q))
    return out


# ---------------------------------------------------------------- V9 shared drawn cells kidney_A / kidney_B
def v9():
    a = set(np.load(R3 / "cache/kidney_A/extraction_cells.npy").tolist()); b = set(np.load(R3 / "cache/kidney_B/extraction_cells.npy").tolist())
    return dict(shared=len(a & b), e10_hardcoded=json.loads((E10 / "partial/kidney_mech.json").read_text())["drawn_cells_shared"])


if __name__ == "__main__":
    which = sys.argv[1:] or ["v9", "v6", "v1", "v2", "v4", "v5", "v8", "v3", "v7"]
    res = {}
    for w in which:
        res[w] = ck(w, globals()[w])
    allf = HERE / "verify_results.json"
    prev = json.loads(allf.read_text()) if allf.exists() else {}
    prev.update(res); allf.write_text(json.dumps(prev, indent=1, default=float))
    log("done " + " ".join(which))
