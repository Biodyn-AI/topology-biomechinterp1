"""E4c verifier part 5 (this verifier): independent recomputation of headline numbers from raw inputs.

Independent of e4c_lib's features, standardisation, solver, inner split, alpha selection and metric. Only the OUTER fold
assignment is taken from e4c_lib.make_folds (it must match so that numbers are comparable; its disjointness was checked
independently in v1_folds.csv). Everything else is re-implemented here:
  - V-centred embedding from the raw local .npy and cache obs_idx; static embedding re-read from the subproject_53 npz.
  - own PCA-16 on training-fold genes; own bilinear features; own standardisation on training pairs.
  - own inner split: training genes dealt at random into 4 groups (verifier seed); inner-train = pairs with both genes in
    groups 1-3, inner-validation = pairs with both genes in group 0 (pair-disjoint, like C3).
  - alpha grid {1e-4,1e-3,1e-2,1e-1} on mean-log-loss + alpha/2||w||^2 (sklearn C = 1/(n alpha)); chosen by own
    within-TF AUROC on the inner validation pairs; refit on the full outer training fold.
  - own within-TF Mann-Whitney AUROC (scipy rankdata), TFs with >= 3 positives and >= 3 non-edges.
Units:
  A  lung scGPT bilinear L2, C3, 5 repeats x 5 folds   (implementer: 0.5968)
  B  lung scGPT target-only L2, C3, 25 folds            (implementer: 0.5459; A-B = +0.0509 [0.0259, 0.0759])
  C  lung static bilinear, C3, 25 folds                 (implementer: 0.5813)
  Cn lung static target-only, C3, 25 folds              (implementer: 0.5550; C-Cn = +0.0263 [0.0028, 0.0497])
  D  lung_external scGPT concat L11, C3, 25 folds + ordered direction test (implementer: 0.594 forward, 453 pairs)
  N  lung scGPT bilinear L2, C3 repeat 0, rewired draws 1-10 (implementer's alpha per column; implementer per-draw means)
usage: python v5_independent.py A|B|C|Cn|D|N   -> verify/v5_<unit>.csv
"""
import json, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import rankdata
from sklearn.linear_model import LogisticRegression

HERE = Path(__file__).resolve().parent; E4C = HERE.parent; R3 = E4C.parent
sys.path.insert(0, str(E4C))
import e4c_lib as L   # make_folds ONLY
ALPHAS = [1e-4, 1e-3, 1e-2, 1e-1]
unit = sys.argv[1]
CFG = {"A": ("lung", "scgpt", "bilinear", 2), "B": ("lung", "scgpt", "tg_only", 2), "C": ("lung", "static", "bilinear", -1),
       "Cn": ("lung", "static", "tg_only", -1), "D": ("lung_external", "scgpt", "concat", 11), "N": ("lung", "scgpt", "bilinear", 2)}
ctx, src, feat, layer = CFG[unit]
d = R3 / "cache" / ctx
genes = json.loads((d / "obs_genes.json").read_text()); oi = np.load(d / "obs_idx.npy")
z = np.load(d / "h3_pairs.npz"); tf, tg, y = z["tf"].astype(int), z["tg"].astype(int), z["y"].astype(bool)
G = len(genes)
if src == "scgpt":
    emb = np.load(Path.home() / f"biomechinterp-local/subproject41_r3/emb/{ctx}.npy", mmap_mode="r")
    E = np.asarray(emb[layer], np.float64)[oi]; assert (np.linalg.norm(E, axis=1) > 0).all(); E -= E.mean(0); extra = None
else:
    zz = np.load(R3.parents[1] / "subproject_53_scgpt_gpl_replication/embeddings/scgpt_gene_embeddings.npz")
    pos = {str(g): i for i, g in enumerate(zz["gene_names"])}; M = zz["embeddings"]
    have = np.array([g in pos for g in genes]); E = np.zeros((G, M.shape[1]))
    E[have] = M[[pos[g] for g in np.array(genes)[have]]]; E[have] -= E[have].mean(0)
    h3g = np.unique(np.r_[tf, tg]); extra = None if have[h3g].all() else (~have).astype(float)[:, None]
    print("static coverage of h3 genes", have[h3g].mean(), flush=True)


def wtf_auc(t, yy, s):
    num = den = 0.0
    for u in np.unique(t):
        m = t == u; a = yy[m]; n1 = a.sum(); n0 = m.sum() - n1
        if n1 < 3 or n0 < 3: continue
        r = rankdata(s[m]); num += r[a].sum() - n1 * (n1 + 1) / 2; den += n1 * n0
    return num / den if den else np.nan


def features(trg):
    if feat == "bilinear":
        mu = E[trg].mean(0); _, _, Vt = np.linalg.svd(E[trg] - mu, full_matrices=False); P = (E - mu) @ Vt[:16].T
        F = (P[tf][:, :, None] * P[tg][:, None, :]).reshape(tf.size, 256)
    elif feat == "tg_only": F = E[tg]
    elif feat == "concat": F = np.hstack([E[tf], E[tg]])
    if extra is not None:
        F = np.hstack([F, extra[tg]]) if feat == "tg_only" else np.hstack([F, extra[tf], extra[tg]])
    return F


def fit(X, yy, a):
    return LogisticRegression(C=1.0 / (X.shape[0] * a), max_iter=3000, tol=1e-7).fit(X, yy)


rows, drows = [], []
if unit == "D":
    te_ = np.load(d / "edges_trrust.npy"); tset = {(int(a), int(b)) for a, b in te_}
    lam = np.load(d / "node_labels.npz")["is_tf_lambert"]
    DIR = np.array([(a, b) for a, b in sorted(tset) if lam[a] and lam[b] and (b, a) not in tset])
    print("direction pairs", len(DIR), flush=True)
if unit == "N":
    zr = np.load(E4C / "rewired" / f"{ctx}.npz"); Yr = np.unpackbits(zr["Y"], axis=1, count=int(zr["n"])).astype(bool)
    ck = pd.read_csv(E4C / "partial" / "m1" / ctx / f"scgpt_bilinear_L{layer}_C3_r0_null.csv")
reps = [0] if unit == "N" else range(5)
CTXOBJ = L.Ctx(ctx)
assert np.array_equal(CTXOBJ.tf, tf) and np.array_equal(CTXOBJ.tg, tg) and np.array_equal(CTXOBJ.y, y)
out_fn = HERE / f"v5_{unit}.csv"
done = pd.read_csv(out_fn) if out_fn.exists() else pd.DataFrame()
rows = done.to_dict("records")
for rep in reps:
    for k, (tr, te) in enumerate(L.make_folds(CTXOBJ, "C3", rep)):
        if len(done) and ((done["repeat"] == rep) & (done.fold == k)).any(): continue
        t0 = time.time()
        trg = np.unique(np.r_[tf[tr], tg[tr]])
        F = features(trg); m_ = F[tr].mean(0); s_ = F[tr].std(0); s_[s_ == 0] = 1
        Xtr = (F[tr] - m_) / s_; Xte = (F[te] - m_) / s_
        cols = [(0, y)] if unit != "N" else [(j, (y if j == 0 else Yr[j - 1])) for j in range(0, 11)]
        for j, lab in cols:
            if unit == "N":
                alpha = float(ck[(ck.fold == k) & (ck.label_col == j)].alpha.values[0]); inner = []
            else:
                # own inner split: training genes dealt randomly into 4 groups (verifier seed)
                rng = np.random.default_rng([777, rep, k]); grp = dict(zip(trg.tolist(), rng.permutation(trg.size) % 4))
                a_ = np.array([grp[int(g)] for g in tf[tr]]); b_ = np.array([grp[int(g)] for g in tg[tr]])
                itr = np.flatnonzero((a_ != 0) & (b_ != 0)); iva = np.flatnonzero((a_ == 0) & (b_ == 0))
                F2 = F[tr]; m2 = F2[itr].mean(0); s2 = F2[itr].std(0); s2[s2 == 0] = 1
                Xi = (F2[itr] - m2) / s2; Xv = (F2[iva] - m2) / s2
                inner = [wtf_auc(tf[tr][iva], lab[tr][iva], fit(Xi, lab[tr][itr], a).decision_function(Xv)) for a in ALPHAS]
                alpha = ALPHAS[int(np.nanargmax(inner))]
            clf = fit(Xtr, lab[tr], alpha); s = clf.decision_function(Xte)
            rows.append(dict(context=ctx, source=src, feature=feat, layer=layer, repeat=rep, fold=k, label_col=j, alpha=alpha,
                             inner=json.dumps([round(float(x), 4) for x in inner]), auc=wtf_auc(tf[te], lab[te], s),
                             n_train=tr.size, n_test=te.size))
            if unit == "D":
                ok = ~np.isin(DIR[:, 0], trg) & ~np.isin(DIR[:, 1], trg); a2, b2 = DIR[ok, 0], DIR[ok, 1]
                Ft = np.hstack([E[a2], E[b2]]); Fr = np.hstack([E[b2], E[a2]])
                dd = clf.decision_function((Ft - m_) / s_) - clf.decision_function((Fr - m_) / s_)
                for x, w_, v in zip(a2, b2, dd): drows.append(dict(repeat=rep, fold=k, i=int(x), j=int(w_), diff=float(v)))
        pd.DataFrame(rows).to_csv(out_fn, index=False)
        if unit == "D":
            dfn = HERE / "v5_D_direction.csv"
            prev = pd.read_csv(dfn) if dfn.exists() else pd.DataFrame()
            pd.concat([prev, pd.DataFrame(drows)], ignore_index=True).to_csv(dfn, index=False); drows = []
        print(f"{unit} rep{rep} fold{k} auc0={rows[-len(cols)]['auc']:.4f} alpha={rows[-len(cols)]['alpha']} {time.time()-t0:.0f}s", flush=True)
R = pd.DataFrame(rows)
print(unit, "mean AUROC (label 0):", R[R.label_col == 0].auc.mean())
