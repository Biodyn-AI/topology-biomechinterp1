"""Shared independent helpers for the E3 verifier (no stats3 / e3_lib imports).

Does NOT import stats3 or e3_lib. Re-implements: observed set, V-centring, LOO mean cosine, detection deciles,
stratified Mann-Whitney AUROC (scipy rankdata per stratum), SVD, Lambert labels from humantfs_v1.01.txt,
static embedding symbol matching, logistic LR test. Writes verify/verify_e3.json and verify/verify_e3.csv.
"""
import json, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import rankdata, chi2

HERE = Path(__file__).resolve().parent; E3 = HERE.parent; R3 = E3.parent; REPO = R3.parents[1]
EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"
LAMBERT = REPO / "subproject_61_two_axes_revision/data/humantfs_v1.01.txt"
ESM = REPO / "subproject_60_uce_esm2_control_gpl_replication/embeddings/uce_esm2_human_gene_embeddings.npz"
SCG = REPO / "subproject_53_scgpt_gpl_replication/embeddings/scgpt_gene_embeddings.npz"
t0 = time.time()
def log(*a): print(f"[{time.time()-t0:6.1f}s]", *a, flush=True)

lam = pd.read_csv(LAMBERT, sep="\t")
LAM = set(lam.loc[lam["Is TF?"].astype(str).str.strip() == "Yes", "HGNC symbol"].astype(str))

def mw(score, y):
    r = rankdata(score); n1 = y.sum(); n0 = y.size - n1
    return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

def strat(score, y, st):
    num = den = 0.0
    for s in np.unique(st):
        m = st == s; yy = y[m]
        if yy.any() and (~yy).any():
            w = yy.sum() * (~yy).sum(); num += mw(score[m], yy) * w; den += w
    return num / den

def deciles(x):
    q = np.unique(np.quantile(x, np.linspace(0, 1, 11)))
    return np.clip(np.searchsorted(q, x, side="right") - 1, 0, len(q) - 2)

def loo(E, y):
    n = np.linalg.norm(E, axis=1); En = E / n[:, None]
    v = En[y].sum(0); s = En @ v; k = y.sum()
    return np.where(y, (s - 1.0) / (k - 1), s / k)

def ctx_data(ctx):
    emb = np.load(EMB / f"{ctx}.npy", mmap_mode="r")
    import h5py
    sys.path.insert(0, str(R3)); import contexts as X        # paths only
    with h5py.File(X.CONTEXTS[ctx]["h5ad"], "r") as f:
        genes_all = [g.decode() if isinstance(g, bytes) else str(g) for g in f["var"]["_index"][:]]
    nrm0 = np.linalg.norm(np.asarray(emb[0], float), axis=1)
    obs = np.flatnonzero(nrm0 > 0)
    genes = [genes_all[i] for i in obs]
    cnt = np.load(Path(X.CONTEXTS[ctx]["emb"]).parent / "layer_gene_embedding_count.npy", mmap_mode="r")[0][obs].astype(float)
    y = np.array([g in LAM for g in genes])
    det = deciles(np.log(np.maximum(cnt, 1)))
    # cross-check against the shared cache
    c = R3 / "cache" / ctx
    assert json.loads((c / "obs_genes.json").read_text()) == genes
    assert np.array_equal(np.load(c / "obs_idx.npy"), obs)
    nl = np.load(c / "node_labels.npz")
    assert np.array_equal(nl["is_tf_lambert"].astype(bool), y) and np.array_equal(nl["det_decile"], det)
    return emb, obs, genes, cnt, y, det

