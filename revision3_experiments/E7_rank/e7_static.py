"""E7 descriptive: family-T (top-k vs full) curves for the static gene-embedding models (E6's model list) on the
lung gene set, PPI / NODE / EDGE metrics as in the primary analysis. No selection, no Delta CI (descriptive).

Gene set: E6's Tier B lung set (E6_crossmodel/work/geneset_tierB_lung.npz, read-only) = lung V intersected with
all 11 static vocabularies (5,226 of lung's 5,648 observed genes), so every model is scored on identical genes and
pairs. Its benchmark arrays were built by E6 with the prep_pairs rules restricted to the set: STRING-physical >= 700
pairs, ONE Maslov-Sneppen rewiring (burn-in 10m) with true pairs removed, Lambert TF labels, lung detection deciles,
and H3 target-matched TRRUST-negative candidate pairs.
Representation: the model's gene vectors on the set, mean-centred over the set, thin SVD; U_k S_k; 'full' = all
columns (cosine-equivalent to the centred embedding). Anchor row: scGPT lung contextual embedding (layer mean of 12
layers) on the same genes, SVD fitted on the same set.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

import e7lib as L

S = L.S
R3 = L.R3
GS = R3 / "E6_crossmodel" / "work" / "geneset_tierB_lung.npz"
STATIC = {   # key -> npz path relative to REPO (E6_crossmodel/e6lib.py STATIC, same files)
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


def bench():
    z = np.load(GS)
    genes = [str(g) for g in z["genes"]]; n = len(genes)
    pa = np.r_[z["pi"], z["ni"]]; pb = np.r_[z["pj"], z["nj"]]
    py = np.r_[np.ones(z["pi"].size, bool), np.zeros(z["ni"].size, bool)]
    is_tf = z["lambert"].astype(bool)
    sp = S.StratifiedPairs(z["e_tf"], z["e_tg"], z["e_y"])
    return dict(genes=genes, n=n, pa=pa, pb=pb, py=py, is_tf=is_tf, tf=np.flatnonzero(is_tf),
                strata=z["det_dec"].astype(int), sp=sp)


def us_of(E):
    Xc = np.asarray(E, np.float64); Xc = Xc - Xc.mean(0, keepdims=True)
    U, Sv, Vt = np.linalg.svd(Xc, full_matrices=False)
    US = U * Sv
    sgn = np.sign(US[np.abs(US).argmax(0), np.arange(US.shape[1])]); sgn[sgn == 0] = 1
    return US * sgn, Sv, Xc


def rows_for(key, US, Xc, D, layer="static"):
    d = US.shape[1]; ks = [k for k in L.K_T if k < d] + [d]
    out = []
    # cosine equivalence of full U S and the centred embedding (spot check)
    A = S.l2norm_rows(Xc[:300]); Bm = S.l2norm_rows(US[:300]); chk = float(np.max(np.abs(A @ A.T - Bm @ Bm.T)))
    assert chk < 1e-6, (key, chk)
    for k in ks:
        En = L.unit_rep(US, slice(0, k))
        for t in L.TASKS:
            out.append(dict(model=key, layer=layer, task=t, k=k, is_full=k == d, dim_full=d, value=L.METRIC[t](En, D)))
    return out


def main():
    out = L.HERE / "static_curves.csv"
    D = bench(); rows = []; meta = {}
    pos = None
    for key, rel in STATIC.items():
        t0 = time.time()
        z = np.load(L.X.REPO / rel); g = np.array([str(x) for x in z["gene_names"]]); E = np.asarray(z["embeddings"], np.float32)
        nz = np.abs(E).sum(1) > 0; g, E = g[nz], E[nz]
        idx = {x: i for i, x in enumerate(g)}
        miss = [x for x in D["genes"] if x not in idx]; assert not miss, (key, len(miss))
        Eg = E[[idx[x] for x in D["genes"]]]
        US, Sv, Xc = us_of(Eg)
        rows += rows_for(key, US, Xc, D)
        meta[key] = dict(npz=rel, sha256=L.file_sha(L.X.REPO / rel), d=int(E.shape[1]), n_genes=D["n"],
                         seconds=time.time() - t0)
        print(key, E.shape, f"{time.time()-t0:.0f}s", flush=True)
        pd.DataFrame(rows).to_csv(out, index=False)
    # anchor: scGPT lung contextual, each layer on the same genes, then layer mean
    t0 = time.time()
    lg = json.loads((L.CACHE / "lung" / "obs_genes.json").read_text()); oi = np.load(L.CACHE / "lung" / "obs_idx.npy")
    p = {x: i for i, x in enumerate(lg)}; rid = oi[[p[x] for x in D["genes"]]]
    emb = L.load_emb("lung"); anc = []
    for l in range(L.N_LAYERS):
        US, Sv, Xc = us_of(np.asarray(emb[l], np.float64)[rid])
        anc += rows_for("scgpt_lung_contextual", US, Xc, D, layer=str(l))
    a = pd.DataFrame(anc); rows += anc
    m = a.groupby(["model", "task", "k", "is_full", "dim_full"], sort=False).value.mean().reset_index(); m["layer"] = "mean"
    df = pd.concat([pd.DataFrame(rows), m], ignore_index=True)
    full = df[df.is_full].set_index(["model", "layer", "task"]).value
    df["full_value"] = [full[(mm, ll, tt)] for mm, ll, tt in zip(df.model, df.layer, df.task)]
    df["diff_vs_full"] = df.value - df.full_value
    df.to_csv(out, index=False)
    meta["scgpt_lung_contextual"] = dict(seconds=time.time() - t0, note="12 layers, SVD fitted on the Tier B lung genes")
    meta["gene_set"] = dict(path=str(GS), sha256=L.file_sha(GS), n=D["n"], n_string=int(D["py"].sum()),
                            n_nonpairs=int((~D["py"]).sum()), n_lambert=int(D["is_tf"].sum()),
                            h3_pairs=int(D["sp"].y.size), h3_pos=int(D["sp"].y.sum()))
    (L.WORK / "static_meta.json").write_text(json.dumps(meta, indent=1))
    print("static done", flush=True)


if __name__ == "__main__":
    main()
