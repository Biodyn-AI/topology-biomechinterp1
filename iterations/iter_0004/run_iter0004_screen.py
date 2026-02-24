"""
iter_0004 - Multi-Hypothesis Screen: Intrinsic Dim, GO-CC, Cross-Layer CKA

Hypotheses:
  H01_intrinsic_dim : Intrinsic dimensionality (TwoNN) per layer — structured changes across depth.
  H02_go_cc         : Within-GO-slim-term kNN CC elevation vs feature-shuffle null (biological anchor).
  H03_layer_cka     : Cross-layer centered kernel alignment (CKA) — layer similarity structure.

Data: scGPT lung embeddings, cycle1 (seed42), 12 layers x 4803 genes x 512 dims.
Environment: subproject40-topology
"""

import numpy as np
import json
import csv
import os
import sys
import re
import urllib.request
from pathlib import Path
from scipy.spatial.distance import cdist

ITER_DIR = Path(__file__).parent
SUBP38_OUTPUTS = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_38_geometric_residual_stream_interpretability"
    "/implementation/outputs"
)
EMB_PATH = SUBP38_OUTPUTS / "cycle1_main" / "layer_gene_embeddings.npy"

RNG = np.random.default_rng(2024)
N_NULL_REPS = 15
K_NN = 10
PCA_DIMS = 20
N_GENES_SAMPLE = 400   # genes to subsample for speed


# ---------------------------------------------------------------------------
# Shared utilities
# ---------------------------------------------------------------------------

def pca_reduce(X, n_components):
    X = X - X.mean(axis=0, keepdims=True)
    _, S, Vt = np.linalg.svd(X, full_matrices=False)
    return X @ Vt[:n_components].T


def build_knn_adj(X, k):
    D = cdist(X, X)
    np.fill_diagonal(D, np.inf)
    n = X.shape[0]
    adj = np.zeros((n, n), dtype=bool)
    for i in range(n):
        nn = np.argpartition(D[i], k)[:k]
        adj[i, nn] = True
        adj[nn, i] = True
    return adj


def clustering_coefficient(adj):
    n = adj.shape[0]
    cc = []
    for i in range(n):
        nbrs = np.where(adj[i])[0]
        d = len(nbrs)
        if d < 2:
            cc.append(0.0)
            continue
        sub = adj[np.ix_(nbrs, nbrs)]
        edges = sub.sum() / 2
        cc.append(float(edges) / (d * (d - 1) / 2))
    return float(np.mean(cc))


def feature_shuffle(X, rng):
    Xs = X.copy()
    for j in range(X.shape[1]):
        Xs[:, j] = rng.permutation(X[:, j])
    return Xs


# ---------------------------------------------------------------------------
# H01: TwoNN Intrinsic Dimensionality per layer
# ---------------------------------------------------------------------------

def twonn_id(X):
    """TwoNN intrinsic dimension estimator (Facco et al. 2017)."""
    n = X.shape[0]
    D = cdist(X, X)
    np.fill_diagonal(D, np.inf)
    sorted_d = np.sort(D, axis=1)
    r1 = sorted_d[:, 0]
    r2 = sorted_d[:, 1]
    # mu = r2 / r1, filter mu > 1
    mu = r2 / np.maximum(r1, 1e-12)
    mu = mu[mu > 1]
    if len(mu) < 10:
        return np.nan
    # empirical cumulative
    mu_sorted = np.sort(mu)
    n_mu = len(mu_sorted)
    F = np.arange(1, n_mu + 1) / n_mu
    # linear regression in log space: log(1-F) = -d * log(mu)
    log_mu = np.log(mu_sorted)
    log_1mF = np.log(np.maximum(1 - F, 1e-9))
    # weighted OLS
    valid = log_1mF > -20
    if valid.sum() < 5:
        return np.nan
    lm = log_mu[valid]
    lf = log_1mF[valid]
    d_hat = -np.dot(lm, lf) / np.dot(lm, lm)
    return float(d_hat)


def run_h01_intrinsic_dim(embeddings):
    """Compute TwoNN ID for each layer on gene-subspace PCA projection."""
    n_layers = embeddings.shape[0]
    n_genes = embeddings.shape[1]

    # subsample genes
    idx = RNG.choice(n_genes, N_GENES_SAMPLE, replace=False)

    rows = []
    id_per_layer = []
    for layer in range(n_layers):
        X = embeddings[layer][idx]
        # PCA reduce
        Xpca = pca_reduce(X, PCA_DIMS)
        id_val = twonn_id(Xpca)
        id_per_layer.append(id_val)
        rows.append({"layer": layer, "intrinsic_dim_twonn": id_val})
        print(f"  Layer {layer:2d}: TwoNN ID = {id_val:.3f}", flush=True)

    # Null: shuffle each layer independently
    null_ids = []
    for rep in range(N_NULL_REPS):
        rep_ids = []
        for layer in range(n_layers):
            X = embeddings[layer][idx]
            Xs = feature_shuffle(X, RNG)
            Xpca = pca_reduce(Xs, PCA_DIMS)
            rep_ids.append(twonn_id(Xpca))
        null_ids.append(rep_ids)

    null_ids = np.array(null_ids)  # (n_reps, n_layers)
    null_mean = null_ids.mean(axis=0)
    null_std = null_ids.std(axis=0)

    id_arr = np.array(id_per_layer)
    z_scores = (id_arr - null_mean) / np.maximum(null_std, 1e-9)

    for i, row in enumerate(rows):
        row["null_mean"] = float(null_mean[i])
        row["null_std"] = float(null_std[i])
        row["z_score"] = float(z_scores[i])

    return rows, id_per_layer, null_ids


# ---------------------------------------------------------------------------
# H02: GO slim CC elevation within GO term neighbourhoods
# ---------------------------------------------------------------------------

GO_SLIM_TERMS = {
    # Selected high-coverage GO slim molecular/biological terms (human genes)
    "cell_cycle": ["CDK1","CDK2","CDK4","CDK6","CCNA2","CCNB1","CCND1","CCND3","CCNE1","CCNE2",
                   "CDC20","CDC25A","CDC25B","CDC25C","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7",
                   "PCNA","RB1","E2F1","E2F2","E2F3","CDKN1A","CDKN1B","CDKN2A","CDKN2B","TP53"],
    "dna_repair": ["BRCA1","BRCA2","RAD51","RAD52","RAD54L","ATM","ATR","CHEK1","CHEK2",
                   "PARP1","PARP2","LIG4","XRCC4","XRCC5","XRCC6","MLH1","MSH2","MSH6","PMS2",
                   "ERCC1","ERCC2","ERCC3","XPC","XPA","FANCD2","FANCA","FANCB","FANCC"],
    "apoptosis": ["BCL2","BCL2L1","MCL1","BAX","BAD","BID","CASP3","CASP7","CASP8","CASP9",
                  "CASP6","APAF1","CYCS","TP53","FAS","FASLG","TNFSF10","XIAP","DIABLO","BIK"],
    "immune_response": ["IL6","IL10","IL12A","TNF","IFNG","IL2","IL4","CXCL10","CCL2","CCL5",
                        "STAT1","STAT3","STAT4","STAT6","JAK1","JAK2","TLR4","TLR9","MYD88","NFKB1"],
    "transcription": ["TP53","MYC","MYCN","BRCA1","E2F1","SP1","AP1","CREB1","NFKB1","STAT1",
                      "CTNNB1","TCF7L2","HIF1A","EGR1","NR3C1","RELA","FOS","JUN","ATF3","KLF4"],
    "metabolism": ["GAPDH","PKM","LDHA","G6PD","PGAM1","ENO1","ALDOA","PGK1","TPI1","IDH1",
                   "IDH2","ACLY","FASN","SCD","HMGCR","FDFT1","SQLE","DHCR7","ACAT1","HADHB"],
    "signaling": ["EGFR","ERBB2","FGFR1","PDGFRA","MET","ALK","RET","KRAS","BRAF","PIK3CA",
                  "AKT1","MTOR","MAPK1","MAP2K1","JAK2","STAT3","PTEN","NF1","TSC1","TSC2"],
    "translation": ["EIF4E","EIF4A1","EIF4B","EIF4G1","EIF2A","EIF3A","RPL5","RPL11","RPL23",
                    "RPS6","RPS3","MTOR","RPS6KB1","EIF2AK3","EIF2AK4","PABPC1","YTHDF1","HNRNPA1"],
}


def get_gene_list_from_embeddings(emb_path):
    """Try to load gene list from associated metadata; fallback to synthetic index."""
    meta_candidates = [
        emb_path.parent / "gene_list.txt",
        emb_path.parent / "gene_names.txt",
        emb_path.parent / "cycle1_edge_dataset.tsv",
    ]
    for p in meta_candidates:
        if p.exists() and p.suffix == ".txt":
            with open(p) as f:
                genes = [l.strip() for l in f if l.strip()]
            print(f"  Loaded {len(genes)} genes from {p}", flush=True)
            return genes
    # Try to parse the edge dataset for gene names
    edge_p = emb_path.parent / "cycle1_edge_dataset.tsv"
    if edge_p.exists():
        genes_seen = set()
        genes_ordered = []
        with open(edge_p) as f:
            reader = csv.DictReader(f, delimiter="\t")
            cols = reader.fieldnames or []
            gene_cols = [c for c in cols if "gene" in c.lower() or "symbol" in c.lower()]
            for row in reader:
                for col in gene_cols:
                    g = row.get(col, "").strip()
                    if g and g not in genes_seen:
                        genes_seen.add(g)
                        genes_ordered.append(g)
        if genes_ordered:
            print(f"  Parsed {len(genes_ordered)} unique genes from {edge_p}", flush=True)
            return genes_ordered
    print("  No gene list found; using synthetic gene-index labels", flush=True)
    return None


def run_h02_go_cc(embeddings, gene_list):
    """Within-GO-term CC elevation vs feature-shuffle null per layer."""
    n_layers = embeddings.shape[0]

    # Map gene names -> indices in embedding matrix
    term_results = []
    global_rows = []

    for term_name, term_genes in GO_SLIM_TERMS.items():
        if gene_list is None:
            continue
        # Find indices
        gl_upper = {g.upper(): i for i, g in enumerate(gene_list)}
        idxs = [gl_upper[g.upper()] for g in term_genes if g.upper() in gl_upper]
        if len(idxs) < 8:
            print(f"  GO term {term_name}: only {len(idxs)} genes found, skipping", flush=True)
            continue
        print(f"  GO term {term_name}: {len(idxs)} genes matched", flush=True)

        for layer in range(n_layers):
            X = embeddings[layer][idxs]
            Xpca = pca_reduce(X, min(PCA_DIMS, X.shape[0] - 1, X.shape[1]))
            k = min(K_NN, len(idxs) - 1)
            adj = build_knn_adj(Xpca, k)
            cc_obs = clustering_coefficient(adj)

            # Null
            null_ccs = []
            for _ in range(N_NULL_REPS):
                Xs = feature_shuffle(X, RNG)
                Xsp = pca_reduce(Xs, min(PCA_DIMS, X.shape[0] - 1, X.shape[1]))
                adj_n = build_knn_adj(Xsp, k)
                null_ccs.append(clustering_coefficient(adj_n))

            null_arr = np.array(null_ccs)
            z = (cc_obs - null_arr.mean()) / max(null_arr.std(), 1e-9)
            global_rows.append({
                "term": term_name, "layer": layer, "n_genes": len(idxs),
                "cc_obs": cc_obs, "null_mean": null_arr.mean(),
                "null_std": null_arr.std(), "z_score": z
            })

    return global_rows


# ---------------------------------------------------------------------------
# H03: Cross-layer CKA (centered kernel alignment)
# ---------------------------------------------------------------------------

def linear_cka(X, Y):
    """Linear CKA between two (n_samples x d) matrices."""
    n = X.shape[0]
    # Center
    Xc = X - X.mean(0)
    Yc = Y - Y.mean(0)
    # Gram matrices
    Kx = Xc @ Xc.T
    Ky = Yc @ Yc.T
    num = np.linalg.norm(Ky @ Kx, 'fro') ** 2
    denom = np.linalg.norm(Kx @ Kx, 'fro') * np.linalg.norm(Ky @ Ky, 'fro')
    if denom < 1e-12:
        return np.nan
    return float(num / denom)


def run_h03_cross_layer_cka(embeddings):
    """Compute pairwise CKA between all 12 scGPT layers on gene subspace."""
    n_layers = embeddings.shape[0]
    n_genes = embeddings.shape[1]

    # subsample genes, use same subset for all layers
    idx = RNG.choice(n_genes, N_GENES_SAMPLE, replace=False)

    # PCA reduce each layer independently
    Xpca_list = []
    for layer in range(n_layers):
        X = embeddings[layer][idx]
        Xpca_list.append(pca_reduce(X, PCA_DIMS))

    cka_matrix = np.zeros((n_layers, n_layers))
    for i in range(n_layers):
        for j in range(n_layers):
            if i <= j:
                cka_matrix[i, j] = linear_cka(Xpca_list[i], Xpca_list[j])
            else:
                cka_matrix[i, j] = cka_matrix[j, i]

    # Null: feature-shuffle layer 0 vs all layers
    null_cka_vs0 = []
    for _ in range(N_NULL_REPS):
        Xs = feature_shuffle(embeddings[0][idx], RNG)
        Xsp = pca_reduce(Xs, PCA_DIMS)
        vals = [linear_cka(Xsp, Xpca_list[j]) for j in range(n_layers)]
        null_cka_vs0.append(vals)

    null_arr = np.array(null_cka_vs0)  # (n_reps, n_layers)

    rows = []
    for i in range(n_layers):
        for j in range(n_layers):
            rows.append({"layer_i": i, "layer_j": j, "cka": float(cka_matrix[i, j])})

    # adjacent layer CKA (diagonal-1)
    adj_ckas = [cka_matrix[i, i+1] for i in range(n_layers - 1)]

    return rows, cka_matrix, null_arr, adj_ckas


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Loading embeddings...", flush=True)
    embeddings = np.load(EMB_PATH)
    print(f"  Shape: {embeddings.shape}", flush=True)

    # Gene list
    gene_list = get_gene_list_from_embeddings(EMB_PATH)
    n_found = sum(1 for g in sum(GO_SLIM_TERMS.values(), [])
                  if gene_list and g.upper() in {x.upper() for x in gene_list}) if gene_list else 0
    print(f"  GO genes found in embedding gene list: {n_found}", flush=True)

    # ---- H01: Intrinsic Dimensionality ----
    print("\n=== H01: TwoNN Intrinsic Dimensionality per layer ===", flush=True)
    h01_rows, id_per_layer, null_ids = run_h01_intrinsic_dim(embeddings)

    id_arr = np.array(id_per_layer)
    z_arr = np.array([r["z_score"] for r in h01_rows])
    print(f"\n  ID range: {id_arr.min():.2f} – {id_arr.max():.2f}")
    print(f"  Z-score range: {z_arr.min():.2f} – {z_arr.max():.2f}")
    print(f"  Mean Z: {z_arr.mean():.2f}")
    print(f"  Layers with |z| > 2: {(np.abs(z_arr) > 2).sum()}/{len(z_arr)}")

    h01_csv = ITER_DIR / "h01_intrinsic_dim_per_layer.csv"
    with open(h01_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["layer", "intrinsic_dim_twonn", "null_mean", "null_std", "z_score"])
        writer.writeheader()
        writer.writerows(h01_rows)
    print(f"  Saved: {h01_csv}")

    # ---- H02: GO-CC ----
    print("\n=== H02: GO slim CC elevation vs feature-shuffle null ===", flush=True)
    h02_rows = run_h02_go_cc(embeddings, gene_list)

    if h02_rows:
        z_vals = np.array([r["z_score"] for r in h02_rows])
        sig_pos = (np.array([r["z_score"] for r in h02_rows]) > 1.96).sum()
        print(f"  Total layer-term tests: {len(h02_rows)}")
        print(f"  Significant positive (z>1.96): {sig_pos}/{len(h02_rows)}")
        print(f"  Mean z: {z_vals.mean():.3f}, Median z: {np.median(z_vals):.3f}")

        h02_csv = ITER_DIR / "h02_go_cc_by_term_layer.csv"
        with open(h02_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["term", "layer", "n_genes", "cc_obs", "null_mean", "null_std", "z_score"])
            writer.writeheader()
            writer.writerows(h02_rows)
        print(f"  Saved: {h02_csv}")

        # Summary by term
        from collections import defaultdict
        by_term = defaultdict(list)
        for r in h02_rows:
            by_term[r["term"]].append(r["z_score"])
        h02_summary = [{"term": t, "mean_z": np.mean(zs), "n_layers_sig": sum(z > 1.96 for z in zs), "n_layers": len(zs)}
                       for t, zs in by_term.items()]
        h02_sum_csv = ITER_DIR / "h02_go_cc_term_summary.csv"
        with open(h02_sum_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["term", "mean_z", "n_layers_sig", "n_layers"])
            writer.writeheader()
            writer.writerows(h02_summary)
        print(f"  Saved: {h02_sum_csv}")
        for s in sorted(h02_summary, key=lambda x: -x["mean_z"]):
            print(f"    {s['term']:25s}: mean_z={s['mean_z']:+.2f}, sig_layers={s['n_layers_sig']}/{s['n_layers']}")
    else:
        print("  No GO terms matched — gene list not available. H02 blocked.", flush=True)
        h02_rows = []

    # ---- H03: Cross-layer CKA ----
    print("\n=== H03: Cross-layer linear CKA ===", flush=True)
    h03_rows, cka_matrix, null_cka, adj_ckas = run_h03_cross_layer_cka(embeddings)

    print(f"  Adjacent-layer CKA (diagonal): {[f'{v:.3f}' for v in adj_ckas]}")
    print(f"  Mean adjacent CKA: {np.mean(adj_ckas):.3f}")
    null_diag = null_cka.mean(axis=0)
    print(f"  Null CKA vs layer0 (mean over {N_NULL_REPS} reps): {null_diag.mean():.4f}")
    obs_vs0 = np.array([cka_matrix[0, j] for j in range(embeddings.shape[0])])
    z_cka = (obs_vs0 - null_diag) / np.maximum(null_cka.std(axis=0), 1e-9)
    print(f"  Z-score (observed vs null) for layer0 row: {[f'{v:.1f}' for v in z_cka]}")

    h03_csv = ITER_DIR / "h03_cross_layer_cka_matrix.csv"
    with open(h03_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["layer_i", "layer_j", "cka"])
        writer.writeheader()
        writer.writerows(h03_rows)
    print(f"  Saved: {h03_csv}")

    np.save(ITER_DIR / "h03_cka_matrix.npy", cka_matrix)
    np.save(ITER_DIR / "h03_cka_null.npy", null_cka)
    print(f"  Saved: h03_cka_matrix.npy, h03_cka_null.npy")

    # ---- Summary JSON ----
    h01_status = "tested"
    h01_direction = "positive" if z_arr.mean() > 2 else ("negative" if z_arr.mean() < -2 else "inconclusive")
    h01_decision = "promising" if z_arr.mean() > 2 else ("negative" if z_arr.mean() < -2 else "inconclusive")

    h02_z_all = np.array([r["z_score"] for r in h02_rows]) if h02_rows else np.array([])
    h02_status = "tested" if h02_rows else "blocked"
    h02_direction = ("positive" if len(h02_z_all) > 0 and h02_z_all.mean() > 1.5 else
                     ("negative" if len(h02_z_all) > 0 and h02_z_all.mean() < -1.5 else "inconclusive"))
    h02_decision = ("promising" if h02_direction == "positive" else
                    ("negative" if h02_direction == "negative" else "inconclusive"))

    # CKA: adjacent layers should be much higher than null
    mean_adj = float(np.mean(adj_ckas))
    z_cka_mean = float(z_cka[1:].mean())  # skip self-comparison
    h03_direction = "positive" if mean_adj > 0.5 and z_cka_mean > 2 else "inconclusive"
    h03_decision = "promising" if h03_direction == "positive" else "inconclusive"

    screen = {
        "iteration": "iter_0004",
        "hypotheses": [
            {
                "id": "H01",
                "name": "TwoNN intrinsic dimensionality per layer",
                "family": "intrinsic_dimensionality",
                "split_regime": "other",
                "novelty_type": "new_family",
                "lineage": "none",
                "method": (f"TwoNN estimator on PCA-{PCA_DIMS} projections of {N_GENES_SAMPLE} "
                           f"sampled scGPT gene embeddings, 12 layers, seed42; "
                           f"{N_NULL_REPS} feature-shuffle null replicates."),
                "status": h01_status,
                "primary_metric": "mean_z_score_vs_null",
                "result_value": f"mean_z={z_arr.mean():.2f}, id_range=[{id_arr.min():.1f},{id_arr.max():.1f}]",
                "result_direction": h01_direction,
                "artifact_paths": ["h01_intrinsic_dim_per_layer.csv"],
                "decision": h01_decision,
                "next_action": ("If positive: test ID monotonicity and compare to Geneformer; "
                                "if inconclusive: try MADA or local PCA estimator"),
                "retired": False
            },
            {
                "id": "H02",
                "name": "Within-GO-slim-term kNN CC elevation",
                "family": "module_structure",
                "split_regime": "other",
                "novelty_type": "new_method",
                "lineage": "iter_0002_H01",
                "method": (f"kNN (k={K_NN}) CC on PCA-{PCA_DIMS} embeddings of GO slim gene subsets "
                           f"({len(GO_SLIM_TERMS)} terms), vs {N_NULL_REPS} feature-shuffle replicates; "
                           f"12 layers, seed42."),
                "status": h02_status,
                "primary_metric": "mean_z_score_within_go_terms",
                "result_value": (f"mean_z={h02_z_all.mean():.3f}, "
                                 f"sig_tests={int((h02_z_all > 1.96).sum())}/{len(h02_z_all)}"
                                 if len(h02_z_all) > 0 else "blocked: no gene list"),
                "result_direction": h02_direction if len(h02_z_all) > 0 else "inconclusive",
                "artifact_paths": (["h02_go_cc_by_term_layer.csv", "h02_go_cc_term_summary.csv"]
                                   if h02_rows else []),
                "decision": h02_decision if len(h02_z_all) > 0 else "inconclusive",
                "next_action": ("If positive: extend to full GO annotations; "
                                "if gene list missing: extract gene list from edge_dataset and rerun"),
                "retired": False
            },
            {
                "id": "H03",
                "name": "Cross-layer linear CKA alignment",
                "family": "cross_model_alignment",
                "split_regime": "other",
                "novelty_type": "new_family",
                "lineage": "none",
                "method": (f"Linear CKA on PCA-{PCA_DIMS} projections of {N_GENES_SAMPLE} "
                           f"sampled gene embeddings; all 12x12 layer pairs, seed42; "
                           f"{N_NULL_REPS} feature-shuffle null replicates for layer0 row."),
                "status": "tested",
                "primary_metric": "mean_adjacent_layer_cka",
                "result_value": f"mean_adj_cka={mean_adj:.3f}, mean_z_vs_null={z_cka_mean:.2f}",
                "result_direction": h03_direction,
                "artifact_paths": ["h03_cross_layer_cka_matrix.csv", "h03_cka_matrix.npy", "h03_cka_null.npy"],
                "decision": h03_decision,
                "next_action": ("If positive: compare CKA gradient across depth as signal for representation drift; "
                                "compare to Geneformer layer CKA if embeddings become available"),
                "retired": False
            }
        ]
    }

    screen_json = ITER_DIR / "executor_hypothesis_screen.json"
    with open(screen_json, "w") as f:
        json.dump(screen, f, indent=2)
    print(f"\nSaved: {screen_json}", flush=True)

    print("\n=== DONE ===", flush=True)
    return screen, h01_rows, h02_rows, h03_rows


if __name__ == "__main__":
    main()
