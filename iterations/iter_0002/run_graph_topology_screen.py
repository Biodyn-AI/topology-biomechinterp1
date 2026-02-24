"""
iter_0001 - Graph Topology Surrogates Hypothesis Screen (TEST MODE)

Hypothesis: kNN graph topology metrics (clustering coefficient, transitivity,
assortativity by expression rank) computed on scGPT layer embeddings show
significant structure versus feature-shuffle null.

Family: graph_topology (materially new family vs H01-H12 which focused on
persistent homology and rewiring nulls)

Data: scGPT lung embeddings from subproject_38 cycle1 outputs (3 seeds).
Protocol: PCA(20) -> kNN graph (k=10) -> graph topology metrics -> feature-shuffle null (15 replicates).
"""

import numpy as np
import json
import csv
import os
import sys
from pathlib import Path
from scipy.spatial.distance import cdist
from scipy.sparse import csr_matrix

ITER_DIR = Path(__file__).parent
SUBP38_OUTPUTS = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/subproject_38_geometric_residual_stream_interpretability/implementation/outputs")

SEEDS = {
    "seed42": SUBP38_OUTPUTS / "cycle1_main" / "layer_gene_embeddings.npy",
    "seed43": SUBP38_OUTPUTS / "cycle1_seed43" / "layer_gene_embeddings.npy",
    "seed44": SUBP38_OUTPUTS / "cycle1_seed44" / "layer_gene_embeddings.npy",
}

RNG_SEED = 2024
N_GENES_SAMPLE = 300      # genes per test
K_NN = 10                 # kNN graph k
N_NULL_REPS = 15          # feature-shuffle null replicates
N_LAYERS_MAX = 12         # maximum layers to test (cap for speed in test mode)
PCA_DIMS = 20

# ---- utilities ---------------------------------------------------------------

def pca_reduce(X, n_components):
    """Simple PCA via SVD."""
    X = X - X.mean(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    return U[:, :n_components] * S[:n_components]


def build_knn_graph_adj(X, k):
    """Build undirected kNN adjacency as dense boolean matrix (fast for small n)."""
    n = X.shape[0]
    D = cdist(X, X)
    np.fill_diagonal(D, np.inf)
    adj = np.zeros((n, n), dtype=bool)
    for i in range(n):
        nn = np.argpartition(D[i], k)[:k]
        adj[i, nn] = True
        adj[nn, i] = True  # symmetrize
    return adj


def clustering_coefficient(adj):
    """Mean local clustering coefficient over all nodes."""
    n = adj.shape[0]
    deg = adj.sum(axis=1).astype(float)
    cc_values = []
    for i in range(n):
        d = int(deg[i])
        if d < 2:
            cc_values.append(0.0)
            continue
        neighbors = np.where(adj[i])[0]
        # count edges among neighbors
        sub = adj[np.ix_(neighbors, neighbors)]
        edges_among = sub.sum() / 2  # undirected
        possible = d * (d - 1) / 2
        cc_values.append(float(edges_among) / possible)
    return float(np.mean(cc_values))


def transitivity(adj):
    """Global transitivity (fraction of closed triangles)."""
    n = adj.shape[0]
    adj_f = adj.astype(float)
    # triangles: trace(A^3)/6 for undirected
    A2 = adj_f @ adj_f
    triangles = float(np.trace(adj_f @ A2)) / 6.0
    # possible: sum_i deg_i*(deg_i-1)/2
    deg = adj_f.sum(axis=1)
    wedges = float(np.sum(deg * (deg - 1)) / 2.0)
    if wedges < 1e-9:
        return 0.0
    return 3.0 * triangles / wedges


def graph_metrics(X_pca, k):
    adj = build_knn_graph_adj(X_pca, k)
    cc = clustering_coefficient(adj)
    tr = transitivity(adj)
    return {"clustering_coeff": cc, "transitivity": tr}


def feature_shuffle(X, rng):
    """Shuffle each feature (column) independently."""
    X_sh = X.copy()
    for j in range(X.shape[1]):
        X_sh[:, j] = rng.permutation(X_sh[:, j])
    return X_sh


# ---- main loop ---------------------------------------------------------------

def run_screen():
    rng = np.random.default_rng(RNG_SEED)
    results_by_seed = []

    for seed_name, emb_path in SEEDS.items():
        if not emb_path.exists():
            print(f"[WARN] Missing: {emb_path} — skipping seed {seed_name}", file=sys.stderr)
            continue
        emb = np.load(str(emb_path))  # shape: (n_layers, n_genes, d)
        print(f"[INFO] Loaded {seed_name}: shape={emb.shape}", flush=True)

        n_layers = min(emb.shape[0], N_LAYERS_MAX)
        n_genes = emb.shape[1]

        # sample gene subset for speed
        gene_idx = rng.choice(n_genes, size=min(N_GENES_SAMPLE, n_genes), replace=False)

        for layer_idx in range(n_layers):
            X_layer = emb[layer_idx][gene_idx]  # (N_GENES_SAMPLE, d)

            # PCA
            X_pca = pca_reduce(X_layer, PCA_DIMS)

            # observed metrics
            obs = graph_metrics(X_pca, K_NN)

            # null replicates
            null_cc = []
            null_tr = []
            for _ in range(N_NULL_REPS):
                X_sh = feature_shuffle(X_layer, rng)
                X_sh_pca = pca_reduce(X_sh, PCA_DIMS)
                nm = graph_metrics(X_sh_pca, K_NN)
                null_cc.append(nm["clustering_coeff"])
                null_tr.append(nm["transitivity"])

            null_cc = np.array(null_cc)
            null_tr = np.array(null_tr)

            # effect sizes
            cc_delta = obs["clustering_coeff"] - float(null_cc.mean())
            tr_delta = obs["transitivity"] - float(null_tr.mean())
            cc_z = cc_delta / (float(null_cc.std()) + 1e-12)
            tr_z = tr_delta / (float(null_tr.std()) + 1e-12)

            # empirical p-values (one-tailed: observed > null)
            cc_p = float((null_cc >= obs["clustering_coeff"]).sum() + 1) / (N_NULL_REPS + 1)
            tr_p = float((null_tr >= obs["transitivity"]).sum() + 1) / (N_NULL_REPS + 1)

            row = {
                "seed": seed_name,
                "layer": layer_idx,
                "obs_cc": round(obs["clustering_coeff"], 6),
                "obs_tr": round(obs["transitivity"], 6),
                "null_cc_mean": round(float(null_cc.mean()), 6),
                "null_tr_mean": round(float(null_tr.mean()), 6),
                "null_cc_std": round(float(null_cc.std()), 6),
                "null_tr_std": round(float(null_tr.std()), 6),
                "cc_delta": round(cc_delta, 6),
                "tr_delta": round(tr_delta, 6),
                "cc_z": round(cc_z, 4),
                "tr_z": round(tr_z, 4),
                "cc_p_empirical": round(cc_p, 4),
                "tr_p_empirical": round(tr_p, 4),
            }
            results_by_seed.append(row)
            print(f"  layer={layer_idx:2d}  cc_delta={cc_delta:+.4f} (z={cc_z:+.2f}, p={cc_p:.3f})  "
                  f"tr_delta={tr_delta:+.4f} (z={tr_z:+.2f}, p={tr_p:.3f})", flush=True)

    if not results_by_seed:
        print("[ERROR] No results produced — all seeds missing.", file=sys.stderr)
        sys.exit(1)

    # save per-seed-layer csv
    out_csv = ITER_DIR / "graph_topology_knn_by_seed_layer.csv"
    fieldnames = list(results_by_seed[0].keys())
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results_by_seed)
    print(f"\n[SAVED] {out_csv}", flush=True)

    # aggregate summary by layer
    all_rows = results_by_seed
    layer_set = sorted(set(r["layer"] for r in all_rows))
    layer_summary = []
    for l in layer_set:
        rows_l = [r for r in all_rows if r["layer"] == l]
        cc_deltas = [r["cc_delta"] for r in rows_l]
        tr_deltas = [r["tr_delta"] for r in rows_l]
        cc_ps = [r["cc_p_empirical"] for r in rows_l]
        tr_ps = [r["tr_p_empirical"] for r in rows_l]
        n_cc_sig = sum(p < 0.05 for p in cc_ps)
        n_tr_sig = sum(p < 0.05 for p in tr_ps)
        layer_summary.append({
            "layer": l,
            "n_seeds": len(rows_l),
            "mean_cc_delta": round(float(np.mean(cc_deltas)), 6),
            "mean_tr_delta": round(float(np.mean(tr_deltas)), 6),
            "mean_cc_z": round(float(np.mean([r["cc_z"] for r in rows_l])), 4),
            "mean_tr_z": round(float(np.mean([r["tr_z"] for r in rows_l])), 4),
            "n_seeds_cc_sig": n_cc_sig,
            "n_seeds_tr_sig": n_tr_sig,
        })

    out_layer_csv = ITER_DIR / "graph_topology_knn_layer_summary.csv"
    with open(out_layer_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(layer_summary[0].keys()))
        writer.writeheader()
        writer.writerows(layer_summary)
    print(f"[SAVED] {out_layer_csv}", flush=True)

    # aggregate stats
    all_cc_deltas = [r["cc_delta"] for r in all_rows]
    all_tr_deltas = [r["tr_delta"] for r in all_rows]
    n_total = len(all_rows)
    n_cc_sig_total = sum(r["cc_p_empirical"] < 0.05 for r in all_rows)
    n_tr_sig_total = sum(r["tr_p_empirical"] < 0.05 for r in all_rows)

    summary = {
        "iteration": "iter_0001",
        "domain": "lung_scgpt",
        "n_seeds": len(SEEDS),
        "n_genes_sampled": N_GENES_SAMPLE,
        "k_nn": K_NN,
        "n_null_reps": N_NULL_REPS,
        "pca_dims": PCA_DIMS,
        "n_layer_tests": n_total,
        "mean_cc_delta": round(float(np.mean(all_cc_deltas)), 6),
        "mean_tr_delta": round(float(np.mean(all_tr_deltas)), 6),
        "n_cc_significant_p05": n_cc_sig_total,
        "n_tr_significant_p05": n_tr_sig_total,
        "frac_cc_sig": round(n_cc_sig_total / n_total, 4),
        "frac_tr_sig": round(n_tr_sig_total / n_total, 4),
        "verdict_cc": "positive" if n_cc_sig_total / n_total > 0.5 else ("mixed" if n_cc_sig_total > 0 else "negative"),
        "verdict_tr": "positive" if n_tr_sig_total / n_total > 0.5 else ("mixed" if n_tr_sig_total > 0 else "negative"),
    }

    out_summary = ITER_DIR / "graph_topology_knn_summary.json"
    with open(out_summary, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[SAVED] {out_summary}", flush=True)

    print("\n=== SUMMARY ===")
    for k, v in summary.items():
        print(f"  {k}: {v}")

    return summary


if __name__ == "__main__":
    run_screen()
