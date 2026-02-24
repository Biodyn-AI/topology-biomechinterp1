"""
H02b: TRRUST co-regulatory clustering in scGPT embeddings.
Test: genes regulated by the same TF should have higher kNN CC (or lower mean
pairwise distance) within their co-target group, compared to feature-shuffle null.

Uses the 209 genes mapped in cycle1_edge_dataset.tsv with their embedding indices.
"""
import numpy as np
import csv
import json
from pathlib import Path
from scipy.spatial.distance import cdist
from collections import defaultdict

ITER_DIR = Path(__file__).parent
SUBP38_OUTPUTS = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_38_geometric_residual_stream_interpretability"
    "/implementation/outputs"
)
EMB_PATH = SUBP38_OUTPUTS / "cycle1_main" / "layer_gene_embeddings.npy"
EDGE_PATH = SUBP38_OUTPUTS / "cycle1_main" / "cycle1_edge_dataset.tsv"

RNG = np.random.default_rng(2024)
N_NULL_REPS = 15
K_NN = 6  # smaller k for small groups
PCA_DIMS = 20
MIN_GROUP_SIZE = 5  # minimum co-targets per TF


def pca_reduce(X, n_components):
    X = X - X.mean(0)
    _, S, Vt = np.linalg.svd(X, full_matrices=False)
    return X @ Vt[:n_components].T


def mean_pairwise_dist(X):
    D = cdist(X, X)
    n = X.shape[0]
    return D[np.triu_indices(n, k=1)].mean()


def feature_shuffle(X, rng):
    Xs = X.copy()
    for j in range(X.shape[1]):
        Xs[:, j] = rng.permutation(X[:, j])
    return Xs


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
    cc = []
    for i in range(adj.shape[0]):
        nbrs = np.where(adj[i])[0]
        d = len(nbrs)
        if d < 2:
            cc.append(0.0)
            continue
        sub = adj[np.ix_(nbrs, nbrs)]
        cc.append(float(sub.sum() / 2) / (d * (d - 1) / 2))
    return float(np.mean(cc))


def main():
    print("Loading embeddings...", flush=True)
    embeddings = np.load(EMB_PATH)
    n_layers = embeddings.shape[0]

    # Load gene -> idx mapping
    gene_idx = {}
    tf_targets = defaultdict(set)
    with open(EDGE_PATH) as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            gene_idx[row['source']] = int(row['source_idx'])
            gene_idx[row['target']] = int(row['target_idx'])
            tf_targets[row['source']].add(row['target'])

    print(f"  Unique genes: {len(gene_idx)}", flush=True)
    print(f"  TF groups: {len(tf_targets)}", flush=True)

    # Filter TFs with enough targets
    valid_tfs = {tf: targets for tf, targets in tf_targets.items()
                 if len(targets) >= MIN_GROUP_SIZE and
                 all(t in gene_idx for t in targets)}
    print(f"  TFs with >={MIN_GROUP_SIZE} mapped targets: {len(valid_tfs)}", flush=True)

    rows = []
    for tf, targets in sorted(valid_tfs.items()):
        target_list = sorted(targets)
        idxs = [gene_idx[t] for t in target_list]

        for layer in range(n_layers):
            X = embeddings[layer][idxs]  # (n_targets, 512)
            n_comp = min(PCA_DIMS, X.shape[0] - 1, X.shape[1])
            Xpca = pca_reduce(X, n_comp)
            k = min(K_NN, len(idxs) - 1)

            # Observed: mean pairwise distance
            obs_dist = mean_pairwise_dist(Xpca)

            # CC (if enough genes)
            obs_cc = None
            if len(idxs) >= 4:
                adj = build_knn_adj(Xpca, k)
                obs_cc = clustering_coefficient(adj)

            # Null: feature shuffle on full layer, then recompute for these genes
            null_dists = []
            null_ccs = []
            for _ in range(N_NULL_REPS):
                Xs = feature_shuffle(embeddings[layer], RNG)
                Xsp = pca_reduce(Xs[idxs], n_comp)
                null_dists.append(mean_pairwise_dist(Xsp))
                if obs_cc is not None:
                    adj_n = build_knn_adj(Xsp, k)
                    null_ccs.append(clustering_coefficient(adj_n))

            null_d = np.array(null_dists)
            z_dist = (obs_dist - null_d.mean()) / max(null_d.std(), 1e-9)
            # Negative z_dist means co-targets are closer in real embedding (positive signal)

            z_cc = None
            if null_ccs:
                null_cc_arr = np.array(null_ccs)
                z_cc = (obs_cc - null_cc_arr.mean()) / max(null_cc_arr.std(), 1e-9)

            rows.append({
                "tf": tf,
                "n_targets": len(idxs),
                "layer": layer,
                "obs_dist": obs_dist,
                "null_dist_mean": null_d.mean(),
                "z_dist": z_dist,
                "obs_cc": obs_cc if obs_cc is not None else np.nan,
                "null_cc_mean": np.mean(null_ccs) if null_ccs else np.nan,
                "z_cc": z_cc if z_cc is not None else np.nan,
            })

    print(f"\n  Total tests: {len(rows)}", flush=True)
    z_dists = np.array([r["z_dist"] for r in rows])
    z_ccs = np.array([r["z_cc"] for r in rows if not np.isnan(r["z_cc"])])

    # Closer in real = negative z_dist
    sig_closer = (z_dists < -1.96).sum()
    print(f"  Sig closer (z_dist < -1.96): {sig_closer}/{len(z_dists)}", flush=True)
    print(f"  Mean z_dist: {z_dists.mean():.3f}", flush=True)
    if len(z_ccs):
        sig_cc_pos = (z_ccs > 1.96).sum()
        print(f"  Sig CC positive (z_cc > 1.96): {sig_cc_pos}/{len(z_ccs)}", flush=True)
        print(f"  Mean z_cc: {z_ccs.mean():.3f}", flush=True)

    # Save
    csv_path = ITER_DIR / "h02b_trrust_cotarget_by_tf_layer.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    print(f"  Saved: {csv_path}", flush=True)

    # Summary by TF
    from collections import defaultdict as dd
    tf_z_dists = dd(list)
    tf_z_ccs = dd(list)
    for r in rows:
        tf_z_dists[r["tf"]].append(r["z_dist"])
        if not np.isnan(r["z_cc"]):
            tf_z_ccs[r["tf"]].append(r["z_cc"])

    summary_rows = []
    for tf in sorted(tf_z_dists):
        zd = tf_z_dists[tf]
        zc = tf_z_ccs.get(tf, [])
        summary_rows.append({
            "tf": tf,
            "n_targets": valid_tfs[tf].__len__(),
            "mean_z_dist": np.mean(zd),
            "n_sig_closer": sum(z < -1.96 for z in zd),
            "mean_z_cc": np.mean(zc) if zc else np.nan,
            "n_sig_cc_pos": sum(z > 1.96 for z in zc) if zc else 0,
        })

    sum_csv = ITER_DIR / "h02b_trrust_tf_summary.csv"
    with open(sum_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"  Saved: {sum_csv}", flush=True)

    print("\nTop TFs by co-target clustering (most negative z_dist):", flush=True)
    for s in sorted(summary_rows, key=lambda x: x["mean_z_dist"])[:10]:
        print(f"  {s['tf']:15s}: z_dist={s['mean_z_dist']:+.2f}, n_targets={s['n_targets']}, sig_closer={s['n_sig_closer']}/12")

    # Update hypothesis screen JSON
    screen_path = ITER_DIR / "executor_hypothesis_screen.json"
    with open(screen_path) as f:
        screen = json.load(f)

    h02b_direction = "positive" if z_dists.mean() < -1.5 else ("negative" if z_dists.mean() > 1.5 else "inconclusive")
    h02b_decision = "promising" if h02b_direction == "positive" else ("negative" if h02b_direction == "negative" else "inconclusive")

    # Update H02 entry
    for h in screen["hypotheses"]:
        if h["id"] == "H02":
            h["name"] = "TRRUST co-target kNN distance/CC elevation"
            h["method"] = (f"Mean pairwise distance + kNN CC for co-targets of each TF "
                           f"({len(valid_tfs)} TFs, >={MIN_GROUP_SIZE} mapped targets) "
                           f"vs {N_NULL_REPS} feature-shuffle replicates; 12 layers, seed42.")
            h["status"] = "tested"
            h["primary_metric"] = "mean_z_dist"
            h["result_value"] = (f"mean_z_dist={z_dists.mean():.3f}, "
                                 f"sig_closer={sig_closer}/{len(z_dists)}, "
                                 f"mean_z_cc={z_ccs.mean():.3f}")
            h["result_direction"] = h02b_direction
            h["artifact_paths"] = ["h02b_trrust_cotarget_by_tf_layer.csv", "h02b_trrust_tf_summary.csv"]
            h["decision"] = h02b_decision
            h["next_action"] = ("If positive: test across all 3 seeds and layers; "
                                "extend to GO term co-membership distance analysis")
            break

    with open(screen_path, "w") as f:
        json.dump(screen, f, indent=2)
    print(f"  Updated: {screen_path}", flush=True)

    print("\n=== H02b DONE ===", flush=True)


if __name__ == "__main__":
    main()
