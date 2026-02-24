"""
H02b (fast version): TRRUST co-regulatory clustering in scGPT embeddings.
Null: randomly sample same-sized gene groups from the 209 named genes.
This is faster than full-layer shuffle because we only work on small subsets.
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
N_NULL_REPS = 20
MIN_GROUP_SIZE = 5


def pca_reduce(X, n_components):
    X = X - X.mean(0)
    _, S, Vt = np.linalg.svd(X, full_matrices=False)
    return X @ Vt[:n_components].T


def mean_pairwise_dist(X):
    D = cdist(X, X)
    n = X.shape[0]
    if n < 2:
        return np.nan
    return D[np.triu_indices(n, k=1)].mean()


def main():
    print("Loading embeddings...", flush=True)
    embeddings = np.load(EMB_PATH)
    n_layers = embeddings.shape[0]

    # Load gene mapping
    gene_idx = {}
    tf_targets = defaultdict(set)
    with open(EDGE_PATH) as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            gene_idx[row['source']] = int(row['source_idx'])
            gene_idx[row['target']] = int(row['target_idx'])
            tf_targets[row['source']].add(row['target'])

    # All mapped gene indices (pool for null sampling)
    all_gene_idxs = sorted(gene_idx.values())
    print(f"  Gene pool for null: {len(all_gene_idxs)} genes", flush=True)

    valid_tfs = {tf: sorted(targets) for tf, targets in tf_targets.items()
                 if len(targets) >= MIN_GROUP_SIZE and
                 all(t in gene_idx for t in targets)}
    print(f"  Valid TFs: {len(valid_tfs)}", flush=True)

    rows = []
    for tf, targets in sorted(valid_tfs.items()):
        idxs = [gene_idx[t] for t in targets]
        n_tgt = len(idxs)
        print(f"  TF {tf}: {n_tgt} targets", flush=True)

        for layer in range(n_layers):
            X_all = embeddings[layer]  # (4803, 512)
            X_tgt = X_all[idxs]        # (n_tgt, 512)
            n_comp = min(20, n_tgt - 1)

            # Observed: standardize within group and reduce
            Xpca = pca_reduce(X_tgt, n_comp)
            obs_dist = mean_pairwise_dist(Xpca)

            # Null: random same-size groups from the pool
            null_dists = []
            for _ in range(N_NULL_REPS):
                rand_idxs = RNG.choice(all_gene_idxs, n_tgt, replace=False).tolist()
                X_rand = X_all[rand_idxs]
                Xrp = pca_reduce(X_rand, n_comp)
                null_dists.append(mean_pairwise_dist(Xrp))

            null_arr = np.array(null_dists)
            z_dist = (obs_dist - null_arr.mean()) / max(null_arr.std(), 1e-9)

            rows.append({
                "tf": tf, "n_targets": n_tgt, "layer": layer,
                "obs_dist": obs_dist,
                "null_dist_mean": float(null_arr.mean()),
                "null_dist_std": float(null_arr.std()),
                "z_dist": float(z_dist),
            })

    print(f"\n  Total tests: {len(rows)}", flush=True)
    z_dists = np.array([r["z_dist"] for r in rows])

    # Negative z = co-targets closer than random groups (positive signal)
    sig_closer = (z_dists < -1.96).sum()
    print(f"  Mean z_dist: {z_dists.mean():.3f}", flush=True)
    print(f"  Sig closer (z_dist < -1.96): {sig_closer}/{len(z_dists)}", flush=True)
    print(f"  Sig farther (z_dist > +1.96): {(z_dists > 1.96).sum()}/{len(z_dists)}", flush=True)

    # Save
    csv_path = ITER_DIR / "h02b_trrust_cotarget_by_tf_layer.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Saved: {csv_path}", flush=True)

    # Summary by TF
    tf_zd = defaultdict(list)
    for r in rows:
        tf_zd[r["tf"]].append(r["z_dist"])

    summary = []
    for tf in sorted(tf_zd):
        zd = tf_zd[tf]
        summary.append({
            "tf": tf,
            "n_targets": valid_tfs[tf].__len__(),
            "mean_z_dist": float(np.mean(zd)),
            "n_sig_closer": int(sum(z < -1.96 for z in zd)),
            "n_layers": len(zd),
        })

    sum_csv = ITER_DIR / "h02b_trrust_tf_summary.csv"
    with open(sum_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        writer.writeheader()
        writer.writerows(summary)
    print(f"  Saved: {sum_csv}", flush=True)

    print("\nTop TFs by co-target clustering:", flush=True)
    for s in sorted(summary, key=lambda x: x["mean_z_dist"])[:10]:
        print(f"  {s['tf']:15s}: z_dist={s['mean_z_dist']:+.2f}, sig_closer={s['n_sig_closer']}/12")

    # Update hypothesis screen
    screen_path = ITER_DIR / "executor_hypothesis_screen.json"
    with open(screen_path) as f:
        screen = json.load(f)

    h02b_direction = ("positive" if z_dists.mean() < -1.5 else
                      ("negative" if z_dists.mean() > 1.5 else "inconclusive"))
    h02b_decision = ("promising" if h02b_direction == "positive" else
                     ("negative" if h02b_direction == "negative" else "inconclusive"))

    for h in screen["hypotheses"]:
        if h["id"] == "H02":
            h["name"] = "TRRUST co-target distance clustering (random-group null)"
            h["method"] = (f"Mean pairwise dist for co-targets of {len(valid_tfs)} TFs "
                           f"(>={MIN_GROUP_SIZE} mapped targets each) vs {N_NULL_REPS} "
                           f"random-same-size groups from gene pool; 12 layers, seed42.")
            h["status"] = "tested"
            h["primary_metric"] = "mean_z_dist_vs_random_groups"
            h["result_value"] = (f"mean_z_dist={z_dists.mean():.3f}, "
                                 f"sig_closer={sig_closer}/{len(z_dists)}")
            h["result_direction"] = h02b_direction
            h["artifact_paths"] = ["h02b_trrust_cotarget_by_tf_layer.csv", "h02b_trrust_tf_summary.csv"]
            h["decision"] = h02b_decision
            h["next_action"] = ("If positive: test across multiple seeds and compare to GO-term null; "
                                "if negative: try cosine similarity instead of Euclidean dist")
            break

    with open(screen_path, "w") as f:
        json.dump(screen, f, indent=2)
    print(f"  Updated: {screen_path}", flush=True)
    print("\n=== H02b DONE ===", flush=True)


if __name__ == "__main__":
    main()
