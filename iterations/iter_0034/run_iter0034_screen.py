"""
iter_0034 Multi-Hypothesis Screen

H01 (module_structure / new_method): Layer-wise B-cell community purity emergence
    For each layer L0-L11: recompute k=10 kNN greedy-modularity community partition.
    Identify PC1-negative community (contains most B-cell markers).
    Compute Fisher OR for B-cell marker enrichment.
    Spearman rho(layer, OR) tests progressive B-cell community crystallization.

H02 (intrinsic_dimensionality / new_method): PC2/PC3 cell-type axes
    At L11, compute PC2 and PC3 directions from 195 in-vocab genes.
    Test T-cell (n=15 in vocab) at PC2 positive pole.
    Test Myeloid (n=6 in vocab) at PC3 positive pole.
    Mann-Whitney AUROC vs null.

H03 (null_sensitivity / new_method): Expression-level confound test for B-cell PC1 signal
    Test whether B-cell PC1 signal (negative pole) is driven by lower mean expression.
    Null: shuffle gene labels, compute random 10-gene set AUROC distribution (n=1000 bootstrap).
    Control: test correlation of PC1 score with L2 norm (proxy for expression).
    If B-cell signal persists after controlling for L2 norm, structural geometry is confirmed.
"""

import numpy as np
import json
import csv
from pathlib import Path

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)

def jdump(obj, f, **kw):
    json.dump(obj, f, cls=NpEncoder, **kw)
from scipy.stats import spearmanr, mannwhitneyu, fisher_exact, pearsonr
from scipy.sparse import csr_matrix
from sklearn.neighbors import NearestNeighbors
import warnings
warnings.filterwarnings("ignore")

PROJECT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
               "/subproject_41_claude_topology_hypothesis_screening_autoloop")
ITER_DIR = PROJECT / "iterations" / "iter_0034"
ITER_DIR.mkdir(parents=True, exist_ok=True)

CYCLE1 = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
              "/subproject_38_geometric_residual_stream_interpretability"
              "/implementation/outputs/cycle1_main")

rng = np.random.default_rng(42)

# ─── Load embeddings ──────────────────────────────────────────────────────────
print("Loading embeddings...", flush=True)
emb = np.load(CYCLE1 / "layer_gene_embeddings.npy")   # [12, 4803, 512]
N_LAYERS, N_GENES_TOTAL, N_DIM = emb.shape
print(f"  Shape: {emb.shape}", flush=True)

with open(CYCLE1 / "gene_list.txt") as f:
    vocab_genes = [line.strip() for line in f if line.strip()]
gene_to_emb_idx = {g: i for i, g in enumerate(vocab_genes)}

EDGES_FILE = CYCLE1 / "cycle1_edge_dataset.tsv"
named_genes_set = set()
with open(EDGES_FILE) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        named_genes_set.add(row["source"])
        named_genes_set.add(row["target"])

# Filter to in-vocab, drop OOV
invocab_named = sorted([g for g in named_genes_set if g in gene_to_emb_idx])
L0_NORMS = np.linalg.norm(emb[0, [gene_to_emb_idx[g] for g in invocab_named], :], axis=1)
invocab_named = [g for g, n in zip(invocab_named, L0_NORMS) if n > 1e-8]
print(f"  In-vocab OOV-filtered: {len(invocab_named)}", flush=True)

gene_idxs_195 = np.array([gene_to_emb_idx[g] for g in invocab_named])
emb_195 = emb[:, gene_idxs_195, :]  # [12, 195, 512]
gene_idx_195 = {g: i for i, g in enumerate(invocab_named)}

# B-cell markers (9 confirmed in-vocab from iter_0033)
BCELL_MARKERS = ['CD19', 'MS4A1', 'CD79A', 'PRDM1', 'IRF4', 'CXCR4', 'SELL', 'BLK', 'SPIB']
bcell_in_195 = [g for g in BCELL_MARKERS if g in gene_idx_195]
print(f"  B-cell markers in-vocab: {len(bcell_in_195)}: {bcell_in_195}", flush=True)

TCELL_MARKERS = ['CD3G', 'CD8A', 'CD8B', 'PRF1', 'LCK', 'IFNG', 'IL2RA', 'FOXP3',
                 'EOMES', 'TBX21', 'GATA3', 'RUNX3', 'LEF1', 'CCR7']
tcell_in_195 = [g for g in TCELL_MARKERS if g in gene_idx_195]
print(f"  T-cell markers in-vocab: {len(tcell_in_195)}", flush=True)

MYELOID_MARKERS = ['ITGAX', 'TLR2', 'IL1B', 'IL6', 'CCL3', 'CXCL8', 'CD14', 'LYZ',
                   'S100A8', 'S100A9', 'CST3', 'FCGR3A']
myeloid_in_195 = [g for g in MYELOID_MARKERS if g in gene_idx_195]
print(f"  Myeloid markers in-vocab: {len(myeloid_in_195)}", flush=True)

# ─── Utilities ────────────────────────────────────────────────────────────────

def build_knn_adj(X, k=10):
    """Build symmetric k-NN adjacency (sparse)."""
    n = X.shape[0]
    nbrs = NearestNeighbors(n_neighbors=k+1, metric='euclidean').fit(X)
    distances, indices = nbrs.kneighbors(X)
    rows, cols = [], []
    for i in range(n):
        for j in indices[i, 1:]:
            rows.extend([i, j])
            cols.extend([j, i])
    data = np.ones(len(rows))
    A = csr_matrix((data, (rows, cols)), shape=(n, n))
    A.data[:] = 1.0
    return A

def greedy_modularity_communities(adj_dense):
    """Simple greedy modularity maximization via networkx."""
    import networkx as nx
    G = nx.from_numpy_array(adj_dense)
    communities = list(nx.algorithms.community.greedy_modularity_communities(G))
    modularity = nx.algorithms.community.modularity(G, communities)
    n = adj_dense.shape[0]
    labels = np.zeros(n, dtype=int)
    for ci, comm in enumerate(communities):
        for node in comm:
            labels[node] = ci
    return labels, modularity, communities

def fisher_enrichment_in_community(bcell_idxs, community_labels, target_community):
    """Fisher test: B-cell markers in target community vs background."""
    n_total = len(community_labels)
    in_comm = (community_labels == target_community)
    n_comm = in_comm.sum()
    # B-cell in community, B-cell NOT in community
    b_in = sum(1 for i in bcell_idxs if community_labels[i] == target_community)
    b_out = len(bcell_idxs) - b_in
    nb_in = n_comm - b_in
    nb_out = (n_total - n_comm) - b_out
    table = [[b_in, b_out], [nb_in, nb_out]]
    odds_ratio, p = fisher_exact(table, alternative='greater')
    return {"b_in": b_in, "b_out": b_out, "nb_in": nb_in, "nb_out": nb_out,
            "odds_ratio": float(odds_ratio), "fisher_p": float(p)}

# ═══════════════════════════════════════════════════════════════════════════════
# H01: Layer-wise B-cell community purity emergence
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== H01: Layer-wise B-cell community purity emergence ===", flush=True)

bcell_idxs_195 = np.array([gene_idx_195[g] for g in bcell_in_195])

layer_community_results = []
for layer in range(12):
    X_lay = emb_195[layer]
    A_lay = build_knn_adj(X_lay, k=10)
    labels_lay, mod_lay, comms_lay = greedy_modularity_communities(A_lay.toarray())

    # Find PC1-negative community
    X_c = X_lay - X_lay.mean(axis=0)
    _, _, Vt = np.linalg.svd(X_c, full_matrices=False)
    pc1_scores = X_c @ Vt[0]

    # Identify which community has negative mean PC1
    comm_pc1_means = {}
    for ci in range(len(comms_lay)):
        mask = (labels_lay == ci)
        comm_pc1_means[ci] = pc1_scores[mask].mean()

    # Community with most negative PC1 mean
    neg_comm = min(comm_pc1_means, key=comm_pc1_means.get)

    # Fisher test for B-cell in PC1-negative community
    fisher_res = fisher_enrichment_in_community(bcell_idxs_195, labels_lay, neg_comm)

    # Also test B-cell in any community (best-match)
    best_b_in = 0
    best_comm_id = 0
    for ci in range(len(comms_lay)):
        b_in = sum(1 for i in bcell_idxs_195 if labels_lay[i] == ci)
        if b_in > best_b_in:
            best_b_in = b_in
            best_comm_id = ci

    fisher_best = fisher_enrichment_in_community(bcell_idxs_195, labels_lay, best_comm_id)

    print(f"  L{layer}: n_comms={len(comms_lay)}, mod={mod_lay:.3f}, "
          f"neg_comm={neg_comm} (PC1_mean={comm_pc1_means[neg_comm]:.3f}), "
          f"OR_neg={fisher_res['odds_ratio']:.2f} p={fisher_res['fisher_p']:.4f}, "
          f"best_comm OR={fisher_best['odds_ratio']:.2f} p={fisher_best['fisher_p']:.4f}", flush=True)

    layer_community_results.append({
        "layer": layer,
        "n_communities": len(comms_lay),
        "modularity": float(mod_lay),
        "pc1_neg_community": int(neg_comm),
        "pc1_neg_community_mean": float(comm_pc1_means[neg_comm]),
        "fisher_neg_comm": fisher_res,
        "best_bcell_community": int(best_comm_id),
        "fisher_best_comm": fisher_best,
    })

# Spearman rho: layer vs OR (PC1-neg community)
layers = [r["layer"] for r in layer_community_results]
ors_neg = [r["fisher_neg_comm"]["odds_ratio"] for r in layer_community_results]
ors_best = [r["fisher_best_comm"]["odds_ratio"] for r in layer_community_results]

# Cap infinity (when b_out=0)
ors_neg_capped = [min(x, 100.0) for x in ors_neg]
ors_best_capped = [min(x, 100.0) for x in ors_best]

rho_neg, p_neg = spearmanr(layers, ors_neg_capped)
rho_best, p_best = spearmanr(layers, ors_best_capped)

print(f"\n  Spearman rho(layer, OR_neg_comm) = {rho_neg:.3f}, p={p_neg:.4f}", flush=True)
print(f"  Spearman rho(layer, OR_best_comm) = {rho_best:.3f}, p={p_best:.4f}", flush=True)

h01_result = {
    "hypothesis": "layer_wise_bcell_community_purity_emergence",
    "bcell_markers": bcell_in_195,
    "per_layer_results": layer_community_results,
    "spearman_rho_layer_vs_OR_neg_comm": float(rho_neg),
    "spearman_p_layer_vs_OR_neg_comm": float(p_neg),
    "spearman_rho_layer_vs_OR_best_comm": float(rho_best),
    "spearman_p_layer_vs_OR_best_comm": float(p_best),
    "L11_OR_neg_comm": float(ors_neg[-1]),
    "L11_OR_best_comm": float(ors_best[-1]),
    "summary": "Spearman rho tests progressive B-cell community crystallization across layers"
}

with open(ITER_DIR / "h01_layer_bcell_community_purity.json", "w") as f:
    jdump(h01_result, f, indent=2)
print("  Saved h01_layer_bcell_community_purity.json", flush=True)

# ═══════════════════════════════════════════════════════════════════════════════
# H02: PC2/PC3 cell-type axes at L11
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== H02: PC2/PC3 cell-type axes at L11 ===", flush=True)

X_L11 = emb_195[11]
X_L11_c = X_L11 - X_L11.mean(axis=0)
U, S, Vt = np.linalg.svd(X_L11_c, full_matrices=False)

var_exp = (S**2) / (S**2).sum()
print(f"  PC1 var: {var_exp[0]:.3f}, PC2: {var_exp[1]:.3f}, PC3: {var_exp[2]:.3f}", flush=True)

pc_scores = X_L11_c @ Vt[:5].T  # [195, 5] — first 5 PCs

def test_pc_enrichment(marker_list, pc_idx, direction, all_genes, gene_idx_map, pc_scores_all):
    """Mann-Whitney AUROC for marker set vs background on a given PC."""
    in_vocab = [g for g in marker_list if g in gene_idx_map]
    if len(in_vocab) < 2:
        return {"n": len(in_vocab), "auroc": None, "p": None, "note": "too_few"}
    set_idxs = np.array([gene_idx_map[g] for g in in_vocab])
    set_scores = pc_scores_all[set_idxs, pc_idx]
    mask = np.ones(len(all_genes), dtype=bool)
    mask[set_idxs] = False
    bg_scores = pc_scores_all[mask, pc_idx]
    alt = 'greater' if direction == 'positive' else 'less'
    stat, p = mannwhitneyu(set_scores, bg_scores, alternative=alt)
    auroc = stat / (len(set_scores) * len(bg_scores))
    return {
        "genes": in_vocab,
        "n": len(in_vocab),
        "mean_score": float(set_scores.mean()),
        "bg_mean": float(bg_scores.mean()),
        "auroc": float(auroc),
        "p": float(p),
        "direction": direction
    }

pc_axes_results = {}

# PC1: B-cell at negative pole (verify)
pc1_bcell = test_pc_enrichment(BCELL_MARKERS, 0, 'negative', invocab_named, gene_idx_195, pc_scores)
pc1_tcell = test_pc_enrichment(TCELL_MARKERS, 0, 'positive', invocab_named, gene_idx_195, pc_scores)
pc1_myeloid = test_pc_enrichment(MYELOID_MARKERS, 0, 'positive', invocab_named, gene_idx_195, pc_scores)
print(f"  PC1: B-cell neg pole AUROC={pc1_bcell['auroc']:.3f} p={pc1_bcell['p']:.4f}", flush=True)
print(f"  PC1: T-cell pos pole AUROC={pc1_tcell['auroc']:.3f} p={pc1_tcell['p']:.4f}", flush=True)
print(f"  PC1: Myeloid pos pole AUROC={pc1_myeloid['auroc']:.3f} p={pc1_myeloid['p']:.4f}", flush=True)

# PC2: test T-cell at positive pole, B-cell at negative pole
pc2_bcell_neg = test_pc_enrichment(BCELL_MARKERS, 1, 'negative', invocab_named, gene_idx_195, pc_scores)
pc2_bcell_pos = test_pc_enrichment(BCELL_MARKERS, 1, 'positive', invocab_named, gene_idx_195, pc_scores)
pc2_tcell_pos = test_pc_enrichment(TCELL_MARKERS, 1, 'positive', invocab_named, gene_idx_195, pc_scores)
pc2_tcell_neg = test_pc_enrichment(TCELL_MARKERS, 1, 'negative', invocab_named, gene_idx_195, pc_scores)
pc2_myeloid_pos = test_pc_enrichment(MYELOID_MARKERS, 1, 'positive', invocab_named, gene_idx_195, pc_scores)
pc2_myeloid_neg = test_pc_enrichment(MYELOID_MARKERS, 1, 'negative', invocab_named, gene_idx_195, pc_scores)
print(f"  PC2: B-cell neg AUROC={pc2_bcell_neg['auroc']:.3f} p={pc2_bcell_neg['p']:.4f}", flush=True)
print(f"  PC2: B-cell pos AUROC={pc2_bcell_pos['auroc']:.3f} p={pc2_bcell_pos['p']:.4f}", flush=True)
print(f"  PC2: T-cell pos AUROC={pc2_tcell_pos['auroc']:.3f} p={pc2_tcell_pos['p']:.4f}", flush=True)
print(f"  PC2: T-cell neg AUROC={pc2_tcell_neg['auroc']:.3f} p={pc2_tcell_neg['p']:.4f}", flush=True)
print(f"  PC2: Myeloid pos AUROC={pc2_myeloid_pos['auroc']:.3f} p={pc2_myeloid_pos['p']:.4f}", flush=True)
print(f"  PC2: Myeloid neg AUROC={pc2_myeloid_neg['auroc']:.3f} p={pc2_myeloid_neg['p']:.4f}", flush=True)

# PC3: test Myeloid at positive/negative
pc3_bcell_neg = test_pc_enrichment(BCELL_MARKERS, 2, 'negative', invocab_named, gene_idx_195, pc_scores)
pc3_tcell_pos = test_pc_enrichment(TCELL_MARKERS, 2, 'positive', invocab_named, gene_idx_195, pc_scores)
pc3_myeloid_pos = test_pc_enrichment(MYELOID_MARKERS, 2, 'positive', invocab_named, gene_idx_195, pc_scores)
pc3_myeloid_neg = test_pc_enrichment(MYELOID_MARKERS, 2, 'negative', invocab_named, gene_idx_195, pc_scores)
print(f"  PC3: B-cell neg AUROC={pc3_bcell_neg['auroc']:.3f} p={pc3_bcell_neg['p']:.4f}", flush=True)
print(f"  PC3: T-cell pos AUROC={pc3_tcell_pos['auroc']:.3f} p={pc3_tcell_pos['p']:.4f}", flush=True)
print(f"  PC3: Myeloid pos AUROC={pc3_myeloid_pos['auroc']:.3f} p={pc3_myeloid_pos['p']:.4f}", flush=True)
print(f"  PC3: Myeloid neg AUROC={pc3_myeloid_neg['auroc']:.3f} p={pc3_myeloid_neg['p']:.4f}", flush=True)

h02_result = {
    "hypothesis": "pc2_pc3_cell_type_axes",
    "layer": 11,
    "variance_explained": [float(v) for v in var_exp[:5]],
    "PC1": {
        "bcell_negative": pc1_bcell,
        "tcell_positive": pc1_tcell,
        "myeloid_positive": pc1_myeloid,
    },
    "PC2": {
        "bcell_negative": pc2_bcell_neg,
        "bcell_positive": pc2_bcell_pos,
        "tcell_positive": pc2_tcell_pos,
        "tcell_negative": pc2_tcell_neg,
        "myeloid_positive": pc2_myeloid_pos,
        "myeloid_negative": pc2_myeloid_neg,
    },
    "PC3": {
        "bcell_negative": pc3_bcell_neg,
        "tcell_positive": pc3_tcell_pos,
        "myeloid_positive": pc3_myeloid_pos,
        "myeloid_negative": pc3_myeloid_neg,
    }
}

with open(ITER_DIR / "h02_pc2_pc3_axes.json", "w") as f:
    jdump(h02_result, f, indent=2)
print("  Saved h02_pc2_pc3_axes.json", flush=True)

# ═══════════════════════════════════════════════════════════════════════════════
# H03: Expression-level confound test for B-cell PC1 signal
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== H03: Expression-level confound test for B-cell PC1 signal ===", flush=True)

# Get PC1 scores at L11
pc1_L11 = pc_scores[:, 0]

# L2 norms as proxy for expression level
l2_norms_L11 = np.linalg.norm(X_L11, axis=1)  # [195]

# Correlation PC1 score vs L2 norm
rho_pc1_l2, p_pc1_l2 = pearsonr(pc1_L11, l2_norms_L11)
print(f"  Pearson r(PC1, L2_norm) at L11 = {rho_pc1_l2:.3f}, p={p_pc1_l2:.4f}", flush=True)

# B-cell markers: PC1 scores and L2 norms
bcell_pc1_scores = pc1_L11[bcell_idxs_195]
bcell_l2_norms = l2_norms_L11[bcell_idxs_195]
all_l2_norms = l2_norms_L11

# Test: are B-cell markers at lower L2 norm (expression proxy) than background?
# If yes, PC1 B-cell signal might be confounded by expression
bcell_l2_vs_bg = mannwhitneyu(bcell_l2_norms, all_l2_norms, alternative='less')
auroc_l2 = bcell_l2_vs_bg[0] / (len(bcell_l2_norms) * len(all_l2_norms))
print(f"  B-cell L2 norms: mean={bcell_l2_norms.mean():.3f}, bg={all_l2_norms.mean():.3f}, "
      f"AUROC={auroc_l2:.3f}, p={bcell_l2_vs_bg[1]:.4f}", flush=True)

# Residual analysis: regress PC1 on L2 norm, test if B-cell markers are still negative
# PC1_resid = PC1 - (a * L2_norm + b)
from numpy.polynomial import polynomial as P
coeffs = np.polyfit(l2_norms_L11, pc1_L11, 1)
pc1_resid = pc1_L11 - (coeffs[0] * l2_norms_L11 + coeffs[1])

bcell_resid = pc1_resid[bcell_idxs_195]
bg_mask = np.ones(195, dtype=bool)
bg_mask[bcell_idxs_195] = False
bg_resid = pc1_resid[bg_mask]

stat_resid, p_resid = mannwhitneyu(bcell_resid, bg_resid, alternative='less')
auroc_resid = stat_resid / (len(bcell_resid) * len(bg_resid))
print(f"  B-cell PC1 residual (after L2 regression): AUROC={auroc_resid:.3f}, p={p_resid:.4f}", flush=True)
print(f"  B-cell mean residual={bcell_resid.mean():.3f}, bg={bg_resid.mean():.3f}", flush=True)

# Bootstrap null: 1000 random 9-gene sets, compute their PC1 AUROC (negative pole)
n_bcell = len(bcell_idxs_195)
n_bootstrap = 1000
null_aurocs = []
for _ in range(n_bootstrap):
    null_idxs = rng.choice(195, size=n_bcell, replace=False)
    null_scores = pc1_L11[null_idxs]
    null_bg = pc1_L11[np.setdiff1d(np.arange(195), null_idxs)]
    stat_null, _ = mannwhitneyu(null_scores, null_bg, alternative='less')
    null_aurocs.append(stat_null / (n_bcell * len(null_bg)))

null_aurocs = np.array(null_aurocs)
observed_auroc_bcell = float(mannwhitneyu(bcell_pc1_scores, pc1_L11[bg_mask], alternative='less')[0]) / (n_bcell * len(pc1_L11[bg_mask]))
empirical_p = (null_aurocs <= observed_auroc_bcell).mean()
z_score = (observed_auroc_bcell - null_aurocs.mean()) / null_aurocs.std()

print(f"\n  Bootstrap null: mean null AUROC={null_aurocs.mean():.4f}, std={null_aurocs.std():.4f}", flush=True)
print(f"  Observed B-cell PC1-neg AUROC = {observed_auroc_bcell:.4f}", flush=True)
print(f"  Empirical p (AUROC <= observed) = {empirical_p:.4f}", flush=True)
print(f"  Z-score vs null = {z_score:.3f}", flush=True)

h03_result = {
    "hypothesis": "expression_confound_test_bcell_pc1",
    "layer": 11,
    "pearson_r_pc1_vs_l2norm": float(rho_pc1_l2),
    "pearson_p": float(p_pc1_l2),
    "bcell_l2_norm_mean": float(bcell_l2_norms.mean()),
    "bg_l2_norm_mean": float(all_l2_norms.mean()),
    "bcell_l2_auroc_less": float(auroc_l2),
    "bcell_l2_p": float(bcell_l2_vs_bg[1]),
    "bcell_pc1_resid_auroc": float(auroc_resid),
    "bcell_pc1_resid_p": float(p_resid),
    "bcell_resid_mean": float(bcell_resid.mean()),
    "bg_resid_mean": float(bg_resid.mean()),
    "bootstrap_null_mean_auroc": float(null_aurocs.mean()),
    "bootstrap_null_std_auroc": float(null_aurocs.std()),
    "observed_bcell_auroc": float(observed_auroc_bcell),
    "bootstrap_empirical_p": float(empirical_p),
    "bootstrap_z_score": float(z_score),
    "n_bootstrap": n_bootstrap,
    "interpretation": (
        "If bcell_l2_auroc_less is near 0.5 (B-cell L2 norms NOT lower), "
        "AND bcell_pc1_resid_auroc is significant, then B-cell PC1 signal is structural (not expression-driven)."
    )
}

with open(ITER_DIR / "h03_expression_confound_test.json", "w") as f:
    jdump(h03_result, f, indent=2)
print("  Saved h03_expression_confound_test.json", flush=True)

print("\n=== All iter_0034 experiments complete ===", flush=True)
print(f"Output dir: {ITER_DIR}", flush=True)
