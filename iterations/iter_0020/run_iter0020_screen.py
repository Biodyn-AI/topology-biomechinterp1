"""
iter_0020 Multi-Hypothesis Screen

H01 (null_sensitivity, refinement): Shuffle null for multi-axis composite + magnitude quintile
    Permute gene labels 1000x to build null distribution for the 2.18x enrichment at 3-axis co-polarity.
    Within 3-axis stratum, split by composite magnitude quintile, test monotonic enrichment.
    This validates whether the 2.18x claim survives permutation testing.

H02 (graph_topology, new_method): Attention-SVD joint ROC predictor for TF-target pairs
    Build a joint predictor: attention score (TF-target) + max multi-axis co-polarity count.
    Compute AUROC for TRRUST (activation+repression) vs null using:
      (a) attention only
      (b) SVD multi-axis composite only
      (c) joint (attention + SVD composite)
    Test whether joint model improves over single features.

H03 (persistent_homology, new_family): Persistent homology H0 on gene embedding cloud
    For 209 named genes at each of 12 scGPT layers: compute H0 Betti curves from Vietoris-Rips.
    Compare persistence lifetime distributions between STRING-connected and non-connected pairs.
    Also compute connected component count at distance threshold.
    Compare to random-permuted-embedding null (same distances, shuffled labels).
"""

import numpy as np
import json
import csv
from pathlib import Path
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.metrics import roc_auc_score
import warnings
warnings.filterwarnings("ignore")

# ─── Paths ────────────────────────────────────────────────────────────────────
PROJECT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
               "/subproject_41_claude_topology_hypothesis_screening_autoloop")
ITER_DIR = PROJECT / "iterations" / "iter_0020"
ITER15_DIR = PROJECT / "iterations" / "iter_0015"
ITER19_DIR = PROJECT / "iterations" / "iter_0019"
CYCLE1 = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
              "/subproject_38_geometric_residual_stream_interpretability"
              "/implementation/outputs/cycle1_main")

EMB_PATH = CYCLE1 / "layer_gene_embeddings.npy"
EDGE_PATH = CYCLE1 / "cycle1_edge_dataset.tsv"
STRING_API_CACHE = ITER15_DIR / "string_ppi_score04_cache.json"
TRRUST_PATH = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
                   "/single_cell_mechinterp/external/networks/trrust_human.tsv")
ATT_PATH = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
                "/single_cell_mechinterp/outputs/invariant_causal_edges/lung"
                "/attention_scores.npy")

ITER_DIR.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(42)

print("Loading scGPT embeddings ...", flush=True)
emb = np.load(EMB_PATH)  # [12, 4803, 512]
N_LAYERS, N_GENES_TOTAL, N_DIM = emb.shape
print(f"  emb shape: {emb.shape}", flush=True)

print("Loading gene list ...", flush=True)
gene_list_path = CYCLE1 / "gene_list.txt"
with open(gene_list_path) as f:
    vocab_genes = [line.strip() for line in f]
gene_to_emb_idx = {g: i for i, g in enumerate(vocab_genes)}
print(f"  Vocab size: {len(vocab_genes)}", flush=True)

print("Loading named genes from edge dataset ...", flush=True)
named_gene_set = set()
with open(EDGE_PATH) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        named_gene_set.add(row['source'])
        named_gene_set.add(row['target'])

named_genes = sorted(g for g in named_gene_set if g in gene_to_emb_idx)
named_idx = np.array([gene_to_emb_idx[g] for g in named_genes])
N_NAMED = len(named_genes)
gene_to_named_idx = {g: i for i, g in enumerate(named_genes)}
print(f"  Named genes with embeddings: {N_NAMED}", flush=True)

print("Loading STRING cache ...", flush=True)
with open(STRING_API_CACHE) as f:
    string_cache = json.load(f)
string_set = set()
# Handle both formats: list-of-dicts under "pairs" key, or flat dict
if isinstance(string_cache, dict) and "pairs" in string_cache:
    for item in string_cache["pairs"]:
        if item["score"] >= 0.4:
            a, b = item["g1"], item["g2"]
            string_set.add((min(a, b), max(a, b)))
elif isinstance(string_cache, dict):
    for k, v in string_cache.items():
        if isinstance(v, (int, float)) and v >= 0.4:
            a, b = k.split("__")
            string_set.add((min(a, b), max(a, b)))
print(f"  STRING edges (score>=0.4): {len(string_set)}", flush=True)

# Build pair labels for all named-gene pairs
print("Building pair label arrays ...", flush=True)
n = N_NAMED
all_pairs = []
string_labels = []
for i in range(n):
    for j in range(i+1, n):
        gi, gj = named_genes[i], named_genes[j]
        key = (min(gi, gj), max(gi, gj))
        all_pairs.append((i, j))
        string_labels.append(1 if key in string_set else 0)
all_pairs = np.array(all_pairs)
string_labels = np.array(string_labels)
print(f"  Total pairs: {len(all_pairs)}, STRING positives: {string_labels.sum()}", flush=True)
string_baseline = string_labels.mean()
print(f"  STRING baseline rate: {string_baseline:.4f}", flush=True)

# ─── Helper: compute multi-axis co-polarity count for all pairs ───────────────
def compute_copole_features(emb_named, axes=(1, 2, 3)):
    """
    emb_named: [N_named, 512]
    For each axis in axes, SVD -> check co-polarity (same sign of projection).
    Returns:
      copole_count: [N_pairs] int (0-len(axes))
      composite_mag: [N_pairs] float (sum of |proj_i| * |proj_j| across co-polar axes)
    """
    centered = emb_named - emb_named.mean(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(centered, full_matrices=False)
    projs = []
    for ax in axes:
        p = U[:, ax]  # [N_named]
        projs.append(p)

    n = len(emb_named)
    i_idx = all_pairs[:, 0]
    j_idx = all_pairs[:, 1]

    copole_count = np.zeros(len(all_pairs), dtype=int)
    composite_mag = np.zeros(len(all_pairs), dtype=float)
    for p in projs:
        pi, pj = p[i_idx], p[j_idx]
        same_sign = (pi * pj) > 0
        copole_count += same_sign.astype(int)
        composite_mag += np.where(same_sign, np.abs(pi) * np.abs(pj), 0.0)
    return copole_count, composite_mag

# ─── H01: Shuffle null for multi-axis composite ───────────────────────────────
print("\n=== H01: Shuffle null for multi-axis composite ===", flush=True)

# Compute observed enrichment per count stratum across all layers
observed_enrichments = {}  # count -> list over layers
for layer in range(N_LAYERS):
    emb_named = emb[layer][named_idx]
    copole_count, _ = compute_copole_features(emb_named)
    for cnt in range(4):
        mask = copole_count == cnt
        if mask.sum() > 0:
            enrich = string_labels[mask].mean() / string_baseline
            observed_enrichments.setdefault(cnt, []).append(enrich)

obs_enrich_mean = {cnt: np.mean(vals) for cnt, vals in observed_enrichments.items()}
print(f"  Observed enrichments: {obs_enrich_mean}", flush=True)

# Permutation null: shuffle gene labels 1000 times at a representative layer (layer 8)
N_PERM = 1000
perm_enrich_3 = []
layer_rep = 8
emb_named_rep = emb[layer_rep][named_idx].copy()

print(f"  Running {N_PERM} permutations at layer {layer_rep} ...", flush=True)
for perm_i in range(N_PERM):
    perm_idx = rng.permutation(N_NAMED)
    emb_perm = emb_named_rep[perm_idx]
    copole_count_perm, _ = compute_copole_features(emb_perm)
    mask = copole_count_perm == 3
    if mask.sum() > 0:
        enrich_perm = string_labels[mask].mean() / string_baseline
    else:
        enrich_perm = 1.0
    perm_enrich_3.append(enrich_perm)

perm_enrich_3 = np.array(perm_enrich_3)
# Observed at layer 8
emb_named_rep2 = emb[layer_rep][named_idx]
copole_count_obs, composite_mag_obs = compute_copole_features(emb_named_rep2)
mask3 = copole_count_obs == 3
obs_at_layer8 = string_labels[mask3].mean() / string_baseline if mask3.sum() > 0 else 1.0
perm_p = (perm_enrich_3 >= obs_at_layer8).mean()
perm_mean = perm_enrich_3.mean()
perm_std = perm_enrich_3.std()
z_score = (obs_at_layer8 - perm_mean) / (perm_std + 1e-12)

print(f"  Layer 8 observed enrichment at count=3: {obs_at_layer8:.3f}", flush=True)
print(f"  Perm null: mean={perm_mean:.3f}, std={perm_std:.3f}", flush=True)
print(f"  z={z_score:.2f}, p={perm_p:.4f}", flush=True)

# Within 3-axis stratum: magnitude quintile stratification
print("  Quintile stratification within 3-axis stratum ...", flush=True)
mask3_obs = copole_count_obs == 3
mags3 = composite_mag_obs[mask3_obs]
labels3 = string_labels[mask3_obs]
quintile_enrichments = []
if mask3_obs.sum() >= 50:
    qtiles = np.percentile(mags3, [0, 20, 40, 60, 80, 100])
    for qi in range(5):
        q_mask = (mags3 >= qtiles[qi]) & (mags3 < qtiles[qi+1])
        if qi == 4:
            q_mask = mags3 >= qtiles[qi]
        if q_mask.sum() > 0:
            enrich_q = labels3[q_mask].mean() / string_baseline
            quintile_enrichments.append({"quintile": qi+1, "n": int(q_mask.sum()),
                                          "enrichment": float(enrich_q)})
            print(f"    Q{qi+1}: n={q_mask.sum()}, enrichment={enrich_q:.3f}x", flush=True)

# Spearman between quintile rank and enrichment
if len(quintile_enrichments) >= 4:
    ranks = [x["quintile"] for x in quintile_enrichments]
    enr_vals = [x["enrichment"] for x in quintile_enrichments]
    sp_r, sp_p = spearmanr(ranks, enr_vals)
    print(f"  Quintile Spearman r={sp_r:.3f}, p={sp_p:.4f}", flush=True)
else:
    sp_r, sp_p = 0.0, 1.0

h01_result = {
    "obs_enrichment_count3_layer8": float(obs_at_layer8),
    "obs_enrichment_by_count_mean_across_layers": {str(k): float(v) for k, v in obs_enrich_mean.items()},
    "perm_null_mean": float(perm_mean),
    "perm_null_std": float(perm_std),
    "perm_z_score": float(z_score),
    "perm_p_value": float(perm_p),
    "n_permutations": N_PERM,
    "quintile_enrichments": quintile_enrichments,
    "quintile_spearman_r": float(sp_r),
    "quintile_spearman_p": float(sp_p),
    "n_pairs_count3": int(mask3.sum())
}
with open(ITER_DIR / "h01_shuffle_null_composite.json", "w") as f:
    json.dump(h01_result, f, indent=2)
print("  H01 saved.", flush=True)

# ─── H02: Attention-SVD Joint ROC ────────────────────────────────────────────
print("\n=== H02: Attention-SVD joint ROC predictor ===", flush=True)

print("  Loading attention scores ...", flush=True)
att = np.load(ATT_PATH)  # [8181, 8181]
att_size = att.shape[0]

# Load gene names from h5ad
import h5py
H5AD_PATH = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
                 "/single_cell_mechinterp/outputs/invariant_causal_edges/lung"
                 "/processed.h5ad")

print("  Loading h5ad gene names ...", flush=True)
with h5py.File(H5AD_PATH, "r") as f:
    var_names = f["var"]["_index"][()]
    if isinstance(var_names[0], bytes):
        var_names = [v.decode() for v in var_names]
    else:
        var_names = list(var_names)
print(f"  h5ad genes: {len(var_names)}", flush=True)

att_gene_to_idx = {g: i for i, g in enumerate(var_names)}

# Build TRRUST TF-target pairs
print("  Loading TRRUST ...", flush=True)
trrust_pairs = set()
with open(TRRUST_PATH) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            tf, target = parts[0], parts[1]
            trrust_pairs.add((tf, target))
            trrust_pairs.add((target, tf))  # symmetric

print(f"  TRRUST total TF-target pairs (symmetric): {len(trrust_pairs)}", flush=True)

# For ROC: build feature vectors for all pairs of named genes present in attention matrix
att_named = []
att_named_pairs = []
for g in named_genes:
    if g in att_gene_to_idx:
        att_named.append(att_gene_to_idx[g])
    else:
        att_named.append(None)

# Compute multi-axis co-polarity features at best layer (layer 8)
emb_l8 = emb[8][named_idx]
copole_cnt_l8, comp_mag_l8 = compute_copole_features(emb_l8)

# Build feature matrix
print("  Building ROC feature matrix ...", flush=True)
n_pairs = len(all_pairs)
feat_att = np.zeros(n_pairs, dtype=float)
feat_svd = copole_cnt_l8.astype(float)
feat_mag = comp_mag_l8

for k, (i, j) in enumerate(all_pairs):
    ai = att_named[i]
    aj = att_named[j]
    if ai is not None and aj is not None and ai < att_size and aj < att_size:
        # Symmetric attention
        feat_att[k] = 0.5 * (att[ai, aj] + att[aj, ai])

# Normalize features
feat_att_norm = (feat_att - feat_att.mean()) / (feat_att.std() + 1e-12)
feat_svd_norm = (feat_svd - feat_svd.mean()) / (feat_svd.std() + 1e-12)
feat_mag_norm = (feat_mag - feat_mag.mean()) / (feat_mag.std() + 1e-12)

# TRRUST labels: is this pair in TRRUST?
trrust_labels = np.zeros(n_pairs, dtype=int)
for k, (i, j) in enumerate(all_pairs):
    gi, gj = named_genes[i], named_genes[j]
    if (gi, gj) in trrust_pairs or (gj, gi) in trrust_pairs:
        trrust_labels[k] = 1

print(f"  TRRUST positives: {trrust_labels.sum()}", flush=True)

# Only pairs with non-zero attention (att_named present)
att_valid = feat_att > 0
n_valid = att_valid.sum()
print(f"  Pairs with valid attention scores: {n_valid}", flush=True)

def safe_auroc(y_true, y_score):
    if y_true.sum() == 0 or y_true.sum() == len(y_true):
        return 0.5
    return roc_auc_score(y_true, y_score)

# ROC for STRING prediction
auc_att_str = safe_auroc(string_labels[att_valid], feat_att_norm[att_valid])
auc_svd_str = safe_auroc(string_labels[att_valid], feat_svd_norm[att_valid])
auc_mag_str = safe_auroc(string_labels[att_valid], feat_mag_norm[att_valid])
auc_joint_str = safe_auroc(string_labels[att_valid], feat_att_norm[att_valid] + feat_svd_norm[att_valid])
auc_joint_mag_str = safe_auroc(string_labels[att_valid], feat_att_norm[att_valid] + feat_mag_norm[att_valid])

print(f"  STRING AUROC - att: {auc_att_str:.3f}, svd_cnt: {auc_svd_str:.3f}, mag: {auc_mag_str:.3f}", flush=True)
print(f"  STRING AUROC - joint(att+svd): {auc_joint_str:.3f}, joint(att+mag): {auc_joint_mag_str:.3f}", flush=True)

# ROC for TRRUST prediction
if trrust_labels[att_valid].sum() >= 10:
    auc_att_tr = safe_auroc(trrust_labels[att_valid], feat_att_norm[att_valid])
    auc_svd_tr = safe_auroc(trrust_labels[att_valid], feat_svd_norm[att_valid])
    auc_joint_tr = safe_auroc(trrust_labels[att_valid], feat_att_norm[att_valid] + feat_svd_norm[att_valid])
    print(f"  TRRUST AUROC - att: {auc_att_tr:.3f}, svd: {auc_svd_tr:.3f}, joint: {auc_joint_tr:.3f}", flush=True)
else:
    auc_att_tr = auc_svd_tr = auc_joint_tr = None
    print(f"  TRRUST: insufficient positives ({trrust_labels[att_valid].sum()}) in valid pairs", flush=True)

h02_result = {
    "n_pairs_total": n_pairs,
    "n_pairs_with_attention": int(n_valid),
    "n_string_positives": int(string_labels.sum()),
    "n_trrust_positives": int(trrust_labels.sum()),
    "n_trrust_positives_att_valid": int(trrust_labels[att_valid].sum()),
    "string_auroc": {
        "attention_only": float(auc_att_str),
        "svd_count_only": float(auc_svd_str),
        "magnitude_only": float(auc_mag_str),
        "joint_att_svd": float(auc_joint_str),
        "joint_att_mag": float(auc_joint_mag_str)
    },
    "trrust_auroc": {
        "attention_only": float(auc_att_tr) if auc_att_tr else None,
        "svd_count_only": float(auc_svd_tr) if auc_svd_tr else None,
        "joint_att_svd": float(auc_joint_tr) if auc_joint_tr else None
    }
}
with open(ITER_DIR / "h02_joint_roc.json", "w") as f:
    json.dump(h02_result, f, indent=2)
print("  H02 saved.", flush=True)

# ─── H03: Persistent Homology H0 on gene embedding cloud ─────────────────────
print("\n=== H03: Persistent Homology H0 on gene embedding cloud ===", flush=True)

# Try to import ripser; if unavailable use fallback (connected component counting)
try:
    from ripser import ripser
    HAS_RIPSER = True
    print("  ripser available.", flush=True)
except ImportError:
    HAS_RIPSER = False
    print("  ripser NOT available — using fallback connected component analysis.", flush=True)

def pairwise_distances(X):
    """Compute pairwise Euclidean distance matrix."""
    norms = np.sum(X**2, axis=1)
    D = norms[:, None] + norms[None, :] - 2 * X @ X.T
    D = np.sqrt(np.maximum(D, 0))
    return D

def connected_components_at_threshold(D, eps):
    """Count connected components in epsilon-neighborhood graph."""
    n = len(D)
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry
    for i in range(n):
        for j in range(i+1, n):
            if D[i, j] <= eps:
                union(i, j)
    return len(set(find(i) for i in range(n)))

h03_per_layer = []

for layer in range(N_LAYERS):
    emb_named_l = emb[layer][named_idx]
    # Normalize to unit sphere for comparability
    emb_norm = emb_named_l / (np.linalg.norm(emb_named_l, axis=1, keepdims=True) + 1e-12)
    D = pairwise_distances(emb_norm)

    if HAS_RIPSER:
        # Run H0+H1 persistent homology
        result = ripser(D, maxdim=1, distance_matrix=True)
        dgms = result['dgms']
        # H0 lifetime stats
        h0 = dgms[0]
        h0_finite = h0[h0[:, 1] < np.inf]
        h0_lifetimes = h0_finite[:, 1] - h0_finite[:, 0]
        h0_n_components_at_birth = len(h0)
        h0_max_lifetime = float(h0_lifetimes.max()) if len(h0_lifetimes) > 0 else 0.0
        h0_mean_lifetime = float(h0_lifetimes.mean()) if len(h0_lifetimes) > 0 else 0.0
        # H1 (loops)
        h1 = dgms[1] if len(dgms) > 1 else np.empty((0, 2))
        h1_finite = h1[h1[:, 1] < np.inf] if len(h1) > 0 else np.empty((0, 2))
        h1_n = len(h1_finite)
        h1_mean_lifetime = float((h1_finite[:, 1] - h1_finite[:, 0]).mean()) if h1_n > 0 else 0.0
        method_used = "ripser"
    else:
        # Fallback: threshold sweep for connected components
        d_vals = D[np.triu_indices_from(D, k=1)]
        thresholds = np.percentile(d_vals, [10, 20, 30, 40, 50])
        cc_counts = [connected_components_at_threshold(D, t) for t in thresholds]
        # Approximate H0 max lifetime = threshold where N_CC drops to 1 (or min threshold)
        # For a simple metric: median distance at which N_CC = 10
        eps_at_10 = None
        for t, cc in zip(thresholds, cc_counts):
            if cc <= 10:
                eps_at_10 = float(t)
                break
        h0_max_lifetime = float(thresholds[-1] - thresholds[0])
        h0_mean_lifetime = float(np.diff(thresholds).mean())
        h0_n_components_at_birth = N_NAMED
        h1_n = 0
        h1_mean_lifetime = 0.0
        cc_at_median = cc_counts[2] if len(cc_counts) > 2 else N_NAMED
        method_used = "connected_components_fallback"

    # Compare STRING pairs vs non-STRING pairs by distance distribution
    i_idx = all_pairs[:, 0]
    j_idx = all_pairs[:, 1]
    dists = D[i_idx, j_idx]
    d_string = dists[string_labels == 1]
    d_nonstring = dists[string_labels == 0]
    mw_stat, mw_p = mannwhitneyu(d_string, d_nonstring, alternative='less')
    d_effect = (d_string.mean() - d_nonstring.mean()) / (d_nonstring.std() + 1e-12)

    # Null: shuffle embedding rows
    emb_perm = emb_norm[rng.permutation(N_NAMED)]
    D_perm = pairwise_distances(emb_perm)
    dists_perm = D_perm[i_idx, j_idx]
    d_string_perm = dists_perm[string_labels == 1]
    d_nonstring_perm = dists_perm[string_labels == 0]
    mw_stat_perm, mw_p_perm = mannwhitneyu(d_string_perm, d_nonstring_perm, alternative='less')
    d_effect_perm = (d_string_perm.mean() - d_nonstring_perm.mean()) / (d_nonstring_perm.std() + 1e-12)

    layer_result = {
        "layer": layer,
        "method": method_used,
        "h0_n_features": int(h0_n_components_at_birth) if HAS_RIPSER else N_NAMED,
        "h0_max_lifetime": float(h0_max_lifetime),
        "h0_mean_lifetime": float(h0_mean_lifetime),
        "h1_n_features": int(h1_n) if HAS_RIPSER else 0,
        "h1_mean_lifetime": float(h1_mean_lifetime),
        "distance_string_vs_nonstring_mw_p": float(mw_p),
        "distance_string_mean": float(d_string.mean()),
        "distance_nonstring_mean": float(d_nonstring.mean()),
        "distance_effect_size": float(d_effect),
        "perm_null_mw_p": float(mw_p_perm),
        "perm_null_effect_size": float(d_effect_perm),
    }
    if not HAS_RIPSER:
        layer_result["cc_at_median_threshold"] = int(cc_at_median)
    h03_per_layer.append(layer_result)
    print(f"  Layer {layer:2d}: dist_effect={d_effect:.3f}, mw_p={mw_p:.2e}, "
          f"perm_effect={d_effect_perm:.3f}, perm_p={mw_p_perm:.2e}", flush=True)

# Summary
effects = [r["distance_effect_size"] for r in h03_per_layer]
perm_effects = [r["perm_null_effect_size"] for r in h03_per_layer]
print(f"\n  Mean distance effect (observed): {np.mean(effects):.3f}")
print(f"  Mean distance effect (perm null): {np.mean(perm_effects):.3f}")
print(f"  Observed negative across all layers: {all(e < 0 for e in effects)}")

h03_result = {
    "has_ripser": HAS_RIPSER,
    "n_named_genes": N_NAMED,
    "per_layer": h03_per_layer,
    "summary": {
        "mean_distance_effect_observed": float(np.mean(effects)),
        "mean_distance_effect_perm_null": float(np.mean(perm_effects)),
        "all_layers_negative_effect": bool(all(e < 0 for e in effects)),
        "n_layers_significant_p01": int(sum(r["distance_string_vs_nonstring_mw_p"] < 0.01 for r in h03_per_layer)),
        "n_layers_perm_significant_p01": int(sum(r["perm_null_mw_p"] < 0.01 for r in h03_per_layer))
    }
}
with open(ITER_DIR / "h03_persistent_homology.json", "w") as f:
    json.dump(h03_result, f, indent=2)
print("  H03 saved.", flush=True)

print("\n=== All hypotheses complete ===", flush=True)
print(f"Artifacts in: {ITER_DIR}", flush=True)
