"""
iter_0021 Multi-Hypothesis Screen

H01 (manifold_distance, refinement of H03 iter_0020):
    STRING confidence quintile × pairwise distance gradient.
    Does proximity in unit-sphere space SCALE with STRING confidence?
    Tests: split STRING pairs into 5 quintiles by score;
    compute mean distance per quintile at each layer;
    Spearman(quintile_rank, mean_distance) across layers.
    Also compute AUROC for STRING binary (score>=0.4) as baseline vs quantitative ranking.

H02 (manifold_distance, new_method):
    TRRUST pairs distance test — specificity control.
    Compute pairwise distances for TRRUST pairs (TF-target) vs non-TRRUST pairs.
    Test: do TRRUST regulation pairs show the same proximity effect as PPI (STRING) pairs?
    If TRRUST pairs are NOT closer (p>0.05), this confirms that geometry encodes PPI not regulation.
    Compare direction and magnitude to STRING H01 effect.

H03 (persistent_homology, new_method):
    H1 Betti curve trajectory across layers + bootstrap CIs for composite enrichment.
    Load iter_0020 H03 data (per-layer H1 lifetime mean).
    Test: is H1 mean lifetime non-random across layers (trend test)?
    Bootstrap CIs (N=200) for the multi-axis co-polarity enrichment at count=3 across layers.
    Trajectory: does H1 increase/decrease monotonically with layer depth?
"""

import numpy as np
import json
import csv
from pathlib import Path
from scipy.stats import mannwhitneyu, spearmanr, kendalltau
import warnings
warnings.filterwarnings("ignore")

# ─── Paths ────────────────────────────────────────────────────────────────────
PROJECT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
               "/subproject_41_claude_topology_hypothesis_screening_autoloop")
ITER_DIR = PROJECT / "iterations" / "iter_0021"
ITER15_DIR = PROJECT / "iterations" / "iter_0015"
ITER20_DIR = PROJECT / "iterations" / "iter_0020"
CYCLE1 = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
              "/subproject_38_geometric_residual_stream_interpretability"
              "/implementation/outputs/cycle1_main")

EMB_PATH = CYCLE1 / "layer_gene_embeddings.npy"
EDGE_PATH = CYCLE1 / "cycle1_edge_dataset.tsv"
STRING_API_CACHE = ITER15_DIR / "string_ppi_score04_cache.json"
TRRUST_PATH = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
                   "/single_cell_mechinterp/external/networks/trrust_human.tsv")

ITER_DIR.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(42)

# ─── Load embeddings ──────────────────────────────────────────────────────────
print("Loading scGPT embeddings ...", flush=True)
emb = np.load(EMB_PATH)   # [12, 4803, 512]
N_LAYERS, N_GENES_TOTAL, N_DIM = emb.shape
print(f"  emb shape: {emb.shape}", flush=True)

# ─── Load gene list ───────────────────────────────────────────────────────────
print("Loading gene list ...", flush=True)
gene_list_path = CYCLE1 / "gene_list.txt"
with open(gene_list_path) as f:
    vocab_genes = [line.strip() for line in f]
gene_to_emb_idx = {g: i for i, g in enumerate(vocab_genes)}

# ─── Named genes from edge dataset ───────────────────────────────────────────
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

# ─── Build all pairs ──────────────────────────────────────────────────────────
print("Building pair arrays ...", flush=True)
all_pairs = []
for i in range(N_NAMED):
    for j in range(i+1, N_NAMED):
        all_pairs.append((i, j))
all_pairs = np.array(all_pairs)
N_PAIRS = len(all_pairs)
print(f"  Total pairs: {N_PAIRS}", flush=True)

# ─── Load STRING with scores ──────────────────────────────────────────────────
print("Loading STRING cache (with scores) ...", flush=True)
with open(STRING_API_CACHE) as f:
    string_cache = json.load(f)

string_pair_scores = {}
if isinstance(string_cache, dict) and "pairs" in string_cache:
    for item in string_cache["pairs"]:
        a, b = item["g1"], item["g2"]
        score = item["score"]
        string_pair_scores[(min(a, b), max(a, b))] = score
elif isinstance(string_cache, dict):
    for k, v in string_cache.items():
        if isinstance(v, (int, float)):
            a, b = k.split("__")
            string_pair_scores[(min(a, b), max(a, b))] = v

print(f"  STRING edges with scores: {len(string_pair_scores)}", flush=True)

# Pair-level string scores
string_scores_arr = np.zeros(N_PAIRS, dtype=float)
string_labels_arr = np.zeros(N_PAIRS, dtype=int)
for k, (i, j) in enumerate(all_pairs):
    gi, gj = named_genes[i], named_genes[j]
    key = (min(gi, gj), max(gi, gj))
    if key in string_pair_scores:
        string_scores_arr[k] = string_pair_scores[key]
        if string_pair_scores[key] >= 0.4:
            string_labels_arr[k] = 1

print(f"  STRING positives (>=0.4): {string_labels_arr.sum()}", flush=True)

# ─── Load TRRUST ──────────────────────────────────────────────────────────────
print("Loading TRRUST ...", flush=True)
trrust_pairs = set()
with open(TRRUST_PATH) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            tf, tgt = parts[0].strip().upper(), parts[1].strip().upper()
            trrust_pairs.add((tf, tgt))
            trrust_pairs.add((tgt, tf))

trrust_labels_arr = np.zeros(N_PAIRS, dtype=int)
for k, (i, j) in enumerate(all_pairs):
    gi, gj = named_genes[i].upper(), named_genes[j].upper()
    if (gi, gj) in trrust_pairs or (gj, gi) in trrust_pairs:
        trrust_labels_arr[k] = 1

print(f"  TRRUST positives: {trrust_labels_arr.sum()}", flush=True)

# ─── Distance helper ──────────────────────────────────────────────────────────
def compute_pairwise_distances(emb_layer, named_idx):
    """L2-normalize then compute pairwise Euclidean distances for named genes."""
    E = emb_layer[named_idx]  # [N_NAMED, 512]
    norms = np.linalg.norm(E, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-12)
    E_norm = E / norms
    # Pairwise distances
    i_idx = all_pairs[:, 0]
    j_idx = all_pairs[:, 1]
    diff = E_norm[i_idx] - E_norm[j_idx]
    dists = np.sqrt(np.sum(diff**2, axis=1))
    return dists

# ─── H01: STRING Confidence Quintile × Distance Gradient ─────────────────────
print("\n=== H01: STRING quintile x distance gradient ===", flush=True)

# Only STRING-positive pairs (score >= 0.4)
str_mask = string_labels_arr == 1
str_pair_scores = string_scores_arr[str_mask]
N_STRING = str_mask.sum()
print(f"  STRING pairs for analysis: {N_STRING}", flush=True)

# Quintile cutoffs
quintile_edges = np.percentile(str_pair_scores, [0, 20, 40, 60, 80, 100])
print(f"  Quintile edges: {quintile_edges}", flush=True)

h01_per_layer = []
spearman_per_layer = []

for layer in range(N_LAYERS):
    dists = compute_pairwise_distances(emb[layer], named_idx)
    str_dists = dists[str_mask]

    # Assign quintile bins
    quintile_bins = np.digitize(str_pair_scores, quintile_edges[1:-1])  # 0-indexed 0..4

    quintile_means = []
    quintile_ns = []
    for q in range(5):
        mask_q = quintile_bins == q
        if mask_q.sum() > 0:
            quintile_means.append(float(str_dists[mask_q].mean()))
            quintile_ns.append(int(mask_q.sum()))
        else:
            quintile_means.append(None)
            quintile_ns.append(0)

    # Spearman between quintile rank and mean distance
    valid = [i for i, m in enumerate(quintile_means) if m is not None]
    if len(valid) >= 3:
        rho, pval = spearmanr(valid, [quintile_means[i] for i in valid])
    else:
        rho, pval = 0.0, 1.0

    # Also: Mann-Whitney low (Q1=q0) vs high (Q4=q4) STRING confidence
    low_mask = quintile_bins == 0
    high_mask = quintile_bins == 4
    if low_mask.sum() > 0 and high_mask.sum() > 0:
        mw_stat, mw_p = mannwhitneyu(str_dists[high_mask], str_dists[low_mask], alternative='less')
        high_mean = float(str_dists[high_mask].mean())
        low_mean = float(str_dists[low_mask].mean())
        effect = high_mean - low_mean
    else:
        mw_p, effect = 1.0, 0.0
        high_mean, low_mean = None, None

    h01_per_layer.append({
        "layer": layer,
        "quintile_means": quintile_means,
        "quintile_ns": quintile_ns,
        "spearman_rho": float(rho),
        "spearman_p": float(pval),
        "high_vs_low_mw_p": float(mw_p),
        "high_confidence_mean_dist": high_mean,
        "low_confidence_mean_dist": low_mean,
        "effect_high_minus_low": float(effect)
    })
    spearman_per_layer.append(float(rho))

    print(f"  Layer {layer:2d}: Spearman rho={rho:.3f} p={pval:.4f}, "
          f"Q4 mean dist={high_mean:.4f} Q0 mean dist={low_mean:.4f} "
          f"effect={effect:.4f}", flush=True)

# Summary
mean_rho = float(np.mean(spearman_per_layer))
n_neg_rho = int(sum(r < 0 for r in spearman_per_layer))
n_sig = int(sum(h["spearman_p"] < 0.05 for h in h01_per_layer))
print(f"\n  Mean Spearman rho across layers: {mean_rho:.4f}", flush=True)
print(f"  Layers with negative rho (high conf = closer): {n_neg_rho}/12", flush=True)
print(f"  Layers with p<0.05 Spearman: {n_sig}/12", flush=True)

h01_result = {
    "hypothesis": "H01",
    "n_string_pairs": int(N_STRING),
    "quintile_edges": quintile_edges.tolist(),
    "per_layer": h01_per_layer,
    "summary": {
        "mean_spearman_rho": mean_rho,
        "n_layers_negative_rho": n_neg_rho,
        "n_layers_spearman_p05": n_sig,
        "direction": "negative rho = high confidence pairs are closer (expected direction)"
    }
}
with open(ITER_DIR / "h01_string_quintile_distance.json", "w") as f:
    json.dump(h01_result, f, indent=2)
print("  H01 saved.", flush=True)

# ─── H02: TRRUST Pairs Distance — Specificity Control ────────────────────────
print("\n=== H02: TRRUST pairs distance specificity control ===", flush=True)

N_TRRUST = trrust_labels_arr.sum()
print(f"  TRRUST pairs: {N_TRRUST}", flush=True)
print(f"  Non-TRRUST pairs: {(~trrust_labels_arr.astype(bool)).sum()}", flush=True)

h02_per_layer = []
trrust_effects = []
string_effects = []

for layer in range(N_LAYERS):
    dists = compute_pairwise_distances(emb[layer], named_idx)

    # TRRUST analysis
    tr_dists = dists[trrust_labels_arr == 1]
    non_tr_dists = dists[trrust_labels_arr == 0]

    if len(tr_dists) > 0 and len(non_tr_dists) > 0:
        mw_stat, mw_p = mannwhitneyu(tr_dists, non_tr_dists, alternative='less')
        tr_effect = float(tr_dists.mean() - non_tr_dists.mean())
        tr_mean = float(tr_dists.mean())
        nontr_mean = float(non_tr_dists.mean())
    else:
        mw_p, tr_effect, tr_mean, nontr_mean = 1.0, 0.0, 0.0, 0.0

    # STRING analysis for comparison (same layer)
    str_dists_layer = dists[string_labels_arr == 1]
    nonstr_dists_layer = dists[string_labels_arr == 0]
    if len(str_dists_layer) > 0 and len(nonstr_dists_layer) > 0:
        str_effect = float(str_dists_layer.mean() - nonstr_dists_layer.mean())
        str_mw_stat, str_mw_p = mannwhitneyu(str_dists_layer, nonstr_dists_layer, alternative='less')
    else:
        str_effect, str_mw_p = 0.0, 1.0

    h02_per_layer.append({
        "layer": layer,
        "trrust_mean_dist": tr_mean,
        "non_trrust_mean_dist": nontr_mean,
        "trrust_effect_size": tr_effect,
        "trrust_mw_p": float(mw_p),
        "string_effect_size": str_effect,
        "string_mw_p": float(str_mw_p),
    })
    trrust_effects.append(tr_effect)
    string_effects.append(str_effect)

    print(f"  Layer {layer:2d}: TRRUST effect={tr_effect:.4f} p={mw_p:.4f} | "
          f"STRING effect={str_effect:.4f} p={str_mw_p:.4f}", flush=True)

mean_trrust_effect = float(np.mean(trrust_effects))
mean_string_effect = float(np.mean(string_effects))
n_trrust_sig = int(sum(h["trrust_mw_p"] < 0.05 for h in h02_per_layer))
n_trrust_neg = int(sum(e < 0 for e in trrust_effects))
print(f"\n  Mean TRRUST effect: {mean_trrust_effect:.4f}, n_negative: {n_trrust_neg}/12", flush=True)
print(f"  Mean STRING effect: {mean_string_effect:.4f}", flush=True)
print(f"  TRRUST layers significant (p<0.05): {n_trrust_sig}/12", flush=True)

h02_result = {
    "hypothesis": "H02",
    "n_trrust_pairs": int(N_TRRUST),
    "per_layer": h02_per_layer,
    "summary": {
        "mean_trrust_effect": mean_trrust_effect,
        "mean_string_effect": mean_string_effect,
        "n_trrust_negative_effect": n_trrust_neg,
        "n_trrust_sig_p05": n_trrust_sig,
        "interpretation": "If TRRUST effect ≈ 0 while STRING effect << 0, geometry is PPI-specific not TF-regulation"
    }
}
with open(ITER_DIR / "h02_trrust_distance_control.json", "w") as f:
    json.dump(h02_result, f, indent=2)
print("  H02 saved.", flush=True)

# ─── H03: H1 Persistence Layer Trajectory + Bootstrap CIs ────────────────────
print("\n=== H03: H1 persistence trajectory + bootstrap CIs ===", flush=True)

# Load iter_0020 H03 per-layer data
with open(ITER20_DIR / "h03_persistent_homology.json") as f:
    h03_prev = json.load(f)

per_layer_prev = h03_prev["per_layer"]
h1_lifetimes = [r["h1_mean_lifetime"] for r in per_layer_prev]
h0_lifetimes = [r["h0_mean_lifetime"] for r in per_layer_prev]
layers_seq = list(range(len(h1_lifetimes)))

print(f"  H1 mean lifetimes by layer: {[f'{v:.4f}' for v in h1_lifetimes]}", flush=True)
print(f"  H0 mean lifetimes by layer: {[f'{v:.4f}' for v in h0_lifetimes]}", flush=True)

# Trend test: Spearman(layer, H1_lifetime)
rho_h1, p_h1 = spearmanr(layers_seq, h1_lifetimes)
rho_h0, p_h0 = spearmanr(layers_seq, h0_lifetimes)
print(f"  H1 Spearman(layer, lifetime): rho={rho_h1:.3f}, p={p_h1:.4f}", flush=True)
print(f"  H0 Spearman(layer, lifetime): rho={rho_h0:.3f}, p={p_h0:.4f}", flush=True)

# Bootstrap CIs for multi-axis co-polarity enrichment at count=3 (layer 8)
print("\n  Bootstrap CIs for co-polarity enrichment (count=3, layer=8) ...", flush=True)

# Compute multi-axis co-polarity at layer 8
emb_l8 = emb[8][named_idx]  # [N_NAMED, 512]
centered = emb_l8 - emb_l8.mean(axis=0, keepdims=True)
U, S, Vt = np.linalg.svd(centered, full_matrices=False)

# Axes SV2, SV3, SV4 (index 1, 2, 3)
copole_count = np.zeros(N_PAIRS, dtype=int)
for ax in [1, 2, 3]:
    proj = U[:, ax]
    i_idx = all_pairs[:, 0]
    j_idx = all_pairs[:, 1]
    same_sign = (proj[i_idx] * proj[j_idx]) > 0
    copole_count += same_sign.astype(int)

# Global enrichment: fraction STRING in count==3 vs overall
mask_c3 = copole_count == 3
n_c3 = mask_c3.sum()
n_str_c3 = string_labels_arr[mask_c3].sum()
base_rate = string_labels_arr.mean()
obs_rate_c3 = float(n_str_c3 / n_c3) if n_c3 > 0 else 0.0
obs_enrichment = obs_rate_c3 / base_rate if base_rate > 0 else 0.0
print(f"  Observed: N_c3={n_c3}, string_rate={obs_rate_c3:.4f}, "
      f"baseline={base_rate:.4f}, enrichment={obs_enrichment:.4f}", flush=True)

# Bootstrap CI on enrichment
N_BOOT = 500
boot_enrichments = []
for _ in range(N_BOOT):
    boot_idx = rng.integers(0, N_PAIRS, size=N_PAIRS)
    boot_labels = string_labels_arr[boot_idx]
    boot_counts = copole_count[boot_idx]

    b_mask_c3 = boot_counts == 3
    b_n_c3 = b_mask_c3.sum()
    b_base = boot_labels.mean()
    if b_n_c3 > 0 and b_base > 0:
        b_rate_c3 = boot_labels[b_mask_c3].sum() / b_n_c3
        boot_enrichments.append(float(b_rate_c3 / b_base))

boot_enrichments = np.array(boot_enrichments)
ci_lo, ci_hi = np.percentile(boot_enrichments, [2.5, 97.5])
boot_mean = float(boot_enrichments.mean())
print(f"  Bootstrap mean={boot_mean:.4f}, 95% CI=[{ci_lo:.4f}, {ci_hi:.4f}]", flush=True)
print(f"  CI excludes 1.0 (null): {ci_lo > 1.0}", flush=True)

# Layer-stratified enrichment
layer_enrichments = []
for layer in range(N_LAYERS):
    emb_l = emb[layer][named_idx]
    c_l = emb_l - emb_l.mean(axis=0, keepdims=True)
    U_l, _, _ = np.linalg.svd(c_l, full_matrices=False)
    cnt = np.zeros(N_PAIRS, dtype=int)
    for ax in [1, 2, 3]:
        p_l = U_l[:, ax]
        same = (p_l[all_pairs[:, 0]] * p_l[all_pairs[:, 1]]) > 0
        cnt += same.astype(int)

    m3 = cnt == 3
    n3 = m3.sum()
    if n3 > 0 and base_rate > 0:
        r3 = string_labels_arr[m3].sum() / n3
        enr = float(r3 / base_rate)
    else:
        enr = 1.0
    layer_enrichments.append(enr)
    print(f"  Layer {layer:2d}: N_c3={n3}, enrichment={enr:.4f}", flush=True)

rho_enr, p_enr = spearmanr(layers_seq, layer_enrichments)
print(f"\n  Enrichment vs layer Spearman: rho={rho_enr:.3f}, p={p_enr:.4f}", flush=True)

h03_result = {
    "hypothesis": "H03",
    "h1_trajectory": {
        "h1_mean_lifetimes_by_layer": h1_lifetimes,
        "h0_mean_lifetimes_by_layer": h0_lifetimes,
        "h1_spearman_rho": float(rho_h1),
        "h1_spearman_p": float(p_h1),
        "h0_spearman_rho": float(rho_h0),
        "h0_spearman_p": float(p_h0),
    },
    "bootstrap_enrichment": {
        "n_bootstrap": N_BOOT,
        "observed_enrichment": obs_enrichment,
        "boot_mean": boot_mean,
        "boot_ci_lo": float(ci_lo),
        "boot_ci_hi": float(ci_hi),
        "ci_excludes_null": bool(ci_lo > 1.0),
        "layer": 8
    },
    "layer_stratified_enrichment": layer_enrichments,
    "enrichment_vs_layer_spearman_rho": float(rho_enr),
    "enrichment_vs_layer_spearman_p": float(p_enr),
}
with open(ITER_DIR / "h03_h1_trajectory_bootstrap.json", "w") as f:
    json.dump(h03_result, f, indent=2)
print("  H03 saved.", flush=True)

print("\n=== All iter_0021 hypotheses complete ===", flush=True)
print(f"Artifacts in: {ITER_DIR}", flush=True)
