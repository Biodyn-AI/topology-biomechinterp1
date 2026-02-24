"""
iter_0030 Multi-Hypothesis Screen

H01 (module_structure / new_method): GO + TRRUST enrichment of 14 diverging genes
    The 14 outlier genes (FOS, JUNB, TNF, PTGS2, HLA-A, HLA-DPB1, KLF6, LDHA,
    LGALS1, NCAM1, NCOA3, NR4A3, PAX5, TBXAS1) have increasing centroid distance
    over layers. Test if they are enriched for specific GO BP terms or TRRUST TF
    regulatory hubs relative to the 195 converging genes.
    novelty_type: new_method (biological annotation of topological split)

H02 (topology_stability / new_method): Layer-of-Bifurcation Sliding Window
    For the 14 diverging vs 195 converging genes, compute inter-group distance
    (mean L2 between group centroids) at each of 12 layers. Also compute within-
    group variance. Test if the two groups are already separated at L0 or if
    divergence increases monotonically. Find "bifurcation layer" = layer where
    inter/intra distance ratio peaks.
    novelty_type: new_method (mechanistic anatomy of split)

H03 (graph_topology / new_method): Continuous slope-vs-STRING-degree
    Each of 209 named genes gets a "trajectory slope" = Spearman rho(layer,
    dist_to_centroid). Test if STRING degree (n_edges) correlates with slope.
    High-degree STRING hubs should have more stable centroid proximity (negative
    slope). Also test STRING edge-weight sum (total connectivity).
    novelty_type: new_method (continuous structure-function relationship)
"""

import numpy as np
import json
import csv
from pathlib import Path
from scipy.stats import spearmanr, mannwhitneyu, pearsonr
from sklearn.neighbors import NearestNeighbors
from collections import defaultdict
import warnings
warnings.filterwarnings("ignore")

PROJECT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
               "/subproject_41_claude_topology_hypothesis_screening_autoloop")
ITER_DIR = PROJECT / "iterations" / "iter_0030"
ITER_DIR.mkdir(parents=True, exist_ok=True)

CYCLE1 = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
              "/subproject_38_geometric_residual_stream_interpretability"
              "/implementation/outputs/cycle1_main")
ITER15_DIR = PROJECT / "iterations" / "iter_0015"
ITER16_DIR = PROJECT / "iterations" / "iter_0016"
STRING_CACHE = ITER15_DIR / "string_ppi_score04_cache.json"

rng = np.random.default_rng(42)

# ─── Load embeddings ──────────────────────────────────────────────────────────
print("Loading embeddings...", flush=True)
emb = np.load(CYCLE1 / "layer_gene_embeddings.npy")   # [12, 4803, 512]
N_LAYERS, N_GENES_TOTAL, N_DIM = emb.shape
print(f"  Shape: {emb.shape}", flush=True)

with open(CYCLE1 / "gene_list.txt") as f:
    vocab_genes = [line.strip() for line in f if line.strip()]
gene_to_emb_idx = {g: i for i, g in enumerate(vocab_genes)}

EDGE_PATH = CYCLE1 / "cycle1_edge_dataset.tsv"
named_gene_set = set()
with open(EDGE_PATH) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        named_gene_set.add(row['source'])
        named_gene_set.add(row['target'])

named_genes = np.array(sorted(g for g in named_gene_set if g in gene_to_emb_idx))
named_idx = np.array([gene_to_emb_idx[g] for g in named_genes])
emb_named = emb[:, named_idx, :]   # [12, N_NAMED, 512]
N_NAMED = len(named_idx)
gene_to_named_idx = {g: i for i, g in enumerate(named_genes)}
print(f"  Named genes: {N_NAMED}", flush=True)

# ─── The 14 diverging genes from iter_0029 ────────────────────────────────────
DIVERGING_14 = ["FOS", "JUNB", "TNF", "PTGS2", "HLA-A", "HLA-DPB1",
                "KLF6", "LDHA", "LGALS1", "NCAM1", "NCOA3", "NR4A3",
                "PAX5", "TBXAS1"]
diverging_set = set(DIVERGING_14)

# Filter to those in named_genes
div_in_named = [g for g in DIVERGING_14 if g in gene_to_named_idx]
conv_in_named = [g for g in named_genes if g not in diverging_set]
print(f"  Diverging: {len(div_in_named)}, Converging: {len(conv_in_named)}", flush=True)

div_idx_named = np.array([gene_to_named_idx[g] for g in div_in_named])
conv_idx_named = np.array([gene_to_named_idx[g] for g in conv_in_named])

# ─── Compute per-gene distance-to-centroid trajectories ──────────────────────
print("Computing distance-to-centroid trajectories...", flush=True)
dist_matrix = np.zeros((N_LAYERS, N_NAMED))
for l in range(N_LAYERS):
    centroid = emb_named[l].mean(axis=0)
    dist_matrix[l] = np.linalg.norm(emb_named[l] - centroid, axis=1)

# Per-gene trajectory slope: Spearman rho(layer_index, dist)
layers_arr = np.arange(N_LAYERS, dtype=float)
gene_slopes = np.zeros(N_NAMED)
for i in range(N_NAMED):
    rho, _ = spearmanr(layers_arr, dist_matrix[:, i])
    gene_slopes[i] = rho

# ─── H01: GO + TRRUST enrichment of 14 diverging genes ───────────────────────
print("\n=== H01: Biological annotation of 14 diverging genes ===", flush=True)

# Load TRRUST data
TRRUST_PATH = PROJECT / "iterations" / "iter_0023" / "trrust_rawdata.human.tsv"
trrust_available = TRRUST_PATH.exists()
print(f"  TRRUST available: {trrust_available}", flush=True)

# Load GO data from iter_0023 or iter_0024
GO_CACHE = PROJECT / "iterations" / "iter_0024" / "go_bp_cache.json"
go_available = GO_CACHE.exists()
print(f"  GO cache available: {go_available}", flush=True)

h01_result = {
    "diverging_genes": div_in_named,
    "n_diverging": len(div_in_named),
    "n_converging": len(conv_in_named),
}

# Slope statistics for div vs conv
div_slopes = gene_slopes[div_idx_named]
conv_slopes = gene_slopes[conv_idx_named]
mw_stat, mw_p = mannwhitneyu(div_slopes, conv_slopes, alternative='greater')
auroc = mw_stat / (len(div_slopes) * len(conv_slopes))
h01_result["slope_div_mean"] = float(np.mean(div_slopes))
h01_result["slope_conv_mean"] = float(np.mean(conv_slopes))
h01_result["slope_mw_auroc"] = float(auroc)
h01_result["slope_mw_p"] = float(mw_p)
print(f"  Div slopes mean: {np.mean(div_slopes):.4f}, Conv: {np.mean(conv_slopes):.4f}", flush=True)
print(f"  Mann-Whitney AUROC (div > conv): {auroc:.4f}, p={mw_p:.2e}", flush=True)

# TRRUST: how many of the 14 are TF targets vs regulators
if trrust_available:
    trrust_tfs = set()
    trrust_targets = set()
    with open(TRRUST_PATH) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                trrust_tfs.add(parts[0])
                for tgt in parts[1].split(','):
                    trrust_targets.add(tgt.strip())

    div_as_tf = [g for g in div_in_named if g in trrust_tfs]
    div_as_target = [g for g in div_in_named if g in trrust_targets]
    conv_as_tf = [g for g in conv_in_named if g in trrust_tfs]
    conv_as_target = [g for g in conv_in_named if g in trrust_targets]

    h01_result["trrust_div_tfs"] = div_as_tf
    h01_result["trrust_div_targets"] = div_as_target
    h01_result["trrust_div_tf_frac"] = len(div_as_tf) / len(div_in_named) if div_in_named else 0
    h01_result["trrust_div_target_frac"] = len(div_as_target) / len(div_in_named) if div_in_named else 0
    h01_result["trrust_conv_tf_frac"] = len(conv_as_tf) / len(conv_in_named) if conv_in_named else 0
    h01_result["trrust_conv_target_frac"] = len(conv_as_target) / len(conv_in_named) if conv_in_named else 0

    print(f"  TRRUST: div TF frac={h01_result['trrust_div_tf_frac']:.3f}, conv TF frac={h01_result['trrust_conv_tf_frac']:.3f}", flush=True)
    print(f"  TRRUST: div target frac={h01_result['trrust_div_target_frac']:.3f}, conv target frac={h01_result['trrust_conv_target_frac']:.3f}", flush=True)

# GO terms for diverging genes using mygene
print("  Fetching GO terms for diverging genes via mygene...", flush=True)
try:
    import mygene
    mg = mygene.MyGeneInfo()
    # Query diverging genes
    div_results = mg.querymany(div_in_named, scopes='symbol', fields='go.BP', species='human')

    div_go_terms = defaultdict(list)
    for r in div_results:
        if 'go' in r and 'BP' in r['go']:
            gene = r.get('query', '')
            bps = r['go']['BP']
            if isinstance(bps, dict):
                bps = [bps]
            for bp in bps:
                div_go_terms[gene].append(bp.get('term', ''))

    # Collect all unique GO terms and count
    from collections import Counter
    div_go_flat = [t for terms in div_go_terms.values() for t in terms]
    div_go_counts = Counter(div_go_flat)

    # Top GO terms for diverging genes
    top_div_go = div_go_counts.most_common(20)
    h01_result["div_top_go_bp"] = [{"term": t, "count": c} for t, c in top_div_go]

    # Stress/inflammatory keyword scan
    stress_keywords = ['stress', 'inflam', 'immune', 'innate', 'response to', 'cytokine',
                       'apoptosis', 'death', 'oxidative', 'hypoxia', 'NF-kB']
    stress_hits = [(t, c) for t, c in top_div_go if any(kw.lower() in t.lower() for kw in stress_keywords)]
    h01_result["div_stress_go_hits"] = [{"term": t, "count": c} for t, c in stress_hits]
    h01_result["div_go_n_unique_terms"] = len(div_go_counts)
    h01_result["div_go_n_genes_with_terms"] = len(div_go_terms)

    print(f"  GO terms fetched. Unique terms: {len(div_go_counts)}", flush=True)
    print(f"  Top GO terms: {top_div_go[:5]}", flush=True)
    print(f"  Stress/immune GO hits: {stress_hits[:5]}", flush=True)

except Exception as e:
    print(f"  mygene error: {e}", flush=True)
    h01_result["mygene_error"] = str(e)

# Save H01
h01_path = ITER_DIR / "h01_diverging_annotation.json"
with open(h01_path, 'w') as f:
    json.dump(h01_result, f, indent=2)
print(f"  Saved: {h01_path}", flush=True)

# ─── H02: Layer-of-Bifurcation Sliding Window ────────────────────────────────
print("\n=== H02: Layer-of-Bifurcation Anatomy ===", flush=True)

h02_result = {
    "n_diverging": len(div_in_named),
    "n_converging": len(conv_in_named),
    "by_layer": [],
}

intergroup_dist = np.zeros(N_LAYERS)
intragroup_div = np.zeros(N_LAYERS)
intragroup_conv = np.zeros(N_LAYERS)
ratio_by_layer = np.zeros(N_LAYERS)

for l in range(N_LAYERS):
    div_emb = emb_named[l, div_idx_named, :]    # [14, 512]
    conv_emb = emb_named[l, conv_idx_named, :]  # [195, 512]

    div_centroid = div_emb.mean(axis=0)
    conv_centroid = conv_emb.mean(axis=0)

    # Inter-group distance = distance between centroids
    intergroup = float(np.linalg.norm(div_centroid - conv_centroid))

    # Intra-group variance = mean distance to own centroid
    intra_div = float(np.mean(np.linalg.norm(div_emb - div_centroid, axis=1)))
    intra_conv = float(np.mean(np.linalg.norm(conv_emb - conv_centroid, axis=1)))

    intergroup_dist[l] = intergroup
    intragroup_div[l] = intra_div
    intragroup_conv[l] = intra_conv

    # Ratio: inter-group / mean(intra_div, intra_conv)
    mean_intra = (intra_div + intra_conv) / 2
    ratio = intergroup / mean_intra if mean_intra > 0 else 0
    ratio_by_layer[l] = ratio

    h02_result["by_layer"].append({
        "layer": l,
        "intergroup_dist": intergroup,
        "intra_div": intra_div,
        "intra_conv": intra_conv,
        "ratio": float(ratio),
    })
    print(f"  L{l:02d}: inter={intergroup:.3f}, intra_div={intra_div:.3f}, intra_conv={intra_conv:.3f}, ratio={ratio:.4f}", flush=True)

# Bifurcation layer = where ratio peaks
bifurcation_layer = int(np.argmax(ratio_by_layer))
h02_result["bifurcation_layer"] = bifurcation_layer
h02_result["max_ratio"] = float(ratio_by_layer[bifurcation_layer])
h02_result["l0_ratio"] = float(ratio_by_layer[0])
h02_result["l11_ratio"] = float(ratio_by_layer[11])

# Spearman rho between layer and ratio
rho_ratio, p_ratio = spearmanr(layers_arr, ratio_by_layer)
h02_result["rho_layer_vs_ratio"] = float(rho_ratio)
h02_result["p_layer_vs_ratio"] = float(p_ratio)

# Spearman rho between layer and intergroup distance
rho_inter, p_inter = spearmanr(layers_arr, intergroup_dist)
h02_result["rho_layer_vs_intergroup"] = float(rho_inter)
h02_result["p_layer_vs_intergroup"] = float(p_inter)

print(f"\n  Bifurcation layer: L{bifurcation_layer} (max ratio={ratio_by_layer[bifurcation_layer]:.4f})", flush=True)
print(f"  L0 ratio: {ratio_by_layer[0]:.4f}, L11 ratio: {ratio_by_layer[11]:.4f}", flush=True)
print(f"  Spearman rho(layer, intergroup): {rho_inter:.4f}, p={p_inter:.2e}", flush=True)
print(f"  Spearman rho(layer, ratio): {rho_ratio:.4f}, p={p_ratio:.2e}", flush=True)

h02_path = ITER_DIR / "h02_bifurcation_anatomy.json"
with open(h02_path, 'w') as f:
    json.dump(h02_result, f, indent=2)
print(f"  Saved: {h02_path}", flush=True)

# ─── H03: Continuous slope-vs-STRING-degree ──────────────────────────────────
print("\n=== H03: Trajectory slope vs STRING degree ===", flush=True)

# Load STRING cache
h03_result = {}
if STRING_CACHE.exists():
    with open(STRING_CACHE) as f:
        string_data = json.load(f)

    # Build degree + weighted degree for each named gene
    string_degree = defaultdict(int)
    string_weighted_degree = defaultdict(float)

    for pair_key, score in string_data.items():
        g1, g2 = pair_key.split('__')
        if g1 in gene_to_named_idx and g2 in gene_to_named_idx:
            string_degree[g1] += 1
            string_degree[g2] += 1
            string_weighted_degree[g1] += score
            string_weighted_degree[g2] += score

    # Build arrays aligned to named_genes
    degree_arr = np.array([string_degree[g] for g in named_genes], dtype=float)
    wdegree_arr = np.array([string_weighted_degree[g] for g in named_genes], dtype=float)

    # Genes with at least 1 STRING edge
    has_edge = degree_arr > 0
    n_with_edge = int(has_edge.sum())
    print(f"  Named genes with STRING edge: {n_with_edge}/{N_NAMED}", flush=True)

    # Spearman rho: slope vs degree (all genes)
    rho_deg_all, p_deg_all = spearmanr(degree_arr, gene_slopes)
    rho_wdeg_all, p_wdeg_all = spearmanr(wdegree_arr, gene_slopes)

    # Spearman rho: slope vs degree (only genes with edges)
    rho_deg_sub, p_deg_sub = spearmanr(degree_arr[has_edge], gene_slopes[has_edge])
    rho_wdeg_sub, p_wdeg_sub = spearmanr(wdegree_arr[has_edge], gene_slopes[has_edge])

    print(f"  rho(degree, slope) all: {rho_deg_all:.4f}, p={p_deg_all:.2e}", flush=True)
    print(f"  rho(wdegree, slope) all: {rho_wdeg_all:.4f}, p={p_wdeg_all:.2e}", flush=True)
    print(f"  rho(degree, slope) subset: {rho_deg_sub:.4f}, p={p_deg_sub:.2e}", flush=True)

    # Top vs bottom degree quartile slope comparison
    q75 = np.percentile(degree_arr[has_edge], 75)
    q25 = np.percentile(degree_arr[has_edge], 25)
    high_degree = has_edge & (degree_arr >= q75)
    low_degree = has_edge & (degree_arr <= q25) & (degree_arr > 0)

    slopes_high = gene_slopes[high_degree]
    slopes_low = gene_slopes[low_degree]

    mw_stat2, mw_p2 = mannwhitneyu(slopes_high, slopes_low, alternative='less')
    auroc2 = mw_stat2 / (len(slopes_high) * len(slopes_low))

    print(f"  High-degree slopes mean: {np.mean(slopes_high):.4f}, Low-degree: {np.mean(slopes_low):.4f}", flush=True)
    print(f"  MW AUROC (high < low slope): {auroc2:.4f}, p={mw_p2:.2e}", flush=True)

    # Also test: are the 14 diverging genes degree-depleted?
    div_degrees = degree_arr[div_idx_named]
    conv_degrees = degree_arr[conv_idx_named]
    mw_stat3, mw_p3 = mannwhitneyu(div_degrees, conv_degrees, alternative='less')
    auroc3 = mw_stat3 / (len(div_degrees) * len(conv_degrees))

    print(f"  Div degrees mean: {np.mean(div_degrees):.2f}, Conv: {np.mean(conv_degrees):.2f}", flush=True)
    print(f"  MW AUROC (div < conv degree): {auroc3:.4f}, p={mw_p3:.2e}", flush=True)

    # Per-gene slope and degree
    per_gene_data = []
    for i, g in enumerate(named_genes):
        per_gene_data.append({
            "gene": g,
            "slope": float(gene_slopes[i]),
            "string_degree": int(degree_arr[i]),
            "string_wdegree": float(wdegree_arr[i]),
            "is_diverging": g in diverging_set,
        })

    h03_result = {
        "n_named": N_NAMED,
        "n_with_string_edge": n_with_edge,
        "rho_degree_vs_slope_all": float(rho_deg_all),
        "p_degree_vs_slope_all": float(p_deg_all),
        "rho_wdegree_vs_slope_all": float(rho_wdeg_all),
        "p_wdegree_vs_slope_all": float(p_wdeg_all),
        "rho_degree_vs_slope_subset": float(rho_deg_sub),
        "p_degree_vs_slope_subset": float(p_deg_sub),
        "high_degree_slope_mean": float(np.mean(slopes_high)),
        "low_degree_slope_mean": float(np.mean(slopes_low)),
        "auroc_high_vs_low_degree": float(auroc2),
        "p_high_vs_low_degree": float(mw_p2),
        "div_degree_mean": float(np.mean(div_degrees)),
        "conv_degree_mean": float(np.mean(conv_degrees)),
        "auroc_div_depleted": float(auroc3),
        "p_div_depleted": float(mw_p3),
        "per_gene": per_gene_data,
    }

    h03_path = ITER_DIR / "h03_slope_vs_degree.json"
    with open(h03_path, 'w') as f:
        json.dump(h03_result, f, indent=2)
    print(f"  Saved: {h03_path}", flush=True)
else:
    print(f"  STRING cache not found at {STRING_CACHE}", flush=True)
    h03_result = {"error": "STRING cache not found"}

# ─── Summary printout ─────────────────────────────────────────────────────────
print("\n=== SUMMARY ===", flush=True)
print(f"H01 (module_structure): Div slope mean={h01_result['slope_div_mean']:.4f}, "
      f"Conv slope mean={h01_result['slope_conv_mean']:.4f}, "
      f"AUROC={h01_result['slope_mw_auroc']:.4f}, p={h01_result['slope_mw_p']:.2e}", flush=True)

print(f"H02 (topology_stability): Bifurcation at L{bifurcation_layer}, "
      f"rho(layer,ratio)={h02_result['rho_layer_vs_ratio']:.4f}, "
      f"p={h02_result['p_layer_vs_ratio']:.2e}", flush=True)

if "rho_degree_vs_slope_all" in h03_result:
    print(f"H03 (graph_topology): rho(degree,slope)={h03_result['rho_degree_vs_slope_all']:.4f}, "
          f"p={h03_result['p_degree_vs_slope_all']:.2e}, "
          f"high_deg_slope={h03_result['high_degree_slope_mean']:.4f}, "
          f"low_deg_slope={h03_result['low_degree_slope_mean']:.4f}", flush=True)

print("\nDone.", flush=True)
