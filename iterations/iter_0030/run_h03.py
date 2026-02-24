"""H03 standalone: Continuous slope-vs-STRING-degree"""
import numpy as np
import json
import csv
from pathlib import Path
from scipy.stats import spearmanr, mannwhitneyu
from collections import defaultdict
import warnings
warnings.filterwarnings("ignore")

PROJECT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
               "/subproject_41_claude_topology_hypothesis_screening_autoloop")
ITER_DIR = PROJECT / "iterations" / "iter_0030"
CYCLE1 = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
              "/subproject_38_geometric_residual_stream_interpretability"
              "/implementation/outputs/cycle1_main")
STRING_CACHE = PROJECT / "iterations" / "iter_0015" / "string_ppi_score04_cache.json"

# Load embeddings
emb = np.load(CYCLE1 / "layer_gene_embeddings.npy")
N_LAYERS, N_GENES_TOTAL, N_DIM = emb.shape

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
emb_named = emb[:, named_idx, :]
N_NAMED = len(named_idx)
gene_to_named_idx = {g: i for i, g in enumerate(named_genes)}

DIVERGING_14 = ["FOS", "JUNB", "TNF", "PTGS2", "HLA-A", "HLA-DPB1",
                "KLF6", "LDHA", "LGALS1", "NCAM1", "NCOA3", "NR4A3",
                "PAX5", "TBXAS1"]
diverging_set = set(DIVERGING_14)
div_idx_named = np.array([gene_to_named_idx[g] for g in DIVERGING_14 if g in gene_to_named_idx])
conv_idx_named = np.array([gene_to_named_idx[g] for g in named_genes if g not in diverging_set])

# Compute per-gene trajectory slopes
layers_arr = np.arange(N_LAYERS, dtype=float)
dist_matrix = np.zeros((N_LAYERS, N_NAMED))
for l in range(N_LAYERS):
    centroid = emb_named[l].mean(axis=0)
    dist_matrix[l] = np.linalg.norm(emb_named[l] - centroid, axis=1)

gene_slopes = np.zeros(N_NAMED)
for i in range(N_NAMED):
    rho, _ = spearmanr(layers_arr, dist_matrix[:, i])
    gene_slopes[i] = rho

# Load STRING data
with open(STRING_CACHE) as f:
    string_raw = json.load(f)
pairs = string_raw['pairs']  # list of {g1, g2, score}

string_degree = defaultdict(int)
string_weighted_degree = defaultdict(float)

for p in pairs:
    g1, g2, score = p['g1'], p['g2'], p['score']
    if g1 in gene_to_named_idx and g2 in gene_to_named_idx:
        string_degree[g1] += 1
        string_degree[g2] += 1
        string_weighted_degree[g1] += score
        string_weighted_degree[g2] += score

degree_arr = np.array([string_degree[g] for g in named_genes], dtype=float)
wdegree_arr = np.array([string_weighted_degree[g] for g in named_genes], dtype=float)

has_edge = degree_arr > 0
n_with_edge = int(has_edge.sum())
print(f"Named genes with STRING edge: {n_with_edge}/{N_NAMED}")

# Correlations
rho_deg_all, p_deg_all = spearmanr(degree_arr, gene_slopes)
rho_wdeg_all, p_wdeg_all = spearmanr(wdegree_arr, gene_slopes)
rho_deg_sub, p_deg_sub = spearmanr(degree_arr[has_edge], gene_slopes[has_edge])

print(f"rho(degree, slope) all: {rho_deg_all:.4f}, p={p_deg_all:.2e}")
print(f"rho(wdegree, slope) all: {rho_wdeg_all:.4f}, p={p_wdeg_all:.2e}")
print(f"rho(degree, slope) subset: {rho_deg_sub:.4f}, p={p_deg_sub:.2e}")

# Quartile comparison
q75 = np.percentile(degree_arr[has_edge], 75)
q25 = np.percentile(degree_arr[has_edge], 25)
high_degree = has_edge & (degree_arr >= q75)
low_degree = has_edge & (degree_arr <= q25) & (degree_arr > 0)

slopes_high = gene_slopes[high_degree]
slopes_low = gene_slopes[low_degree]
mw_stat2, mw_p2 = mannwhitneyu(slopes_high, slopes_low, alternative='less')
auroc2 = mw_stat2 / (len(slopes_high) * len(slopes_low))

print(f"High-degree slopes mean: {np.mean(slopes_high):.4f}, Low-degree: {np.mean(slopes_low):.4f}")
print(f"MW AUROC (high < low slope): {auroc2:.4f}, p={mw_p2:.2e}")

# Div vs conv degree
div_degrees = degree_arr[div_idx_named]
conv_degrees = degree_arr[conv_idx_named]
mw_stat3, mw_p3 = mannwhitneyu(div_degrees, conv_degrees, alternative='less')
auroc3 = mw_stat3 / (len(div_degrees) * len(conv_degrees))
print(f"Div degrees mean: {np.mean(div_degrees):.2f}, Conv: {np.mean(conv_degrees):.2f}")
print(f"MW AUROC (div < conv degree): {auroc3:.4f}, p={mw_p3:.2e}")

# Per-gene data
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
print(f"Saved: {h03_path}")
