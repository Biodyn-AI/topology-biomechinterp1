"""
H03 (fixed): scGPT attention co-occurrence geometry for STRING/TRRUST pairs
Using correct attention matrix gene ordering from processed.h5ad var_names.
"""

import numpy as np
import json
import csv
from pathlib import Path
from scipy.stats import mannwhitneyu

ITER_DIR = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0018"
)
ITER15_DIR = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0015"
)
EDGE_PATH = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_38_geometric_residual_stream_interpretability"
    "/implementation/outputs/cycle1_main/cycle1_edge_dataset.tsv"
)
STRING_API_CACHE = ITER15_DIR / "string_ppi_score04_cache.json"
TRRUST_PATH = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/single_cell_mechinterp/external/networks/trrust_human.tsv"
)
ATT_PATH = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/single_cell_mechinterp/outputs/invariant_causal_edges/lung"
    "/attention_scores.npy"
)
H5AD_PATH = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/single_cell_mechinterp/outputs/invariant_causal_edges/lung"
    "/processed.h5ad"
)

ITER_DIR.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(42)

# ─── Load attention matrix with correct gene ordering ─────────────────────────
print("Loading attention matrix with correct gene ordering ...", flush=True)
att = np.load(ATT_PATH)
print(f"  att shape: {att.shape}", flush=True)

import anndata
adata = anndata.read_h5ad(str(H5AD_PATH))
gene_names = adata.var_names.tolist()
att_gene_to_idx = {g: i for i, g in enumerate(gene_names)}
print(f"  Loaded gene->index map: {len(att_gene_to_idx)} genes", flush=True)

# ─── Load named genes (scGPT 209-gene universe) ────────────────────────────────
gene_to_emb_idx = {}
with open(EDGE_PATH) as f:
    rdr = csv.DictReader(f, delimiter="\t")
    for row in rdr:
        src, si = row["source"], row.get("source_idx", "")
        tgt, ti = row["target"], row.get("target_idx", "")
        if src and si and src not in gene_to_emb_idx:
            gene_to_emb_idx[src] = int(si)
        if tgt and ti and tgt not in gene_to_emb_idx:
            gene_to_emb_idx[tgt] = int(ti)
named_genes = sorted(gene_to_emb_idx.keys())
gene_to_pos = {g: i for i, g in enumerate(named_genes)}

# Filter to genes in both scGPT and attention matrix
valid_genes = [g for g in named_genes if g in att_gene_to_idx]
valid_att_idx = np.array([att_gene_to_idx[g] for g in valid_genes])
print(f"  Valid named genes (in att matrix): {len(valid_genes)}", flush=True)

# ─── Load STRING pairs ─────────────────────────────────────────────────────────
with open(STRING_API_CACHE) as f:
    string_cache = json.load(f)

string_pairs = set()
string_score_map = {}
pair_list_raw = string_cache.get("pairs", string_cache) if isinstance(string_cache, dict) else string_cache
for entry in pair_list_raw:
    a = entry.get("g1") or entry.get("preferredName_A")
    b = entry.get("g2") or entry.get("preferredName_B")
    sc = entry.get("score", 0)
    if a and b and a in gene_to_pos and b in gene_to_pos:
        key = (min(a, b), max(a, b))
        string_pairs.add(key)
        string_score_map[key] = sc
print(f"  STRING pairs: {len(string_pairs)}", flush=True)

# ─── Load TRRUST ───────────────────────────────────────────────────────────────
trrust_activation = set()
with open(TRRUST_PATH) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        tf, tgt, reg = parts[0], parts[1], parts[2]
        if tf in gene_to_pos and tgt in gene_to_pos:
            key = (min(tf, tgt), max(tf, tgt))
            if "Activation" in reg:
                trrust_activation.add(key)
print(f"  TRRUST act: {len(trrust_activation)}", flush=True)

# ─── Build pair sets with attention indices ────────────────────────────────────
valid_set = set(valid_genes)

string_att_pairs = []
for (a, b) in string_pairs:
    if a in valid_set and b in valid_set:
        ai, bi = att_gene_to_idx[a], att_gene_to_idx[b]
        string_att_pairs.append((a, b, ai, bi, string_score_map[(min(a,b), max(a,b))]))
print(f"  STRING pairs with valid att idx: {len(string_att_pairs)}", flush=True)

trrust_att_pairs = []
for (a, b) in trrust_activation:
    if a in valid_set and b in valid_set:
        ai, bi = att_gene_to_idx[a], att_gene_to_idx[b]
        trrust_att_pairs.append((a, b, ai, bi))
print(f"  TRRUST act pairs with valid att idx: {len(trrust_att_pairs)}", flush=True)

# ─── Compute attention values ──────────────────────────────────────────────────
string_att_vals = np.array([
    (att[ai, bi] + att[bi, ai]) / 2
    for a, b, ai, bi, sc in string_att_pairs
])
trrust_att_vals = np.array([
    (att[ai, bi] + att[bi, ai]) / 2
    for a, b, ai, bi in trrust_att_pairs
])

# Random background: all pairwise within valid_genes universe
N_RAND_BG = 5000
rng_bg = np.random.default_rng(42)
all_valid_pairs = [(valid_genes[i], valid_genes[j], valid_att_idx[i], valid_att_idx[j])
                   for i in range(len(valid_genes))
                   for j in range(i+1, len(valid_genes))]
print(f"  All valid pairs: {len(all_valid_pairs)}", flush=True)

rand_bg_vals = np.array([
    (att[ai, bi] + att[bi, ai]) / 2
    for a, b, ai, bi in all_valid_pairs
])
print(f"  Random background att values computed: {len(rand_bg_vals)}", flush=True)

# Filter to just non-STRING pairs for cleaner background
non_string_bg = np.array([
    (att[ai, bi] + att[bi, ai]) / 2
    for a, b, ai, bi in all_valid_pairs
    if (min(a, b), max(a, b)) not in string_pairs
])
print(f"  Non-STRING background size: {len(non_string_bg)}", flush=True)

# Stats
mw_stat_s, mw_p_s = mannwhitneyu(string_att_vals, non_string_bg, alternative="greater")
mw_stat_t, mw_p_t = mannwhitneyu(trrust_att_vals, non_string_bg, alternative="greater")

null_mean = non_string_bg.mean()
null_std = non_string_bg.std()
z_string = (string_att_vals.mean() - null_mean) / max(null_std, 1e-9)
z_trrust = (trrust_att_vals.mean() - null_mean) / max(null_std, 1e-9)

print(f"\n  STRING: n={len(string_att_vals)}, mean={string_att_vals.mean():.6f}, "
      f"null_mean={null_mean:.6f}, z={z_string:.3f}, MW_p={mw_p_s:.4e}", flush=True)
print(f"  TRRUST act: n={len(trrust_att_vals)}, mean={trrust_att_vals.mean():.6f}, "
      f"z={z_trrust:.3f}, MW_p={mw_p_t:.4e}", flush=True)

# STRING quintile gradient (confidence score vs mean attention)
from scipy.stats import spearmanr
import csv as csv_mod

string_scores_arr = np.array([sc for a, b, ai, bi, sc in string_att_pairs])
string_att_arr = string_att_vals

quintile_edges = np.quantile(string_scores_arr, [0, 0.2, 0.4, 0.6, 0.8, 1.0])
quintile_means = []
for qi in range(5):
    mask = (string_scores_arr >= quintile_edges[qi]) & (string_scores_arr < quintile_edges[qi+1]) if qi < 4 else (string_scores_arr >= quintile_edges[qi])
    q_vals = string_att_arr[mask]
    quintile_means.append(float(q_vals.mean()) if len(q_vals) > 0 else 0.0)
    print(f"  Q{qi+1} (score [{quintile_edges[qi]:.3f},{quintile_edges[qi+1]:.3f}]): n={mask.sum()}, mean_att={quintile_means[-1]:.6f}", flush=True)

rho, p_rho = spearmanr([1,2,3,4,5], quintile_means)
print(f"  Quintile gradient Spearman rho={rho:.3f}, p={p_rho:.4f}", flush=True)

h03_output = {
    "hypothesis": "H03_attention_copresence_geometry",
    "description": (
        "scGPT aggregated attention scores (8181x8181 from lung processed.h5ad). "
        "STRING/TRRUST pairs vs non-STRING background within 209-gene universe. "
        "Also STRING confidence quintile gradient on attention."
    ),
    "n_valid_named_genes": len(valid_genes),
    "n_string_att_pairs": len(string_att_pairs),
    "n_trrust_att_pairs": len(trrust_att_pairs),
    "n_non_string_bg": len(non_string_bg),
    "string_mean_att": float(string_att_vals.mean()),
    "null_mean_att": float(null_mean),
    "null_std_att": float(null_std),
    "z_score_string": float(z_string),
    "mw_p_string": float(mw_p_s),
    "z_score_trrust": float(z_trrust),
    "mw_p_trrust": float(mw_p_t),
    "quintile_means": quintile_means,
    "quintile_edges": quintile_edges.tolist(),
    "quintile_gradient_spearman_rho": float(rho),
    "quintile_gradient_spearman_p": float(p_rho),
    "conclusion": (
        "positive" if (z_string > 1.96 and mw_p_s < 0.05) else
        "negative" if z_string < 0 else
        "inconclusive"
    ),
}
with open(ITER_DIR / "h03_attention_geometry.json", "w") as f:
    json.dump(h03_output, f, indent=2)
print(f"\nSaved: h03_attention_geometry.json", flush=True)

# Update consolidated results
with open(ITER_DIR / "iter0018_results.json") as f:
    results = json.load(f)
results["h03_attention_geometry"] = h03_output
with open(ITER_DIR / "iter0018_results.json", "w") as f:
    json.dump(results, f, indent=2)
print("Updated iter0018_results.json", flush=True)
print("Done.", flush=True)
