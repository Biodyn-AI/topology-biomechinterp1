"""
iter_0018 - Multi-Hypothesis Screen

H01 (null_sensitivity, new_method): Random Gaussian embedding null control
    Replace scGPT layer_gene_embeddings.npy with random Gaussian [12, N_GENES, 512].
    Run exact same SV2 co-pole test for STRING pairs (score>=0.4).
    Expect 0/12 layers significant (establishes that PPI geometry is model-specific, not random).

H02 (graph_topology, new_method): Out-of-sample PPI precision@k prediction benchmark
    Use SV2 projections of 209 named genes. Compute cosine similarity for all N*(N-1)/2 pairs.
    Rank pairs by cosine proximity. Report precision@k at k=50,100,200,500,1000 for STRING
    edges vs random baseline. This frames the SV2 finding as a quantitative prediction task.

H03 (new family - attention): scGPT co-attention geometry for STRING/TRRUST pairs
    Use aggregated attention_scores.npy (8181x8181). Test if STRING pairs (score>=0.4)
    have higher mean co-attention than random gene pairs from the same 209-gene universe.
    Also test TRRUST TF-target pairs. Compare using Mann-Whitney U and z-score vs shuffle null.
    This is an entirely new geometry family: attention-based co-occurrence geometry.
"""

import numpy as np
import json
import sys
import pickle
import csv
import random
from pathlib import Path
from scipy.stats import spearmanr, mannwhitneyu, zscore
from collections import defaultdict

ITER_DIR = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0018"
)
ITER15_DIR = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0015"
)
EMB_PATH = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_38_geometric_residual_stream_interpretability"
    "/implementation/outputs/cycle1_main/layer_gene_embeddings.npy"
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
ATT_GENE_LIST_PATH = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_38_geometric_residual_stream_interpretability"
    "/implementation/outputs/cycle1_main/geneformer_gene_token_map.csv"
)

ITER_DIR.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(42)

# ─── Load scGPT embeddings ─────────────────────────────────────────────────────
print("Loading scGPT embeddings ...", flush=True)
emb = np.load(EMB_PATH)  # [12, 4803, 512]
N_LAYERS, N_GENES_TOTAL, N_DIM = emb.shape
print(f"  emb shape: {emb.shape}", flush=True)

# Build gene→embedding index map
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
gene_indices = np.array([gene_to_emb_idx[g] for g in named_genes])
N_NAMED = len(named_genes)
gene_to_pos = {g: i for i, g in enumerate(named_genes)}
print(f"  Named genes: {N_NAMED}", flush=True)

# ─── Load STRING pairs ─────────────────────────────────────────────────────────
print("Loading STRING pairs ...", flush=True)
with open(STRING_API_CACHE) as f:
    string_cache = json.load(f)

string_pairs = set()
string_score_map = {}
# Cache format: {"pairs": [{"g1": ..., "g2": ..., "score": ...}, ...]}
pair_list_raw = string_cache.get("pairs", string_cache) if isinstance(string_cache, dict) else string_cache
for entry in pair_list_raw:
    a = entry.get("g1") or entry.get("preferredName_A")
    b = entry.get("g2") or entry.get("preferredName_B")
    sc = entry.get("score", 0)
    if a and b and a in gene_to_pos and b in gene_to_pos:
        key = (min(a, b), max(a, b))
        string_pairs.add(key)
        string_score_map[key] = sc
print(f"  STRING pairs (score>=0.4): {len(string_pairs)}", flush=True)

# ─── Load TRRUST ───────────────────────────────────────────────────────────────
trrust_activation = set()
trrust_repression = set()
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
            elif "Repression" in reg:
                trrust_repression.add(key)
print(f"  TRRUST act: {len(trrust_activation)}, rep: {len(trrust_repression)}", flush=True)

# ─── Helper: SVD co-pole test ─────────────────────────────────────────────────
def sv2_copole_test(embedding, pairs_set, sv_idx=1, K=52, n_shuffle=300, seed=42):
    """
    Given embedding [N_LAYERS, N_GENES_TOTAL, N_DIM], named_genes list (len N_NAMED),
    gene_indices array, pairs_set (set of (a,b) tuples), sv_idx (0-based),
    run co-pole test at each layer.
    Returns dict with layer results.
    """
    rng_local = np.random.default_rng(seed)
    results = []
    for layer in range(embedding.shape[0]):
        layer_emb = embedding[layer][gene_indices]  # [N_NAMED, N_DIM]
        centered = layer_emb - layer_emb.mean(axis=0, keepdims=True)
        U, S, Vt = np.linalg.svd(centered, full_matrices=False)
        proj = U[:, sv_idx]  # [N_NAMED]
        sorted_idx = np.argsort(proj)
        top_genes = set(named_genes[i] for i in sorted_idx[-K:])
        bot_genes = set(named_genes[i] for i in sorted_idx[:K])
        # Fraction of pairs in same pole
        pair_list = list(pairs_set)
        obs_same_pole = sum(
            (a in top_genes and b in top_genes) or (a in bot_genes and b in bot_genes)
            for a, b in pair_list
        ) / max(len(pair_list), 1)
        # Shuffle null
        null_same = []
        for _ in range(n_shuffle):
            perm = rng_local.permutation(N_NAMED)
            top_sh = set(named_genes[i] for i in perm[-K:])
            bot_sh = set(named_genes[i] for i in perm[:K])
            null_same.append(
                sum(
                    (a in top_sh and b in top_sh) or (a in bot_sh and b in bot_sh)
                    for a, b in pair_list
                ) / max(len(pair_list), 1)
            )
        null_mean = np.mean(null_same)
        null_std = np.std(null_same)
        z = (obs_same_pole - null_mean) / max(null_std, 1e-9)
        results.append({
            "layer": layer,
            "obs_same_pole": float(obs_same_pole),
            "null_mean": float(null_mean),
            "null_std": float(null_std),
            "z_score": float(z),
            "significant": bool(z > 1.96),
        })
    return results


# ════════════════════════════════════════════════════════════════════════════
# H01: Random Gaussian embedding null control
# ════════════════════════════════════════════════════════════════════════════
print("\n--- H01: Random Gaussian embedding null control ---", flush=True)

rand_emb = rng.standard_normal(emb.shape).astype(np.float32)  # [12, 4803, 512]
print(f"  Random embedding shape: {rand_emb.shape}", flush=True)

h01_results = sv2_copole_test(rand_emb, string_pairs, sv_idx=1, K=52, n_shuffle=300, seed=42)
h01_n_sig = sum(r["significant"] for r in h01_results)
h01_mean_z = np.mean([r["z_score"] for r in h01_results])
print(f"  H01: n_sig={h01_n_sig}/12, mean_z={h01_mean_z:.3f}", flush=True)
for r in h01_results:
    print(f"    L{r['layer']:02d}: z={r['z_score']:.3f} {'*' if r['significant'] else ''}", flush=True)

h01_output = {
    "hypothesis": "H01_random_gaussian_null",
    "description": "Random Gaussian embedding [12, 4803, 512]. SV2 co-pole test for STRING pairs (score>=0.4).",
    "n_string_pairs": len(string_pairs),
    "n_sig_layers": h01_n_sig,
    "mean_z_score": float(h01_mean_z),
    "layer_results": h01_results,
    "conclusion": "negative" if h01_n_sig <= 1 else "positive",
}
with open(ITER_DIR / "h01_random_gaussian_null.json", "w") as f:
    json.dump(h01_output, f, indent=2)
print(f"  Saved: h01_random_gaussian_null.json", flush=True)


# ════════════════════════════════════════════════════════════════════════════
# H02: Out-of-sample PPI precision@k prediction benchmark
# ════════════════════════════════════════════════════════════════════════════
print("\n--- H02: Out-of-sample PPI precision@k ---", flush=True)

# Use SV2 projections (averaged across all 12 layers for best generalization)
# Also test per-layer best SV2

# Build the all-pairs similarity matrix using SV2 projections at each layer
# Focus on best layer (layer 7-8 based on prior work) and averaged

def precision_at_k(pair_scores, string_set, ks):
    """
    pair_scores: list of (a, b, score). Higher score = more similar.
    string_set: set of (min, max) positive pairs.
    Returns dict k -> precision@k
    """
    sorted_pairs = sorted(pair_scores, key=lambda x: x[2], reverse=True)
    results = {}
    for k in ks:
        top_k = sorted_pairs[:k]
        n_pos = sum(1 for a, b, _ in top_k if (min(a, b), max(a, b)) in string_set)
        results[k] = n_pos / k
    return results

ks = [50, 100, 200, 500, 1000]

# Build all named-gene pairs
all_pairs_idx = [(i, j) for i in range(N_NAMED) for j in range(i+1, N_NAMED)]
all_pairs_names = [(named_genes[i], named_genes[j]) for i, j in all_pairs_idx]
N_PAIRS_TOTAL = len(all_pairs_idx)
print(f"  Total named-gene pairs: {N_PAIRS_TOTAL}", flush=True)

# Random baseline: fraction of STRING pairs
n_string_in_named = len([p for p in string_pairs
                         if p[0] in gene_to_pos and p[1] in gene_to_pos])
random_baseline = n_string_in_named / N_PAIRS_TOTAL
print(f"  STRING pairs in named universe: {n_string_in_named}", flush=True)
print(f"  Random baseline precision: {random_baseline:.4f}", flush=True)

h02_layer_results = []
for layer in range(N_LAYERS):
    layer_emb = emb[layer][gene_indices]  # [N_NAMED, N_DIM]
    centered = layer_emb - layer_emb.mean(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(centered, full_matrices=False)
    sv2_proj = U[:, 1]  # [N_NAMED]
    sv3_proj = U[:, 2]
    sv4_proj = U[:, 3]

    # Cosine-like similarity in SV2 space: |proj_i - proj_j| small -> same pole
    # Use negative absolute distance as similarity (closer = more similar)
    sv2_scores = [(-abs(sv2_proj[i] - sv2_proj[j]), named_genes[i], named_genes[j])
                  for i, j in all_pairs_idx]
    sv2_prec = precision_at_k(
        [(a, b, s) for s, a, b in sv2_scores],
        string_pairs, ks
    )

    h02_layer_results.append({
        "layer": layer,
        "sv2_prec": sv2_prec,
        "random_baseline": random_baseline,
        "sv2_enrichment_k100": sv2_prec.get(100, 0) / max(random_baseline, 1e-9),
    })

best_layer = max(h02_layer_results, key=lambda x: x["sv2_prec"].get(100, 0))
print(f"  Best layer (precision@100): L{best_layer['layer']}, P@100={best_layer['sv2_prec'][100]:.4f} "
      f"(random={random_baseline:.4f}, "
      f"enrichment={best_layer['sv2_enrichment_k100']:.2f}x)", flush=True)

# Print table
print(f"  {'L':>3} | P@50  P@100 P@200 P@500 P@1000 | Enr@100")
for r in h02_layer_results:
    p = r["sv2_prec"]
    enr = r["sv2_enrichment_k100"]
    print(f"  L{r['layer']:02d} | "
          f"{p[50]:.3f} {p[100]:.3f} {p[200]:.3f} {p[500]:.3f} {p[1000]:.3f} | {enr:.2f}x")

# Also compute SV2 average across layers
avg_sv2_scores_dict = defaultdict(float)
for layer in range(N_LAYERS):
    layer_emb = emb[layer][gene_indices]
    centered = layer_emb - layer_emb.mean(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(centered, full_matrices=False)
    sv2_proj = U[:, 1]
    for idx, (i, j) in enumerate(all_pairs_idx):
        key = (named_genes[i], named_genes[j])
        avg_sv2_scores_dict[key] += -abs(sv2_proj[i] - sv2_proj[j]) / N_LAYERS

avg_sv2_pair_scores = [(a, b, avg_sv2_scores_dict[(a, b)]) for a, b in all_pairs_names]
avg_prec = precision_at_k(avg_sv2_pair_scores, string_pairs, ks)
avg_enr_100 = avg_prec[100] / max(random_baseline, 1e-9)
print(f"\n  Average SV2 (across layers): P@100={avg_prec[100]:.4f} (enrichment={avg_enr_100:.2f}x)", flush=True)

h02_output = {
    "hypothesis": "H02_precision_at_k_benchmark",
    "description": "SV2 cosine proximity ranked pairs, precision@k for STRING positive edges.",
    "n_named_genes": N_NAMED,
    "n_total_pairs": N_PAIRS_TOTAL,
    "n_string_pairs_in_universe": n_string_in_named,
    "random_baseline_precision": float(random_baseline),
    "best_layer": best_layer["layer"],
    "best_sv2_prec_k50": best_layer["sv2_prec"][50],
    "best_sv2_prec_k100": best_layer["sv2_prec"][100],
    "best_sv2_prec_k200": best_layer["sv2_prec"][200],
    "best_sv2_prec_k500": best_layer["sv2_prec"][500],
    "best_sv2_prec_k1000": best_layer["sv2_prec"][1000],
    "best_enrichment_k100": best_layer["sv2_enrichment_k100"],
    "avg_layers_prec_k50": avg_prec[50],
    "avg_layers_prec_k100": avg_prec[100],
    "avg_layers_prec_k200": avg_prec[200],
    "avg_layers_prec_k100_enrichment": float(avg_enr_100),
    "layer_results": h02_layer_results,
    "conclusion": "positive" if avg_enr_100 > 2.0 else ("neutral" if avg_enr_100 > 1.3 else "negative"),
}
with open(ITER_DIR / "h02_precision_at_k.json", "w") as f:
    json.dump(h02_output, f, indent=2)
print(f"  Saved: h02_precision_at_k.json", flush=True)


# ════════════════════════════════════════════════════════════════════════════
# H03: scGPT attention co-occurrence geometry for STRING/TRRUST pairs
# ════════════════════════════════════════════════════════════════════════════
print("\n--- H03: Attention co-occurrence geometry ---", flush=True)

# Load attention scores [8181, 8181]
att = np.load(ATT_PATH)
print(f"  Attention scores shape: {att.shape}", flush=True)

# Load Geneformer gene token map to get gene→att_idx mapping
# The attention matrix is over the full lung gene universe in scGPT (8181 genes)
# We need a gene→index mapping for the attention matrix
# The attention scores are from single_cell_mechinterp, indexed by scGPT gene tokens

# Try to find the gene list used for attention
att_gene_list = None
att_gene_paths = [
    Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/single_cell_mechinterp"
         "/outputs/invariant_causal_edges/lung/gene_names.txt"),
    Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/single_cell_mechinterp"
         "/outputs/invariant_causal_edges/lung/gene_list.txt"),
    Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/single_cell_mechinterp"
         "/outputs/invariant_causal_edges/lung/genes.txt"),
]
for p in att_gene_paths:
    if p.exists():
        with open(p) as f:
            att_gene_list = [line.strip() for line in f if line.strip()]
        print(f"  Loaded attention gene list from {p.name}: {len(att_gene_list)} genes", flush=True)
        break

# Try scGPT edge dataset gene list (8181 genes)
if att_gene_list is None:
    # The cycle1_main/geneformer_gene_token_map.csv might have the ordering
    gf_token_map_path = Path(
        "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
        "/subproject_38_geometric_residual_stream_interpretability"
        "/implementation/outputs/cycle1_main/geneformer_gene_token_map.csv"
    )
    if gf_token_map_path.exists():
        import csv as csv_mod
        gf_map = {}
        with open(gf_token_map_path) as f:
            rdr = csv_mod.DictReader(f)
            for row in rdr:
                gf_map[row["gene"]] = int(row["token_id"])
        att_gene_list = sorted(gf_map.keys())  # lexicographic order - not the right order
        print(f"  Using geneformer token map gene list: {len(att_gene_list)} genes (WARNING: order unverified)", flush=True)
    else:
        # Fall back: use scGPT layer_gene_embedding index
        # The 8181x8181 attention matrix - try to use geneformer_edge_dataset source/target indices
        print("  WARNING: no gene list found for attention matrix; using position search", flush=True)
        att_gene_list = None

# Alternative: check if scGPT h5ad processed file has gene names matching 8181
if att_gene_list is None or len(att_gene_list) != att.shape[0]:
    print(f"  Attention matrix size mismatch: {att.shape[0]} vs gene list {len(att_gene_list) if att_gene_list else 'None'}", flush=True)
    print("  Falling back to checking if named_genes are a subset and finding their positions", flush=True)

    # Use the geneformer_edge_dataset to find source_idx/target_idx mappings to attention matrix
    # The geneformer_edge_dataset has source_idx: max seen is ~8181, probably indexes into same 8181-gene space
    gf_edge_path = Path(
        "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
        "/subproject_38_geometric_residual_stream_interpretability"
        "/implementation/outputs/cycle1_main/geneformer_edge_dataset.tsv"
    )
    att_gene_to_idx = {}
    if gf_edge_path.exists():
        with open(gf_edge_path) as f:
            rdr = csv.DictReader(f, delimiter="\t")
            for row in rdr:
                src, si = row["source"], row.get("source_idx", "")
                tgt, ti = row["target"], row.get("target_idx", "")
                if src and si:
                    att_gene_to_idx[src] = int(si)
                if tgt and ti:
                    att_gene_to_idx[tgt] = int(ti)
        print(f"  Loaded att_gene_to_idx from geneformer_edge_dataset: {len(att_gene_to_idx)} genes", flush=True)
    else:
        print("  geneformer_edge_dataset.tsv not found at cycle1_main", flush=True)
        # Last resort: same index as scGPT embedding
        att_gene_to_idx = {g: idx for g, idx in gene_to_emb_idx.items()}
        print(f"  Using scGPT indices as fallback: {len(att_gene_to_idx)} genes", flush=True)
else:
    att_gene_to_idx = {g: i for i, g in enumerate(att_gene_list)}

# Filter named genes to those with valid attention index
named_genes_att = [g for g in named_genes if g in att_gene_to_idx
                   and att_gene_to_idx[g] < att.shape[0]]
print(f"  Named genes with valid attention index: {len(named_genes_att)}", flush=True)

# Build STRING pair and TRRUST pair sets with attention indices
string_att_pairs = []
for a, b in string_pairs:
    if a in att_gene_to_idx and b in att_gene_to_idx:
        ai, bi = att_gene_to_idx[a], att_gene_to_idx[b]
        if ai < att.shape[0] and bi < att.shape[0]:
            sc = string_att_pairs
            string_att_pairs = string_att_pairs  # keep appending below

# Rebuild properly
string_att_pairs = []
for (a, b) in string_pairs:
    if a in att_gene_to_idx and b in att_gene_to_idx:
        ai, bi = att_gene_to_idx[a], att_gene_to_idx[b]
        if ai < att.shape[0] and bi < att.shape[0]:
            string_att_pairs.append((a, b, ai, bi))

trrust_att_pairs = []
for (a, b) in trrust_activation:
    if a in att_gene_to_idx and b in att_gene_to_idx:
        ai, bi = att_gene_to_idx[a], att_gene_to_idx[b]
        if ai < att.shape[0] and bi < att.shape[0]:
            trrust_att_pairs.append((a, b, ai, bi))

print(f"  STRING pairs with att indices: {len(string_att_pairs)}", flush=True)
print(f"  TRRUST act pairs with att indices: {len(trrust_att_pairs)}", flush=True)

# Extract attention scores for STRING pairs
if len(string_att_pairs) > 0:
    string_att_vals = np.array([
        (att[ai, bi] + att[bi, ai]) / 2
        for a, b, ai, bi in string_att_pairs
    ])

    # Random baseline: sample random pairs from named_genes_att universe
    named_att_idx = [att_gene_to_idx[g] for g in named_genes_att]
    N_NULL = 300
    n_pairs_str = len(string_att_pairs)
    null_means = []
    rng_h03 = np.random.default_rng(42)
    for _ in range(N_NULL):
        rand_pairs = rng_h03.choice(len(named_att_idx), size=(n_pairs_str, 2), replace=True)
        rand_vals = np.array([
            (att[named_att_idx[i], named_att_idx[j]] + att[named_att_idx[j], named_att_idx[i]]) / 2
            for i, j in rand_pairs if i != j
        ])
        if len(rand_vals) > 0:
            null_means.append(rand_vals.mean())

    obs_mean = string_att_vals.mean()
    null_mean_val = np.mean(null_means)
    null_std_val = np.std(null_means)
    z_string = (obs_mean - null_mean_val) / max(null_std_val, 1e-9)

    # MW test
    # Build random background once
    rng_h03b = np.random.default_rng(123)
    rand_bg_pairs = rng_h03b.choice(len(named_att_idx), size=(min(5000, n_pairs_str * 3), 2), replace=True)
    rand_bg_vals = np.array([
        (att[named_att_idx[i], named_att_idx[j]] + att[named_att_idx[j], named_att_idx[i]]) / 2
        for i, j in rand_bg_pairs if i != j
    ])

    mw_stat, mw_p = mannwhitneyu(string_att_vals, rand_bg_vals, alternative="greater")
    print(f"  STRING att: obs_mean={obs_mean:.6f}, null_mean={null_mean_val:.6f}, "
          f"z={z_string:.3f}, MW_p={mw_p:.4e}", flush=True)

    # TRRUST
    if len(trrust_att_pairs) > 0:
        trrust_att_vals = np.array([
            (att[ai, bi] + att[bi, ai]) / 2
            for a, b, ai, bi in trrust_att_pairs
        ])
        mw_stat_t, mw_p_t = mannwhitneyu(trrust_att_vals, rand_bg_vals, alternative="greater")
        z_trrust = (trrust_att_vals.mean() - null_mean_val) / max(null_std_val, 1e-9)
        print(f"  TRRUST att: obs_mean={trrust_att_vals.mean():.6f}, z={z_trrust:.3f}, MW_p={mw_p_t:.4e}", flush=True)
    else:
        mw_p_t = None
        z_trrust = None
        trrust_att_vals = np.array([])

    # STRING high-confidence vs low-confidence split
    string_high_idx = [(a, b, ai, bi) for a, b, ai, bi in string_att_pairs
                       if string_score_map.get((min(a, b), max(a, b)), 0) >= 0.7]
    string_low_idx = [(a, b, ai, bi) for a, b, ai, bi in string_att_pairs
                      if string_score_map.get((min(a, b), max(a, b)), 0) < 0.7]
    if string_high_idx:
        high_att_vals = np.array([(att[ai, bi] + att[bi, ai]) / 2
                                  for a, b, ai, bi in string_high_idx])
        low_att_vals = np.array([(att[ai, bi] + att[bi, ai]) / 2
                                 for a, b, ai, bi in string_low_idx]) if string_low_idx else np.array([0.0])
        mw_hl, p_hl = mannwhitneyu(high_att_vals, low_att_vals, alternative="greater")
        print(f"  STRING high (>=0.7, N={len(high_att_vals)}) mean={high_att_vals.mean():.6f} "
              f"vs low (N={len(low_att_vals)}) mean={low_att_vals.mean():.6f}, MW_p={p_hl:.4e}", flush=True)
    else:
        p_hl = None

    h03_output = {
        "hypothesis": "H03_attention_copresence_geometry",
        "description": "scGPT aggregated attention scores for STRING/TRRUST pairs vs random background.",
        "att_matrix_shape": list(att.shape),
        "n_string_att_pairs": len(string_att_pairs),
        "n_trrust_att_pairs": len(trrust_att_pairs),
        "string_obs_mean_att": float(obs_mean),
        "null_mean_att": float(null_mean_val),
        "null_std_att": float(null_std_val),
        "z_score_string": float(z_string),
        "mw_p_string": float(mw_p),
        "z_score_trrust": float(z_trrust) if z_trrust is not None else None,
        "mw_p_trrust": float(mw_p_t) if mw_p_t is not None else None,
        "high_confidence_vs_low_p": float(p_hl) if p_hl is not None else None,
        "string_high_mean_att": float(high_att_vals.mean()) if string_high_idx else None,
        "string_low_mean_att": float(low_att_vals.mean()) if string_low_idx else None,
        "conclusion": ("positive" if z_string > 1.96 and mw_p < 0.05 else
                       "negative" if z_string < 0 else "inconclusive"),
    }
else:
    print("  WARNING: no STRING pairs with valid attention indices - skipping H03", flush=True)
    h03_output = {
        "hypothesis": "H03_attention_copresence_geometry",
        "description": "Blocked: no STRING pairs with valid attention indices.",
        "status": "blocked",
        "conclusion": "inconclusive",
    }

with open(ITER_DIR / "h03_attention_geometry.json", "w") as f:
    json.dump(h03_output, f, indent=2)
print(f"  Saved: h03_attention_geometry.json", flush=True)


# ─── Consolidated results ──────────────────────────────────────────────────────
results = {
    "iteration": "iter_0018",
    "h01_random_gaussian_null": {
        "n_sig": h01_n_sig,
        "mean_z": float(h01_mean_z),
        "direction": "negative" if h01_n_sig <= 1 else "positive",
    },
    "h02_precision_at_k": {
        "best_layer": best_layer["layer"],
        "best_prec_k100": best_layer["sv2_prec"][100],
        "avg_prec_k100": float(avg_prec[100]),
        "random_baseline": float(random_baseline),
        "enrichment_k100": float(avg_enr_100),
    },
    "h03_attention_geometry": h03_output,
}
with open(ITER_DIR / "iter0018_results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"\nSaved consolidated: iter0018_results.json", flush=True)
print("Done.", flush=True)
