"""
iter_0019 Multi-Hypothesis Screen

H01 (cross_model_alignment, new_family): Geneformer static token embedding geometry vs scGPT SV2 co-poles
    Load Geneformer model embedding weights for 209 named genes.
    Compute SVD of Geneformer static embeddings; test if STRING pairs are more co-polar.
    Cross-model correlation: compare Geneformer SV1 projections vs scGPT SV2 projections.

H02 (graph_topology, new_method): Multi-axis composite co-pole P@k
    Use SV2+SV3+SV4 projections. For each gene pair: count number of axes where co-polar.
    Test: P@k at k=100 stratified by co-polarity count (0,1,2,3).
    Expect monotonic increase with co-polarity count.

H03 (null_sensitivity/attention, refinement): TRRUST repression vs activation scGPT attention
    Extend iter_0018 H03: activation pairs showed 2x co-attention (p=9.9e-9).
    Test repression pairs (N=64): same or different co-attention pattern?
    Compare activation vs repression vs STRING vs null.
"""

import numpy as np
import json
import csv
import pickle
from pathlib import Path
from scipy.stats import spearmanr, mannwhitneyu, pearsonr
from collections import defaultdict

# ─── Paths ────────────────────────────────────────────────────────────────────
PROJECT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
               "/subproject_41_claude_topology_hypothesis_screening_autoloop")
ITER_DIR = PROJECT / "iterations" / "iter_0019"
ITER15_DIR = PROJECT / "iterations" / "iter_0015"
CYCLE1 = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
              "/subproject_38_geometric_residual_stream_interpretability"
              "/implementation/outputs/cycle1_main")
CYCLE12 = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
               "/subproject_38_geometric_residual_stream_interpretability"
               "/implementation/outputs/cycle12_geneformer_lung_bootstrap")

EMB_PATH = CYCLE1 / "layer_gene_embeddings.npy"
EDGE_PATH = CYCLE1 / "cycle1_edge_dataset.tsv"
STRING_API_CACHE = ITER15_DIR / "string_ppi_score04_cache.json"
TRRUST_PATH = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
                   "/single_cell_mechinterp/external/networks/trrust_human.tsv")
ATT_PATH = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
                "/single_cell_mechinterp/outputs/invariant_causal_edges/lung"
                "/attention_scores.npy")
H5AD_PATH = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
                 "/single_cell_mechinterp/outputs/invariant_causal_edges/lung"
                 "/processed.h5ad")
GF_TOKEN_MAP = CYCLE12 / "geneformer_gene_token_map.csv"
GF_MODEL_CACHE = Path("/Users/ihorkendiukhov/.cache/huggingface/hub"
                      "/models--ctheodoris--Geneformer/snapshots"
                      "/05fcbeb8a27d49e0a7a4349152202ee2c1cbfd28")

ITER_DIR.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(42)

# ─── Load scGPT embeddings ────────────────────────────────────────────────────
print("Loading scGPT embeddings ...", flush=True)
emb = np.load(EMB_PATH)  # [12, 4803, 512]
N_LAYERS, N_GENES_TOTAL, N_DIM = emb.shape
print(f"  emb shape: {emb.shape}", flush=True)

# ─── Load gene list ───────────────────────────────────────────────────────────
print("Loading gene list ...", flush=True)
gene_list_path = CYCLE1 / "gene_list.txt"
with open(gene_list_path) as f:
    vocab_genes = [line.strip() for line in f]
vocab_gene_set = set(vocab_genes)
gene_to_emb_idx = {g: i for i, g in enumerate(vocab_genes)}
print(f"  Vocab size: {len(vocab_genes)}", flush=True)

# ─── Load named genes from edge dataset ───────────────────────────────────────
print("Loading named genes from edge dataset ...", flush=True)
named_gene_set = set()
with open(EDGE_PATH) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        named_gene_set.add(row['source'])
        named_gene_set.add(row['target'])
named_gene_set = named_gene_set & vocab_gene_set
named_genes = sorted(named_gene_set)
gene_indices = np.array([gene_to_emb_idx[g] for g in named_genes])
N_NAMED = len(named_genes)
gene_to_pos = {g: i for i, g in enumerate(named_genes)}
print(f"  Named genes: {N_NAMED}", flush=True)

# ─── Load STRING pairs ────────────────────────────────────────────────────────
print("Loading STRING pairs ...", flush=True)
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
print(f"  STRING pairs (score>=0.4): {len(string_pairs)}", flush=True)

# ─── Load TRRUST pairs ────────────────────────────────────────────────────────
print("Loading TRRUST pairs ...", flush=True)
trrust_activation = set()
trrust_repression = set()
with open(TRRUST_PATH) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 3:
            tf, target, direction = parts[0], parts[1], parts[2]
            if tf in gene_to_pos and target in gene_to_pos:
                key = (min(tf, target), max(tf, target))
                if 'Activation' in direction:
                    trrust_activation.add(key)
                if 'Repression' in direction:
                    trrust_repression.add(key)
print(f"  TRRUST activation pairs: {len(trrust_activation)}", flush=True)
print(f"  TRRUST repression pairs: {len(trrust_repression)}", flush=True)

# All named-gene pair universe
all_named_pairs = set()
for i in range(N_NAMED):
    for j in range(i+1, N_NAMED):
        all_named_pairs.add((named_genes[i], named_genes[j]))
print(f"  Total named-gene pairs: {len(all_named_pairs)}", flush=True)

# ─── Extract named gene embeddings ───────────────────────────────────────────
named_emb = emb[:, gene_indices, :]  # [12, N_NAMED, 512]

print("\n" + "="*60, flush=True)
print("=== H01: Geneformer Static Token Embedding Cross-Model Test ===", flush=True)
print("="*60, flush=True)

# Load Geneformer gene-to-token map
gf_gene_to_token = {}
with open(GF_TOKEN_MAP) as f:
    reader = csv.DictReader(f)
    for row in reader:
        gf_gene_to_token[row['gene']] = int(row['token_id'])
print(f"  Geneformer token map: {len(gf_gene_to_token)} genes", flush=True)

# Find 209 named genes that have Geneformer token IDs
gf_named_genes = [g for g in named_genes if g in gf_gene_to_token]
gf_token_ids = [gf_gene_to_token[g] for g in gf_named_genes]
print(f"  Named genes with GF tokens: {len(gf_named_genes)}", flush=True)

# Load Geneformer model embedding weights
print("  Loading Geneformer model ...", flush=True)
try:
    import torch
    from safetensors.torch import load_file
    # Try main model first, then V2
    safetensor_path = GF_MODEL_CACHE / "model.safetensors"
    if not safetensor_path.exists():
        safetensor_path = GF_MODEL_CACHE / "Geneformer-V2-316M" / "model.safetensors"

    model_weights = load_file(str(safetensor_path))
    # Find embedding weight key
    emb_keys = [k for k in model_weights.keys() if 'embed' in k.lower() and 'weight' in k.lower()]
    print(f"  Embedding keys: {emb_keys}", flush=True)

    # Use word embeddings (bert.embeddings.word_embeddings.weight)
    # Explicitly prefer word_embeddings key
    gf_emb_key = None
    for k in model_weights.keys():
        if 'word_embeddings.weight' in k:
            gf_emb_key = k
            break
    if gf_emb_key is None:
        for k in emb_keys:
            if 'word' in k.lower():
                gf_emb_key = k
                break
    if gf_emb_key is None and emb_keys:
        gf_emb_key = emb_keys[0]

    print(f"  Using key: {gf_emb_key}", flush=True)
    gf_word_emb = model_weights[gf_emb_key].numpy()  # [vocab_size, hidden_dim]
    print(f"  GF embedding shape: {gf_word_emb.shape}", flush=True)

    # Extract named gene embeddings
    gf_token_ids_arr = np.array(gf_token_ids)
    valid_mask = gf_token_ids_arr < gf_word_emb.shape[0]
    gf_named_genes_valid = [g for g, v in zip(gf_named_genes, valid_mask) if v]
    gf_token_ids_valid = gf_token_ids_arr[valid_mask]
    gf_gene_emb = gf_word_emb[gf_token_ids_valid]  # [N_valid, hidden_dim]
    print(f"  Valid named gene embeddings: {len(gf_named_genes_valid)}", flush=True)
    gf_emb_loaded = True

except Exception as e:
    print(f"  ERROR loading GF model: {e}", flush=True)
    gf_emb_loaded = False

if gf_emb_loaded:
    # Mean-center GF embeddings
    gf_gene_emb_centered = gf_gene_emb - gf_gene_emb.mean(axis=0, keepdims=True)

    # SVD of GF embeddings
    U_gf, S_gf, Vt_gf = np.linalg.svd(gf_gene_emb_centered, full_matrices=False)
    gf_proj = U_gf * S_gf[np.newaxis, :]  # [N_genes, n_components]
    print(f"  GF SVD singular values (top 5): {S_gf[:5]}", flush=True)

    K = len(gf_named_genes_valid) // 4  # Top/bottom K for poles

    # Build GF gene position map
    gf_gene_to_local = {g: i for i, g in enumerate(gf_named_genes_valid)}

    # Test co-pole structure in GF for STRING pairs
    n_string_gf = 0
    string_copole_counts = []
    string_scores_gf = []

    non_string_pairs = []
    all_gf_pairs = []
    pair_labels = []

    for g1, g2 in all_named_pairs:
        if g1 not in gf_gene_to_local or g2 not in gf_gene_to_local:
            continue
        i1, i2 = gf_gene_to_local[g1], gf_gene_to_local[g2]
        key = (min(g1, g2), max(g1, g2))
        is_string = key in string_pairs

        # Compute cosine similarity of GF static embeddings
        e1 = gf_gene_emb_centered[i1]
        e2 = gf_gene_emb_centered[i2]
        n1, n2 = np.linalg.norm(e1), np.linalg.norm(e2)
        cos = np.dot(e1, e2) / (n1 * n2 + 1e-10)

        all_gf_pairs.append(cos)
        pair_labels.append(1 if is_string else 0)

    all_gf_pairs = np.array(all_gf_pairs)
    pair_labels = np.array(pair_labels)

    gf_string_sims = all_gf_pairs[pair_labels == 1]
    gf_nonstring_sims = all_gf_pairs[pair_labels == 0]

    mw_gf = mannwhitneyu(gf_string_sims, gf_nonstring_sims, alternative='greater')
    z_gf = (gf_string_sims.mean() - gf_nonstring_sims.mean()) / (gf_nonstring_sims.std() + 1e-10)
    print(f"\n  GF static cosine for STRING pairs vs non-STRING:", flush=True)
    print(f"    N_string_valid={len(gf_string_sims)}, N_nonstring={len(gf_nonstring_sims)}", flush=True)
    print(f"    mean_string={gf_string_sims.mean():.4f}, mean_nonstring={gf_nonstring_sims.mean():.4f}", flush=True)
    print(f"    z={z_gf:.3f}, MW_p={mw_gf.pvalue:.3e}", flush=True)

    # Cross-model: scGPT SV2 at each layer vs GF SV1 projections
    print(f"\n  Cross-model correlation: scGPT SV2 vs GF SV1-4 projections", flush=True)
    cross_model_results = []

    # Get shared genes (in both scGPT and GF)
    shared_genes = [g for g in named_genes if g in gf_gene_to_local]
    scgpt_shared_idx = [gene_to_pos[g] for g in shared_genes]
    gf_shared_idx = [gf_gene_to_local[g] for g in shared_genes]

    for layer in range(N_LAYERS):
        # scGPT SVD at this layer
        layer_emb = named_emb[layer, :, :]  # [N_NAMED, 512]
        layer_centered = layer_emb - layer_emb.mean(axis=0, keepdims=True)
        U_sc, S_sc, Vt_sc = np.linalg.svd(layer_centered, full_matrices=False)
        scgpt_sv2 = (U_sc * S_sc[np.newaxis, :])[scgpt_shared_idx, 1]  # SV2 (index 1)

        # GF SV1
        gf_sv1 = gf_proj[gf_shared_idx, 0]

        rho, _ = spearmanr(np.abs(scgpt_sv2), np.abs(gf_sv1))
        rho_signed, _ = spearmanr(scgpt_sv2, gf_sv1)
        cross_model_results.append({
            'layer': layer,
            'spearman_abs': float(rho),
            'spearman_signed': float(rho_signed)
        })
        if layer in [0, 5, 11]:
            print(f"    L{layer:02d}: spearman_abs={rho:.3f}, spearman_signed={rho_signed:.3f}", flush=True)

    mean_spearman_abs = np.mean([r['spearman_abs'] for r in cross_model_results])
    print(f"\n  Mean spearman(abs) across layers: {mean_spearman_abs:.3f}", flush=True)

    h01_result = {
        "n_gf_genes": len(gf_named_genes_valid),
        "n_string_gf_valid": int(len(gf_string_sims)),
        "gf_string_mean_cosine": float(gf_string_sims.mean()),
        "gf_nonstring_mean_cosine": float(gf_nonstring_sims.mean()),
        "gf_z_score": float(z_gf),
        "gf_mw_pvalue": float(mw_gf.pvalue),
        "cross_model_mean_spearman_abs": float(mean_spearman_abs),
        "cross_model_by_layer": cross_model_results
    }
else:
    h01_result = {"error": "Geneformer model loading failed"}

# Save H01
h01_path = ITER_DIR / "h01_geneformer_cross_model.json"
with open(h01_path, 'w') as f:
    json.dump(h01_result, f, indent=2)
print(f"  H01 saved to {h01_path}", flush=True)

print("\n" + "="*60, flush=True)
print("=== H02: Multi-Axis Composite Co-Pole P@k ===", flush=True)
print("="*60, flush=True)

K = N_NAMED // 4  # 52 top/bottom per pole
N_SHUFFLE = 300
axes_to_use = [1, 2, 3]  # SV2, SV3, SV4 (0-indexed)

# Compute P@k by co-polarity count, averaged across layers
h02_layer_results = []

for layer in range(N_LAYERS):
    layer_emb = named_emb[layer, :, :]
    layer_centered = layer_emb - layer_emb.mean(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(layer_centered, full_matrices=False)
    projections = U * S[np.newaxis, :]  # [N_NAMED, n_components]

    # For each axis: get top/bottom K poles
    poles = {}
    for ax in axes_to_use:
        proj_ax = projections[:, ax]
        top_k_idx = set(np.argsort(proj_ax)[-K:])
        bot_k_idx = set(np.argsort(proj_ax)[:K])
        poles[ax] = (top_k_idx, bot_k_idx)

    # For each pair: compute composite co-polarity count
    pair_copole_count = {}
    for g1, g2 in all_named_pairs:
        i1, i2 = gene_to_pos[g1], gene_to_pos[g2]
        count = 0
        for ax in axes_to_use:
            top_k, bot_k = poles[ax]
            if (i1 in top_k and i2 in top_k) or (i1 in bot_k and i2 in bot_k):
                count += 1
        key = (min(g1, g2), max(g1, g2))
        pair_copole_count[key] = count

    # P@k by co-polarity count: for each count level, compute fraction that are STRING edges
    layer_result = {'layer': layer, 'by_count': {}}
    for ct in range(4):  # 0, 1, 2, 3
        pairs_at_count = [k for k, v in pair_copole_count.items() if v == ct]
        if len(pairs_at_count) == 0:
            continue
        n_string = sum(1 for p in pairs_at_count if p in string_pairs)
        frac = n_string / len(pairs_at_count)
        layer_result['by_count'][ct] = {
            'n_pairs': len(pairs_at_count),
            'n_string': n_string,
            'fraction_string': float(frac)
        }
    h02_layer_results.append(layer_result)

# Average across layers
count_summary = defaultdict(lambda: {'n_pairs_sum': 0, 'n_string_sum': 0, 'frac_list': []})
for lr in h02_layer_results:
    for ct, v in lr['by_count'].items():
        count_summary[ct]['n_pairs_sum'] += v['n_pairs']
        count_summary[ct]['n_string_sum'] += v['n_string']
        count_summary[ct]['frac_list'].append(v['fraction_string'])

random_baseline = len(string_pairs) / len(all_named_pairs)
print(f"  Random baseline: {random_baseline:.4f}", flush=True)
print(f"  Results by co-polarity count (averaged across layers):", flush=True)
for ct in sorted(count_summary.keys()):
    fracs = count_summary[ct]['frac_list']
    mean_frac = np.mean(fracs)
    enrichment = mean_frac / random_baseline
    print(f"    Count={ct}: N_pairs={count_summary[ct]['n_pairs_sum']//N_LAYERS:.0f}, "
          f"mean_frac={mean_frac:.4f}, enrichment={enrichment:.2f}x", flush=True)

# Summary: compute enrichment for count=3 vs count=0
frac_ct0 = np.mean(count_summary[0]['frac_list']) if 0 in count_summary else 0
frac_ct1 = np.mean(count_summary[1]['frac_list']) if 1 in count_summary else 0
frac_ct2 = np.mean(count_summary[2]['frac_list']) if 2 in count_summary else 0
frac_ct3 = np.mean(count_summary[3]['frac_list']) if 3 in count_summary else 0

print(f"\n  Enrichment vs baseline:", flush=True)
print(f"    0-axis: {frac_ct0/random_baseline:.2f}x", flush=True)
print(f"    1-axis: {frac_ct1/random_baseline:.2f}x", flush=True)
print(f"    2-axis: {frac_ct2/random_baseline:.2f}x", flush=True)
print(f"    3-axis: {frac_ct3/random_baseline:.2f}x", flush=True)

h02_result = {
    "random_baseline": float(random_baseline),
    "by_count": {
        ct: {
            "mean_frac": float(np.mean(count_summary[ct]['frac_list'])),
            "enrichment": float(np.mean(count_summary[ct]['frac_list']) / random_baseline),
            "n_layers_observed": len(count_summary[ct]['frac_list'])
        }
        for ct in sorted(count_summary.keys())
    },
    "layer_results": h02_layer_results
}

h02_path = ITER_DIR / "h02_multiaxis_composite_pak.json"
with open(h02_path, 'w') as f:
    json.dump(h02_result, f, indent=2)
print(f"\n  H02 saved to {h02_path}", flush=True)


print("\n" + "="*60, flush=True)
print("=== H03: TRRUST Repression vs Activation in scGPT Attention ===", flush=True)
print("="*60, flush=True)

# Load attention matrix
print("  Loading attention matrix ...", flush=True)
att = np.load(ATT_PATH)  # [8181, 8181]
print(f"  att shape: {att.shape}", flush=True)

# Map 209 named genes to attention matrix position
import anndata as ad
adata = ad.read_h5ad(H5AD_PATH)
att_gene_list = list(adata.var_names)
att_gene_to_idx = {g: i for i, g in enumerate(att_gene_list)}

# Compute symmetric attention
att_sym = (att + att.T) / 2

# Build null distribution (all named-gene pairs)
print("  Computing null distribution ...", flush=True)
null_atts = []
for g1, g2 in list(all_named_pairs)[:5000]:  # Sample 5000 for speed
    i1 = att_gene_to_idx.get(g1)
    i2 = att_gene_to_idx.get(g2)
    if i1 is not None and i2 is not None:
        null_atts.append(float(att_sym[i1, i2]))
null_atts = np.array(null_atts)
print(f"  Null size: {len(null_atts)}, mean: {null_atts.mean():.6f}", flush=True)

# STRING pairs attention
string_atts = []
for g1, g2 in string_pairs:
    i1 = att_gene_to_idx.get(g1)
    i2 = att_gene_to_idx.get(g2)
    if i1 is not None and i2 is not None:
        string_atts.append(float(att_sym[i1, i2]))
string_atts = np.array(string_atts)

mw_str = mannwhitneyu(string_atts, null_atts, alternative='greater')
print(f"\n  STRING: N={len(string_atts)}, mean={string_atts.mean():.6f}, MW_p={mw_str.pvalue:.3e}", flush=True)

# TRRUST activation attention
act_atts = []
for g1, g2 in trrust_activation:
    i1 = att_gene_to_idx.get(g1)
    i2 = att_gene_to_idx.get(g2)
    if i1 is not None and i2 is not None:
        act_atts.append(float(att_sym[i1, i2]))
act_atts = np.array(act_atts)

mw_act = mannwhitneyu(act_atts, null_atts, alternative='greater')
enrichment_act = act_atts.mean() / null_atts.mean()
print(f"  TRRUST Activation: N={len(act_atts)}, mean={act_atts.mean():.6f}, "
      f"enrichment={enrichment_act:.2f}x, MW_p={mw_act.pvalue:.3e}", flush=True)

# TRRUST repression attention
rep_atts = []
for g1, g2 in trrust_repression:
    i1 = att_gene_to_idx.get(g1)
    i2 = att_gene_to_idx.get(g2)
    if i1 is not None and i2 is not None:
        rep_atts.append(float(att_sym[i1, i2]))
rep_atts = np.array(rep_atts)

mw_rep = mannwhitneyu(rep_atts, null_atts, alternative='greater')
enrichment_rep = rep_atts.mean() / null_atts.mean()
print(f"  TRRUST Repression: N={len(rep_atts)}, mean={rep_atts.mean():.6f}, "
      f"enrichment={enrichment_rep:.2f}x, MW_p={mw_rep.pvalue:.3e}", flush=True)

# Activation vs Repression comparison
if len(act_atts) > 0 and len(rep_atts) > 0:
    mw_act_vs_rep = mannwhitneyu(act_atts, rep_atts, alternative='greater')
    print(f"\n  Activation vs Repression: MW_p={mw_act_vs_rep.pvalue:.3e}", flush=True)
else:
    mw_act_vs_rep = None
    print(f"  Activation vs Repression: insufficient data", flush=True)

# Pairs in both activation AND STRING
overlap_act_string = trrust_activation & string_pairs
overlap_rep_string = trrust_repression & string_pairs
print(f"\n  Activation ∩ STRING: {len(overlap_act_string)}", flush=True)
print(f"  Repression ∩ STRING: {len(overlap_rep_string)}", flush=True)

# Attention for activation-exclusive pairs (not in STRING)
act_exclusive = trrust_activation - string_pairs
act_excl_atts = []
for g1, g2 in act_exclusive:
    i1 = att_gene_to_idx.get(g1)
    i2 = att_gene_to_idx.get(g2)
    if i1 is not None and i2 is not None:
        act_excl_atts.append(float(att_sym[i1, i2]))
act_excl_atts = np.array(act_excl_atts)
if len(act_excl_atts) > 5:
    mw_excl = mannwhitneyu(act_excl_atts, null_atts, alternative='greater')
    print(f"\n  Activation-exclusive (not in STRING): N={len(act_excl_atts)}, "
          f"mean={act_excl_atts.mean():.6f}, enrichment={act_excl_atts.mean()/null_atts.mean():.2f}x, "
          f"MW_p={mw_excl.pvalue:.3e}", flush=True)
else:
    mw_excl = None
    print(f"  Activation-exclusive: too few pairs ({len(act_excl_atts)})", flush=True)

h03_result = {
    "null_size": int(len(null_atts)),
    "null_mean": float(null_atts.mean()),
    "string": {
        "n": int(len(string_atts)),
        "mean": float(string_atts.mean()),
        "enrichment": float(string_atts.mean() / null_atts.mean()),
        "mw_pvalue": float(mw_str.pvalue)
    },
    "trrust_activation": {
        "n": int(len(act_atts)),
        "mean": float(act_atts.mean()),
        "enrichment": float(enrichment_act),
        "mw_pvalue": float(mw_act.pvalue)
    },
    "trrust_repression": {
        "n": int(len(rep_atts)),
        "mean": float(rep_atts.mean() if len(rep_atts) > 0 else 0),
        "enrichment": float(enrichment_rep if len(rep_atts) > 0 else 0),
        "mw_pvalue": float(mw_rep.pvalue)
    },
    "activation_vs_repression_mw_pvalue": float(mw_act_vs_rep.pvalue) if mw_act_vs_rep else None,
    "activation_exclusive": {
        "n": int(len(act_excl_atts)),
        "mean": float(act_excl_atts.mean() if len(act_excl_atts) > 0 else 0),
        "enrichment": float(act_excl_atts.mean() / null_atts.mean() if len(act_excl_atts) > 0 else 0),
        "mw_pvalue": float(mw_excl.pvalue) if mw_excl else None
    }
}

h03_path = ITER_DIR / "h03_trrust_repression_vs_activation_attention.json"
with open(h03_path, 'w') as f:
    json.dump(h03_result, f, indent=2)
print(f"\n  H03 saved to {h03_path}", flush=True)

print("\n" + "="*60, flush=True)
print("=== SUMMARY ===", flush=True)
print("="*60, flush=True)
print(f"H01 (GF cross-model): GF static cosine z={h01_result.get('gf_z_score', 'N/A'):.3f}, "
      f"MW_p={h01_result.get('gf_mw_pvalue', 'N/A'):.3e}, "
      f"cross-model spearman_abs={h01_result.get('cross_model_mean_spearman_abs', 'N/A'):.3f}", flush=True)
print(f"H02 (multi-axis P@k): enrichment ct=3: {frac_ct3/random_baseline:.2f}x", flush=True)
print(f"H03 (attention): activation {enrichment_act:.2f}x (p={mw_act.pvalue:.2e}), "
      f"repression {enrichment_rep:.2f}x (p={mw_rep.pvalue:.2e})", flush=True)

print("\nAll done.", flush=True)
