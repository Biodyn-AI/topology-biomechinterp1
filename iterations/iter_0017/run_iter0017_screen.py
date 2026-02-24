"""
iter_0017 - Multi-Hypothesis Screen

H-A (SV4 confidence quintile gradient):
    Replicate SV2/SV3 quintile sweep on SV4 axis. Split 3092 STRING pairs (score>=0.4)
    into 5 quintiles. Compute mean SV4 co-pole z-score per quintile. Report Spearman rho.

H-B (Axis-dominant pair GO composition):
    Use axis-dominant pair labels from iter_0016 H02. For each axis-dominant gene set
    (SV2/SV3/SV4-dominant pairs), run Fisher's exact GO enrichment vs. full named gene
    background. Report top 5 GO terms per axis-dominant set + inter-subset Jaccard.

H-E (TRRUST signed regulation in SV3/SV4):
    Extend TRRUST activation/repression copole test to SV3 and SV4 axes.
    Compare layer-significant count to SV2 baseline (activation 12/12 layers).
"""

import numpy as np
import json
import sys
import pickle
from pathlib import Path
from collections import defaultdict
from scipy.stats import spearmanr, fisher_exact, mannwhitneyu
import csv

ITER_DIR = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0017"
)
ITER15_DIR = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0015"
)
ITER16_DIR = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0016"
)
EMB_PATH = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_38_geometric_residual_stream_interpretability"
    "/implementation/outputs/cycle1_main/layer_gene_embeddings.npy"
)
GENE2GO_PKL = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/single_cell_mechinterp/data/perturb/gene2go_all.pkl"
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

ITER_DIR.mkdir(parents=True, exist_ok=True)

# ─── Load embeddings ───────────────────────────────────────────────────────────
print("Loading embeddings ...", flush=True)
emb = np.load(EMB_PATH)  # [12, 4803, 512]
N_LAYERS, N_GENES, N_DIM = emb.shape
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
N_named = len(named_genes)
gene_to_local = {g: i for i, g in enumerate(named_genes)}
print(f"  Named genes: {N_named}", flush=True)

emb_named = emb[:, gene_indices, :]
emb_mc = emb_named - emb_named.mean(axis=1, keepdims=True)

# ─── Load STRING cache ─────────────────────────────────────────────────────────
print("Loading STRING cache ...", flush=True)
with open(STRING_API_CACHE) as f:
    string_cache = json.load(f)

raw_pairs = string_cache.get("pairs", string_cache) if isinstance(string_cache, dict) else string_cache
string_pairs_04 = []
for entry in raw_pairs:
    ga = entry.get("g1") or entry.get("gene_a")
    gb = entry.get("g2") or entry.get("gene_b")
    sc = entry.get("score", 0.0)
    if ga and gb and ga in gene_to_local and gb in gene_to_local:
        string_pairs_04.append((ga, gb, sc))
print(f"  STRING pairs (score>=0.4): {len(string_pairs_04)}", flush=True)

# ─── SVD helper ───────────────────────────────────────────────────────────────
def get_sv_projections(layer, sv_idx, K=52):
    """Return (sv_proj [N_named], top_genes_set, bot_genes_set) for a given layer/sv_idx."""
    U, S, Vt = np.linalg.svd(emb_mc[layer], full_matrices=False)
    # emb_mc[layer]: [N_named, N_DIM] → U: [N_named, min], Vt: [min, N_DIM]
    sv_proj = U[:, sv_idx] * S[sv_idx]
    sorted_idx = np.argsort(sv_proj)
    top_genes = set(named_genes[i] for i in sorted_idx[-K:])
    bot_genes = set(named_genes[i] for i in sorted_idx[:K])
    return sv_proj, top_genes, bot_genes

# ─── H-A: SV4 confidence quintile gradient ────────────────────────────────────
print("\n=== H-A: SV4 confidence quintile gradient ===", flush=True)

SV4_IDX = 3  # 0-indexed: SV1=0, SV2=1, SV3=2, SV4=3
K = 52
N_QUINTILES = 5
FOCUS_LAYERS = [7, 8, 11]  # 0-indexed layers for quick test
ALL_LAYERS = list(range(12))

# Compute co-pole z-scores for all STRING pairs at each layer for SV4
scores = np.array([sc for _, _, sc in string_pairs_04])
quintile_edges = np.percentile(scores, np.linspace(0, 100, N_QUINTILES + 1))
quintile_labels = np.digitize(scores, quintile_edges[1:-1])  # 0..4

# For each pair, compute SV4 co-pole z at focus layers
# z = (mean_copole_indicator - baseline_rate) / se_baseline
ha_layer_results = {}
for layer in ALL_LAYERS:
    sv_proj, top_genes, bot_genes = get_sv_projections(layer, SV4_IDX, K=K)
    pole_genes = top_genes | bot_genes
    copole_flags = []
    for ga, gb, sc in string_pairs_04:
        in_top = (ga in top_genes and gb in top_genes)
        in_bot = (ga in bot_genes and gb in bot_genes)
        copole_flags.append(1 if (in_top or in_bot) else 0)
    copole_flags = np.array(copole_flags)

    # Compute copole rate per quintile
    q_rates = []
    for q in range(N_QUINTILES):
        mask = quintile_labels == q
        rate = copole_flags[mask].mean() if mask.sum() > 0 else 0.0
        q_rates.append(float(rate))

    # Baseline rate
    baseline = copole_flags.mean()
    baseline_std = np.sqrt(baseline * (1 - baseline) / len(copole_flags))

    # z-scores per quintile
    q_z = [(r - baseline) / (baseline_std + 1e-9) for r in q_rates]
    ha_layer_results[layer] = {
        "quintile_rates": q_rates,
        "quintile_z": q_z,
        "baseline_rate": float(baseline)
    }

# Compute Spearman rho across quintiles using mean z across all layers
quintile_midpoints = []
for q in range(N_QUINTILES):
    mask = quintile_labels == q
    midpoint = float(scores[mask].mean()) if mask.sum() > 0 else 0.0
    quintile_midpoints.append(midpoint)

mean_z_per_quintile = [
    float(np.mean([ha_layer_results[l]["quintile_z"][q] for l in ALL_LAYERS]))
    for q in range(N_QUINTILES)
]

rho_ha, p_ha = spearmanr(quintile_midpoints, mean_z_per_quintile)
q5_mean_z_ha = mean_z_per_quintile[-1]

print(f"  SV4 Spearman rho={rho_ha:.3f}, p={p_ha:.4f}", flush=True)
print(f"  Quintile mean z (all layers): {[round(z,3) for z in mean_z_per_quintile]}", flush=True)
print(f"  Q5 mean_z={q5_mean_z_ha:.3f}", flush=True)

ha_out = {
    "sv_idx": SV4_IDX,
    "K": K,
    "n_quintiles": N_QUINTILES,
    "n_pairs": len(string_pairs_04),
    "quintile_midpoints": quintile_midpoints,
    "mean_z_per_quintile_all_layers": mean_z_per_quintile,
    "spearman_rho": round(rho_ha, 4),
    "spearman_p": p_ha,
    "Q5_mean_z": round(q5_mean_z_ha, 3),
    "layer_results": ha_layer_results
}
with open(ITER_DIR / "h_a_sv4_confidence_gradient.json", "w") as f:
    json.dump(ha_out, f, indent=2)
print("  Saved h_a_sv4_confidence_gradient.json", flush=True)

# ─── Load GO annotations ──────────────────────────────────────────────────────
print("\nLoading GO annotations ...", flush=True)
with open(GENE2GO_PKL, "rb") as f:
    gene2go_raw = pickle.load(f)

# Build GO → gene set (among named genes) and gene → GO set
go_to_genes = defaultdict(set)
gene_to_go = defaultdict(set)
for gene, go_terms in gene2go_raw.items():
    if gene in gene_to_local:
        for go in go_terms:
            go_to_genes[go].add(gene)
            gene_to_go[gene].add(go)

# Filter GO terms: size 3–50 among named genes
valid_go_terms = [go for go, genes in go_to_genes.items() if 3 <= len(genes) <= 50]
print(f"  Valid GO terms (3-50 named genes): {len(valid_go_terms)}", flush=True)

# ─── H-B: Axis-dominant pair GO composition ────────────────────────────────────
print("\n=== H-B: Axis-dominant pair GO composition ===", flush=True)

# Load axis-dominant pair info from iter_0016
with open(ITER16_DIR / "h02_axis_independence.json") as f:
    h02_data = json.load(f)

# Recompute axis-dominant pair labels (using layers 7,8,11 as representative)
# A pair is "SV-k dominant" if it has more co-pole layers on axis k than others
# Method: use layers 7,8,11 to count co-pole hits per axis
FOCUS_LAYERS_B = [7, 8, 11]
axis_sv_indices = {"SV2": 1, "SV3": 2, "SV4": 3}

# Precompute poles for each axis at focus layers
axis_poles = {}
for ax_name, sv_idx in axis_sv_indices.items():
    axis_poles[ax_name] = {}
    for layer in FOCUS_LAYERS_B:
        _, top_g, bot_g = get_sv_projections(layer, sv_idx, K=K)
        axis_poles[ax_name][layer] = (top_g, bot_g)

# For each STRING pair, count co-pole hits per axis
pair_axis_counts = []
for ga, gb, sc in string_pairs_04:
    counts = {}
    for ax_name in axis_sv_indices:
        cnt = 0
        for layer in FOCUS_LAYERS_B:
            top_g, bot_g = axis_poles[ax_name][layer]
            if (ga in top_g and gb in top_g) or (ga in bot_g and gb in bot_g):
                cnt += 1
        counts[ax_name] = cnt
    pair_axis_counts.append(counts)

# Assign dominant axis (must have >= 1 co-pole hit, more than any other axis)
axis_dominant_gene_sets = {"SV2": set(), "SV3": set(), "SV4": set(), "none": set()}
axis_dominant_pair_counts = {"SV2": 0, "SV3": 0, "SV4": 0, "none": 0}
for i, (ga, gb, sc) in enumerate(string_pairs_04):
    counts = pair_axis_counts[i]
    max_count = max(counts.values())
    if max_count == 0:
        dominant = "none"
    else:
        # All axes with max count
        top_axes = [ax for ax, c in counts.items() if c == max_count]
        dominant = top_axes[0] if len(top_axes) == 1 else "none"

    axis_dominant_pair_counts[dominant] += 1
    if dominant != "none":
        axis_dominant_gene_sets[dominant].add(ga)
        axis_dominant_gene_sets[dominant].add(gb)

print(f"  Axis-dominant pair counts: {axis_dominant_pair_counts}", flush=True)
for ax_name in ["SV2", "SV3", "SV4"]:
    print(f"  {ax_name}-dominant gene set size: {len(axis_dominant_gene_sets[ax_name])}", flush=True)

# Fisher's exact GO enrichment per axis-dominant gene set vs. background
N_bg = N_named
def fisher_go_enrichment(gene_set, go_to_genes, valid_go_terms, N_bg):
    """Run Fisher's exact for each GO term; return sorted results."""
    results = []
    n_set = len(gene_set)
    for go in valid_go_terms:
        go_genes = go_to_genes[go]
        N_go = len(go_genes)
        k = len(gene_set & go_genes)
        if k == 0:
            continue
        # 2x2 contingency: [k, n_set-k; N_go-k, N_bg-n_set-N_go+k]
        a = k
        b = n_set - k
        c = N_go - k
        d = N_bg - n_set - N_go + k
        if d < 0:
            continue
        _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        results.append({"go": go, "p": p, "k": k, "N_go": N_go, "n_set": n_set})
    results.sort(key=lambda x: x["p"])
    return results

hb_axis_results = {}
for ax_name in ["SV2", "SV3", "SV4"]:
    gene_set = axis_dominant_gene_sets[ax_name]
    if len(gene_set) < 3:
        hb_axis_results[ax_name] = {"error": "too few genes"}
        continue
    enrich = fisher_go_enrichment(gene_set, go_to_genes, valid_go_terms, N_bg)
    top5 = enrich[:5]
    hb_axis_results[ax_name] = {
        "n_genes": len(gene_set),
        "n_go_tested": len([x for x in enrich]),
        "top5": top5
    }
    print(f"  {ax_name}: top GO = {top5[0]['go'] if top5 else 'none'} (p={top5[0]['p']:.4f} k={top5[0]['k']}/{top5[0]['N_go']})" if top5 else f"  {ax_name}: no enrichment", flush=True)

# Inter-subset GO Jaccard of top-20 GO hits
def go_jaccard(results_a, results_b, top_k=20):
    set_a = set(r["go"] for r in results_a[:top_k])
    set_b = set(r["go"] for r in results_b[:top_k])
    inter = len(set_a & set_b)
    union = len(set_a | set_b)
    return inter / union if union > 0 else 0.0

# Need to store full enrichment for Jaccard
hb_full_enrich = {}
for ax_name in ["SV2", "SV3", "SV4"]:
    gene_set = axis_dominant_gene_sets[ax_name]
    hb_full_enrich[ax_name] = fisher_go_enrichment(gene_set, go_to_genes, valid_go_terms, N_bg) if len(gene_set) >= 3 else []

jaccard_sv2_sv3 = go_jaccard(hb_full_enrich["SV2"], hb_full_enrich["SV3"])
jaccard_sv2_sv4 = go_jaccard(hb_full_enrich["SV2"], hb_full_enrich["SV4"])
jaccard_sv3_sv4 = go_jaccard(hb_full_enrich["SV3"], hb_full_enrich["SV4"])
print(f"  GO Jaccard SV2–SV3={jaccard_sv2_sv3:.3f}, SV2–SV4={jaccard_sv2_sv4:.3f}, SV3–SV4={jaccard_sv3_sv4:.3f}", flush=True)

hb_out = {
    "axis_dominant_pair_counts": axis_dominant_pair_counts,
    "axis_gene_set_sizes": {ax: len(axis_dominant_gene_sets[ax]) for ax in ["SV2","SV3","SV4"]},
    "axis_go_enrichment": hb_axis_results,
    "go_top20_jaccard": {
        "SV2_SV3": round(jaccard_sv2_sv3, 4),
        "SV2_SV4": round(jaccard_sv2_sv4, 4),
        "SV3_SV4": round(jaccard_sv3_sv4, 4)
    }
}
with open(ITER_DIR / "h_b_axis_go_composition.json", "w") as f:
    json.dump(hb_out, f, indent=2)
print("  Saved h_b_axis_go_composition.json", flush=True)

# ─── H-E: TRRUST signed regulation in SV3/SV4 ────────────────────────────────
print("\n=== H-E: TRRUST signed regulation in SV3/SV4 ===", flush=True)

# Load TRRUST pairs (activation and repression)
act_pairs = []  # (tf, target) for activation
rep_pairs = []  # (tf, target) for repression
with open(TRRUST_PATH) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        tf, target, reg_type = parts[0], parts[1], parts[2]
        if tf in gene_to_local and target in gene_to_local:
            if "Activation" in reg_type:
                act_pairs.append((tf, target))
            if "Repression" in reg_type:
                rep_pairs.append((tf, target))

act_pairs = list(set(act_pairs))
rep_pairs = list(set(rep_pairs))
print(f"  TRRUST activation pairs (both in named): {len(act_pairs)}", flush=True)
print(f"  TRRUST repression pairs (both in named): {len(rep_pairs)}", flush=True)

# For each axis (SV3, SV4) at each layer:
# - compute SV projections
# - for activation: test if TF and target tend to co-pole (same end)
# - for repression: test if TF and target tend to anti-pole (opposite ends)
# Use Mann-Whitney U: compare |SV_proj(tf) - SV_proj(target)| for act vs rep
# (repression pairs should have larger separation if anti-pole)

he_results = {}
for ax_name, sv_idx in [("SV3", 2), ("SV4", 3)]:
    layer_results = {}
    for layer in ALL_LAYERS:
        sv_proj, top_genes, bot_genes = get_sv_projections(layer, sv_idx, K=K)
        gene_proj = {g: float(sv_proj[gene_to_local[g]]) for g in named_genes}

        # Signed projection differences (TF proj - target proj)
        act_diffs = [gene_proj[tf] - gene_proj[tgt] for tf, tgt in act_pairs]
        rep_diffs = [gene_proj[tf] - gene_proj[tgt] for tf, tgt in rep_pairs]

        # Co-pole indicator: sign(tf_proj) == sign(target_proj)
        act_copole = [1 if np.sign(gene_proj[tf]) == np.sign(gene_proj[tgt]) else 0 for tf, tgt in act_pairs]
        rep_copole = [1 if np.sign(gene_proj[tf]) == np.sign(gene_proj[tgt]) else 0 for tf, tgt in rep_pairs]

        act_copole_rate = np.mean(act_copole)
        rep_copole_rate = np.mean(rep_copole)

        # Test: abs separation
        act_abs = [abs(d) for d in act_diffs]
        rep_abs = [abs(d) for d in rep_diffs]
        stat, p_mwu = mannwhitneyu(act_abs, rep_abs, alternative="less")

        # Expected for regulation: activation = co-pole (smaller abs diff), repression = anti-pole (larger abs diff)
        layer_results[layer] = {
            "act_copole_rate": round(act_copole_rate, 4),
            "rep_copole_rate": round(rep_copole_rate, 4),
            "act_minus_rep_copole": round(act_copole_rate - rep_copole_rate, 4),
            "mwu_stat": float(stat),
            "mwu_p": float(p_mwu),
            "significant": p_mwu < 0.05
        }

    n_sig = sum(1 for r in layer_results.values() if r["significant"])
    mean_act_cop = np.mean([r["act_copole_rate"] for r in layer_results.values()])
    mean_rep_cop = np.mean([r["rep_copole_rate"] for r in layer_results.values()])

    print(f"  {ax_name}: n_layers_sig={n_sig}/12, mean_act_copole={mean_act_cop:.3f}, mean_rep_copole={mean_rep_cop:.3f}", flush=True)
    he_results[ax_name] = {
        "sv_idx": sv_idx,
        "n_act_pairs": len(act_pairs),
        "n_rep_pairs": len(rep_pairs),
        "n_layers_significant": n_sig,
        "mean_act_copole_rate": round(mean_act_cop, 4),
        "mean_rep_copole_rate": round(mean_rep_cop, 4),
        "mean_act_minus_rep": round(mean_act_cop - mean_rep_cop, 4),
        "layer_results": layer_results
    }

with open(ITER_DIR / "h_e_trrust_sv3_sv4_signed.json", "w") as f:
    json.dump(he_results, f, indent=2)
print("  Saved h_e_trrust_sv3_sv4_signed.json", flush=True)

# ─── Combined results ─────────────────────────────────────────────────────────
results = {
    "iter": "iter_0017",
    "H_A_sv4_confidence_gradient": {
        "spearman_rho": round(rho_ha, 4),
        "spearman_p": p_ha,
        "quintile_mean_z": [round(z, 3) for z in mean_z_per_quintile],
        "Q5_mean_z": round(q5_mean_z_ha, 3)
    },
    "H_B_axis_go_composition": {
        "axis_dominant_pair_counts": axis_dominant_pair_counts,
        "axis_gene_set_sizes": {ax: len(axis_dominant_gene_sets[ax]) for ax in ["SV2","SV3","SV4"]},
        "go_jaccard_top20": {
            "SV2_SV3": round(jaccard_sv2_sv3, 4),
            "SV2_SV4": round(jaccard_sv2_sv4, 4),
            "SV3_SV4": round(jaccard_sv3_sv4, 4)
        },
        "top_go_per_axis": {
            ax: (hb_full_enrich[ax][0]["go"] if hb_full_enrich[ax] else "none")
            for ax in ["SV2", "SV3", "SV4"]
        }
    },
    "H_E_trrust_sv3_sv4_signed": {
        ax: {
            "n_layers_sig": he_results[ax]["n_layers_significant"],
            "mean_act_copole": he_results[ax]["mean_act_copole_rate"],
            "mean_rep_copole": he_results[ax]["mean_rep_copole_rate"],
            "mean_delta": he_results[ax]["mean_act_minus_rep"]
        }
        for ax in ["SV3", "SV4"]
    }
}

with open(ITER_DIR / "iter0017_results.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nSaved iter0017_results.json", flush=True)
print("\n=== DONE ===", flush=True)
print(json.dumps(results, indent=2), flush=True)
