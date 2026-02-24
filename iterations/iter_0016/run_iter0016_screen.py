"""
iter_0016 - Multi-Hypothesis Screen

H01 (SV3 STRING confidence gradient): Replicate the SV2 rho=1.000 finding on SV3 axis.
    Split 3092 STRING pairs into 5 quintiles by score. Compute mean SV3 co-pole z-score
    per quintile. Report Spearman rho of quintile midpoint vs mean z.

H02 (SV2/SV3/SV4 axis independence): For each axis, extract top-K=52 and bottom-K=52 genes.
    Compute pairwise Jaccard overlap of pole gene sets across axes. Then cluster STRING pairs
    by which axis they co-localize on (sig co-pole layers). Test if different axes capture
    distinct sets of STRING partners (orthogonal biology test).

H03 (SV4 GO biology profile): For each of 12 layers, extract SV4 top-K=52 and bottom-K=52
    genes. For each GO term (size 3-50 among 209 named genes), compute enrichment in each pole
    via Fisher's exact test. Report top enriched GO terms per pole, compare to SV2/SV3 profiles.
"""

import numpy as np
import json
import sys
import pickle
from pathlib import Path
from collections import defaultdict
from scipy.stats import spearmanr, fisher_exact
from scipy.special import comb
import csv

ITER_DIR = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0016"
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

# Slice embeddings to named genes
# emb_named: [12, N_named, 512]
emb_named = emb[:, gene_indices, :]

# Mean-center per layer
emb_mc = emb_named - emb_named.mean(axis=1, keepdims=True)

# ─── Load STRING cache (score >= 0.4 pairs among named genes) ──────────────────
print("Loading STRING cache ...", flush=True)
with open(STRING_API_CACHE) as f:
    string_cache = json.load(f)

# string_cache: dict with "pairs" list of {"g1":..., "g2":..., "score":...}
raw_pairs = string_cache.get("pairs", string_cache) if isinstance(string_cache, dict) else string_cache
string_pairs_04 = []
for entry in raw_pairs:
    ga = entry.get("g1") or entry.get("gene_a")
    gb = entry.get("g2") or entry.get("gene_b")
    sc = entry.get("score", 0.0)
    if ga and gb and ga in gene_to_local and gb in gene_to_local:
        string_pairs_04.append((ga, gb, sc))
print(f"  STRING pairs (score>=0.4, in named set): {len(string_pairs_04)}", flush=True)

# Sort for determinism
string_pairs_04.sort(key=lambda x: (x[0], x[1]))

# ─── SVD helper ───────────────────────────────────────────────────────────────
def compute_svd_projections(layer_idx, n_components=6):
    """Return projections [N_named, n_components] after SVD of mean-centered emb at layer_idx."""
    X = emb_mc[layer_idx]  # [N_named, 512]
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    # U: [N_named, min(N_named,512)], projections = U * S
    projs = U[:, :n_components] * S[:n_components]
    return projs  # [N_named, n_components]

# Pre-compute SVD projections for all layers
print("Pre-computing SVD projections for all layers...", flush=True)
all_projs = []  # list of [N_named, 6] for each layer
for l in range(N_LAYERS):
    all_projs.append(compute_svd_projections(l, n_components=6))
print("  SVD done.", flush=True)

K = 52  # top/bottom pole size
N_NULL = 300

def copole_zscore(proj_vals, pairs, k=K, n_null=N_NULL, rng_seed=42):
    """
    proj_vals: 1D array [N_named] of projections for one SV at one layer.
    pairs: list of (ga, gb) tuples (local indices or gene names with gene_to_local).
    Returns z-score of co-pole rate vs shuffle null.
    """
    rng = np.random.default_rng(rng_seed)
    sorted_idx = np.argsort(proj_vals)
    top_set = set(sorted_idx[-k:].tolist())
    bot_set = set(sorted_idx[:k].tolist())

    def copole_rate(p_list):
        cnt = 0
        for ga, gb in p_list:
            ia, ib = gene_to_local[ga], gene_to_local[gb]
            if (ia in top_set and ib in top_set) or (ia in bot_set and ib in bot_set):
                cnt += 1
        return cnt / len(p_list) if p_list else 0.0

    obs = copole_rate([(ga, gb) for ga, gb, _ in pairs])

    # Null: shuffle gene labels
    null_rates = []
    genes_arr = np.array(named_genes)
    for _ in range(n_null):
        perm = rng.permutation(N_named)
        shuffled = {g: perm[i] for i, g in enumerate(named_genes)}
        cnt = 0
        for ga, gb, _ in pairs:
            ia, ib = shuffled[ga], shuffled[gb]
            if (ia in top_set and ib in top_set) or (ia in bot_set and ib in bot_set):
                cnt += 1
        null_rates.append(cnt / len(pairs))

    null_mean = np.mean(null_rates)
    null_std = np.std(null_rates)
    z = (obs - null_mean) / null_std if null_std > 0 else 0.0
    emp_p = np.mean([r >= obs for r in null_rates])
    return z, emp_p, obs, null_mean, null_std


# ═══════════════════════════════════════════════════════════════════════════════
# H01: SV3 STRING confidence gradient
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== H01: SV3 confidence gradient ===", flush=True)

# Sort pairs by score, split into 5 quintiles
scores_arr = np.array([sc for _, _, sc in string_pairs_04])
quintile_edges = np.percentile(scores_arr, [0, 20, 40, 60, 80, 100])
print(f"  Quintile edges: {[round(x, 3) for x in quintile_edges]}", flush=True)

# Assign quintile to each pair
quintile_labels = np.digitize(scores_arr, quintile_edges[1:-1])  # 0-4

SV_IDX = 2  # SV3 (0-indexed)

# For each quintile, compute mean z across 12 layers
h01_quintile_data = []
for q in range(5):
    q_mask = quintile_labels == q
    q_pairs = [string_pairs_04[i] for i in range(len(string_pairs_04)) if q_mask[i]]
    q_scores = scores_arr[q_mask]
    q_mid = (quintile_edges[q] + quintile_edges[q + 1]) / 2

    layer_z = []
    for l in range(N_LAYERS):
        proj = all_projs[l][:, SV_IDX]
        z, ep, obs, nm, ns = copole_zscore(proj, q_pairs, k=K, n_null=N_NULL, rng_seed=42 + l)
        layer_z.append(z)

    mean_z = float(np.mean(layer_z))
    print(f"  Q{q+1} n={len(q_pairs)} score=[{quintile_edges[q]:.3f},{quintile_edges[q+1]:.3f}] mid={q_mid:.3f} mean_z={mean_z:.3f}", flush=True)
    h01_quintile_data.append({
        "quintile": q + 1,
        "n_pairs": len(q_pairs),
        "score_lo": float(quintile_edges[q]),
        "score_hi": float(quintile_edges[q + 1]),
        "score_mid": float(q_mid),
        "layer_z": [round(z, 4) for z in layer_z],
        "mean_z": round(mean_z, 4),
    })

quintile_mids = [d["score_mid"] for d in h01_quintile_data]
quintile_zs = [d["mean_z"] for d in h01_quintile_data]
rho_sv3, p_sv3 = spearmanr(quintile_mids, quintile_zs)
print(f"  SV3 gradient rho={rho_sv3:.4f}, p={p_sv3:.3e}", flush=True)

h01_result = {
    "sv_index": SV_IDX,
    "sv_name": "SV3",
    "quintile_data": h01_quintile_data,
    "spearman_rho": round(rho_sv3, 6),
    "spearman_p": float(p_sv3),
    "quintile_mean_z": [round(z, 3) for z in quintile_zs],
}
with open(ITER_DIR / "h01_sv3_confidence_gradient.json", "w") as f:
    json.dump(h01_result, f, indent=2)
print("  Saved h01_sv3_confidence_gradient.json", flush=True)


# ═══════════════════════════════════════════════════════════════════════════════
# H02: SV2/SV3/SV4 axis independence test
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== H02: Axis independence test (SV2/SV3/SV4) ===", flush=True)

# For each axis, compute pole gene sets at each layer
# Then compute average Jaccard overlap across layers between axes
# Also: for each STRING pair, determine on which axes it is co-pole (sig z > 1.96)
# across 12 layers → axis assignment

SV_AXES = [1, 2, 3]  # SV2, SV3, SV4 (0-indexed)
SV_NAMES = ["SV2", "SV3", "SV4"]

# Step 1: Pairwise pole Jaccard overlap across axes (layer-averaged)
axis_top_sets = {sv: [] for sv in SV_AXES}  # list over layers of frozenset
axis_bot_sets = {sv: [] for sv in SV_AXES}

for sv in SV_AXES:
    for l in range(N_LAYERS):
        proj = all_projs[l][:, sv]
        sorted_idx = np.argsort(proj)
        axis_top_sets[sv].append(frozenset(sorted_idx[-K:].tolist()))
        axis_bot_sets[sv].append(frozenset(sorted_idx[:K].tolist()))

def jaccard(a, b):
    return len(a & b) / len(a | b) if (a | b) else 0.0

h02_jaccard = {}
for i, sv_a in enumerate(SV_AXES):
    for j, sv_b in enumerate(SV_AXES):
        if j <= i:
            continue
        # Compute layer-averaged Jaccard for top and bottom poles separately
        top_jaccards = [jaccard(axis_top_sets[sv_a][l], axis_top_sets[sv_b][l]) for l in range(N_LAYERS)]
        bot_jaccards = [jaccard(axis_bot_sets[sv_a][l], axis_bot_sets[sv_b][l]) for l in range(N_LAYERS)]
        cross_jaccards_ab = [jaccard(axis_top_sets[sv_a][l], axis_bot_sets[sv_b][l]) for l in range(N_LAYERS)]
        cross_jaccards_ba = [jaccard(axis_top_sets[sv_b][l], axis_bot_sets[sv_a][l]) for l in range(N_LAYERS)]
        key = f"{SV_NAMES[i]}_vs_{SV_NAMES[j]}"
        h02_jaccard[key] = {
            "top_top_mean_jaccard": round(float(np.mean(top_jaccards)), 4),
            "bot_bot_mean_jaccard": round(float(np.mean(bot_jaccards)), 4),
            "top_bot_mean_jaccard": round(float(np.mean(cross_jaccards_ab)), 4),
            "bot_top_mean_jaccard": round(float(np.mean(cross_jaccards_ba)), 4),
        }
        print(f"  {key}: top-top={h02_jaccard[key]['top_top_mean_jaccard']:.4f}, bot-bot={h02_jaccard[key]['bot_bot_mean_jaccard']:.4f}", flush=True)

# Step 2: For each STRING pair, compute per-axis mean co-pole z across 12 layers
# Then assign pair to axis with highest mean_z (if > 1.96)
print("  Computing per-pair per-axis copole z...", flush=True)

# For efficiency, compute per-layer co-pole indicator for each pair x axis
# Then average over layers
pair_axis_z = np.zeros((len(string_pairs_04), len(SV_AXES)))  # [N_pairs, 3]

for ax_i, sv in enumerate(SV_AXES):
    layer_z_per_pair = []
    for l in range(N_LAYERS):
        proj = all_projs[l][:, sv]
        sorted_idx = np.argsort(proj)
        top_set = set(sorted_idx[-K:].tolist())
        bot_set = set(sorted_idx[:K].tolist())
        z_l = []
        for ga, gb, sc in string_pairs_04:
            ia, ib = gene_to_local[ga], gene_to_local[gb]
            copole = int((ia in top_set and ib in top_set) or (ia in bot_set and ib in bot_set))
            z_l.append(copole)
        layer_z_per_pair.append(z_l)
    # Average over layers (as rate, not z-score per pair)
    layer_z_arr = np.array(layer_z_per_pair).T  # [N_pairs, 12]
    pair_axis_z[:, ax_i] = layer_z_arr.mean(axis=1)

# Pair assignment: argmax axis (only if max > 0)
pair_assignments = np.argmax(pair_axis_z, axis=1)  # 0,1,2 = SV2,SV3,SV4
# How many pairs dominated by each axis?
axis_counts = [int((pair_assignments == ax_i).sum()) for ax_i in range(3)]
print(f"  Pair dominant axis counts: {dict(zip(SV_NAMES, axis_counts))}", flush=True)

# Also compute correlation of axis copole rates across pairs (independence test)
corr_12 = float(np.corrcoef(pair_axis_z[:, 0], pair_axis_z[:, 1])[0, 1])
corr_13 = float(np.corrcoef(pair_axis_z[:, 0], pair_axis_z[:, 2])[0, 1])
corr_23 = float(np.corrcoef(pair_axis_z[:, 1], pair_axis_z[:, 2])[0, 1])
print(f"  Inter-axis correlations: SV2-SV3={corr_12:.4f}, SV2-SV4={corr_13:.4f}, SV3-SV4={corr_23:.4f}", flush=True)

h02_result = {
    "axes": SV_NAMES,
    "K": K,
    "n_pairs": len(string_pairs_04),
    "jaccard_overlaps": h02_jaccard,
    "axis_dominant_pair_counts": dict(zip(SV_NAMES, axis_counts)),
    "inter_axis_copole_rate_correlations": {
        "SV2_SV3": round(corr_12, 4),
        "SV2_SV4": round(corr_13, 4),
        "SV3_SV4": round(corr_23, 4),
    },
}
with open(ITER_DIR / "h02_axis_independence.json", "w") as f:
    json.dump(h02_result, f, indent=2)
print("  Saved h02_axis_independence.json", flush=True)


# ═══════════════════════════════════════════════════════════════════════════════
# H03: SV4 GO biology profile
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== H03: SV4 GO biology profile ===", flush=True)

# Load gene2go
with open(GENE2GO_PKL, "rb") as f:
    gene2go_all = pickle.load(f)

print(f"  gene2go_all type: {type(gene2go_all)}", flush=True)
# Determine structure
if isinstance(gene2go_all, dict):
    sample_key = next(iter(gene2go_all))
    print(f"  Sample key: {sample_key!r}, value type: {type(gene2go_all[sample_key])}", flush=True)
    if isinstance(gene2go_all[sample_key], (set, list, frozenset)):
        gene2go = {g: set(terms) for g, terms in gene2go_all.items() if g in gene_to_local}
        print(f"  Using direct gene→GO mapping.", flush=True)
    else:
        # Might be nested or different format
        print(f"  Value sample: {str(gene2go_all[sample_key])[:200]}", flush=True)
        gene2go = {}

# Build GO→gene mapping restricted to named genes
go2genes = defaultdict(set)
for g, terms in gene2go.items():
    for t in terms:
        go2genes[t].add(g)

# Filter GO terms with 3-50 named genes
valid_go = {t: gs for t, gs in go2genes.items() if 3 <= len(gs) <= 50}
print(f"  Valid GO terms (3-50 named genes): {len(valid_go)}", flush=True)

# For SV4 at focus layers (7, 8, 11), run GO enrichment on each pole
SV4_IDX = 3  # SV4 (0-indexed, so index 3)
FOCUS_LAYERS = [7, 8, 11]

h03_results = {}
for l in FOCUS_LAYERS:
    proj = all_projs[l][:, SV4_IDX]
    sorted_idx = np.argsort(proj)
    top_genes = set(named_genes[i] for i in sorted_idx[-K:])
    bot_genes = set(named_genes[i] for i in sorted_idx[:K])

    enrichments = []
    for go_term, go_genes in valid_go.items():
        for pole_name, pole_genes in [("top", top_genes), ("bottom", bot_genes)]:
            overlap = len(pole_genes & go_genes)
            if overlap == 0:
                continue
            # Fisher's exact test
            a = overlap  # in pole AND in GO
            b = len(pole_genes) - overlap  # in pole, not in GO
            c = len(go_genes) - overlap  # not in pole, in GO
            d = N_named - len(pole_genes) - c  # neither
            if d < 0:
                continue
            _, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
            enrichments.append({
                "go_term": go_term,
                "pole": pole_name,
                "n_go": len(go_genes),
                "overlap": overlap,
                "pval": float(pval),
                "fold_enrichment": round((overlap / len(pole_genes)) / (len(go_genes) / N_named), 3) if len(go_genes) > 0 else 0,
                "genes": sorted(pole_genes & go_genes),
            })

    # Sort by pval, take top 20
    enrichments.sort(key=lambda x: x["pval"])
    top20 = enrichments[:20]
    print(f"  Layer {l}: top GO hit = {top20[0]['go_term']} (pole={top20[0]['pole']}, p={top20[0]['pval']:.3e}, genes={top20[0]['genes'][:5]})", flush=True)

    h03_results[f"layer_{l}"] = {
        "layer": l,
        "sv": "SV4",
        "top_K_genes": sorted(top_genes),
        "bottom_K_genes": sorted(bot_genes),
        "top_go_enrichments": top20,
    }

with open(ITER_DIR / "h03_sv4_go_biology.json", "w") as f:
    json.dump(h03_results, f, indent=2)
print("  Saved h03_sv4_go_biology.json", flush=True)


# ═══════════════════════════════════════════════════════════════════════════════
# Combined results summary
# ═══════════════════════════════════════════════════════════════════════════════
combined = {
    "iter": "iter_0016",
    "H01_sv3_confidence_gradient": {
        "spearman_rho": h01_result["spearman_rho"],
        "spearman_p": h01_result["spearman_p"],
        "quintile_mean_z": h01_result["quintile_mean_z"],
    },
    "H02_axis_independence": {
        "inter_axis_correlations": h02_result["inter_axis_copole_rate_correlations"],
        "axis_dominant_pair_counts": h02_result["axis_dominant_pair_counts"],
        "jaccard_overlaps_summary": {
            k: v["top_top_mean_jaccard"] for k, v in h02_result["jaccard_overlaps"].items()
        },
    },
    "H03_sv4_go_biology": {
        f"layer_{l}": {
            "top_go": h03_results[f"layer_{l}"]["top_go_enrichments"][0]["go_term"],
            "top_pval": h03_results[f"layer_{l}"]["top_go_enrichments"][0]["pval"],
            "top_pole": h03_results[f"layer_{l}"]["top_go_enrichments"][0]["pole"],
            "overlap_genes": h03_results[f"layer_{l}"]["top_go_enrichments"][0]["genes"][:5],
        }
        for l in FOCUS_LAYERS
    },
}
with open(ITER_DIR / "iter0016_results.json", "w") as f:
    json.dump(combined, f, indent=2)
print("\nSaved iter0016_results.json", flush=True)
print("=== iter_0016 done ===", flush=True)
