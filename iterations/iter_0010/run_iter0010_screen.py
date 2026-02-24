"""
iter_0010 - Multi-Hypothesis Screen

H01: SV2/SV3 12-layer × 9-compartment scan with null control (N=200 shuffles per cell).
     Do higher SVD components (SV2, SV3) also encode GO compartment structure?
     Which poles and which layers? Novel: extends SV1-only scan to SV2/SV3.

H02: TRRUST TF-target co-pole test in SVD projection space (layer 11).
     Test whether TF-target pairs co-localize in SV1 poles more than random gene pairs.
     Null: 1000 random same-size pairs. Novel: first regulatory graph topology → SVD geometry link.

H03: SV1/SV2/SV3 variance-explained + annotation density confounder check.
     Track variance explained by SV1, SV2, SV3 across 12 layers.
     Also check whether enrichment results are driven by annotation density (n_ann) confound.
     Validation of prior claims.

Data: layer_gene_embeddings.npy [12, 4803, 512], cycle1_edge_dataset.tsv (209 named genes)
TRRUST: trrust_human.tsv (TF, target, mode, pmid)
GO annotations from local pkl cache.
"""

import numpy as np
import json
import csv
import pickle
import sys
from pathlib import Path
from scipy.stats import fisher_exact, spearmanr, mannwhitneyu
from collections import defaultdict

ITER_DIR = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0010"
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
GENE2GO_PKL = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/single_cell_mechinterp/data/perturb/gene2go_all.pkl"
)
TRRUST_PATH = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/single_cell_mechinterp/external/networks/trrust_human.tsv"
)

RNG = np.random.default_rng(42)
N_SHUFFLE_H1 = 200    # per layer-compartment-SV cell
N_SHUFFLE_H2 = 1000   # for TF-target co-pole null
TOP_K = 52            # ~25% of 209 genes

ITER_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Load shared data
# ---------------------------------------------------------------------------
print("Loading embeddings...", flush=True)
emb = np.load(EMB_PATH)   # [12, 4803, 512]
print(f"  Shape: {emb.shape}", flush=True)

gene2idx = {}
with open(EDGE_PATH) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        gene2idx[row['source'].strip().upper()] = int(row['source_idx'])
        gene2idx[row['target'].strip().upper()] = int(row['target_idx'])
gene_list = sorted(gene2idx.keys())
gene_arr = np.array(gene_list)
n_genes = len(gene_list)
print(f"  Named genes: {n_genes}", flush=True)

# Build index array for named genes in full embedding vocab
named_idx = np.array([gene2idx[g] for g in gene_list])

print("Loading GO annotations...", flush=True)
with open(GENE2GO_PKL, "rb") as f:
    gene2go_raw = pickle.load(f)
gene_set = set(gene_list)
term2genes = defaultdict(set)
for gene, terms in gene2go_raw.items():
    g = gene.strip().upper()
    if g not in gene_set:
        continue
    for t in terms:
        if isinstance(t, dict):
            go_id = t.get("GO_ID", t.get("go_id", ""))
        else:
            go_id = str(t)
        if go_id.startswith("GO:"):
            term2genes[go_id].add(g)

# ---------------------------------------------------------------------------
# GO compartment definitions (same as iter_0009 + 1 new: lysosome)
# ---------------------------------------------------------------------------
COMPARTMENTS = {
    "mitochondrion":       "GO:0005739",
    "ER_lumen":            "GO:0005788",
    "extracellular_vesicle": "GO:0070062",
    "cytoskeleton":        "GO:0005856",
    "plasma_membrane":     "GO:0005886",
    "nucleus":             "GO:0005634",
    "secreted":            "GO:0005615",
    "cytoplasm":           "GO:0005737",
    "lysosome":            "GO:0005764",   # new in iter_0010
}

universe_genes = set(gene_list)

def fisher_or_p(top_set, term_genes, universe_genes):
    top = set(top_set) & universe_genes
    ann = term_genes & universe_genes
    a = len(top & ann)
    b = len(ann - top)
    c = len(top - ann)
    d = len(universe_genes - top - ann)
    if a == 0:
        return 1.0, 0.0, a
    _, p = fisher_exact([[a, b], [c, d]], alternative='greater')
    or_ = (a * d) / (b * c) if (b * c) > 0 else np.inf
    return p, or_, a

def svd_top_bottom(emb_matrix, named_idx, top_k, sv_idx):
    """Return top-k and bottom-k gene labels for SVD component sv_idx."""
    mat = emb_matrix[named_idx, :]   # [n_genes, 512]
    mat_c = mat - mat.mean(axis=0)
    U, S, Vt = np.linalg.svd(mat_c, full_matrices=False)
    proj = U[:, sv_idx] * S[sv_idx]   # signed projection
    order = np.argsort(proj)
    bottom_genes = set(gene_arr[order[:top_k]])
    top_genes = set(gene_arr[order[-top_k:]])
    return top_genes, bottom_genes, S

# ============================================================
# H01: SV2/SV3 12-layer × 9-compartment scan
# ============================================================
print("\n=== H01: SV2/SV3 12-layer × 9-compartment scan ===", flush=True)

h01_results = {}  # key: (sv_name, layer, comp_name) → dict

for sv_idx, sv_name in [(1, "SV2"), (2, "SV3")]:
    h01_results[sv_name] = {}
    for layer in range(12):
        h01_results[sv_name][layer] = {}
        top_genes, bot_genes, S = svd_top_bottom(emb[layer], named_idx, TOP_K, sv_idx)

        for comp_name, go_id in COMPARTMENTS.items():
            term_genes = term2genes.get(go_id, set()) & universe_genes
            n_ann = len(term_genes)
            if n_ann < 3:
                h01_results[sv_name][layer][comp_name] = {"skip": True, "n_ann": n_ann}
                continue

            # Test both poles
            best_p, best_or, best_pole, best_a = 1.0, 0.0, "none", 0
            for pole_name, pole_set in [("top", top_genes), ("bot", bot_genes)]:
                p, or_, a = fisher_or_p(pole_set, term_genes, universe_genes)
                if or_ > best_or:
                    best_p, best_or, best_pole, best_a = p, or_, pole_name, a

            # Empirical null: gene-label shuffle
            null_ors = []
            for _ in range(N_SHUFFLE_H1):
                shuffled = gene_arr.copy()
                RNG.shuffle(shuffled)
                shuf_top = set(shuffled[:TOP_K])
                shuf_bot = set(shuffled[-TOP_K:])
                for pole_set in [shuf_top, shuf_bot]:
                    _, or_n, _ = fisher_or_p(pole_set, term_genes, universe_genes)
                    null_ors.append(or_n)
            null_ors = np.array(null_ors)
            emp_p = float(np.mean(null_ors >= best_or))

            h01_results[sv_name][layer][comp_name] = {
                "or": round(best_or, 3),
                "p": round(best_p, 6),
                "emp_p": round(emp_p, 4),
                "n_ann": n_ann,
                "pole": best_pole,
                "a": best_a,
            }
        print(f"  {sv_name} layer {layer} done", flush=True)

# Save H01 results
with open(ITER_DIR / "h01_sv2sv3_layer_compartment_scan.json", "w") as f:
    json.dump(h01_results, f, indent=2)

# Find significant hits (emp_p < 0.05)
h01_hits = []
for sv_name in ["SV2", "SV3"]:
    for layer in range(12):
        for comp_name, res in h01_results[sv_name][layer].items():
            if res.get("skip"):
                continue
            if res["emp_p"] < 0.05:
                h01_hits.append({
                    "sv": sv_name, "layer": layer, "comp": comp_name,
                    "OR": res["or"], "p": res["p"], "emp_p": res["emp_p"],
                    "n_ann": res["n_ann"], "pole": res["pole"]
                })
h01_hits.sort(key=lambda x: x["emp_p"])
print(f"  Significant hits (emp_p<0.05): {len(h01_hits)}", flush=True)
for h in h01_hits[:10]:
    print(f"    {h['sv']} L{h['layer']} {h['comp']}: OR={h['OR']}, emp_p={h['emp_p']}, pole={h['pole']}", flush=True)

# ============================================================
# H02: TRRUST TF-target co-pole test at layer 11
# ============================================================
print("\n=== H02: TRRUST TF-target co-pole test ===", flush=True)

# Load TRRUST
tf_target_pairs = []
with open(TRRUST_PATH) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            tf = parts[0].strip().upper()
            tgt = parts[1].strip().upper()
            tf_target_pairs.append((tf, tgt))

# Filter to genes in our named set
named_set = set(gene_list)
valid_pairs = [(tf, tgt) for tf, tgt in tf_target_pairs
               if tf in named_set and tgt in named_set and tf != tgt]
print(f"  TRRUST pairs in our gene set: {len(valid_pairs)}", flush=True)

if len(valid_pairs) >= 5:
    # Use layer 11 embeddings
    layer = 11
    mat = emb[layer][named_idx, :]
    mat_c = mat - mat.mean(axis=0)
    U, S, Vt = np.linalg.svd(mat_c, full_matrices=False)

    gene2rank = {g: int(np.argsort(U[:, 0])[i]) for i, g in enumerate(gene_arr)}
    # Rank in [0, n_genes-1]; co-pole = both in top-K or both in bottom-K
    sv1_order = np.argsort(U[:, 0])
    top_set = set(gene_arr[sv1_order[-TOP_K:]])
    bot_set = set(gene_arr[sv1_order[:TOP_K]])

    def is_copole(g1, g2):
        return (g1 in top_set and g2 in top_set) or (g1 in bot_set and g2 in bot_set)

    # Observed co-pole rate
    obs_copole = sum(1 for tf, tgt in valid_pairs if is_copole(tf, tgt))
    obs_rate = obs_copole / len(valid_pairs)

    # Null: random same-size pairs (from named genes)
    gene_arr_local = np.array(gene_list)
    null_rates = []
    for _ in range(N_SHUFFLE_H2):
        idx1 = RNG.integers(0, n_genes, len(valid_pairs))
        idx2 = RNG.integers(0, n_genes, len(valid_pairs))
        pairs_null = [(gene_arr_local[i], gene_arr_local[j]) for i, j in zip(idx1, idx2)
                      if gene_arr_local[i] != gene_arr_local[j]]
        if len(pairs_null) == 0:
            null_rates.append(0.0)
            continue
        null_copole = sum(1 for g1, g2 in pairs_null if is_copole(g1, g2))
        null_rates.append(null_copole / len(pairs_null))
    null_rates = np.array(null_rates)
    emp_p_h2 = float(np.mean(null_rates >= obs_rate))
    null_mean = float(null_rates.mean())
    null_std = float(null_rates.std())

    # Also test SV2 co-pole
    sv2_order = np.argsort(U[:, 1])
    top2_set = set(gene_arr[sv2_order[-TOP_K:]])
    bot2_set = set(gene_arr[sv2_order[:TOP_K]])

    def is_copole2(g1, g2):
        return (g1 in top2_set and g2 in top2_set) or (g1 in bot2_set and g2 in bot2_set)

    obs_copole2 = sum(1 for tf, tgt in valid_pairs if is_copole2(tf, tgt))
    obs_rate2 = obs_copole2 / len(valid_pairs)
    null_rates2 = []
    for _ in range(N_SHUFFLE_H2):
        idx1 = RNG.integers(0, n_genes, len(valid_pairs))
        idx2 = RNG.integers(0, n_genes, len(valid_pairs))
        pairs_null = [(gene_arr_local[i], gene_arr_local[j]) for i, j in zip(idx1, idx2)
                      if gene_arr_local[i] != gene_arr_local[j]]
        if len(pairs_null) == 0:
            null_rates2.append(0.0)
            continue
        null_copole2 = sum(1 for g1, g2 in pairs_null if is_copole2(g1, g2))
        null_rates2.append(null_copole2 / len(pairs_null))
    null_rates2 = np.array(null_rates2)
    emp_p_h2_sv2 = float(np.mean(null_rates2 >= obs_rate2))

    h02_result = {
        "n_valid_pairs": len(valid_pairs),
        "layer": layer,
        "SV1": {
            "obs_copole_rate": round(obs_rate, 4),
            "obs_copole_count": obs_copole,
            "null_mean_rate": round(null_mean, 4),
            "null_std_rate": round(null_std, 4),
            "emp_p": round(emp_p_h2, 4),
        },
        "SV2": {
            "obs_copole_rate": round(obs_rate2, 4),
            "obs_copole_count": obs_copole2,
            "null_mean_rate": round(float(null_rates2.mean()), 4),
            "null_std_rate": round(float(null_rates2.std()), 4),
            "emp_p": round(emp_p_h2_sv2, 4),
        }
    }
    np.save(ITER_DIR / "h02_trrust_copole_null_sv1.npy", null_rates)
    np.save(ITER_DIR / "h02_trrust_copole_null_sv2.npy", null_rates2)
    print(f"  H02 SV1: obs={obs_rate:.4f}, null_mean={null_mean:.4f}, emp_p={emp_p_h2:.4f}", flush=True)
    print(f"  H02 SV2: obs={obs_rate2:.4f}, null_mean={float(null_rates2.mean()):.4f}, emp_p={emp_p_h2_sv2:.4f}", flush=True)
else:
    h02_result = {"error": "insufficient TRRUST pairs in gene set", "n_valid_pairs": len(valid_pairs)}
    print(f"  WARNING: only {len(valid_pairs)} valid pairs, skipping", flush=True)

with open(ITER_DIR / "h02_trrust_copole_result.json", "w") as f:
    json.dump(h02_result, f, indent=2)

# ============================================================
# H03: Spectral ratio profile + annotation density confounder
# ============================================================
print("\n=== H03: Spectral ratio profile + annotation density confounder ===", flush=True)

h03_spectral = []
for layer in range(12):
    mat = emb[layer][named_idx, :]
    mat_c = mat - mat.mean(axis=0)
    U, S, Vt = np.linalg.svd(mat_c, full_matrices=False)
    total_var = float(np.sum(S**2))
    sv1_var = float(S[0]**2) / total_var
    sv2_var = float(S[1]**2) / total_var
    sv3_var = float(S[2]**2) / total_var
    sv1_sv2_ratio = float(S[0] / S[1]) if S[1] > 0 else np.inf
    h03_spectral.append({
        "layer": layer,
        "SV1_var_frac": round(sv1_var, 4),
        "SV2_var_frac": round(sv2_var, 4),
        "SV3_var_frac": round(sv3_var, 4),
        "SV1_SV2_ratio": round(sv1_sv2_ratio, 3),
        "SV1_magnitude": round(float(S[0]), 3),
        "SV2_magnitude": round(float(S[1]), 3),
        "SV3_magnitude": round(float(S[2]), 3),
    })

# Annotation density confounder: does enrichment track n_ann?
# Use iter_0009 H02 results as reference (ER_lumen, mito, EV across layers)
# Test: Spearman r between n_ann and OR for each significant compartment-SV combo from H01
ann_density_check = {}
for sv_name in ["SV2", "SV3"]:
    layers_arr = list(range(12))
    for comp_name in COMPARTMENTS:
        ors = []
        n_anns = []
        for layer in layers_arr:
            res = h01_results[sv_name][layer].get(comp_name, {})
            if res.get("skip"):
                continue
            ors.append(res["or"])
            n_anns.append(res["n_ann"])
        if len(ors) >= 5:
            r, p = spearmanr(n_anns, ors)
            ann_density_check[f"{sv_name}_{comp_name}"] = {
                "spearman_r": round(float(r), 3),
                "p": round(float(p), 4),
                "n": len(ors)
            }

h03_result = {
    "spectral_profile": h03_spectral,
    "annotation_density_confounder": ann_density_check
}

# Save CSV for spectral profile
with open(ITER_DIR / "h03_spectral_profile.csv", "w", newline='') as f:
    writer = csv.DictWriter(f, fieldnames=["layer","SV1_var_frac","SV2_var_frac","SV3_var_frac","SV1_SV2_ratio","SV1_magnitude","SV2_magnitude","SV3_magnitude"])
    writer.writeheader()
    writer.writerows(h03_spectral)

with open(ITER_DIR / "h03_annotation_density_confounder.json", "w") as f:
    json.dump(ann_density_check, f, indent=2)

print("  Spectral profile (SV1/SV2/SV3 var fractions):", flush=True)
for row in h03_spectral:
    print(f"    L{row['layer']}: SV1={row['SV1_var_frac']:.3f} SV2={row['SV2_var_frac']:.3f} SV3={row['SV3_var_frac']:.3f} ratio={row['SV1_SV2_ratio']:.2f}", flush=True)

# Confounder check summary
sig_conf = [(k, v) for k, v in ann_density_check.items() if v["p"] < 0.05]
print(f"  Ann density confounders (n_ann ~ OR, p<0.05): {len(sig_conf)}", flush=True)
for k, v in sig_conf[:5]:
    print(f"    {k}: r={v['spearman_r']}, p={v['p']}", flush=True)

# ============================================================
# Consolidate and save master results
# ============================================================
master_results = {
    "iteration": "iter_0010",
    "H01_sv2sv3_compartment_scan": {
        "n_sig_hits_emp05": len(h01_hits),
        "top_hits": h01_hits[:10],
    },
    "H02_trrust_copole": h02_result,
    "H03_spectral_profile": {
        "spectral_profile": h03_spectral,
        "n_ann_confounders_p05": len(sig_conf),
    }
}

with open(ITER_DIR / "iter0010_results.json", "w") as f:
    json.dump(master_results, f, indent=2)

print("\n=== iter_0010 DONE ===", flush=True)
print(f"  H01 sig hits: {len(h01_hits)}", flush=True)
if "n_valid_pairs" in h02_result:
    print(f"  H02 SV1 emp_p: {h02_result.get('SV1', {}).get('emp_p', 'N/A')}", flush=True)
print(f"  H03 spectral layers: {len(h03_spectral)}", flush=True)
