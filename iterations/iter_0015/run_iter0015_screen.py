"""
iter_0015 - Multi-Hypothesis Screen

H01 (SV4/SV5 PPI co-pole test): Extend the SV2/SV3 PPI co-pole analysis to SV4 and SV5.
    Tests whether PPI geometry extends beyond the first 3 singular vectors.
    Zero new infrastructure - natural extension of proven pipeline.
    Uses STRING score>=0.4 (N=3364 pairs) for adequate power.

H02 (STRING confidence gradient): Test if PPI co-pole enrichment (z-score) correlates with
    STRING interaction confidence score (continuous). Bin pairs by score quintile, measure
    mean z-score per quintile. Tests whether geometry continuously encodes interaction strength.

H03 (GO co-annotation vs PPI): For each gene pair, compute GO Jaccard similarity (# shared terms
    / # union terms using BP+CC). Split pairs into GO-only (high GO Jaccard, low STRING score),
    PPI-only (low GO Jaccard, high STRING score), and both. Test SV2 co-pole enrichment per group.
    Discriminates functional vs physical proximity as driver of SV2 geometry.

Data: layer_gene_embeddings.npy [12, 4803, 512], gene2go_all.pkl, string_ppi_named_genes.json
"""

import numpy as np
import json
import sys
import pickle
from pathlib import Path
from collections import defaultdict
from scipy.stats import spearmanr, pearsonr

ITER_DIR = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0015"
)
ITER12_DIR = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0012"
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
STRING_JSON_04 = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0014"
) / "h03_hub_degree_control_04.json"  # Has STRING edges score>=0.4 baked in? Check.
EDGE_PATH = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_38_geometric_residual_stream_interpretability"
    "/implementation/outputs/cycle1_main/cycle1_edge_dataset.tsv"
)
STRING_API_CACHE = ITER_DIR / "string_ppi_score04_cache.json"

ITER_DIR.mkdir(parents=True, exist_ok=True)

# ─────────────────────────────────────────────
# Load embeddings and gene list
# ─────────────────────────────────────────────
print("Loading embeddings ...", flush=True)
emb = np.load(EMB_PATH)  # [12, 4803, 512]
N_LAYERS, N_GENES, N_DIM = emb.shape
print(f"  emb shape: {emb.shape}", flush=True)

# gene_set built below via edge dataset with idx

# Build gene→index map from embeddings file (need to load separately)
# Use source_idx/target_idx from edge dataset (actual embedding row indices)
import csv
gene_to_emb_idx = {}
with open(EDGE_PATH) as f:
    rdr = csv.DictReader(f, delimiter="\t")
    for row in rdr:
        src, si = row["source"], row.get("source_idx","")
        tgt, ti = row["target"], row.get("target_idx","")
        if src and si and src not in gene_to_emb_idx:
            gene_to_emb_idx[src] = int(si)
        if tgt and ti and tgt not in gene_to_emb_idx:
            gene_to_emb_idx[tgt] = int(ti)

named_genes = sorted(gene_to_emb_idx.keys())
named_idx = np.array([gene_to_emb_idx[g] for g in named_genes])
print(f"  named genes: {len(named_genes)}, idx range [{named_idx.min()},{named_idx.max()}]", flush=True)

# ─────────────────────────────────────────────
# Load STRING pairs (score>=0.4) — fetch or reconstruct
# ─────────────────────────────────────────────
print("Loading STRING pairs ...", flush=True)

# Try to load from iter_0012 string json (score>=0.7) first, then fetch broader set
STRING_JSON_07 = ITER12_DIR / "string_ppi_named_genes.json"
with open(STRING_JSON_07) as f:
    string07 = json.load(f)

# string07 has pairs with score>=0.7. We need score>=0.4 for the gradient test.
# Check if we have a cached version from iter_0014
# iter_0014 h03 ran hub-degree control with score>=0.4 — it may store pairs in layers dict
# Reconstruct by fetching from STRING API
if STRING_API_CACHE.exists():
    with open(STRING_API_CACHE) as f:
        string04_data = json.load(f)
    all_pairs_raw = string04_data["pairs"]
    print(f"  Loaded {len(all_pairs_raw)} pairs from cache", flush=True)
else:
    print("  Fetching STRING score>=0.4 from API ...", flush=True)
    import urllib.request, urllib.parse
    genes_str = "%0d".join(named_genes[:200])  # API limit ~200 at a time
    # Batch into chunks
    chunk_size = 200
    all_pairs_raw = []
    seen_pairs = set()
    for ci in range(0, len(named_genes), chunk_size):
        chunk = named_genes[ci:ci+chunk_size]
        genes_str = "%0d".join(chunk)
        url = (
            "https://string-db.org/api/json/network?"
            f"identifiers={urllib.parse.quote(genes_str)}&"
            "species=9606&"
            "required_score=400&"
            "caller_identity=topology_screen_iter0015"
        )
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                data = json.loads(resp.read())
            for e in data:
                g1, g2 = e.get("preferredName_A",""), e.get("preferredName_B","")
                score = e.get("score", 0)
                key = tuple(sorted([g1,g2]))
                if key not in seen_pairs and g1 in gene2idx and g2 in gene2idx:
                    seen_pairs.add(key)
                    all_pairs_raw.append({"g1": g1, "g2": g2, "score": score})
        except Exception as ex:
            print(f"  chunk {ci}: fetch failed: {ex}", flush=True)
    print(f"  Fetched {len(all_pairs_raw)} unique pairs", flush=True)
    with open(STRING_API_CACHE, "w") as f:
        json.dump({"pairs": all_pairs_raw}, f)

# Filter to pairs where both genes are in named_genes
named_set = set(named_genes)
pairs_filtered = [(p["g1"], p["g2"], p["score"]) for p in all_pairs_raw
                  if p["g1"] in named_set and p["g2"] in named_set]
print(f"  Filtered pairs (both named): {len(pairs_filtered)}", flush=True)

# ─────────────────────────────────────────────
# SVD helper (reuse across hypotheses)
# ─────────────────────────────────────────────
def compute_svd_projections(layer_idx, sv_indices):
    """Return dict sv→ projections array [N_named_genes]"""
    X = emb[layer_idx, named_idx, :]  # [N_named, 512]
    X_c = X - X.mean(axis=0)
    U, S, Vt = np.linalg.svd(X_c, full_matrices=False)
    # U: [N_named, K]; columns are left singular vectors
    result = {}
    for sv in sv_indices:
        result[sv] = U[:, sv]  # projection of each gene onto sv-th axis
    return result

def copole_rate(proj, g1_list, g2_list, K=52):
    """Fraction of pairs where both genes in same pole (top-K or bottom-K)."""
    order = np.argsort(proj)
    bottom_set = set(order[:K])
    top_set = set(order[-K:])
    gene2i = {g: i for i, g in enumerate(named_genes)}
    copole = 0
    valid = 0
    for g1, g2 in zip(g1_list, g2_list):
        if g1 not in gene2i or g2 not in gene2i:
            continue
        i1, i2 = gene2i[g1], gene2i[g2]
        valid += 1
        if (i1 in top_set and i2 in top_set) or (i1 in bottom_set and i2 in bottom_set):
            copole += 1
    return copole / valid if valid > 0 else np.nan, valid

def shuffle_null(proj, g1_list, g2_list, K=52, N=500, rng=None):
    """Shuffle gene labels N times and compute copole rate each time."""
    if rng is None:
        rng = np.random.default_rng(42)
    null_rates = []
    for _ in range(N):
        perm = rng.permutation(proj)
        r, _ = copole_rate(perm, g1_list, g2_list, K=K)
        null_rates.append(r)
    return np.array(null_rates)

# ─────────────────────────────────────────────
# H01: SV4/SV5 PPI co-pole test
# ─────────────────────────────────────────────
print("\n=== H01: SV4/SV5 PPI co-pole test ===", flush=True)

g1_list = [p[0] for p in pairs_filtered]
g2_list = [p[1] for p in pairs_filtered]
scores_list = [p[2] for p in pairs_filtered]

N_NULL = 300
K = 52

h01_layers = []
rng = np.random.default_rng(0)

for layer in range(N_LAYERS):
    projs = compute_svd_projections(layer, sv_indices=[0,1,2,3,4])
    layer_entry = {"layer": layer, "svd_axes": {}}
    for sv_idx, sv_label in [(3, "SV4"), (4, "SV5")]:
        proj = projs[sv_idx]
        obs, n_valid = copole_rate(proj, g1_list, g2_list, K=K)
        null = shuffle_null(proj, g1_list, g2_list, K=K, N=N_NULL, rng=rng)
        z = (obs - null.mean()) / (null.std() + 1e-12)
        emp_p = (null >= obs).mean()
        layer_entry["svd_axes"][sv_label] = {
            "obs_copole": float(obs),
            "null_mean": float(null.mean()),
            "null_std": float(null.std()),
            "z_score": float(z),
            "emp_p": float(emp_p),
            "n_pairs": int(n_valid)
        }
    h01_layers.append(layer_entry)
    sv4 = layer_entry["svd_axes"]["SV4"]
    sv5 = layer_entry["svd_axes"]["SV5"]
    print(f"  L{layer:02d} SV4 z={sv4['z_score']:.3f} emp_p={sv4['emp_p']:.3f} | "
          f"SV5 z={sv5['z_score']:.3f} emp_p={sv5['emp_p']:.3f}", flush=True)

# Summaries
sv4_zscores = [l["svd_axes"]["SV4"]["z_score"] for l in h01_layers]
sv5_zscores = [l["svd_axes"]["SV5"]["z_score"] for l in h01_layers]
sv4_sig = sum(l["svd_axes"]["SV4"]["emp_p"] < 0.05 for l in h01_layers)
sv5_sig = sum(l["svd_axes"]["SV5"]["emp_p"] < 0.05 for l in h01_layers)

h01_result = {
    "meta": {
        "hypothesis": "SV4/SV5 PPI co-pole test",
        "K": K, "N_pairs": len(pairs_filtered), "N_null": N_NULL
    },
    "summary": {
        "SV4_mean_z": float(np.mean(sv4_zscores)),
        "SV4_sig_layers": int(sv4_sig),
        "SV5_mean_z": float(np.mean(sv5_zscores)),
        "SV5_sig_layers": int(sv5_sig)
    },
    "layers": h01_layers
}
with open(ITER_DIR / "h01_sv45_ppi_copole.json", "w") as f:
    json.dump(h01_result, f, indent=2)
print(f"  SV4 mean_z={np.mean(sv4_zscores):.3f} sig={sv4_sig}/12")
print(f"  SV5 mean_z={np.mean(sv5_zscores):.3f} sig={sv5_sig}/12")

# ─────────────────────────────────────────────
# H02: STRING confidence gradient
# ─────────────────────────────────────────────
print("\n=== H02: STRING confidence gradient ===", flush=True)

# Use score values from pairs_filtered; quantize into 5 bins
scores_arr = np.array(scores_list)
# Quintile edges
q_edges = np.quantile(scores_arr, [0, 0.2, 0.4, 0.6, 0.8, 1.0])
print(f"  Score quintile edges: {q_edges.round(3)}", flush=True)

bins = []
for qi in range(5):
    lo, hi = q_edges[qi], q_edges[qi+1]
    if qi < 4:
        mask = (scores_arr >= lo) & (scores_arr < hi)
    else:
        mask = (scores_arr >= lo) & (scores_arr <= hi)
    bin_pairs = [(g1_list[i], g2_list[i]) for i in range(len(g1_list)) if mask[i]]
    bins.append({"quintile": qi+1, "lo": float(lo), "hi": float(hi),
                 "n_pairs": len(bin_pairs), "g1g2": bin_pairs})
    print(f"  Quintile {qi+1}: [{lo:.3f},{hi:.3f}] N={len(bin_pairs)}", flush=True)

# Focus on layer 7 (middle) and layer 11 (deepest) for speed
# Full cross-layer: compute mean z per quintile at each layer
# but do full 12 layers to get correlation

N_NULL2 = 200
h02_layers = []
rng2 = np.random.default_rng(1)

for layer in range(N_LAYERS):
    proj_sv2 = compute_svd_projections(layer, sv_indices=[1])[1]  # SV2
    quintile_results = []
    for b in bins:
        if b["n_pairs"] < 5:
            quintile_results.append(None)
            continue
        obs, nv = copole_rate(proj_sv2, [p[0] for p in b["g1g2"]], [p[1] for p in b["g1g2"]], K=K)
        null = shuffle_null(proj_sv2, [p[0] for p in b["g1g2"]], [p[1] for p in b["g1g2"]],
                            K=K, N=N_NULL2, rng=rng2)
        z = (obs - null.mean()) / (null.std() + 1e-12)
        emp_p = float((null >= obs).mean())
        quintile_results.append({"q": b["quintile"], "lo": b["lo"], "hi": b["hi"],
                                  "n": nv, "obs": float(obs), "z": float(z), "emp_p": emp_p})
    h02_layers.append({"layer": layer, "quintiles": quintile_results})

# Compute Spearman corr between quintile midpoint score and mean z across 12 layers
q_mids = [((b["lo"] + b["hi"]) / 2) for b in bins]
# mean z per quintile across layers
q_mean_z = []
for qi in range(5):
    zvals = [l["quintiles"][qi]["z"] for l in h02_layers if l["quintiles"][qi] is not None]
    q_mean_z.append(np.mean(zvals) if zvals else np.nan)

rho, pval = spearmanr(q_mids, q_mean_z)
print(f"  Spearman(quintile_mid, mean_z): rho={rho:.3f}, p={pval:.4f}", flush=True)
print(f"  Per-quintile mean_z: {[round(z,3) for z in q_mean_z]}", flush=True)

h02_result = {
    "meta": {"hypothesis": "STRING confidence gradient", "K": K, "N_null": N_NULL2},
    "quintile_edges": list(q_edges.tolist()),
    "quintile_mids": q_mids,
    "quintile_mean_z_across_layers": q_mean_z,
    "spearman_rho": float(rho),
    "spearman_p": float(pval),
    "layers": h02_layers
}
with open(ITER_DIR / "h02_string_confidence_gradient.json", "w") as f:
    json.dump(h02_result, f, indent=2)

# ─────────────────────────────────────────────
# H03: GO co-annotation vs PPI
# ─────────────────────────────────────────────
print("\n=== H03: GO co-annotation vs PPI ===", flush=True)

print("  Loading gene2go ...", flush=True)
with open(GENE2GO_PKL, "rb") as f:
    gene2go_raw = pickle.load(f)

# gene2go_raw may be {gene: set(go_terms)} or similar
# Check structure
sample_key = list(gene2go_raw.keys())[0]
sample_val = list(gene2go_raw.values())[0]
print(f"  gene2go sample key={sample_key!r} val_type={type(sample_val).__name__} len={len(sample_val)}", flush=True)

# Build gene→set(go_terms) normalized to uppercase gene names
gene2go = {}
for k, v in gene2go_raw.items():
    gene_name = k.upper() if isinstance(k, str) else str(k)
    if isinstance(v, (set, list, frozenset)):
        gene2go[gene_name] = set(v)
    elif isinstance(v, dict):
        # might be {go_id: evidence_codes}
        gene2go[gene_name] = set(v.keys())
    else:
        gene2go[gene_name] = set()

def go_jaccard(g1, g2):
    s1 = gene2go.get(g1.upper(), set()) | gene2go.get(g1, set())
    s2 = gene2go.get(g2.upper(), set()) | gene2go.get(g2, set())
    if not s1 or not s2:
        return np.nan
    inter = len(s1 & s2)
    union = len(s1 | s2)
    return inter / union if union > 0 else 0.0

# Compute GO Jaccard for all STRING pairs (score>=0.4)
print("  Computing GO Jaccard for pairs ...", flush=True)
go_jac = np.array([go_jaccard(g1, g2) for g1, g2 in zip(g1_list, g2_list)])
valid_mask = ~np.isnan(go_jac)
print(f"  Valid GO Jaccard: {valid_mask.sum()} / {len(go_jac)}", flush=True)
print(f"  GO Jaccard: mean={np.nanmean(go_jac):.3f} std={np.nanstd(go_jac):.3f}", flush=True)

# Median split on GO Jaccard and STRING score among valid pairs
go_median = np.nanmedian(go_jac[valid_mask])
str_median = np.median(scores_arr[valid_mask])
print(f"  GO Jaccard median: {go_median:.3f}, STRING score median: {str_median:.3f}", flush=True)

# Groups (among valid pairs):
# PPI-only: high STRING score (>=str_median), low GO Jaccard (<go_median)
# GO-only:  low STRING score (<str_median), high GO Jaccard (>=go_median)
# Both:     high STRING, high GO
# Neither:  low STRING, low GO

def group_mask(go_jac, scores_arr, valid_mask, go_thresh, str_thresh, go_hi, str_hi):
    m = valid_mask.copy()
    if go_hi:
        m &= (go_jac >= go_thresh)
    else:
        m &= (go_jac < go_thresh)
    if str_hi:
        m &= (scores_arr >= str_thresh)
    else:
        m &= (scores_arr < str_thresh)
    return m

mask_ppi_only = group_mask(go_jac, scores_arr, valid_mask, go_median, str_median, False, True)
mask_go_only  = group_mask(go_jac, scores_arr, valid_mask, go_median, str_median, True, False)
mask_both     = group_mask(go_jac, scores_arr, valid_mask, go_median, str_median, True, True)
mask_neither  = group_mask(go_jac, scores_arr, valid_mask, go_median, str_median, False, False)

print(f"  PPI-only N={mask_ppi_only.sum()}, GO-only N={mask_go_only.sum()}, "
      f"Both N={mask_both.sum()}, Neither N={mask_neither.sum()}", flush=True)

N_NULL3 = 300
rng3 = np.random.default_rng(2)

def copole_from_mask(m, g1_list, g2_list, proj, K=52):
    mg1 = [g1_list[i] for i in range(len(g1_list)) if m[i]]
    mg2 = [g2_list[i] for i in range(len(g2_list)) if m[i]]
    return mg1, mg2

group_names = ["PPI-only", "GO-only", "Both", "Neither"]
group_masks = [mask_ppi_only, mask_go_only, mask_both, mask_neither]

h03_layers = []
for layer in range(N_LAYERS):
    proj_sv2 = compute_svd_projections(layer, sv_indices=[1])[1]
    layer_entry = {"layer": layer, "groups": {}}
    for gname, gmask in zip(group_names, group_masks):
        g1s, g2s = copole_from_mask(gmask, g1_list, g2_list, proj_sv2)
        if len(g1s) < 5:
            layer_entry["groups"][gname] = {"n": len(g1s), "z": None, "emp_p": None}
            continue
        obs, nv = copole_rate(proj_sv2, g1s, g2s, K=K)
        null = shuffle_null(proj_sv2, g1s, g2s, K=K, N=N_NULL3, rng=rng3)
        z = float((obs - null.mean()) / (null.std() + 1e-12))
        emp_p = float((null >= obs).mean())
        layer_entry["groups"][gname] = {"n": nv, "obs": float(obs), "z": z, "emp_p": emp_p}
    h03_layers.append(layer_entry)

# Summarize
for gname in group_names:
    zvals = [l["groups"][gname]["z"] for l in h03_layers if l["groups"][gname]["z"] is not None]
    sig = sum(l["groups"][gname]["emp_p"] < 0.05 for l in h03_layers
              if l["groups"][gname]["emp_p"] is not None)
    if zvals:
        print(f"  {gname}: mean_z={np.mean(zvals):.3f} sig={sig}/12", flush=True)

h03_result = {
    "meta": {"hypothesis": "GO co-annotation vs PPI driver", "K": K, "N_null": N_NULL3,
             "go_jaccard_median": float(go_median), "string_score_median": float(str_median)},
    "group_sizes": {gn: int(gm.sum()) for gn, gm in zip(group_names, group_masks)},
    "layers": h03_layers
}
with open(ITER_DIR / "h03_go_vs_ppi_driver.json", "w") as f:
    json.dump(h03_result, f, indent=2)

# ─────────────────────────────────────────────
# Compile overall results
# ─────────────────────────────────────────────
results_summary = {
    "iter": "iter_0015",
    "H01_sv45_ppi_copole": h01_result["summary"],
    "H02_confidence_gradient": {
        "spearman_rho": float(rho),
        "spearman_p": float(pval),
        "quintile_mean_z": [round(z,3) for z in q_mean_z]
    },
    "H03_go_vs_ppi_driver": {
        "group_sizes": h03_result["group_sizes"],
        "group_mean_z": {
            gn: round(np.mean([l["groups"][gn]["z"] for l in h03_layers
                              if l["groups"][gn]["z"] is not None]), 3)
            for gn in group_names
        }
    }
}
with open(ITER_DIR / "iter0015_results.json", "w") as f:
    json.dump(results_summary, f, indent=2)

print("\n=== iter_0015 complete ===", flush=True)
print(json.dumps(results_summary, indent=2))
