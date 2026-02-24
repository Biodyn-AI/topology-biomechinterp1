"""
iter_0009 - Multi-Hypothesis Screen

H01: SV2 bottom-pole null — extracellular vesicle (GO:0070062) enrichment vs N=1000 gene-label shuffles.
     Critical missing control from iter_0008 (observed OR=7.99, p=8.7e-9).
H02: 12-layer × 8-compartment systematic scan — gene-label-null-controlled OR for 8 key GO terms across all 12 layers.
     Produces a "layer compartment map" heatmap artifact.
H03: Cross-layer SV1 rank stability — Spearman r between gene SV1 projections at each pair of adjacent layers.
     Also compute SV1 variance explained at each layer. Novel: tests whether secretory axis is layer-stable.

Data: layer_gene_embeddings.npy [12, 4803, 512], cycle1_edge_dataset.tsv (209 named genes)
GO annotations from local pkl cache.
"""

import numpy as np
import json
import csv
import pickle
import sys
from pathlib import Path
from scipy.stats import fisher_exact, spearmanr
from collections import defaultdict

ITER_DIR = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0009"
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

RNG = np.random.default_rng(42)
N_SHUFFLE_H1 = 1000   # null replicates for H01
N_SHUFFLE_H2 = 200    # null replicates for H02 (8 terms × 12 layers = 96 tests; keep feasible)
TOP_K = 52            # top/bottom-52 = ~25% of 209 genes

MIN_TERM = 5
MAX_TERM = 200

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
            aspect = t.get("aspect", t.get("Category", "P"))
        else:
            go_id = str(t)
            aspect = "P"
        # Accept CC terms too (needed for compartment tests)
        if go_id.startswith("GO:"):
            term2genes[go_id].add(g)

def fisher_or_p(top_set, term_genes, universe_genes):
    top = set(top_set) & universe_genes
    ann = term_genes & universe_genes
    a = len(top & ann)
    b = len(ann - top)
    c = len(top - ann)
    d = len(universe_genes - top - ann)
    if a == 0:
        return 1.0, 0.0, a, b, c, d
    _, p = fisher_exact([[a, b], [c, d]], alternative='greater')
    or_ = (a * d) / (b * c) if (b * c) > 0 else np.inf
    return p, or_, a, b, c, d

def gene_svd_projection(layer_emb_full, gene_idx_list, sv_idx=0):
    """Extract named gene sub-matrix, mean-center, SVD, return projection on sv_idx."""
    sub = layer_emb_full[gene_idx_list]   # [n_genes, 512]
    centered = sub - sub.mean(axis=0)
    U, S, Vt = np.linalg.svd(centered, full_matrices=False)
    return U[:, sv_idx] * S[sv_idx], S  # projections and singular values

# Precompute gene index array (positional in gene_list order)
gene_idx_arr = np.array([gene2idx[g] for g in gene_list])  # shape [n_genes]
universe = set(gene_list)

# ============================================================
# H01: SV2 bottom-pole null (GO:0070062 extracellular vesicle)
# ============================================================
print("\n--- H01: SV2 bottom-pole null ---", flush=True)
LAYER = 11
SV_IDX = 1   # SV2 = index 1

sub11 = emb[LAYER][gene_idx_arr]  # [209, 512]
centered11 = sub11 - sub11.mean(axis=0)
U11, S11, Vt11 = np.linalg.svd(centered11, full_matrices=False)
sv2_proj = U11[:, SV_IDX] * S11[SV_IDX]

# Bottom pole: genes with lowest SV2 projection
sorted_bot = gene_arr[np.argsort(sv2_proj)[:TOP_K]]
term_ev = term2genes.get("GO:0070062", set()) & universe  # extracellular vesicle

obs_p_h1, obs_or_h1, a_h1, b_h1, c_h1, d_h1 = fisher_or_p(sorted_bot, term_ev, universe)
print(f"  Observed: OR={obs_or_h1:.3f}, p={obs_p_h1:.2e}, a={a_h1}, b={b_h1}, n_ev={len(term_ev)}", flush=True)

# Null: shuffle gene labels N=1000 times
null_ps_h1 = []
for _ in range(N_SHUFFLE_H1):
    shuffled = RNG.choice(gene_arr, size=TOP_K, replace=False)
    p_null, _, _, _, _, _ = fisher_or_p(shuffled, term_ev, universe)
    null_ps_h1.append(p_null)
null_ps_h1 = np.array(null_ps_h1)
emp_p_h1 = float((null_ps_h1 <= obs_p_h1).mean())
print(f"  Null mean p={null_ps_h1.mean():.4f}, emp_p={emp_p_h1:.4f}", flush=True)

np.save(ITER_DIR / "h01_sv2_bot_null_ps.npy", null_ps_h1)

# Also identify specific EV genes in bottom pole
ev_genes_in_bot = sorted(set(sorted_bot) & term_ev)
print(f"  EV genes in bottom pole: {ev_genes_in_bot}", flush=True)

# ============================================================
# H02: 12-layer × 8-compartment systematic scan
# ============================================================
print("\n--- H02: 12-layer × 8-compartment scan ---", flush=True)

# Key compartment terms (CC and BP mixed, biologically meaningful)
COMPARTMENTS = {
    "mitochondrion":         "GO:0005739",
    "ER_lumen":              "GO:0005788",
    "extracellular_vesicle": "GO:0070062",
    "cytoskeleton":          "GO:0005856",
    "plasma_membrane":       "GO:0005886",
    "nucleus":               "GO:0005634",
    "secreted":              "GO:0005615",
    "ribosome":              "GO:0005840",
}

scan_results = {}  # layer -> term -> {or, p, emp_p}

for layer in range(12):
    print(f"  Layer {layer}...", flush=True, end=" ")
    sub = emb[layer][gene_idx_arr]
    centered = sub - sub.mean(axis=0)
    U, S, Vt = np.linalg.svd(centered, full_matrices=False)
    sv1_proj = U[:, 0] * S[0]
    # top/bottom TOP_K by SV1
    sorted_top_sv1 = gene_arr[np.argsort(sv1_proj)[-TOP_K:]]
    sorted_bot_sv1 = gene_arr[np.argsort(sv1_proj)[:TOP_K]]

    layer_res = {}
    for comp_name, go_term in COMPARTMENTS.items():
        term_g = term2genes.get(go_term, set()) & universe
        if len(term_g) < MIN_TERM:
            layer_res[comp_name] = {"or": None, "p": None, "emp_p": None, "n_ann": len(term_g)}
            continue
        # Test both poles; use min p (better signal)
        p_top, or_top, a_top, *_ = fisher_or_p(sorted_top_sv1, term_g, universe)
        p_bot, or_bot, a_bot, *_ = fisher_or_p(sorted_bot_sv1, term_g, universe)
        if p_top < p_bot:
            obs_p, obs_or, pole = p_top, or_top, "top"
        else:
            obs_p, obs_or, pole = p_bot, or_bot, "bot"

        # Null (N=N_SHUFFLE_H2)
        null_ps = []
        for _ in range(N_SHUFFLE_H2):
            sh = RNG.choice(gene_arr, size=TOP_K, replace=False)
            pn, _, _, _, _, _ = fisher_or_p(sh, term_g, universe)
            null_ps.append(pn)
        null_ps = np.array(null_ps)
        emp_p = float((null_ps <= obs_p).mean())

        layer_res[comp_name] = {
            "or": float(obs_or) if not np.isinf(obs_or) else 999.0,
            "p": float(obs_p),
            "emp_p": emp_p,
            "n_ann": len(term_g),
            "pole": pole,
            "a": int(a_top) if pole == "top" else int(a_bot),
        }
    scan_results[layer] = layer_res
    print("done", flush=True)

# Save scan results
with open(ITER_DIR / "h02_layer_compartment_scan.json", "w") as f:
    json.dump(scan_results, f, indent=2)

# Build a simple OR table CSV
csv_rows = []
for layer, lres in scan_results.items():
    for comp, vals in lres.items():
        csv_rows.append({
            "layer": layer,
            "compartment": comp,
            "go_term": COMPARTMENTS[comp],
            "n_ann": vals.get("n_ann"),
            "OR": vals.get("or"),
            "p": vals.get("p"),
            "emp_p": vals.get("emp_p"),
            "pole": vals.get("pole"),
        })
with open(ITER_DIR / "h02_layer_compartment_map.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["layer","compartment","go_term","n_ann","OR","p","emp_p","pole"])
    w.writeheader()
    w.writerows(csv_rows)
print(f"  Saved {len(csv_rows)} rows to h02_layer_compartment_map.csv", flush=True)

# ============================================================
# H03: Cross-layer SV1 rank stability
# ============================================================
print("\n--- H03: Cross-layer SV1 rank stability ---", flush=True)

sv1_projections = []   # [12, n_genes]
sv1_var_explained = [] # [12]
for layer in range(12):
    sub = emb[layer][gene_idx_arr]
    centered = sub - sub.mean(axis=0)
    U, S, Vt = np.linalg.svd(centered, full_matrices=False)
    proj = U[:, 0] * S[0]
    sv1_projections.append(proj)
    var_exp = (S[0]**2) / (S**2).sum()
    sv1_var_explained.append(float(var_exp))

sv1_projections = np.array(sv1_projections)  # [12, 209]
np.save(ITER_DIR / "h03_sv1_projections_12layers.npy", sv1_projections)

# Spearman correlation matrix between all layer pairs
spearman_mat = np.zeros((12, 12))
for i in range(12):
    for j in range(12):
        r, _ = spearmanr(sv1_projections[i], sv1_projections[j])
        spearman_mat[i, j] = r
np.save(ITER_DIR / "h03_sv1_spearman_mat.npy", spearman_mat)

# Adjacent-layer Spearman
adjacent_r = [float(spearman_mat[i, i+1]) for i in range(11)]
print(f"  Adjacent-layer Spearman r: {[f'{r:.3f}' for r in adjacent_r]}", flush=True)
print(f"  SV1 variance explained by layer: {[f'{v:.3f}' for v in sv1_var_explained]}", flush=True)
print(f"  Min/Max adj r: {min(adjacent_r):.3f} / {max(adjacent_r):.3f}", flush=True)

# Null: random column permutation preserves rank structure → use random-gene-subset control
# Null hypothesis: SV1 rank is random (no stable axis). Simulate by projecting shuffled embeddings.
null_adj_r_list = []
for rep in range(200):
    null_proj_layers = []
    for layer in range(12):
        sub = emb[layer][gene_idx_arr].copy()
        # shuffle genes (rows) independently per layer
        perm = RNG.permutation(n_genes)
        null_proj_layers.append(sub[perm, 0])  # take first raw feature as surrogate
    # compute adjacent r for null
    null_adj = []
    for i in range(11):
        r, _ = spearmanr(null_proj_layers[i], null_proj_layers[i+1])
        null_adj.append(r)
    null_adj_r_list.append(np.mean(null_adj))
null_adj_r_arr = np.array(null_adj_r_list)
emp_p_h3 = float((null_adj_r_arr >= np.mean(adjacent_r)).mean())
print(f"  Null mean adj r: {null_adj_r_arr.mean():.4f}, obs mean adj r: {np.mean(adjacent_r):.4f}, emp_p: {emp_p_h3:.4f}", flush=True)

# Save CSV summary
with open(ITER_DIR / "h03_sv1_stability.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["layer_pair", "spearman_r"])
    for i, r in enumerate(adjacent_r):
        w.writerow([f"L{i}-L{i+1}", f"{r:.6f}"])
    w.writerow(["mean_adjacent_r", f"{np.mean(adjacent_r):.6f}"])
    w.writerow(["null_mean_adj_r", f"{null_adj_r_arr.mean():.6f}"])
    w.writerow(["emp_p", f"{emp_p_h3:.4f}"])

# ============================================================
# Aggregate results
# ============================================================
print("\n--- Aggregating results ---", flush=True)

# Find significant cells in H02 scan (emp_p <= 0.05)
h2_sig = []
for comp, go_term in COMPARTMENTS.items():
    for layer, lres in scan_results.items():
        v = lres.get(comp, {})
        if v.get("emp_p") is not None and v["emp_p"] <= 0.05:
            h2_sig.append((int(layer), comp, v["or"], v["p"], v["emp_p"]))
h2_sig.sort(key=lambda x: x[4])
print(f"  H02: {len(h2_sig)} significant (layer, compartment) cells at emp_p<=0.05:", flush=True)
for layer, comp, or_, p, ep in h2_sig[:10]:
    print(f"    L{layer} {comp}: OR={or_:.2f}, p={p:.2e}, emp_p={ep:.3f}", flush=True)

results = {
    "iteration": "iter_0009",
    "H01_sv2_bot_null": {
        "term": "GO:0070062",
        "label": "extracellular_vesicle",
        "n_annotated": len(term_ev),
        "n_in_bottom_pole": len(ev_genes_in_bot),
        "obs_or": float(obs_or_h1) if not np.isinf(obs_or_h1) else 999.0,
        "obs_p": float(obs_p_h1),
        "null_mean_p": float(null_ps_h1.mean()),
        "empirical_p": emp_p_h1,
        "ev_genes": ev_genes_in_bot,
        "pass": emp_p_h1 < 0.01,
    },
    "H02_layer_compartment_scan": {
        "n_compartments": len(COMPARTMENTS),
        "n_layers": 12,
        "n_sig_cells_emp05": len(h2_sig),
        "top_hits": [{"layer": l, "comp": c, "OR": o, "p": p, "emp_p": ep}
                     for l, c, o, p, ep in h2_sig[:10]],
    },
    "H03_sv1_stability": {
        "adjacent_spearman_r": adjacent_r,
        "mean_adjacent_r": float(np.mean(adjacent_r)),
        "min_adjacent_r": float(min(adjacent_r)),
        "sv1_var_explained": sv1_var_explained,
        "null_mean_adj_r": float(null_adj_r_arr.mean()),
        "empirical_p": emp_p_h3,
        "pass": emp_p_h3 < 0.05 and np.mean(adjacent_r) > 0.5,
    },
}

with open(ITER_DIR / "iter0009_results.json", "w") as f:
    json.dump(results, f, indent=2)
print("Saved iter0009_results.json", flush=True)
print("Done.", flush=True)
