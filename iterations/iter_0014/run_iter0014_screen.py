"""
iter_0014 - Multi-Hypothesis Screen

H01 (GO enrichment of SV3 poles): For each of 12 layers, identify top-K=52 and bottom-K=52 genes
    by SV3 projection. Run Fisher exact test against GO BP+CC annotations. Compare to gene-label
    shuffle null (N=300). Report top enriched GO terms per layer, and whether SV3 encodes a
    distinct biological theme from SV2.
    Novel: SV3 GO biology never previously characterized.

H02 (Layer-depth Spearman trend on SV2/SV3 z-scores): Use iter_0013 H01 results to fit Spearman
    rank correlation between layer index [0-11] and SV2/SV3 PPI co-pole z-score.
    Tests whether PPI geometry strengthens or weakens monotonically with depth.
    Zero new data needed - uses saved JSON.

H03 (Expanded hub-degree control with score>=0.4 STRING edges): Re-run hub vs non-hub SV2 co-pole
    test with expanded STRING (score>=0.4 instead of >=0.7). Expect ~5-10x more non-hub edges.
    Tests whether PPI geometry in SV2 is driven by hubs or reflects general PPI network structure.

Data: layer_gene_embeddings.npy [12, 4803, 512], gene2go_all.pkl, string_ppi_named_genes.json
"""

import numpy as np
import json
import csv
import sys
import pickle
from pathlib import Path
from collections import defaultdict
from scipy.stats import fisher_exact, spearmanr

ITER_DIR = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0014"
)
PREV_ITER_DIR = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0013"
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
EDGE_PATH = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_38_geometric_residual_stream_interpretability"
    "/implementation/outputs/cycle1_main/cycle1_edge_dataset.tsv"
)
GENE2GO_PKL = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/single_cell_mechinterp/data/perturb/gene2go_all.pkl"
)
STRING_JSON_07 = ITER12_DIR / "string_ppi_named_genes.json"  # score>=0.7

ITER_DIR.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(42)
N_SHUFFLE = 300
K_POLE = 52  # top/bottom K for SV poles


# ============================================================
# Data loading
# ============================================================

print("Loading embeddings...", flush=True)
emb = np.load(EMB_PATH)  # [12, 4803, 512]
n_layers, n_genes_total, n_dim = emb.shape
print(f"  Embeddings: {emb.shape}", flush=True)

# Load named gene set
gene_to_idx = {}
with open(EDGE_PATH) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        src = row.get("source", "")
        tgt = row.get("target", "")
        src_idx = row.get("source_idx", "")
        tgt_idx = row.get("target_idx", "")
        if src and src_idx and src not in gene_to_idx:
            gene_to_idx[src] = int(src_idx)
        if tgt and tgt_idx and tgt not in gene_to_idx:
            gene_to_idx[tgt] = int(tgt_idx)

print(f"  gene_to_idx: {len(gene_to_idx)} entries", flush=True)
named_in_emb = sorted(gene_to_idx.items(), key=lambda x: x[0])
named_gene_names = [g for g, _ in named_in_emb]
named_gene_idxs = [i for _, i in named_in_emb]
named_emb = emb[:, named_gene_idxs, :]  # [12, n_named, 512]
n_named = named_emb.shape[1]
named_name_to_idx = {g: i for i, g in enumerate(named_gene_names)}
gene_set = set(named_gene_names)
print(f"  Named: {n_named}", flush=True)

# Load GO annotations
print("Loading GO annotations...", flush=True)
with open(GENE2GO_PKL, "rb") as f:
    gene2go_raw = pickle.load(f)

MIN_TERM, MAX_TERM = 3, 60
term2genes_all = defaultdict(set)
for gene, terms in gene2go_raw.items():
    g = gene.strip().upper()
    if g not in gene_set:
        continue
    for entry in terms:
        if isinstance(entry, tuple):
            go_id, aspect = entry[0], entry[1] if len(entry) > 1 else "U"
        else:
            go_id = str(entry)
            aspect = "U"
        if go_id.startswith("GO:"):
            term2genes_all[go_id].add(g)

term2genes_filtered = {t: gs for t, gs in term2genes_all.items()
                       if MIN_TERM <= len(gs) <= MAX_TERM}
print(f"  GO terms (size {MIN_TERM}-{MAX_TERM}): {len(term2genes_filtered)}", flush=True)


# ============================================================
# SVD helper
# ============================================================

def compute_svd_projections(layer_emb):
    """layer_emb: [n_named, 512]. Returns projections [n_named, 3] for SV1, SV2, SV3."""
    mc = layer_emb - layer_emb.mean(axis=0)
    U, S, Vt = np.linalg.svd(mc, full_matrices=False)
    # project onto top 3 singular vectors
    projs = mc @ Vt[:3].T  # [n_named, 3]
    return projs, S


# ============================================================
# H01: GO enrichment of SV3 poles
# ============================================================
print("\n=== H01: GO enrichment of SV3 poles ===", flush=True)

def fisher_pole_enrichment(pole_genes, universe_genes, term2genes):
    """Run Fisher exact test for all terms. Returns list of (term, OR, p, k_term_in_pole, k_term_total)."""
    universe = set(universe_genes)
    pole_set = set(pole_genes)
    N = len(universe)
    results = []
    for term, tgenes in term2genes.items():
        tg = tgenes & universe
        if len(tg) == 0:
            continue
        a = len(pole_set & tg)
        b = len(pole_set) - a
        c = len(tg) - a
        d = N - a - b - c
        _, p = fisher_exact([[a, b], [c, d]], alternative='greater')
        OR = (a * d) / (b * c + 1e-9) if (b * c) > 0 else float('inf')
        results.append((term, OR, p, a, len(tg)))
    results.sort(key=lambda x: x[2])
    return results

h01_results = {"meta": {"K": K_POLE, "n_named": n_named, "n_go_terms": len(term2genes_filtered), "n_shuffle": N_SHUFFLE}, "layers": []}

for layer in range(n_layers):
    layer_emb = named_emb[layer]
    projs, S = compute_svd_projections(layer_emb)
    sv3_proj = projs[:, 2]  # SV3 projection

    top_idx = np.argsort(sv3_proj)[-K_POLE:]
    bot_idx = np.argsort(sv3_proj)[:K_POLE]
    top_genes = [named_gene_names[i] for i in top_idx]
    bot_genes = [named_gene_names[i] for i in bot_idx]

    # Observed Fisher results for top pole and bottom pole
    top_results = fisher_pole_enrichment(top_genes, named_gene_names, term2genes_filtered)
    bot_results = fisher_pole_enrichment(bot_genes, named_gene_names, term2genes_filtered)

    # Get best 5 terms from each pole (lowest p)
    top5_top = [(t, OR, p, a, n) for t, OR, p, a, n in top_results[:5]]
    top5_bot = [(t, OR, p, a, n) for t, OR, p, a, n in bot_results[:5]]

    # Shuffle null for best term in top pole
    if top_results:
        best_term = top_results[0][0]
        best_obs_p = top_results[0][2]
        null_ps = []
        for _ in range(N_SHUFFLE):
            shuffled = list(RNG.choice(named_gene_names, size=K_POLE, replace=False))
            res = fisher_pole_enrichment(shuffled, named_gene_names, {best_term: term2genes_filtered[best_term]})
            null_ps.append(res[0][2] if res else 1.0)
        emp_p_top = np.mean(np.array(null_ps) <= best_obs_p)
    else:
        emp_p_top = 1.0

    if bot_results:
        best_term_bot = bot_results[0][0]
        best_obs_p_bot = bot_results[0][2]
        null_ps_bot = []
        for _ in range(N_SHUFFLE):
            shuffled = list(RNG.choice(named_gene_names, size=K_POLE, replace=False))
            res = fisher_pole_enrichment(shuffled, named_gene_names, {best_term_bot: term2genes_filtered[best_term_bot]})
            null_ps_bot.append(res[0][2] if res else 1.0)
        emp_p_bot = np.mean(np.array(null_ps_bot) <= best_obs_p_bot)
    else:
        emp_p_bot = 1.0

    layer_data = {
        "layer": layer,
        "sv3_variance_fraction": float(S[2]**2 / (S**2).sum()),
        "top_pole_genes": top_genes[:10],  # first 10 for compactness
        "bot_pole_genes": bot_genes[:10],
        "top_pole_top5_terms": [(t, round(OR, 2), round(p, 6), a, n) for t, OR, p, a, n in top5_top],
        "bot_pole_top5_terms": [(t, round(OR, 2), round(p, 6), a, n) for t, OR, p, a, n in top5_bot],
        "top_best_emp_p": round(emp_p_top, 4),
        "bot_best_emp_p": round(emp_p_bot, 4),
    }
    h01_results["layers"].append(layer_data)
    sig_top = "***" if emp_p_top < 0.05 else ""
    sig_bot = "***" if emp_p_bot < 0.05 else ""
    best_top_str = f"{top5_top[0][0]} OR={top5_top[0][1]:.1f} p={top5_top[0][2]:.4f}" if top5_top else "none"
    best_bot_str = f"{top5_bot[0][0]} OR={top5_bot[0][1]:.1f} p={top5_bot[0][2]:.4f}" if top5_bot else "none"
    print(f"  Layer {layer:2d}: TOP={best_top_str} emp_p={emp_p_top:.3f}{sig_top} | BOT={best_bot_str} emp_p={emp_p_bot:.3f}{sig_bot}", flush=True)

with open(ITER_DIR / "h01_sv3_go_enrichment.json", "w") as f:
    json.dump(h01_results, f, indent=2)
print("  H01 saved.", flush=True)


# ============================================================
# H02: Layer-depth Spearman trend on SV2/SV3 z-scores
# ============================================================
print("\n=== H02: Layer-depth Spearman on SV2/SV3 z-scores ===", flush=True)

# Load iter_0013 H01 data
h01_path = PREV_ITER_DIR / "h01_sv_axis_specificity.json"
with open(h01_path) as f:
    sv_axis_data = json.load(f)

layers_list = [r["layer"] for r in sv_axis_data["results"]]
sv2_z = [r["sv2_z"] for r in sv_axis_data["results"]]
sv3_z = [r["sv3_z"] for r in sv_axis_data["results"]]
sv1_z = [r["sv1_z"] for r in sv_axis_data["results"]]

rho_sv2, p_sv2 = spearmanr(layers_list, sv2_z)
rho_sv3, p_sv3 = spearmanr(layers_list, sv3_z)
rho_sv1, p_sv1 = spearmanr(layers_list, sv1_z)

print(f"  SV1: rho={rho_sv1:.3f} p={p_sv1:.4f}", flush=True)
print(f"  SV2: rho={rho_sv2:.3f} p={p_sv2:.4f}", flush=True)
print(f"  SV3: rho={rho_sv3:.3f} p={p_sv3:.4f}", flush=True)

# Also compute per-axis: peak layer and trough layer
sv2_peak_layer = layers_list[int(np.argmax(sv2_z))]
sv3_peak_layer = layers_list[int(np.argmax(sv3_z))]

h02_results = {
    "meta": {"source": "iter_0013/h01_sv_axis_specificity.json"},
    "sv1": {"z_scores": sv1_z, "spearman_rho": round(rho_sv1, 4), "spearman_p": round(p_sv1, 5)},
    "sv2": {"z_scores": sv2_z, "spearman_rho": round(rho_sv2, 4), "spearman_p": round(p_sv2, 5), "peak_layer": sv2_peak_layer},
    "sv3": {"z_scores": sv3_z, "spearman_rho": round(rho_sv3, 4), "spearman_p": round(p_sv3, 5), "peak_layer": sv3_peak_layer},
    "interpretation": "Positive rho=increasing z with depth. Negative rho=decreasing z with depth."
}

with open(ITER_DIR / "h02_layer_depth_trend.json", "w") as f:
    json.dump(h02_results, f, indent=2)
print("  H02 saved.", flush=True)


# ============================================================
# H03: Expanded hub-degree control with STRING score>=0.4
# ============================================================
print("\n=== H03: Hub-degree control STRING score>=0.4 ===", flush=True)

# Load STRING data at score>=0.7 to know which genes we have
with open(STRING_JSON_07) as f:
    string_data_07 = json.load(f)

named_genes_in_string = set(string_data_07.get("named_genes", []))
print(f"  Named genes in STRING (score>=0.7): {len(named_genes_in_string)}", flush=True)

# Download STRING at score>=0.4
import urllib.request

STRING_9606_URL = "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz"
LOCAL_STRING_GZ = ITER_DIR / "9606.protein.links.v12.0.txt.gz"
LOCAL_STRING_TXT = ITER_DIR / "9606.protein.links.v12.0.txt"

# Check if we already have the protein info mapping
STRING_INFO_URL = "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz"
LOCAL_INFO_GZ = ITER_DIR / "9606.protein.info.v12.0.txt.gz"
LOCAL_INFO_TXT = ITER_DIR / "9606.protein.info.v12.0.txt"

# Check if prior iterations have STRING raw data
PREV_LINKS = ITER12_DIR / "9606.protein.links.v12.0.txt"
PREV_INFO = ITER12_DIR / "9606.protein.info.v12.0.txt"

if PREV_LINKS.exists() and PREV_INFO.exists():
    print(f"  Using cached STRING from iter_0012", flush=True)
    STRING_LINKS_PATH = PREV_LINKS
    STRING_INFO_PATH = PREV_INFO
else:
    # Try to find in iter_0011 or earlier
    for prev_iter in ["iter_0011", "iter_0010"]:
        pp = Path(str(ITER12_DIR).replace("iter_0012", prev_iter))
        if (pp / "9606.protein.links.v12.0.txt").exists():
            STRING_LINKS_PATH = pp / "9606.protein.links.v12.0.txt"
            STRING_INFO_PATH = pp / "9606.protein.info.v12.0.txt"
            print(f"  Using cached STRING from {prev_iter}", flush=True)
            break
    else:
        # Download
        print("  Downloading STRING PPI database (may take time)...", flush=True)
        import gzip, shutil
        try:
            urllib.request.urlretrieve(STRING_9606_URL, LOCAL_STRING_GZ)
            with gzip.open(LOCAL_STRING_GZ, 'rb') as f_in:
                with open(LOCAL_STRING_TXT, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            urllib.request.urlretrieve(STRING_INFO_URL, LOCAL_INFO_GZ)
            with gzip.open(LOCAL_INFO_GZ, 'rb') as f_in:
                with open(LOCAL_INFO_TXT, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            STRING_LINKS_PATH = LOCAL_STRING_TXT
            STRING_INFO_PATH = LOCAL_INFO_TXT
            print("  Downloaded STRING OK.", flush=True)
        except Exception as e:
            print(f"  ERROR downloading STRING: {e}", flush=True)
            STRING_LINKS_PATH = None
            STRING_INFO_PATH = None

if STRING_LINKS_PATH is None or STRING_INFO_PATH is None:
    print("  BLOCKED: STRING files unavailable. Skipping H03.", flush=True)
    h03_results = {"error": "STRING files unavailable", "status": "blocked"}
else:
    # Build protein_id -> gene_name mapping
    print("  Building protein->gene mapping...", flush=True)
    prot2gene = {}
    with open(STRING_INFO_PATH) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                prot_id = parts[0]
                gene_name = parts[1].strip().upper()
                prot2gene[prot_id] = gene_name

    # Load STRING edges at score>=0.4 for our named genes
    print("  Parsing STRING edges (score>=400)...", flush=True)
    SCORE_THRESHOLD = 400
    edges_04 = []
    named_prots = {prot: gene for prot, gene in prot2gene.items() if gene in gene_set}
    named_prots_set = set(named_prots.keys())

    with open(STRING_LINKS_PATH) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            p1, p2, score = parts[0], parts[1], int(parts[2])
            if score >= SCORE_THRESHOLD and p1 in named_prots_set and p2 in named_prots_set:
                g1 = named_prots[p1]
                g2 = named_prots[p2]
                if g1 != g2 and g1 in named_name_to_idx and g2 in named_name_to_idx:
                    edges_04.append((g1, g2, score))

    # Deduplicate
    edge_set = set()
    edges_04_dedup = []
    for g1, g2, s in edges_04:
        key = (min(g1,g2), max(g1,g2))
        if key not in edge_set:
            edge_set.add(key)
            edges_04_dedup.append((g1, g2, s))

    print(f"  STRING edges (score>=400) between named genes: {len(edges_04_dedup)}", flush=True)

    # Compute degree for each named gene in the 0.4 network
    degree = defaultdict(int)
    for g1, g2, _ in edges_04_dedup:
        degree[g1] += 1
        degree[g2] += 1

    degrees_list = [degree.get(g, 0) for g in named_gene_names]
    median_degree = np.median([d for d in degrees_list if d > 0])
    print(f"  Median degree (non-zero genes): {median_degree}", flush=True)

    # Split edges into hub and non-hub
    hub_edges = [(g1, g2) for g1, g2, s in edges_04_dedup if degree[g1] > median_degree and degree[g2] > median_degree]
    nonhub_edges = [(g1, g2) for g1, g2, s in edges_04_dedup if degree[g1] <= median_degree and degree[g2] <= median_degree]
    print(f"  Hub edges: {len(hub_edges)}, Non-hub edges: {len(nonhub_edges)}", flush=True)

    def compute_copole_rate_from_edges(edges, named_names, projs_k, K):
        """For given edges, compute fraction that co-appear in top-K or bottom-K pole."""
        top_set = set(np.argsort(projs_k)[-K:])
        bot_set = set(np.argsort(projs_k)[:K])
        n_copole = 0
        n_valid = 0
        for g1, g2 in edges:
            i1 = named_name_to_idx.get(g1)
            i2 = named_name_to_idx.get(g2)
            if i1 is None or i2 is None:
                continue
            n_valid += 1
            if (i1 in top_set and i2 in top_set) or (i1 in bot_set and i2 in bot_set):
                n_copole += 1
        return n_copole / n_valid if n_valid > 0 else 0.0, n_valid

    def run_copole_test(edges, layer_emb, axis_idx, K, n_shuf=300):
        """Run co-pole enrichment test for given edges, returns obs_rate, z, emp_p."""
        mc = layer_emb - layer_emb.mean(axis=0)
        U, S, Vt = np.linalg.svd(mc, full_matrices=False)
        projs = mc @ Vt[axis_idx]  # [n_named]
        obs_rate, n_valid = compute_copole_rate_from_edges(edges, named_gene_names, projs, K)
        if n_valid == 0:
            return 0, 0, 1.0, 0

        # Shuffle null: randomly permute gene labels
        null_rates = []
        for _ in range(n_shuf):
            perm = RNG.permutation(len(named_gene_names))
            perm_projs = projs[perm]
            null_rate, _ = compute_copole_rate_from_edges(edges, named_gene_names, perm_projs, K)
            null_rates.append(null_rate)
        null_arr = np.array(null_rates)
        z = (obs_rate - null_arr.mean()) / (null_arr.std() + 1e-9)
        emp_p = np.mean(null_arr >= obs_rate)
        return obs_rate, z, emp_p, n_valid

    h03_layer_results = []
    for layer in range(n_layers):
        layer_emb = named_emb[layer]
        hub_obs, hub_z, hub_emp_p, hub_n = run_copole_test(hub_edges, layer_emb, 1, K_POLE, N_SHUFFLE)
        nh_obs, nh_z, nh_emp_p, nh_n = run_copole_test(nonhub_edges, layer_emb, 1, K_POLE, N_SHUFFLE)
        h03_layer_results.append({
            "layer": layer,
            "hub_copole_obs": round(hub_obs, 4),
            "hub_z": round(hub_z, 3),
            "hub_emp_p": round(hub_emp_p, 4),
            "hub_n": hub_n,
            "nonhub_copole_obs": round(nh_obs, 4),
            "nonhub_z": round(nh_z, 3),
            "nonhub_emp_p": round(nh_emp_p, 4),
            "nonhub_n": nh_n,
        })
        sig_hub = "*" if hub_emp_p < 0.05 else ""
        sig_nh = "*" if nh_emp_p < 0.05 else ""
        print(f"  Layer {layer:2d}: hub z={hub_z:.2f}{sig_hub} (N={hub_n}) | nonhub z={nh_z:.2f}{sig_nh} (N={nh_n})", flush=True)

    hub_z_mean = np.mean([r["hub_z"] for r in h03_layer_results])
    nh_z_mean = np.mean([r["nonhub_z"] for r in h03_layer_results])
    hub_n_sig = sum(1 for r in h03_layer_results if r["hub_emp_p"] < 0.05)
    nh_n_sig = sum(1 for r in h03_layer_results if r["nonhub_emp_p"] < 0.05)

    h03_results = {
        "meta": {
            "score_threshold": SCORE_THRESHOLD,
            "K": K_POLE,
            "n_edges_total": len(edges_04_dedup),
            "n_hub_edges": len(hub_edges),
            "n_nonhub_edges": len(nonhub_edges),
            "median_degree": float(median_degree),
            "n_shuffle": N_SHUFFLE,
        },
        "summary": {
            "hub_mean_z": round(hub_z_mean, 3),
            "hub_n_sig_layers": hub_n_sig,
            "nonhub_mean_z": round(nh_z_mean, 3),
            "nonhub_n_sig_layers": nh_n_sig,
        },
        "layers": h03_layer_results
    }

    print(f"  Hub mean z={hub_z_mean:.3f} ({hub_n_sig}/12 sig)", flush=True)
    print(f"  Non-hub mean z={nh_z_mean:.3f} ({nh_n_sig}/12 sig)", flush=True)

with open(ITER_DIR / "h03_hub_degree_control_04.json", "w") as f:
    json.dump(h03_results, f, indent=2)
print("  H03 saved.", flush=True)

print("\n=== ALL HYPOTHESES COMPLETE ===", flush=True)
