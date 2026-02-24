"""
iter_0013 - Multi-Hypothesis Screen

H01 (SV axis specificity): Run STRING PPI co-pole enrichment on SV1, SV2, SV3 for all 12 layers.
    Tests whether SV2 is uniquely or non-uniquely the 'interaction geometry axis'.
    Uses cached STRING PPI data and layer embeddings.

H02 (Repression anti-pole / cross-pole): For TRRUST activation vs repression edges:
    Compute cross-pole rate (gene i in top pole, gene j in bottom pole OR vice versa).
    Test if repression pairs show elevated cross-pole rates vs null (opposite of activation).
    Novel: tests if SV2 encodes *regulatory sign* (act=same, rep=opposite).

H03 (Hub-degree control): Regress out STRING degree from co-pole rate.
    Check if PPI co-pole enrichment is driven by high-degree hubs.
    Compare co-pole enrichment for hub edges vs non-hub edges.

Data: layer_gene_embeddings.npy [12, 4803, 512], string_ppi_named_genes.json, trrust_human.tsv
"""

import numpy as np
import json
import csv
import sys
from pathlib import Path
from collections import defaultdict

ITER_DIR = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0013"
)
PREV_ITER_DIR = Path(
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
TRRUST_PATH = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/single_cell_mechinterp/external/networks/trrust_human.tsv"
)
STRING_JSON = PREV_ITER_DIR / "string_ppi_named_genes.json"

ITER_DIR.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(42)
N_SHUFFLE = 500
K_POLE = 52  # top/bottom K for SV poles


# ============================================================
# Data loading
# ============================================================

print("Loading embeddings...", flush=True)
emb = np.load(EMB_PATH)  # [12, 4803, 512]
n_layers, n_genes_total, n_dim = emb.shape
print(f"  Embeddings: {emb.shape}", flush=True)

# Load named gene set from edge dataset using source/target + source_idx/target_idx columns
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

# Find named gene positions in embedding
named_in_emb = sorted(gene_to_idx.items(), key=lambda x: x[0])  # list of (name, idx)
print(f"  Named genes with embedding positions: {len(named_in_emb)}", flush=True)

# Extract named gene embeddings: [n_layers, n_named, n_dim]
named_gene_names = [g for g, _ in named_in_emb]
named_gene_idxs = [i for _, i in named_in_emb]
named_emb = emb[:, named_gene_idxs, :]  # [12, n_named, 512]
print(f"  Named embedding shape: {named_emb.shape}", flush=True)
n_named = named_emb.shape[1]

# Gene name to named-index mapping
named_name_to_idx = {g: i for i, g in enumerate(named_gene_names)}


# ============================================================
# SVD helper: compute SV projections for all axes
# ============================================================

def compute_sv_projections(layer_emb, n_components=4):
    """Center and SVD; return projections onto SV1..SVn."""
    centered = layer_emb - layer_emb.mean(axis=0)
    U, S, Vt = np.linalg.svd(centered, full_matrices=False)
    # projections[k] = score on (k+1)-th singular vector [n_genes]
    projections = [U[:, k] * S[k] for k in range(n_components)]
    return projections, S


# ============================================================
# Load STRING PPI data
# ============================================================

print("Loading STRING PPI...", flush=True)
with open(STRING_JSON) as f:
    string_edges_raw = json.load(f)

SCORE_THRESH = 0.7
string_pairs = []
for edge in string_edges_raw:
    if float(edge.get("score", 0)) >= SCORE_THRESH:
        gA = edge.get("preferredName_A", "")
        gB = edge.get("preferredName_B", "")
        if gA in named_name_to_idx and gB in named_name_to_idx:
            ia = named_name_to_idx[gA]
            ib = named_name_to_idx[gB]
            if ia != ib:
                string_pairs.append((ia, ib, gA, gB, float(edge["score"])))

print(f"  STRING pairs (score>=0.7, both named): {len(string_pairs)}", flush=True)


# ============================================================
# Load TRRUST data
# ============================================================

print("Loading TRRUST...", flush=True)
act_pairs = []
rep_pairs = []
with open(TRRUST_PATH) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        tf, target, effect = parts[0], parts[1], parts[2]
        if tf in named_name_to_idx and target in named_name_to_idx:
            i, j = named_name_to_idx[tf], named_name_to_idx[target]
            if i != j:
                if "Activation" in effect:
                    act_pairs.append((i, j))
                elif "Repression" in effect:
                    rep_pairs.append((i, j))

print(f"  Activation pairs: {len(act_pairs)}, Repression pairs: {len(rep_pairs)}", flush=True)


# ============================================================
# Utility: co-pole and cross-pole rates
# ============================================================

def copole_rate(pairs, proj, K):
    """Fraction of pairs where both genes in same pole (top-K or bottom-K)."""
    n = len(proj)
    sorted_idx = np.argsort(proj)
    top_set = set(sorted_idx[-K:])
    bot_set = set(sorted_idx[:K])
    count = 0
    for (i, j) in pairs:
        same = (i in top_set and j in top_set) or (i in bot_set and j in bot_set)
        if same:
            count += 1
    return count / len(pairs) if pairs else 0.0


def crosspole_rate(pairs, proj, K):
    """Fraction of pairs where genes in OPPOSITE poles (one top-K, one bottom-K)."""
    n = len(proj)
    sorted_idx = np.argsort(proj)
    top_set = set(sorted_idx[-K:])
    bot_set = set(sorted_idx[:K])
    count = 0
    for (i, j) in pairs:
        cross = (i in top_set and j in bot_set) or (i in bot_set and j in top_set)
        if cross:
            count += 1
    return count / len(pairs) if pairs else 0.0


def null_copole_rates(n_named, n_pairs, proj, K, n_shuffle=500, rng=None):
    if rng is None:
        rng = np.random.default_rng(0)
    sorted_idx = np.argsort(proj)
    top_set = set(sorted_idx[-K:])
    bot_set = set(sorted_idx[:K])
    rates = []
    for _ in range(n_shuffle):
        idxs_a = rng.integers(0, n_named, size=n_pairs)
        idxs_b = rng.integers(0, n_named, size=n_pairs)
        count = sum(
            (a in top_set and b in top_set) or (a in bot_set and b in bot_set)
            for a, b in zip(idxs_a, idxs_b)
        )
        rates.append(count / n_pairs)
    return np.array(rates)


def null_crosspole_rates(n_named, n_pairs, proj, K, n_shuffle=500, rng=None):
    if rng is None:
        rng = np.random.default_rng(0)
    sorted_idx = np.argsort(proj)
    top_set = set(sorted_idx[-K:])
    bot_set = set(sorted_idx[:K])
    rates = []
    for _ in range(n_shuffle):
        idxs_a = rng.integers(0, n_named, size=n_pairs)
        idxs_b = rng.integers(0, n_named, size=n_pairs)
        count = sum(
            (a in top_set and b in bot_set) or (a in bot_set and b in top_set)
            for a, b in zip(idxs_a, idxs_b)
        )
        rates.append(count / n_pairs)
    return np.array(rates)


# ============================================================
# H01: SV axis specificity — STRING PPI co-pole on SV1, SV2, SV3
# ============================================================

print("\n=== H01: SV axis specificity (STRING co-pole on SV1/SV2/SV3) ===", flush=True)

string_idx_pairs = [(i, j) for (i, j, *_) in string_pairs]
n_string = len(string_idx_pairs)
h01_results = []

for layer in range(n_layers):
    projs, S = compute_sv_projections(named_emb[layer], n_components=4)
    row = {"layer": layer}
    for sv_idx in range(3):  # SV1=0, SV2=1, SV3=2
        proj = projs[sv_idx]
        obs_rate = copole_rate(string_idx_pairs, proj, K_POLE)
        null_rates = null_copole_rates(n_named, n_string, proj, K_POLE, N_SHUFFLE, RNG)
        z = (obs_rate - null_rates.mean()) / (null_rates.std() + 1e-9)
        p_emp = (null_rates >= obs_rate).mean()
        row[f"sv{sv_idx+1}_obs"] = round(float(obs_rate), 4)
        row[f"sv{sv_idx+1}_null_mean"] = round(float(null_rates.mean()), 4)
        row[f"sv{sv_idx+1}_z"] = round(float(z), 3)
        row[f"sv{sv_idx+1}_emp_p"] = round(float(p_emp), 4)
    h01_results.append(row)
    print(f"  Layer {layer:2d}: SV1 z={row['sv1_z']:.2f}  SV2 z={row['sv2_z']:.2f}  SV3 z={row['sv3_z']:.2f}", flush=True)

with open(ITER_DIR / "h01_sv_axis_specificity.json", "w") as f:
    json.dump({"meta": {"method": "STRING_copole_SV1_SV2_SV3", "K": K_POLE, "n_string_pairs": n_string,
                        "n_shuffle": N_SHUFFLE}, "results": h01_results}, f, indent=2)
print("  Saved h01_sv_axis_specificity.json", flush=True)

# Summary: count significant layers per SV axis
for sv_idx in range(3):
    sig = sum(1 for r in h01_results if r[f"sv{sv_idx+1}_emp_p"] < 0.05)
    mean_z = np.mean([r[f"sv{sv_idx+1}_z"] for r in h01_results])
    print(f"  SV{sv_idx+1}: {sig}/12 layers significant, mean z={mean_z:.2f}", flush=True)


# ============================================================
# H02: Repression anti-pole — cross-pole rates for act vs rep
# ============================================================

print("\n=== H02: Repression anti-pole (cross-pole + co-pole for act/rep) ===", flush=True)

# Use SV2 (sv_idx=1) as the main axis based on prior evidence
h02_results = []
for layer in range(n_layers):
    projs, S = compute_sv_projections(named_emb[layer], n_components=3)
    proj_sv2 = projs[1]

    # Activation co-pole
    act_cop = copole_rate(act_pairs, proj_sv2, K_POLE)
    null_act_cop = null_copole_rates(n_named, len(act_pairs), proj_sv2, K_POLE, N_SHUFFLE, RNG)
    z_act_cop = (act_cop - null_act_cop.mean()) / (null_act_cop.std() + 1e-9)

    # Repression co-pole
    rep_cop = copole_rate(rep_pairs, proj_sv2, K_POLE)
    null_rep_cop = null_copole_rates(n_named, len(rep_pairs), proj_sv2, K_POLE, N_SHUFFLE, RNG)
    z_rep_cop = (rep_cop - null_rep_cop.mean()) / (null_rep_cop.std() + 1e-9)

    # Repression cross-pole (anti-pole)
    rep_xpole = crosspole_rate(rep_pairs, proj_sv2, K_POLE)
    null_rep_xpole = null_crosspole_rates(n_named, len(rep_pairs), proj_sv2, K_POLE, N_SHUFFLE, RNG)
    z_rep_xpole = (rep_xpole - null_rep_xpole.mean()) / (null_rep_xpole.std() + 1e-9)

    # Activation cross-pole (control)
    act_xpole = crosspole_rate(act_pairs, proj_sv2, K_POLE)
    null_act_xpole = null_crosspole_rates(n_named, len(act_pairs), proj_sv2, K_POLE, N_SHUFFLE, RNG)
    z_act_xpole = (act_xpole - null_act_xpole.mean()) / (null_act_xpole.std() + 1e-9)

    row = {
        "layer": layer,
        "act_copole_obs": round(float(act_cop), 4),
        "act_copole_z": round(float(z_act_cop), 3),
        "rep_copole_obs": round(float(rep_cop), 4),
        "rep_copole_z": round(float(z_rep_cop), 3),
        "rep_xpole_obs": round(float(rep_xpole), 4),
        "rep_xpole_z": round(float(z_rep_xpole), 3),
        "act_xpole_obs": round(float(act_xpole), 4),
        "act_xpole_z": round(float(z_act_xpole), 3),
        "rep_xpole_null_mean": round(float(null_rep_xpole.mean()), 4),
    }
    h02_results.append(row)
    print(f"  Layer {layer:2d}: act_cop z={z_act_cop:.2f}  rep_cop z={z_rep_cop:.2f}  rep_xpole z={z_rep_xpole:.2f}  act_xpole z={z_act_xpole:.2f}", flush=True)

with open(ITER_DIR / "h02_repression_antipole.json", "w") as f:
    json.dump({"meta": {"method": "SV2_copole_crosspole_act_rep", "K": K_POLE,
                        "n_act": len(act_pairs), "n_rep": len(rep_pairs),
                        "n_shuffle": N_SHUFFLE}, "results": h02_results}, f, indent=2)
print("  Saved h02_repression_antipole.json", flush=True)

# Summaries
for key, label in [("act_copole_z", "ACT co-pole"), ("rep_copole_z", "REP co-pole"),
                    ("rep_xpole_z", "REP cross-pole"), ("act_xpole_z", "ACT cross-pole")]:
    vals = [r[key] for r in h02_results]
    print(f"  {label}: mean z={np.mean(vals):.2f}, n_pos_z>1.5={sum(v>1.5 for v in vals)}/12", flush=True)


# ============================================================
# H03: Hub-degree control for STRING PPI co-pole
# ============================================================

print("\n=== H03: Hub-degree control for STRING PPI co-pole ===", flush=True)

# Compute degree for each named gene in STRING PPI
degree = defaultdict(int)
for (i, j, gA, gB, score) in string_pairs:
    degree[i] += 1
    degree[j] += 1

# Compute median degree to split hub vs non-hub edges
all_degrees = [degree.get(i, 0) for (i, j, *_) in string_pairs] + \
              [degree.get(j, 0) for (i, j, *_) in string_pairs]
# An edge is "hub" if either gene has degree > median
med_degree = np.median([degree[g] for g in degree if degree[g] > 0])
print(f"  Median STRING degree (named): {med_degree:.1f}", flush=True)

hub_pairs = [(i, j) for (i, j, gA, gB, score) in string_pairs
             if degree.get(i, 0) > med_degree or degree.get(j, 0) > med_degree]
nonhub_pairs = [(i, j) for (i, j, gA, gB, score) in string_pairs
                if degree.get(i, 0) <= med_degree and degree.get(j, 0) <= med_degree]
print(f"  Hub edges: {len(hub_pairs)}, Non-hub edges: {len(nonhub_pairs)}", flush=True)

h03_results = []
for layer in range(n_layers):
    projs, S = compute_sv_projections(named_emb[layer], n_components=3)
    proj_sv2 = projs[1]

    # All STRING
    all_rate = copole_rate(string_idx_pairs, proj_sv2, K_POLE)
    null_all = null_copole_rates(n_named, len(string_idx_pairs), proj_sv2, K_POLE, 300, RNG)
    z_all = (all_rate - null_all.mean()) / (null_all.std() + 1e-9)

    # Hub
    hub_rate = copole_rate(hub_pairs, proj_sv2, K_POLE) if hub_pairs else 0.0
    null_hub = null_copole_rates(n_named, max(1, len(hub_pairs)), proj_sv2, K_POLE, 300, RNG) if hub_pairs else np.array([0.0])
    z_hub = (hub_rate - null_hub.mean()) / (null_hub.std() + 1e-9) if hub_pairs else 0.0

    # Non-hub
    nhub_rate = copole_rate(nonhub_pairs, proj_sv2, K_POLE) if nonhub_pairs else 0.0
    null_nhub = null_copole_rates(n_named, max(1, len(nonhub_pairs)), proj_sv2, K_POLE, 300, RNG) if nonhub_pairs else np.array([0.0])
    z_nhub = (nhub_rate - null_nhub.mean()) / (null_nhub.std() + 1e-9) if nonhub_pairs else 0.0

    row = {
        "layer": layer,
        "all_copole_obs": round(float(all_rate), 4),
        "all_copole_z": round(float(z_all), 3),
        "hub_copole_obs": round(float(hub_rate), 4),
        "hub_copole_z": round(float(z_hub), 3),
        "nonhub_copole_obs": round(float(nhub_rate), 4),
        "nonhub_copole_z": round(float(z_nhub), 3),
    }
    h03_results.append(row)
    print(f"  Layer {layer:2d}: all z={z_all:.2f}  hub z={z_hub:.2f}  non-hub z={z_nhub:.2f}", flush=True)

with open(ITER_DIR / "h03_hub_degree_control.json", "w") as f:
    json.dump({"meta": {"method": "STRING_copole_hub_vs_nonhub_SV2", "K": K_POLE,
                        "median_degree": float(med_degree),
                        "n_hub_pairs": len(hub_pairs), "n_nonhub_pairs": len(nonhub_pairs),
                        "n_shuffle": 300}, "results": h03_results}, f, indent=2)
print("  Saved h03_hub_degree_control.json", flush=True)

# Summary
for key, label in [("all_copole_z", "ALL"), ("hub_copole_z", "HUB"), ("nonhub_copole_z", "NON-HUB")]:
    vals = [r[key] for r in h03_results]
    sig = sum(1 for v in vals if v > 1.96)
    print(f"  {label}: mean z={np.mean(vals):.2f}, n_sig(z>1.96)={sig}/12", flush=True)


# ============================================================
# Consolidated results JSON
# ============================================================

results_summary = {
    "iteration": "iter_0013",
    "H01_sv_axis_specificity": {
        "sv1_mean_z": round(float(np.mean([r["sv1_z"] for r in h01_results])), 3),
        "sv2_mean_z": round(float(np.mean([r["sv2_z"] for r in h01_results])), 3),
        "sv3_mean_z": round(float(np.mean([r["sv3_z"] for r in h01_results])), 3),
        "sv1_sig_layers": sum(1 for r in h01_results if r["sv1_emp_p"] < 0.05),
        "sv2_sig_layers": sum(1 for r in h01_results if r["sv2_emp_p"] < 0.05),
        "sv3_sig_layers": sum(1 for r in h01_results if r["sv3_emp_p"] < 0.05),
    },
    "H02_repression_antipole": {
        "act_copole_mean_z": round(float(np.mean([r["act_copole_z"] for r in h02_results])), 3),
        "rep_copole_mean_z": round(float(np.mean([r["rep_copole_z"] for r in h02_results])), 3),
        "rep_xpole_mean_z": round(float(np.mean([r["rep_xpole_z"] for r in h02_results])), 3),
        "act_xpole_mean_z": round(float(np.mean([r["act_xpole_z"] for r in h02_results])), 3),
        "rep_xpole_sig_layers": sum(1 for r in h02_results if r["rep_xpole_z"] > 1.96),
    },
    "H03_hub_degree_control": {
        "all_mean_z": round(float(np.mean([r["all_copole_z"] for r in h03_results])), 3),
        "hub_mean_z": round(float(np.mean([r["hub_copole_z"] for r in h03_results])), 3),
        "nonhub_mean_z": round(float(np.mean([r["nonhub_copole_z"] for r in h03_results])), 3),
        "nonhub_sig_layers": sum(1 for r in h03_results if r["nonhub_copole_z"] > 1.96),
        "hub_sig_layers": sum(1 for r in h03_results if r["hub_copole_z"] > 1.96),
        "median_degree": float(med_degree),
    }
}

with open(ITER_DIR / "iter0013_results.json", "w") as f:
    json.dump(results_summary, f, indent=2)

print("\n=== SUMMARY ===", flush=True)
print(json.dumps(results_summary, indent=2), flush=True)
print("\nDone.", flush=True)
