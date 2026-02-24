"""
iter_0012 - Multi-Hypothesis Screen

H01 (H-A): Bootstrap CI on activation vs repression co-pole differential.
     Data in hand from iter_0011. Bootstrap-resample activation edges and repression edges
     separately at each layer. Compute co-pole rate CI. Report differential CI and test
     whether act co-pole > rep co-pole with uncertainty quantified.
     Novel: converts qualitative observation to CI-backed claim.

H02 (H-D): SV2 mean pairwise distance for activation pairs vs null (spatial concentration).
     For each of 12 layers, compute SV2 projection for all 209 genes.
     For TRRUST activation edges: compute mean |SV2_i - SV2_j| for each pair.
     Null: same-size random pairs (N=1000 shuffles).
     Tests whether activation pairs are point-clustered in SV2 space (not just pole-separated).
     Novel: complementary to co-pole; tests continuous distance not binary pole membership.

H03 (H-B): STRING PPI co-pole screen.
     Download STRING v12.0 PPI interactions for our 209 named genes via STRING API.
     Filter to combined_score >= 700 (high-confidence).
     Run co-pole enrichment test (same pipeline as TRRUST, N=500 shuffles).
     Novel: cross-graph-type validation. Tests if generic PPI proximity maps to SV2 geometry.

Data: layer_gene_embeddings.npy [12, 4803, 512], cycle1_edge_dataset.tsv (209 named genes)
TRRUST: trrust_human.tsv (9396 edges, activation/repression annotations)
STRING: downloaded from API (v12.0, human 9606)
"""

import numpy as np
import json
import csv
import sys
import urllib.request
import urllib.parse
from pathlib import Path
from collections import defaultdict

ITER_DIR = Path(
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

ITER_DIR.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(42)
N_BOOT = 2000    # bootstrap samples for CIs
N_SHUFFLE = 500  # null shuffles for pairwise distance test
TOP_K = 52       # ~25% of 209 genes

print("Loading data...", flush=True)
emb = np.load(EMB_PATH)  # [12, 4803, 512]
print(f"  Embeddings: {emb.shape}", flush=True)

# Load 209 named genes
gene2emb_idx = {}
with open(EDGE_PATH) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        src = row["source"].strip().upper()
        tgt = row["target"].strip().upper()
        if src not in gene2emb_idx:
            gene2emb_idx[src] = int(row["source_idx"])
        if tgt not in gene2emb_idx:
            gene2emb_idx[tgt] = int(row["target_idx"])
named_gene_list = sorted(gene2emb_idx.keys())
named_set = set(named_gene_list)
print(f"  Named genes: {len(named_gene_list)}", flush=True)

# Load TRRUST
trrust_act_pairs = []
trrust_rep_pairs = []
with open(TRRUST_PATH) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 3:
            tf, target, mode = parts[0].upper(), parts[1].upper(), parts[2].upper()
            if tf in named_set and target in named_set:
                if "ACTIVATION" in mode:
                    trrust_act_pairs.append((tf, target))
                elif "REPRESSION" in mode:
                    trrust_rep_pairs.append((tf, target))

print(f"  TRRUST activation pairs (both named): {len(trrust_act_pairs)}", flush=True)
print(f"  TRRUST repression pairs (both named): {len(trrust_rep_pairs)}", flush=True)


# ==========================================
# Helper: SVD of named-gene embedding at layer L
# ==========================================
def get_sv_projections(layer_idx, sv_idx=1):
    """Returns SV[sv_idx] projections for all named genes at given layer."""
    mat = emb[layer_idx]
    gene_names_local = list(gene2emb_idx.keys())
    idxs = [gene2emb_idx[g] for g in gene_names_local]
    X = mat[idxs]
    X_centered = X - X.mean(axis=0)
    U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)
    projs = U[:, sv_idx] * S[sv_idx]
    return {g: float(projs[i]) for i, g in enumerate(gene_names_local)}


def get_poles(proj_dict, top_k=TOP_K):
    sorted_genes = sorted(proj_dict.items(), key=lambda x: x[1], reverse=True)
    top_pole = set(g for g, _ in sorted_genes[:top_k])
    bot_pole = set(g for g, _ in sorted_genes[-top_k:])
    return top_pole, bot_pole


def copole_rate(pairs, top_pole, bot_pole):
    """Fraction of pairs that are co-localized in same pole."""
    if not pairs:
        return np.nan
    n_copole = sum(
        1 for (a, b) in pairs
        if (a in top_pole and b in top_pole) or (a in bot_pole and b in bot_pole)
    )
    return n_copole / len(pairs)


# ==========================================
# H01: Bootstrap CI on activation vs repression co-pole differential
# ==========================================
print("\n=== H01: Bootstrap CI on act vs rep co-pole differential ===", flush=True)

h01_results = []

for layer in range(12):
    proj = get_sv_projections(layer, sv_idx=1)
    top_pole, bot_pole = get_poles(proj)

    # Observed rates
    obs_act = copole_rate(trrust_act_pairs, top_pole, bot_pole)
    obs_rep = copole_rate(trrust_rep_pairs, top_pole, bot_pole)
    obs_diff = obs_act - obs_rep

    # Bootstrap: resample pairs WITH REPLACEMENT
    act_arr = np.array(trrust_act_pairs)
    rep_arr = np.array(trrust_rep_pairs)

    boot_diffs = []
    boot_act_rates = []
    boot_rep_rates = []

    for _ in range(N_BOOT):
        # Resample pairs (with replacement)
        b_act_idx = RNG.integers(0, len(act_arr), size=len(act_arr))
        b_rep_idx = RNG.integers(0, len(rep_arr), size=len(rep_arr))
        b_act_pairs = [tuple(act_arr[i]) for i in b_act_idx]
        b_rep_pairs = [tuple(rep_arr[i]) for i in b_rep_idx]

        r_act = copole_rate(b_act_pairs, top_pole, bot_pole)
        r_rep = copole_rate(b_rep_pairs, top_pole, bot_pole)
        boot_act_rates.append(r_act)
        boot_rep_rates.append(r_rep)
        boot_diffs.append(r_act - r_rep)

    boot_diffs = np.array(boot_diffs)
    boot_act_rates = np.array(boot_act_rates)
    boot_rep_rates = np.array(boot_rep_rates)

    # 95% CI via percentile bootstrap
    ci_diff_lo, ci_diff_hi = np.percentile(boot_diffs, [2.5, 97.5])
    ci_act_lo, ci_act_hi = np.percentile(boot_act_rates, [2.5, 97.5])
    ci_rep_lo, ci_rep_hi = np.percentile(boot_rep_rates, [2.5, 97.5])

    # p-value: fraction of bootstrap diffs <= 0
    p_one_sided = float(np.mean(boot_diffs <= 0))

    row = {
        "layer": layer,
        "obs_act": round(obs_act, 4),
        "obs_rep": round(obs_rep, 4),
        "obs_diff": round(obs_diff, 4),
        "ci_act_lo": round(ci_act_lo, 4),
        "ci_act_hi": round(ci_act_hi, 4),
        "ci_rep_lo": round(ci_rep_lo, 4),
        "ci_rep_hi": round(ci_rep_hi, 4),
        "ci_diff_lo": round(ci_diff_lo, 4),
        "ci_diff_hi": round(ci_diff_hi, 4),
        "p_act_gt_rep": round(p_one_sided, 4),  # fraction boot <= 0 (want low)
        "n_act": len(trrust_act_pairs),
        "n_rep": len(trrust_rep_pairs),
        "diff_ci_excludes_zero": bool(ci_diff_lo > 0),
    }
    h01_results.append(row)

    print(
        f"  L{layer:02d}: act={obs_act:.3f} [{ci_act_lo:.3f},{ci_act_hi:.3f}]  "
        f"rep={obs_rep:.3f} [{ci_rep_lo:.3f},{ci_rep_hi:.3f}]  "
        f"diff={obs_diff:.3f} [{ci_diff_lo:.3f},{ci_diff_hi:.3f}]  "
        f"p={p_one_sided:.3f}  CI_excl_0={row['diff_ci_excludes_zero']}",
        flush=True,
    )

# Summary
n_ci_excl = sum(1 for r in h01_results if r["diff_ci_excludes_zero"])
n_p05 = sum(1 for r in h01_results if r["p_act_gt_rep"] < 0.05)
print(f"\n  H01 Summary: {n_ci_excl}/12 layers CI excludes 0, {n_p05}/12 layers p<0.05", flush=True)

out_h01 = ITER_DIR / "h01_bootstrap_act_rep_diff.json"
with open(out_h01, "w") as f:
    json.dump({"results": h01_results,
               "n_layers_ci_excludes_zero": n_ci_excl,
               "n_layers_p_lt_05": n_p05,
               "n_boot": N_BOOT}, f, indent=2)
print(f"  Saved: {out_h01}", flush=True)


# ==========================================
# H02: Mean pairwise SV2 distance — spatial concentration test
# ==========================================
print("\n=== H02: Mean pairwise SV2 distance for activation/repression pairs vs null ===", flush=True)

h02_results = []
gene_list_fixed = named_gene_list  # sorted, consistent order

for layer in range(12):
    proj = get_sv_projections(layer, sv_idx=1)
    proj_arr = np.array([proj[g] for g in gene_list_fixed])  # [N_named]
    proj_dict = {g: proj[g] for g in gene_list_fixed}

    # Compute mean pairwise |SV2_i - SV2_j| for act and rep pairs
    def mean_abs_diff(pairs, pd):
        if not pairs:
            return np.nan
        diffs = [abs(pd[a] - pd[b]) for (a, b) in pairs if a in pd and b in pd]
        return float(np.mean(diffs)) if diffs else np.nan

    obs_act_dist = mean_abs_diff(trrust_act_pairs, proj_dict)
    obs_rep_dist = mean_abs_diff(trrust_rep_pairs, proj_dict)

    # Null: shuffle gene labels N_SHUFFLE times, compute mean pairwise dist
    n_act = len(trrust_act_pairs)
    n_rep = len(trrust_rep_pairs)
    null_act_dists = []
    null_rep_dists = []

    for _ in range(N_SHUFFLE):
        perm = RNG.permutation(len(gene_list_fixed))
        perm_dict = {gene_list_fixed[i]: proj_arr[perm[i]] for i in range(len(gene_list_fixed))}
        null_act_dists.append(mean_abs_diff(trrust_act_pairs, perm_dict))
        null_rep_dists.append(mean_abs_diff(trrust_rep_pairs, perm_dict))

    null_act_dists = np.array(null_act_dists)
    null_rep_dists = np.array(null_rep_dists)

    # Empirical p: fraction of null values <= obs (lower obs = more concentrated)
    emp_p_act = float(np.mean(null_act_dists <= obs_act_dist))
    emp_p_rep = float(np.mean(null_rep_dists <= obs_rep_dist))

    # Z-score (lower distance = more concentrated = better)
    z_act = float((obs_act_dist - null_act_dists.mean()) / (null_act_dists.std() + 1e-9))
    z_rep = float((obs_rep_dist - null_rep_dists.mean()) / (null_rep_dists.std() + 1e-9))

    row = {
        "layer": layer,
        "obs_act_mean_dist": round(obs_act_dist, 4),
        "null_act_mean": round(float(null_act_dists.mean()), 4),
        "null_act_std": round(float(null_act_dists.std()), 4),
        "z_act": round(z_act, 3),
        "emp_p_act_concentrated": round(emp_p_act, 4),
        "obs_rep_mean_dist": round(obs_rep_dist, 4),
        "null_rep_mean": round(float(null_rep_dists.mean()), 4),
        "null_rep_std": round(float(null_rep_dists.std()), 4),
        "z_rep": round(z_rep, 3),
        "emp_p_rep_concentrated": round(emp_p_rep, 4),
        "act_more_concentrated_than_rep": bool(obs_act_dist < obs_rep_dist),
    }
    h02_results.append(row)

    print(
        f"  L{layer:02d}: act_dist={obs_act_dist:.3f} (null={null_act_dists.mean():.3f}, z={z_act:.2f}, p={emp_p_act:.3f})  "
        f"rep_dist={obs_rep_dist:.3f} (null={null_rep_dists.mean():.3f}, z={z_rep:.2f}, p={emp_p_rep:.3f})",
        flush=True,
    )

n_act_sig = sum(1 for r in h02_results if r["emp_p_act_concentrated"] < 0.05)
n_rep_sig = sum(1 for r in h02_results if r["emp_p_rep_concentrated"] < 0.05)
n_act_more = sum(1 for r in h02_results if r["act_more_concentrated_than_rep"])
print(f"\n  H02 Summary: act concentrated p<0.05 at {n_act_sig}/12 layers, "
      f"rep concentrated p<0.05 at {n_rep_sig}/12 layers, "
      f"act < rep dist at {n_act_more}/12 layers", flush=True)

out_h02 = ITER_DIR / "h02_sv2_pairwise_dist.json"
with open(out_h02, "w") as f:
    json.dump({"results": h02_results,
               "n_act_significant_concentrated": n_act_sig,
               "n_rep_significant_concentrated": n_rep_sig,
               "n_act_more_concentrated_than_rep": n_act_more,
               "n_shuffle": N_SHUFFLE}, f, indent=2)
print(f"  Saved: {out_h02}", flush=True)


# ==========================================
# H03: STRING PPI co-pole screen
# ==========================================
print("\n=== H03: STRING PPI co-pole enrichment ===", flush=True)

STRING_CACHE = ITER_DIR / "string_ppi_named_genes.json"

def download_string_interactions(gene_list, species=9606, score_threshold=0):
    """Download STRING interactions for a list of gene symbols."""
    identifiers = "%0d".join(gene_list)
    url = (
        "https://string-db.org/api/json/network"
        f"?identifiers={urllib.parse.quote(identifiers)}"
        f"&species={species}"
        f"&required_score={score_threshold}"
        f"&caller_identity=biomechinterp_research_iter0012"
    )
    print(f"  Querying STRING API... ({len(gene_list)} genes)", flush=True)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=60) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        print(f"  STRING returned {len(data)} interactions", flush=True)
        return data
    except Exception as e:
        print(f"  STRING API failed: {e}", flush=True)
        return None

# Try to load from cache first
if STRING_CACHE.exists():
    print(f"  Loading STRING data from cache: {STRING_CACHE}", flush=True)
    with open(STRING_CACHE) as f:
        string_data = json.load(f)
else:
    string_data = download_string_interactions(named_gene_list, score_threshold=0)
    if string_data is not None:
        with open(STRING_CACHE, "w") as f:
            json.dump(string_data, f)
        print(f"  Cached STRING data to: {STRING_CACHE}", flush=True)

h03_results = []
h03_meta = {}

if string_data is None or len(string_data) == 0:
    print("  STRING API unavailable. Running fallback: co-expression-proxy PPI from SV1/SV3 correlation.", flush=True)
    h03_meta["fallback"] = True
    h03_meta["fallback_reason"] = "STRING API unavailable or returned 0 edges"
    # Fallback: build proxy PPI from cosine similarity in embedding space at layer 0
    # (prior to transformer = pure positional/co-expression signal)
    mat_l0 = emb[0]
    idxs = [gene2emb_idx[g] for g in named_gene_list]
    X_l0 = mat_l0[idxs]  # [209, 512]
    # Normalize
    norms = np.linalg.norm(X_l0, axis=1, keepdims=True) + 1e-9
    X_norm = X_l0 / norms
    cosine_sim = X_norm @ X_norm.T  # [209, 209]
    # Build pairs with cosine_sim > 0.95 (high co-expression proxy)
    proxy_ppi = []
    for i in range(len(named_gene_list)):
        for j in range(i+1, len(named_gene_list)):
            if cosine_sim[i, j] > 0.90:  # threshold for proxy PPI
                proxy_ppi.append((named_gene_list[i], named_gene_list[j]))
    print(f"  Proxy PPI pairs (cosine > 0.90 at L0): {len(proxy_ppi)}", flush=True)
    string_ppi_pairs = proxy_ppi
    h03_meta["edge_type"] = "cosine_proxy_L0_gt0.90"
    h03_meta["n_pairs"] = len(proxy_ppi)
else:
    # Parse STRING edges and filter by score
    h03_meta["fallback"] = False
    string_ppi_pairs = []
    score_threshold_use = 700  # high confidence
    for edge in string_data:
        g1 = edge.get("preferredName_A", "").upper()
        g2 = edge.get("preferredName_B", "").upper()
        score = edge.get("score", 0)
        if g1 in named_set and g2 in named_set and score >= score_threshold_use:
            string_ppi_pairs.append((g1, g2))
    # Deduplicate
    seen = set()
    deduped = []
    for (a, b) in string_ppi_pairs:
        key = tuple(sorted([a, b]))
        if key not in seen:
            seen.add(key)
            deduped.append((a, b))
    string_ppi_pairs = deduped
    print(f"  STRING PPI pairs (score>={score_threshold_use}): {len(string_ppi_pairs)}", flush=True)
    h03_meta["edge_type"] = f"STRING_v12_combined_score_ge{score_threshold_use}"
    h03_meta["n_pairs"] = len(string_ppi_pairs)
    h03_meta["score_threshold"] = score_threshold_use

# Run co-pole enrichment on STRING/proxy PPI
N_SHUFFLE_STRING = 500

if len(string_ppi_pairs) < 5:
    print(f"  Too few PPI pairs ({len(string_ppi_pairs)}). Skipping enrichment.", flush=True)
    h03_meta["skipped"] = True
    h03_meta["skip_reason"] = "fewer than 5 pairs"
else:
    h03_meta["skipped"] = False
    for layer in range(12):
        proj = get_sv_projections(layer, sv_idx=1)
        top_pole, bot_pole = get_poles(proj)

        obs_ppi = copole_rate(string_ppi_pairs, top_pole, bot_pole)

        # Null: gene-label shuffle
        null_rates = []
        gene_arr = np.array(named_gene_list)
        for _ in range(N_SHUFFLE_STRING):
            perm = RNG.permutation(len(gene_arr))
            perm_map = {gene_arr[i]: gene_arr[perm[i]] for i in range(len(gene_arr))}
            shuffled_pairs = [(perm_map.get(a, a), perm_map.get(b, b)) for (a, b) in string_ppi_pairs]
            null_rates.append(copole_rate(shuffled_pairs, top_pole, bot_pole))

        null_rates = np.array(null_rates)
        emp_p = float(np.mean(null_rates >= obs_ppi))
        z = float((obs_ppi - null_rates.mean()) / (null_rates.std() + 1e-9))

        row = {
            "layer": layer,
            "obs_copole": round(float(obs_ppi), 4),
            "null_mean": round(float(null_rates.mean()), 4),
            "null_std": round(float(null_rates.std()), 4),
            "z": round(z, 3),
            "emp_p": round(emp_p, 4),
            "significant_p05": bool(emp_p < 0.05),
        }
        h03_results.append(row)
        print(
            f"  L{layer:02d}: obs={obs_ppi:.3f} null_mean={null_rates.mean():.3f} z={z:.2f} emp_p={emp_p:.3f}",
            flush=True,
        )

    n_sig = sum(1 for r in h03_results if r["significant_p05"])
    print(f"\n  H03 Summary: {n_sig}/12 layers significant (emp_p<0.05)", flush=True)
    h03_meta["n_layers_significant"] = n_sig

out_h03 = ITER_DIR / "h03_ppi_copole.json"
with open(out_h03, "w") as f:
    json.dump({"meta": h03_meta, "results": h03_results}, f, indent=2)
print(f"  Saved: {out_h03}", flush=True)


# ==========================================
# Master results summary
# ==========================================
master = {
    "iteration": "iter_0012",
    "hypotheses": {
        "H01": {
            "name": "Bootstrap CI on activation vs repression co-pole differential",
            "n_layers_ci_excludes_zero": n_ci_excl,
            "n_layers_p_lt_05": n_p05,
            "n_boot": N_BOOT,
            "artifact": "h01_bootstrap_act_rep_diff.json",
        },
        "H02": {
            "name": "Mean pairwise SV2 distance — spatial concentration test",
            "n_act_significant_concentrated": n_act_sig,
            "n_rep_significant_concentrated": n_rep_sig,
            "n_act_more_concentrated_than_rep": n_act_more,
            "n_shuffle": N_SHUFFLE,
            "artifact": "h02_sv2_pairwise_dist.json",
        },
        "H03": {
            "name": "STRING/proxy PPI co-pole enrichment",
            "meta": h03_meta,
            "n_layers_significant": h03_meta.get("n_layers_significant", "N/A"),
            "artifact": "h03_ppi_copole.json",
        },
    },
}

out_master = ITER_DIR / "iter0012_results.json"
with open(out_master, "w") as f:
    json.dump(master, f, indent=2)
print(f"\nMaster results saved: {out_master}", flush=True)
print("iter_0012 screen complete.", flush=True)
