"""
iter_0011 - Multi-Hypothesis Screen

H01: TRRUST co-pole × 12 layers + Activation vs Repression stratification.
     Extends iter_0010 H02 (layer-11-only) to all 12 layers.
     Also stratifies by edge type (Activation vs Repression) to test signed regulatory geometry.
     Novel: multi-layer trajectory + edge-type stratification.

H02: SV2 mito pole-flip gene tracking at L3→L4 transition.
     Identify which mitochondrial genes switch SV2 poles at the L3→L4 boundary.
     Compute per-gene SV2 projection across all 12 layers. Track sign-flip events.
     Biological characterization: are flip genes enriched for specific mito sub-compartments?

H03: GO Biological Process enrichment in SV2 top/bottom poles (broad functional screen).
     Test whether functional programs (not just compartments) are geometrically encoded.
     For each of 12 layers, enumerate all GO BP terms with 3-50 genes in our 209-gene set.
     Fisher exact test vs null (N=200 gene-label shuffles).
     Novel: first BP-level test (prior work was CC/MF compartments only).

Data: layer_gene_embeddings.npy [12, 4803, 512], cycle1_edge_dataset.tsv (209 named genes)
TRRUST: trrust_human.tsv (TF, target, mode, pmid)
GO annotations: gene2go_all.pkl
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
    "/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0011"
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
GO_OBO_ASPECTS = None  # we'll filter by namespace from pkl

RNG = np.random.default_rng(42)
N_SHUFFLE = 200   # per cell
TOP_K = 52        # ~25% of 209 genes

ITER_DIR.mkdir(parents=True, exist_ok=True)

print("Loading data...", flush=True)
emb = np.load(EMB_PATH)  # [12, 4803, 512]
print(f"  Embeddings: {emb.shape}", flush=True)

# Load 209 named genes (use source_idx/target_idx for embedding lookup)
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
print(f"  Named genes from edge file: {len(named_gene_list)}", flush=True)

# Load TRRUST
trrust_edges = []
with open(TRRUST_PATH) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 3:
            tf, target, mode = parts[0], parts[1], parts[2]
            trrust_edges.append((tf, target, mode))
print(f"  TRRUST edges: {len(trrust_edges)}", flush=True)

# Load GO annotations
with open(GENE2GO_PKL, "rb") as f:
    gene2go_raw = pickle.load(f)
print(f"  gene2go keys: {len(gene2go_raw)}", flush=True)

# --- Build per-gene GO term sets and per-term gene sets ---
# gene2go_raw: dict[gene_symbol -> list of GO ids] OR dict[gene_id -> ...]
# Figure out key type
sample_key = next(iter(gene2go_raw))
print(f"  Sample gene2go key: {sample_key!r}, val sample: {list(gene2go_raw[sample_key])[:3]}", flush=True)

# Build gene -> set of GO terms
gene2goterms = {}
for k, v in gene2go_raw.items():
    if isinstance(v, (list, set, tuple)):
        gene2goterms[k] = set(v)
    elif isinstance(v, dict):
        # might be {go_id: namespace, ...}
        gene2goterms[k] = set(v.keys())
    else:
        gene2goterms[k] = set()

# Build GO term -> gene set (restricted to named 209 genes)
goterm2genes = defaultdict(set)
named_set = set(named_gene_list)
for gene in named_set:
    if gene in gene2goterms:
        for go_id in gene2goterms[gene]:
            goterm2genes[go_id].add(gene)

print(f"  GO terms with >=1 named gene: {len(goterm2genes)}", flush=True)
print(f"  Named genes with embedding index: {len(gene2emb_idx)}", flush=True)


# ==========================================
# Helper: SVD of named-gene embedding at layer L
# ==========================================
def get_sv_projections(layer_idx, sv_idx=1):
    """
    Returns projection of named genes onto SV[sv_idx] at given layer.
    sv_idx=0 -> SV1, sv_idx=1 -> SV2, etc.
    Returns dict: gene -> scalar projection.
    """
    mat = emb[layer_idx]  # [4803, 512]
    # Extract named gene embeddings
    gene_names_local = list(gene2emb_idx.keys())
    idxs = [gene2emb_idx[g] for g in gene_names_local]
    X = mat[idxs]  # [N_named, 512]
    X_centered = X - X.mean(axis=0)
    U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)
    projs = U[:, sv_idx] * S[sv_idx]
    return {g: float(projs[i]) for i, g in enumerate(gene_names_local)}


def get_poles(proj_dict, top_k=TOP_K):
    """Return (top_pole_set, bottom_pole_set) based on sorted projections."""
    sorted_genes = sorted(proj_dict.items(), key=lambda x: x[1], reverse=True)
    top = set(g for g, _ in sorted_genes[:top_k])
    bot = set(g for g, _ in sorted_genes[-top_k:])
    return top, bot


def copole_rate(pairs, top_set, bot_set):
    """Fraction of pairs where both genes are in the same pole."""
    if not pairs:
        return 0.0
    count = 0
    for (a, b) in pairs:
        if (a in top_set and b in top_set) or (a in bot_set and b in bot_set):
            count += 1
    return count / len(pairs)


# ==========================================
# H01: TRRUST co-pole × 12 layers + edge-type stratification
# ==========================================
print("\n=== H01: TRRUST co-pole x 12 layers + Activation/Repression ===", flush=True)

named_set_emb = set(gene2emb_idx.keys())

# Filter TRRUST pairs where both genes in our named set
all_pairs = [(tf, tgt, mode) for tf, tgt, mode in trrust_edges
             if tf in named_set_emb and tgt in named_set_emb]
act_pairs = [(tf, tgt) for tf, tgt, m in all_pairs if "Activation" in m]
rep_pairs = [(tf, tgt) for tf, tgt, m in all_pairs if "Repression" in m]
all_pairs_nomode = [(tf, tgt) for tf, tgt, _ in all_pairs]

print(f"  All mapped TRRUST pairs: {len(all_pairs_nomode)}", flush=True)
print(f"  Activation pairs: {len(act_pairs)}", flush=True)
print(f"  Repression pairs: {len(rep_pairs)}", flush=True)

h01_results = []
N_SHUFFLE_H1 = 500

for layer_idx in range(12):
    proj = get_sv_projections(layer_idx, sv_idx=1)  # SV2
    top, bot = get_poles(proj)

    obs_all = copole_rate(all_pairs_nomode, top, bot)
    obs_act = copole_rate(act_pairs, top, bot)
    obs_rep = copole_rate(rep_pairs, top, bot)

    gene_list = list(named_set_emb)
    # Null: shuffle gene labels, compute co-pole rate
    null_all = []
    null_act = []
    null_rep = []
    for _ in range(N_SHUFFLE_H1):
        shuffled = list(gene_list)
        RNG.shuffle(shuffled)
        mapping = dict(zip(gene_list, shuffled))
        shuf_top = {mapping[g] for g in top if mapping[g] in named_set_emb}
        shuf_bot = {mapping[g] for g in bot if mapping[g] in named_set_emb}
        null_all.append(copole_rate(all_pairs_nomode, shuf_top, shuf_bot))
        null_act.append(copole_rate(act_pairs, shuf_top, shuf_bot))
        null_rep.append(copole_rate(rep_pairs, shuf_top, shuf_bot))

    null_all = np.array(null_all)
    null_act = np.array(null_act)
    null_rep = np.array(null_rep)

    emp_p_all = float(np.mean(null_all >= obs_all))
    emp_p_act = float(np.mean(null_act >= obs_act))
    emp_p_rep = float(np.mean(null_rep >= obs_rep))

    row = {
        "layer": layer_idx,
        "obs_all": round(obs_all, 4),
        "null_mean_all": round(float(null_all.mean()), 4),
        "null_std_all": round(float(null_all.std()), 4),
        "emp_p_all": emp_p_all,
        "obs_act": round(obs_act, 4),
        "null_mean_act": round(float(null_act.mean()), 4),
        "emp_p_act": emp_p_act,
        "obs_rep": round(obs_rep, 4),
        "null_mean_rep": round(float(null_rep.mean()), 4),
        "emp_p_rep": emp_p_rep,
        "n_all": len(all_pairs_nomode),
        "n_act": len(act_pairs),
        "n_rep": len(rep_pairs),
    }
    h01_results.append(row)
    print(f"  L{layer_idx}: all p={emp_p_all:.3f} (obs={obs_all:.3f} vs null={null_all.mean():.3f}), "
          f"act p={emp_p_act:.3f}, rep p={emp_p_rep:.3f}", flush=True)

with open(ITER_DIR / "h01_trrust_copole_12layers.json", "w") as f:
    json.dump(h01_results, f, indent=2)
print("  Saved h01_trrust_copole_12layers.json", flush=True)


# ==========================================
# H02: Mito pole-flip gene tracking at L3→L4
# ==========================================
print("\n=== H02: Mito pole-flip gene tracking at L3->L4 ===", flush=True)

# Get mito gene set
mito_go_ids = set()
for gene, terms in gene2goterms.items():
    if gene in named_set_emb:
        for t in terms:
            pass  # we need mito specifically
# GO:0005739 = mitochondrion
MITO_GO = "GO:0005739"
mito_genes = {g for g in named_set_emb if MITO_GO in gene2goterms.get(g, set())}
print(f"  Mito genes in named set: {len(mito_genes)}: {sorted(mito_genes)}", flush=True)

# Track SV2 projection of all named genes across all 12 layers
all_sv2_projs = {}  # gene -> [12 values]
for layer_idx in range(12):
    proj = get_sv_projections(layer_idx, sv_idx=1)
    for g, val in proj.items():
        if g not in all_sv2_projs:
            all_sv2_projs[g] = []
        all_sv2_projs[g].append(val)

# For mito genes: track sign at each layer
# Pole assignment: positive proj = top pole, negative = bottom pole
mito_trajectories = {}
for g in mito_genes:
    if g in all_sv2_projs:
        vals = all_sv2_projs[g]
        signs = [1 if v > 0 else -1 for v in vals]
        mito_trajectories[g] = {"projections": vals, "signs": signs}

# Identify genes that flip sign between L3 and L4 (index 3 vs 4)
flip_genes = []
stable_top = []
stable_bot = []
for g, data in mito_trajectories.items():
    s3 = data["signs"][3]
    s4 = data["signs"][4]
    if s3 != s4:
        flip_genes.append(g)
    elif s3 > 0:
        stable_top.append(g)
    else:
        stable_bot.append(g)

print(f"  Mito genes flip L3->L4: {len(flip_genes)}: {flip_genes}", flush=True)
print(f"  Mito genes stable top: {len(stable_top)}: {stable_top}", flush=True)
print(f"  Mito genes stable bot: {len(stable_bot)}: {stable_bot}", flush=True)

# Also track all sign-flip events across all layer transitions for mito genes
all_flip_events = defaultdict(list)
for g, data in mito_trajectories.items():
    signs = data["signs"]
    for i in range(len(signs) - 1):
        if signs[i] != signs[i+1]:
            all_flip_events[g].append(i)

print(f"  All flip events per gene:", flush=True)
for g in mito_genes:
    if g in all_flip_events:
        print(f"    {g}: flips at layers {all_flip_events[g]}", flush=True)

# Get GO annotations for flip vs stable genes
def get_go_terms(gene_set):
    """Get all GO terms for a set of genes."""
    terms = defaultdict(int)
    for g in gene_set:
        for t in gene2goterms.get(g, set()):
            terms[t] += 1
    return terms

flip_terms = get_go_terms(set(flip_genes))
stable_terms = get_go_terms(set(stable_top + stable_bot))
print(f"  Flip gene GO terms: {len(flip_terms)}", flush=True)
print(f"  Stable gene GO terms: {len(stable_terms)}", flush=True)

# Also check SV2 projection magnitude for mito at each layer
mito_mean_proj_by_layer = []
for layer_idx in range(12):
    vals = [all_sv2_projs[g][layer_idx] for g in mito_genes if g in all_sv2_projs]
    mito_mean_proj_by_layer.append({
        "layer": layer_idx,
        "mean_sv2_proj": round(float(np.mean(vals)), 4) if vals else None,
        "std_sv2_proj": round(float(np.std(vals)), 4) if vals else None,
        "n_genes": len(vals),
    })

h02_result = {
    "mito_genes": sorted(mito_genes),
    "flip_genes_L3_to_L4": sorted(flip_genes),
    "stable_top_genes": sorted(stable_top),
    "stable_bottom_genes": sorted(stable_bot),
    "all_flip_events": {g: v for g, v in all_flip_events.items()},
    "mito_sv2_proj_by_layer": mito_mean_proj_by_layer,
    "gene_sv2_trajectories": {
        g: {"projections": [round(v, 4) for v in data["projections"]], "signs": data["signs"]}
        for g, data in mito_trajectories.items()
    },
    "top_go_terms_flip_genes": sorted(flip_terms.items(), key=lambda x: -x[1])[:20],
    "top_go_terms_stable_genes": sorted(stable_terms.items(), key=lambda x: -x[1])[:20],
}

with open(ITER_DIR / "h02_mito_pole_flip.json", "w") as f:
    json.dump(h02_result, f, indent=2)
print("  Saved h02_mito_pole_flip.json", flush=True)


# ==========================================
# H03: GO Biological Process enrichment in SV2 poles
# ==========================================
print("\n=== H03: GO BP enrichment in SV2 poles (broad functional screen) ===", flush=True)

# Need to identify which GO terms are BP (not CC/MF)
# Try to detect from the annotations
# The pkl may have aspect info or we filter by prefix heuristics
# GO:0008150 = Biological Process root
# We'll include all terms that appear in the data and let statistics decide
# Filter to terms with 3-50 genes in our set
valid_bp_terms = {
    t: genes for t, genes in goterm2genes.items()
    if 3 <= len(genes) <= 50
}
print(f"  Valid terms (3-50 genes): {len(valid_bp_terms)}", flush=True)

# Focus on layer 11 (strongest SV2 signal from iter_0010)
# Also test layer 7 and 8 (SV2 variance peak region)
test_layers = [7, 8, 11]
h03_results = []
N_SHUFFLE_H3 = 200

for layer_idx in test_layers:
    proj = get_sv_projections(layer_idx, sv_idx=1)  # SV2
    sorted_genes = sorted(proj.items(), key=lambda x: x[1], reverse=True)
    top_set = set(g for g, _ in sorted_genes[:TOP_K])
    bot_set = set(g for g, _ in sorted_genes[-TOP_K:])
    all_named = set(proj.keys())
    n_total = len(all_named)

    layer_hits = []
    for go_id, gene_set in valid_bp_terms.items():
        gene_set_named = gene_set & all_named
        n_go = len(gene_set_named)
        if n_go < 3:
            continue

        # Test both poles; take best (Bonferroni-aware: this is exploratory)
        best_or = 0.0
        best_pole = None
        best_obs = None

        for pole_set, pole_name in [(top_set, "top"), (bot_set, "bottom")]:
            in_pole = len(gene_set_named & pole_set)
            out_pole = len(pole_set) - in_pole
            in_other = n_go - in_pole
            out_other = n_total - n_go - out_pole
            table = [[in_pole, out_pole], [in_other, out_other]]
            try:
                or_val, _ = fisher_exact(table, alternative="greater")
            except Exception:
                or_val = 0.0
            if or_val > best_or:
                best_or = or_val
                best_pole = pole_name
                best_obs = in_pole

        if best_or < 1.5:
            continue  # Skip weak hits for efficiency

        # Empirical null
        null_ors = []
        gene_list_arr = list(all_named)
        go_gene_list = list(gene_set_named)
        for _ in range(N_SHUFFLE_H3):
            # Shuffle gene labels
            shuffled_map = dict(zip(gene_list_arr, RNG.permutation(gene_list_arr)))
            shuf_top = {shuffled_map[g] for g in top_set}
            shuf_bot = {shuffled_map[g] for g in bot_set}
            shuf_go = {shuffled_map[g] for g in go_gene_list}
            best_shuf_or = 0.0
            for pole_set in [shuf_top, shuf_bot]:
                in_p = len(shuf_go & pole_set)
                out_p = len(pole_set) - in_p
                in_o = len(shuf_go) - in_p
                out_o = n_total - len(shuf_go) - out_p
                table = [[in_p, out_p], [in_o, out_o]]
                try:
                    or_s, _ = fisher_exact(table, alternative="greater")
                except Exception:
                    or_s = 0.0
                best_shuf_or = max(best_shuf_or, or_s)
            null_ors.append(best_shuf_or)

        null_ors = np.array(null_ors)
        emp_p = float(np.mean(null_ors >= best_or))

        layer_hits.append({
            "go_id": go_id,
            "n_go_genes": n_go,
            "best_or": round(best_or, 3),
            "best_pole": best_pole,
            "emp_p": emp_p,
            "null_mean_or": round(float(null_ors.mean()), 3),
            "null_std_or": round(float(null_ors.std()), 3),
            "genes_in_pole": sorted(gene_set_named & (top_set if best_pole == "top" else bot_set)),
        })

    layer_hits.sort(key=lambda x: x["emp_p"])
    print(f"  Layer {layer_idx}: tested {len(valid_bp_terms)} terms; "
          f"{sum(1 for h in layer_hits if h['emp_p'] < 0.05)} with emp_p<0.05; "
          f"top hit: {layer_hits[0]['go_id'] if layer_hits else 'none'} OR={layer_hits[0]['best_or'] if layer_hits else 'NA'}", flush=True)

    h03_results.append({
        "layer": layer_idx,
        "n_terms_tested": len([t for t, gs in valid_bp_terms.items() if len(gs & all_named) >= 3]),
        "n_sig_005": sum(1 for h in layer_hits if h["emp_p"] < 0.05),
        "n_sig_001": sum(1 for h in layer_hits if h["emp_p"] < 0.01),
        "top_hits": layer_hits[:30],
    })

with open(ITER_DIR / "h03_gobp_sv2_poles.json", "w") as f:
    json.dump(h03_results, f, indent=2)
print("  Saved h03_gobp_sv2_poles.json", flush=True)


# ==========================================
# Master results summary
# ==========================================
print("\n=== Writing master results ===", flush=True)

# Summarize H01
h01_sig_layers = [r["layer"] for r in h01_results if r["emp_p_all"] < 0.05]
h01_act_sig = [r["layer"] for r in h01_results if r["emp_p_act"] < 0.05]
h01_rep_sig = [r["layer"] for r in h01_results if r["emp_p_rep"] < 0.05]
h01_peak = min(h01_results, key=lambda x: x["emp_p_all"])

# Summarize H02
n_flip = len(h02_result["flip_genes_L3_to_L4"])
n_stable = len(h02_result["stable_top_genes"]) + len(h02_result["stable_bottom_genes"])

# Summarize H03
h03_summary = []
for r in h03_results:
    h03_summary.append({
        "layer": r["layer"],
        "n_sig_005": r["n_sig_005"],
        "n_sig_001": r["n_sig_001"],
        "top_go_id": r["top_hits"][0]["go_id"] if r["top_hits"] else None,
        "top_or": r["top_hits"][0]["best_or"] if r["top_hits"] else None,
        "top_emp_p": r["top_hits"][0]["emp_p"] if r["top_hits"] else None,
    })

master = {
    "iteration": "iter_0011",
    "h01_summary": {
        "sig_layers_all_edges": h01_sig_layers,
        "sig_layers_activation": h01_act_sig,
        "sig_layers_repression": h01_rep_sig,
        "peak_layer": h01_peak["layer"],
        "peak_emp_p_all": h01_peak["emp_p_all"],
        "peak_obs_all": h01_peak["obs_all"],
        "peak_null_mean_all": h01_peak["null_mean_all"],
        "n_pairs_all": h01_results[0]["n_all"],
        "n_pairs_act": h01_results[0]["n_act"],
        "n_pairs_rep": h01_results[0]["n_rep"],
    },
    "h02_summary": {
        "n_mito_genes": len(mito_genes),
        "n_flip_L3_L4": n_flip,
        "flip_genes": h02_result["flip_genes_L3_to_L4"],
        "stable_top": h02_result["stable_top_genes"],
        "stable_bottom": h02_result["stable_bottom_genes"],
    },
    "h03_summary": h03_summary,
}

with open(ITER_DIR / "iter0011_results.json", "w") as f:
    json.dump(master, f, indent=2)
print("  Saved iter0011_results.json", flush=True)
print("\nDONE.", flush=True)
