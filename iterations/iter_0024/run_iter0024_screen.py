"""
iter_0024 Multi-Hypothesis Screen

H01 (manifold_distance, new_method: Dorothea confidence stratification):
    Use Dorothea regulon database (confidence A-E) to test if high-confidence
    TF-target pairs (A/B) are closer in scGPT embedding space than low-confidence
    (D/E) pairs. This replicates STRING-score stratification but in a regulatory
    network (independent database, different biological concept).
    Control: label-shuffled Dorothea pairs.

H02 (module_structure, new_method: TF activation hub centrality):
    From TRRUST, identify TFs ranked by number of activation targets vs repression
    targets. Test if high-activation-degree TFs are geometrically central
    (lower mean distance to all named genes) compared to high-repression TFs.
    Prediction: activation TFs should be geometric hubs in the manifold.

H03 (module_structure, new_family: GO ontology BP vs MF vs CC comparison):
    Fetch GO annotations for named genes via mygene (BP, MF, CC).
    Compute pairwise Jaccard similarity for each ontology.
    Test Spearman(Jaccard, -L2_distance) at layer 8 for all three.
    Identify which ontology has strongest geometric embedding signal.
"""

import numpy as np
import json
import csv
from pathlib import Path
from scipy.stats import mannwhitneyu, spearmanr, rankdata
import warnings
warnings.filterwarnings("ignore")

# ─── Paths ────────────────────────────────────────────────────────────────────
PROJECT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
               "/subproject_41_claude_topology_hypothesis_screening_autoloop")
ITER_DIR = PROJECT / "iterations" / "iter_0024"
ITER15_DIR = PROJECT / "iterations" / "iter_0015"
CYCLE1 = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
              "/subproject_38_geometric_residual_stream_interpretability"
              "/implementation/outputs/cycle1_main")
TRRUST_PATH = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
                   "/single_cell_mechinterp/external/networks/trrust_human.tsv")
DOROTHEA_PATH = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
                     "/single_cell_mechinterp/external/networks/dorothea_human.tsv")

ITER_DIR.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(42)

# ─── Load embeddings ──────────────────────────────────────────────────────────
print("Loading scGPT embeddings ...", flush=True)
EMB_PATH = CYCLE1 / "layer_gene_embeddings.npy"
emb = np.load(EMB_PATH)   # [12, 4803, 512]
N_LAYERS, N_GENES_TOTAL, N_DIM = emb.shape
print(f"  emb shape: {emb.shape}", flush=True)

# ─── Load gene list ───────────────────────────────────────────────────────────
with open(CYCLE1 / "gene_list.txt") as f:
    vocab_genes = [line.strip() for line in f]
gene_to_emb_idx = {g: i for i, g in enumerate(vocab_genes) if g}

# ─── Named genes from edge dataset ───────────────────────────────────────────
EDGE_PATH = CYCLE1 / "cycle1_edge_dataset.tsv"
named_gene_set = set()
with open(EDGE_PATH) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        named_gene_set.add(row['source'])
        named_gene_set.add(row['target'])

named_genes = sorted(g for g in named_gene_set if g in gene_to_emb_idx)
named_idx = np.array([gene_to_emb_idx[g] for g in named_genes])
N_NAMED = len(named_genes)
gene_to_named_idx = {g: i for i, g in enumerate(named_genes)}
print(f"  Named genes with embeddings: {N_NAMED}", flush=True)

# ─── Pre-compute L2 distances at each layer for named genes ──────────────────
print("Pre-computing pairwise L2 distances ...", flush=True)
named_emb = emb[:, named_idx, :]  # [12, N_NAMED, 512]

def pairwise_l2(X):
    """X: [N, D] → [N, N] L2 distances"""
    sq = np.sum(X**2, axis=1, keepdims=True)
    return np.sqrt(np.maximum(sq + sq.T - 2 * X @ X.T, 0))

dists = []
for l in range(N_LAYERS):
    dists.append(pairwise_l2(named_emb[l]))
dists = np.stack(dists)  # [12, N_NAMED, N_NAMED]
print(f"  dists shape: {dists.shape}", flush=True)

# ─── Load STRING pairs (for non-pair baseline) ───────────────────────────────
STRING_CACHE = ITER15_DIR / "string_ppi_score04_cache.json"
string_data = json.load(open(STRING_CACHE))
string_pairs_raw = string_data["pairs"]
string_pairs = set()
for p in string_pairs_raw:
    a, b = p["g1"], p["g2"]
    if a in gene_to_named_idx and b in gene_to_named_idx:
        i, j = gene_to_named_idx[a], gene_to_named_idx[b]
        string_pairs.add((min(i,j), max(i,j)))

print(f"  STRING pairs in named set: {len(string_pairs)}", flush=True)


# ═══════════════════════════════════════════════════════════════════════════════
# H01: Dorothea confidence stratification
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== H01: Dorothea confidence stratification ===", flush=True)

dorothea_pairs = {}  # (i,j) -> confidence_level
with open(DOROTHEA_PATH) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        src, tgt, conf = row['source'], row['target'], row['confidence']
        if src in gene_to_named_idx and tgt in gene_to_named_idx:
            i, j = gene_to_named_idx[src], gene_to_named_idx[tgt]
            key = (min(i,j), max(i,j))
            # Keep highest confidence (A > B > C > D > E)
            if key not in dorothea_pairs or conf < dorothea_pairs[key]:
                dorothea_pairs[key] = conf

print(f"  Dorothea pairs in named set: {len(dorothea_pairs)}", flush=True)

# Split by confidence tier
high_conf = [(i,j) for (i,j), c in dorothea_pairs.items() if c in ('A', 'B')]
low_conf  = [(i,j) for (i,j), c in dorothea_pairs.items() if c in ('D', 'E')]
print(f"  High-conf (A/B): {len(high_conf)}, Low-conf (D/E): {len(low_conf)}", flush=True)

# All non-Dorothea pairs (background)
non_dorothea_pairs = []
for ni in range(N_NAMED):
    for nj in range(ni+1, N_NAMED):
        key = (ni, nj)
        if key not in dorothea_pairs and key not in string_pairs:
            non_dorothea_pairs.append(key)
# Subsample for speed
rng.shuffle(non_dorothea_pairs_arr := np.array(non_dorothea_pairs))
n_bg = min(5000, len(non_dorothea_pairs_arr))
bg_pairs = [tuple(p) for p in non_dorothea_pairs_arr[:n_bg]]

h01_per_layer = []
for l in range(N_LAYERS):
    D = dists[l]
    high_d = np.array([D[i,j] for i,j in high_conf])
    low_d  = np.array([D[i,j] for i,j in low_conf])
    bg_d   = np.array([D[i,j] for i,j in bg_pairs])

    # AUROC: high_conf vs background (lower distance = positive)
    if len(high_conf) > 0 and len(bg_d) > 0:
        u_stat, p_high = mannwhitneyu(-high_d, -bg_d, alternative='greater')
        auroc_high = u_stat / (len(high_d) * len(bg_d))
    else:
        auroc_high, p_high = np.nan, np.nan

    # AUROC: low_conf vs background
    if len(low_conf) > 0 and len(bg_d) > 0:
        u_stat, p_low = mannwhitneyu(-low_d, -bg_d, alternative='greater')
        auroc_low = u_stat / (len(low_d) * len(bg_d))
    else:
        auroc_low, p_low = np.nan, np.nan

    # High vs low direct
    if len(high_conf) > 0 and len(low_conf) > 0:
        u_stat, p_hl = mannwhitneyu(-high_d, -low_d, alternative='greater')
        auroc_hl = u_stat / (len(high_d) * len(low_d))
        effect_hl = np.mean(high_d) - np.mean(low_d)
    else:
        auroc_hl, p_hl, effect_hl = np.nan, np.nan, np.nan

    h01_per_layer.append({
        "layer": l,
        "auroc_high_vs_bg": float(auroc_high),
        "p_high_vs_bg": float(p_high),
        "auroc_low_vs_bg": float(auroc_low),
        "p_low_vs_bg": float(p_low),
        "auroc_high_vs_low": float(auroc_hl),
        "p_high_vs_low": float(p_hl),
        "effect_high_vs_low_mean_dist": float(effect_hl),
        "mean_dist_high": float(np.mean(high_d)) if len(high_conf) > 0 else np.nan,
        "mean_dist_low": float(np.mean(low_d)) if len(low_conf) > 0 else np.nan,
    })

# Shuffle control at layer 8
l8 = 8
D8 = dists[l8]
# Shuffle gene labels for high_conf pairs
shuffle_aurocs = []
all_named_idx_arr = np.arange(N_NAMED)
for _ in range(500):
    perm = rng.permutation(N_NAMED)
    shuf_high_d = np.array([D8[perm[i], perm[j]] for i,j in high_conf])
    u, _ = mannwhitneyu(-shuf_high_d, -bg_d, alternative='greater')
    shuffle_aurocs.append(u / (len(high_conf) * len(bg_d)))

real_auroc_l8 = h01_per_layer[l8]['auroc_high_vs_bg']
perm_p = np.mean(np.array(shuffle_aurocs) >= real_auroc_l8)
shuffle_mean = float(np.mean(shuffle_aurocs))
shuffle_std  = float(np.std(shuffle_aurocs))
print(f"  Layer 8 AUROC high-conf vs bg: {real_auroc_l8:.4f}  (shuffle: {shuffle_mean:.4f}±{shuffle_std:.4f}, perm_p={perm_p:.4f})", flush=True)

n_sig_high = sum(1 for r in h01_per_layer if r['p_high_vs_bg'] < 0.05)
n_sig_hl   = sum(1 for r in h01_per_layer if r['p_high_vs_low'] < 0.05)
print(f"  Layers sig (high vs bg): {n_sig_high}/12  (high vs low): {n_sig_hl}/12", flush=True)

h01_result = {
    "hypothesis": "H01_dorothea_confidence_stratification",
    "n_high_conf_AB": len(high_conf),
    "n_low_conf_DE": len(low_conf),
    "n_bg_pairs": len(bg_pairs),
    "n_layers_high_sig_vs_bg": n_sig_high,
    "n_layers_high_vs_low_sig": n_sig_hl,
    "layer8_auroc_high_vs_bg": real_auroc_l8,
    "layer8_perm_p": float(perm_p),
    "shuffle_mean_auroc": shuffle_mean,
    "shuffle_std_auroc": shuffle_std,
    "per_layer": h01_per_layer,
    "summary_direction": "positive" if n_sig_high >= 8 else ("inconclusive" if n_sig_high >= 4 else "negative")
}
with open(ITER_DIR / "h01_dorothea_confidence.json", "w") as f:
    json.dump(h01_result, f, indent=2)
print(f"  Saved h01_dorothea_confidence.json", flush=True)


# ═══════════════════════════════════════════════════════════════════════════════
# H02: TF activation hub centrality
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== H02: TF activation hub centrality ===", flush=True)

# Load TRRUST
trrust_act = {}  # TF -> set of activation targets (among named genes)
trrust_rep = {}  # TF -> set of repression targets (among named genes)
with open(TRRUST_PATH) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 3:
            continue
        tf, target, mode = parts[0], parts[1], parts[2]
        if tf in gene_to_named_idx and target in gene_to_named_idx:
            if 'Activation' in mode:
                trrust_act.setdefault(tf, set()).add(target)
            if 'Repression' in mode:
                trrust_rep.setdefault(tf, set()).add(target)

# Only TFs present in named genes
act_tfs = {tf: len(tgts) for tf, tgts in trrust_act.items() if tf in gene_to_named_idx}
rep_tfs = {tf: len(tgts) for tf, tgts in trrust_rep.items() if tf in gene_to_named_idx}
print(f"  TFs with activation targets: {len(act_tfs)}, repression: {len(rep_tfs)}", flush=True)

# Sort by degree
act_sorted = sorted(act_tfs.items(), key=lambda x: -x[1])
rep_sorted = sorted(rep_tfs.items(), key=lambda x: -x[1])
print(f"  Top-5 activation TFs: {act_sorted[:5]}", flush=True)
print(f"  Top-5 repression TFs: {rep_sorted[:5]}", flush=True)

# For each TF, compute mean L2 to ALL other named genes (= geometric centrality proxy)
# centrality = mean pairwise distance (lower = more central)
def tf_mean_dist_to_all(tf_name, layer):
    if tf_name not in gene_to_named_idx:
        return np.nan
    tf_i = gene_to_named_idx[tf_name]
    D = dists[layer]
    drow = np.concatenate([D[tf_i, :tf_i], D[tf_i, tf_i+1:]])
    return float(np.mean(drow))

# Test at each layer: Spearman(activation_degree, -mean_dist_to_all)
# i.e., do high-activation TFs have lower mean distance (more central)?
h02_per_layer = []
for l in range(N_LAYERS):
    # All TFs present in named genes
    tfs_in_named = [tf for tf in (set(act_tfs) | set(rep_tfs))]

    act_degs = []
    rep_degs = []
    centralities = []
    for tf in tfs_in_named:
        ad = act_tfs.get(tf, 0)
        rd = rep_tfs.get(tf, 0)
        c  = tf_mean_dist_to_all(tf, l)
        act_degs.append(ad)
        rep_degs.append(rd)
        centralities.append(c)

    act_degs = np.array(act_degs)
    rep_degs = np.array(rep_degs)
    centralities = np.array(centralities)

    # Spearman(act_degree, proximity) where proximity = -distance
    sp_act, p_act = spearmanr(act_degs, -centralities)
    sp_rep, p_rep = spearmanr(rep_degs, -centralities)

    # Also: top-K act TFs vs top-K rep TFs — which are more central?
    K = min(15, len(act_sorted), len(rep_sorted))
    top_act_tfs = [tf for tf, _ in act_sorted[:K] if tf in gene_to_named_idx]
    top_rep_tfs = [tf for tf, _ in rep_sorted[:K] if tf in gene_to_named_idx]
    act_cents = np.array([tf_mean_dist_to_all(tf, l) for tf in top_act_tfs])
    rep_cents = np.array([tf_mean_dist_to_all(tf, l) for tf in top_rep_tfs])

    if len(act_cents) > 2 and len(rep_cents) > 2:
        u, p_topk = mannwhitneyu(-act_cents, -rep_cents, alternative='greater')
        auroc_topk = u / (len(act_cents) * len(rep_cents))
        effect_topk = float(np.mean(act_cents) - np.mean(rep_cents))
    else:
        auroc_topk, p_topk, effect_topk = np.nan, np.nan, np.nan

    h02_per_layer.append({
        "layer": l,
        "spearman_act_degree_vs_proximity": float(sp_act),
        "p_act": float(p_act),
        "spearman_rep_degree_vs_proximity": float(sp_rep),
        "p_rep": float(p_rep),
        "auroc_top_act_vs_top_rep": float(auroc_topk),
        "p_topk": float(p_topk),
        "effect_topk_mean_dist_diff": effect_topk,
        "n_tfs_total": len(tfs_in_named),
        "K_topk": K,
    })

n_sig_act = sum(1 for r in h02_per_layer if r['p_act'] < 0.05 and r['spearman_act_degree_vs_proximity'] > 0)
n_sig_topk = sum(1 for r in h02_per_layer if r['p_topk'] < 0.05)
print(f"  Layers where act_degree->proximity sig: {n_sig_act}/12", flush=True)
print(f"  Layers where top-act vs top-rep centrality sig: {n_sig_topk}/12", flush=True)
print(f"  Layer 8: sp_act={h02_per_layer[8]['spearman_act_degree_vs_proximity']:.3f} p={h02_per_layer[8]['p_act']:.4f}", flush=True)
print(f"  Layer 8: AUROC_topk={h02_per_layer[8]['auroc_top_act_vs_top_rep']:.3f} p={h02_per_layer[8]['p_topk']:.4f}", flush=True)

h02_result = {
    "hypothesis": "H02_TF_activation_hub_centrality",
    "n_act_tfs": len(act_tfs),
    "n_rep_tfs": len(rep_tfs),
    "K_topk": min(15, len(act_sorted), len(rep_sorted)),
    "n_layers_act_degree_sig": n_sig_act,
    "n_layers_topk_sig": n_sig_topk,
    "top_act_tfs": act_sorted[:10],
    "top_rep_tfs": rep_sorted[:10],
    "per_layer": h02_per_layer,
    "summary_direction": "positive" if n_sig_act >= 6 or n_sig_topk >= 6 else ("inconclusive" if n_sig_act >= 3 or n_sig_topk >= 3 else "negative")
}
with open(ITER_DIR / "h02_tf_hub_centrality.json", "w") as f:
    json.dump(h02_result, f, indent=2)
print(f"  Saved h02_tf_hub_centrality.json", flush=True)


# ═══════════════════════════════════════════════════════════════════════════════
# H03: GO BP vs MF vs CC comparison
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== H03: GO BP vs MF vs CC comparison ===", flush=True)

try:
    import mygene
    mg = mygene.MyGeneInfo()

    # Fetch GO annotations for all named genes
    print(f"  Querying mygene for {N_NAMED} named genes ...", flush=True)
    results = mg.querymany(named_genes, scopes='symbol', fields='go', species='human',
                           returnall=False, verbose=False)

    # Parse GO annotations per gene per ontology
    go_bp = {}  # gene -> set of GO IDs
    go_mf = {}
    go_cc = {}
    for hit in results:
        if 'notfound' in hit:
            continue
        q = hit.get('query', '')
        go_data = hit.get('go', {})
        # BP
        bp_terms = go_data.get('BP', [])
        if isinstance(bp_terms, dict):
            bp_terms = [bp_terms]
        go_bp[q] = set(t['id'] for t in bp_terms if 'id' in t)
        # MF
        mf_terms = go_data.get('MF', [])
        if isinstance(mf_terms, dict):
            mf_terms = [mf_terms]
        go_mf[q] = set(t['id'] for t in mf_terms if 'id' in t)
        # CC
        cc_terms = go_data.get('CC', [])
        if isinstance(cc_terms, dict):
            cc_terms = [cc_terms]
        go_cc[q] = set(t['id'] for t in cc_terms if 'id' in t)

    n_with_bp = sum(1 for g in named_genes if go_bp.get(g))
    n_with_mf = sum(1 for g in named_genes if go_mf.get(g))
    n_with_cc = sum(1 for g in named_genes if go_cc.get(g))
    print(f"  Genes with BP: {n_with_bp}, MF: {n_with_mf}, CC: {n_with_cc}", flush=True)

    def jaccard(s1, s2):
        if not s1 or not s2:
            return np.nan
        u = len(s1 | s2)
        return len(s1 & s2) / u if u > 0 else 0.0

    # Sample pairs with both annotations at layer 8
    TARGET_LAYER = 8
    D8 = dists[TARGET_LAYER]

    bp_jacc, mf_jacc, cc_jacc, pair_dists = [], [], [], []
    for ni in range(N_NAMED):
        for nj in range(ni+1, N_NAMED):
            g1, g2 = named_genes[ni], named_genes[nj]
            bj = jaccard(go_bp.get(g1, set()), go_bp.get(g2, set()))
            mj = jaccard(go_mf.get(g1, set()), go_mf.get(g2, set()))
            cj = jaccard(go_cc.get(g1, set()), go_cc.get(g2, set()))
            if not (np.isnan(bj) and np.isnan(mj) and np.isnan(cj)):
                bp_jacc.append(bj if not np.isnan(bj) else 0.0)
                mf_jacc.append(mj if not np.isnan(mj) else 0.0)
                cc_jacc.append(cj if not np.isnan(cj) else 0.0)
                pair_dists.append(D8[ni, nj])

    bp_jacc = np.array(bp_jacc)
    mf_jacc = np.array(mf_jacc)
    cc_jacc = np.array(cc_jacc)
    pair_dists = np.array(pair_dists)
    print(f"  Total pairs: {len(pair_dists)}", flush=True)

    # Spearman at layer 8
    sp_bp, p_bp = spearmanr(bp_jacc, -pair_dists)
    sp_mf, p_mf = spearmanr(mf_jacc, -pair_dists)
    sp_cc, p_cc = spearmanr(cc_jacc, -pair_dists)
    print(f"  BP Spearman: {sp_bp:.4f} p={p_bp:.2e}", flush=True)
    print(f"  MF Spearman: {sp_mf:.4f} p={p_mf:.2e}", flush=True)
    print(f"  CC Spearman: {sp_cc:.4f} p={p_cc:.2e}", flush=True)

    # Also at all layers
    h03_per_layer = []
    for l in range(N_LAYERS):
        Dl = dists[l]
        pd_l = np.array([Dl[ni, nj] for ni in range(N_NAMED) for nj in range(ni+1, N_NAMED)
                         if len(bp_jacc) > 0])
        # Recompute pair distances for each layer
        pd_l = []
        for ni in range(N_NAMED):
            for nj in range(ni+1, N_NAMED):
                pd_l.append(Dl[ni, nj])
        pd_l = np.array(pd_l)
        # Need matching pairs — use same order as above
        # Rebuild for layer l
        pd_this = []
        for ni in range(N_NAMED):
            for nj in range(ni+1, N_NAMED):
                pd_this.append(Dl[ni, nj])
        pd_this = np.array(pd_this)

        # Note: bp_jacc etc. were computed above in same order (all pairs in order ni<nj)
        if len(bp_jacc) == len(pd_this):
            sp_bp_l, p_bp_l = spearmanr(bp_jacc, -pd_this)
            sp_mf_l, p_mf_l = spearmanr(mf_jacc, -pd_this)
            sp_cc_l, p_cc_l = spearmanr(cc_jacc, -pd_this)
        else:
            sp_bp_l = sp_mf_l = sp_cc_l = np.nan
            p_bp_l = p_mf_l = p_cc_l = np.nan

        h03_per_layer.append({
            "layer": l,
            "spearman_bp": float(sp_bp_l), "p_bp": float(p_bp_l),
            "spearman_mf": float(sp_mf_l), "p_mf": float(p_mf_l),
            "spearman_cc": float(sp_cc_l), "p_cc": float(p_cc_l),
        })

    n_sig_bp = sum(1 for r in h03_per_layer if r['p_bp'] < 0.05 and r['spearman_bp'] > 0)
    n_sig_mf = sum(1 for r in h03_per_layer if r['p_mf'] < 0.05 and r['spearman_mf'] > 0)
    n_sig_cc = sum(1 for r in h03_per_layer if r['p_cc'] < 0.05 and r['spearman_cc'] > 0)
    print(f"  Layers sig: BP={n_sig_bp}/12, MF={n_sig_mf}/12, CC={n_sig_cc}/12", flush=True)

    best_ontology = max([('BP', n_sig_bp, sp_bp), ('MF', n_sig_mf, sp_mf), ('CC', n_sig_cc, sp_cc)],
                        key=lambda x: (x[1], abs(x[2])))[0]

    h03_result = {
        "hypothesis": "H03_GO_ontology_comparison_BP_MF_CC",
        "n_pairs_with_annotations": int(len(pair_dists)),
        "n_genes_with_bp": n_with_bp, "n_genes_with_mf": n_with_mf, "n_genes_with_cc": n_with_cc,
        "layer8_spearman_bp": float(sp_bp), "layer8_p_bp": float(p_bp),
        "layer8_spearman_mf": float(sp_mf), "layer8_p_mf": float(p_mf),
        "layer8_spearman_cc": float(sp_cc), "layer8_p_cc": float(p_cc),
        "n_layers_sig_bp": n_sig_bp, "n_layers_sig_mf": n_sig_mf, "n_layers_sig_cc": n_sig_cc,
        "best_ontology": best_ontology,
        "per_layer": h03_per_layer,
        "summary_direction": "positive" if max(n_sig_bp, n_sig_mf, n_sig_cc) >= 6 else "inconclusive",
        "go_query_success": True
    }
    go_query_success = True

except Exception as e:
    print(f"  ERROR in H03: {e}", flush=True)
    h03_result = {
        "hypothesis": "H03_GO_ontology_comparison_BP_MF_CC",
        "error": str(e),
        "go_query_success": False,
        "summary_direction": "inconclusive"
    }
    go_query_success = False

with open(ITER_DIR / "h03_go_ontology_comparison.json", "w") as f:
    json.dump(h03_result, f, indent=2)
print(f"  Saved h03_go_ontology_comparison.json", flush=True)


# ═══════════════════════════════════════════════════════════════════════════════
# Compile summary
# ═══════════════════════════════════════════════════════════════════════════════
summary = {
    "iteration": "iter_0024",
    "h01_dorothea": {
        "n_high_conf_AB": len(high_conf),
        "n_low_conf_DE": len(low_conf),
        "n_layers_sig_high_vs_bg": n_sig_high,
        "layer8_auroc_high_vs_bg": real_auroc_l8,
        "layer8_perm_p": float(perm_p),
        "direction": h01_result["summary_direction"]
    },
    "h02_hub": {
        "n_act_tfs": len(act_tfs),
        "n_rep_tfs": len(rep_tfs),
        "n_layers_sig_act_degree": n_sig_act,
        "n_layers_sig_topk": n_sig_topk,
        "layer8_spearman_act": h02_per_layer[8]['spearman_act_degree_vs_proximity'],
        "direction": h02_result["summary_direction"]
    },
    "h03_go": h03_result if not go_query_success else {
        "layer8_sp_bp": float(sp_bp),
        "layer8_sp_mf": float(sp_mf),
        "layer8_sp_cc": float(sp_cc),
        "n_layers_sig_bp": n_sig_bp,
        "n_layers_sig_mf": n_sig_mf,
        "n_layers_sig_cc": n_sig_cc,
        "best_ontology": best_ontology,
        "direction": h03_result["summary_direction"]
    }
}
with open(ITER_DIR / "iter0024_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print("\n=== All done. Summary saved. ===", flush=True)
print(json.dumps(summary, indent=2), flush=True)
