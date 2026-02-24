"""
iter_0023 Multi-Hypothesis Screen

H01 (module_structure, refinement from iter_0022 H02):
    Cell-type marker separation — expanded to 5+ types + contamination control.
    (a) Expand cell-type marker panel to include macrophage, NK, and myeloid if present.
    (b) Contamination control: test if non-marker named genes (random same-sized subsets)
        show comparable separation — if NOT, confirms the signal is marker-specific.
    (c) Cross-type AUROC per cell-type pair to identify which pairs are best separated.

H02 (module_structure, new_method: GO biological process proximity):
    Use mygene to fetch GO Biological Process terms for named genes.
    Compute GO BP Jaccard similarity for all named-gene pairs with shared annotations.
    Test: higher GO BP Jaccard → lower embedding L2 distance (Spearman).
    Compare AUROC for high-GO-Jaccard pairs vs low-GO-Jaccard pairs at each of 12 layers.
    This provides a 3rd independent biological anchor beyond STRING+TRRUST.

H03 (manifold_distance, new_method: TRRUST activation vs repression proximity split):
    Split TRRUST-exclusive pairs (NOT in STRING, N~141) by direction:
        Activation, Repression, Unknown.
    Test if activation vs repression pairs show different embedding proximity.
    If directional semantics are encoded in geometry, we expect a difference.
    Also: fit sigmoid / linear model to the per-layer AUROC curve from iter_0022 H02
    to characterize the functional form of geometric deepening.
"""

import numpy as np
import json
import csv
from pathlib import Path
from scipy.stats import mannwhitneyu, spearmanr, linregress
import warnings
warnings.filterwarnings("ignore")

# ─── Paths ────────────────────────────────────────────────────────────────────
PROJECT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
               "/subproject_41_claude_topology_hypothesis_screening_autoloop")
ITER_DIR = PROJECT / "iterations" / "iter_0023"
ITER15_DIR = PROJECT / "iterations" / "iter_0015"
ITER22_DIR = PROJECT / "iterations" / "iter_0022"
CYCLE1 = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
              "/subproject_38_geometric_residual_stream_interpretability"
              "/implementation/outputs/cycle1_main")
TRRUST_PATH = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
                   "/single_cell_mechinterp/external/networks/trrust_human.tsv")

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

# ─── Load STRING pairs (for reference) ───────────────────────────────────────
STRING_CACHE = ITER15_DIR / "string_ppi_score04_cache.json"
string_data = json.load(open(STRING_CACHE))
string_pairs_raw = string_data["pairs"]
string_pairs = []
for p in string_pairs_raw:
    g1, g2, sc = p["g1"], p["g2"], p["score"]
    if g1 in gene_to_named_idx and g2 in gene_to_named_idx:
        i, j = gene_to_named_idx[g1], gene_to_named_idx[g2]
        if i != j:
            string_pairs.append((i, j, sc))
string_set = {(min(i,j), max(i,j)) for i,j,_ in string_pairs}
N_STRING = len(string_pairs)
print(f"  STRING pairs: {N_STRING}", flush=True)

# Non-STRING sample for baseline
all_named_pairs = [(i,j) for i in range(N_NAMED) for j in range(i+1,N_NAMED)]
non_string_pairs_all = [(i,j) for i,j in all_named_pairs if (i,j) not in string_set]
ns_rng = np.random.default_rng(42)
ns_sel = ns_rng.choice(len(non_string_pairs_all), size=min(N_STRING*3, len(non_string_pairs_all)), replace=False)
non_string_pairs = [non_string_pairs_all[k] for k in ns_sel]
ns_idx1 = np.array([p[0] for p in non_string_pairs])
ns_idx2 = np.array([p[1] for p in non_string_pairs])

# ─── Helper: L2-normalize embeddings at a layer ───────────────────────────────
def get_unit_emb(layer):
    e = emb[layer][named_idx]
    norms = np.linalg.norm(e, axis=1, keepdims=True)
    norms[norms == 0] = 1
    return e / norms


# ═══════════════════════════════════════════════════════════════════════════════
# H01: Cell-type marker separation — expanded + contamination control
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== H01: Cell-type marker separation (expanded + contamination) ===", flush=True)

# Expanded canonical lung cell-type markers (broader panel)
CELL_TYPE_MARKERS_EXPANDED = {
    "T_cell":       ["CD3G", "CD8A", "CD8B", "CD6", "RUNX3", "FOXP3", "PRF1", "GZMB", "CD247"],
    "B_cell":       ["CD19", "MS4A1", "CD79A", "SPIB", "PAX5", "BLNK"],
    "macrophage":   ["ALOX5", "PTPRC", "CD68", "CSF1R", "MRC1", "ITGAM"],
    "NK_cell":      ["NCAM1", "KLRB1", "KLRD1", "NKG7", "GNLY"],
    "fibroblast":   ["DCN", "COL1A1", "VIM", "FAP", "PDGFRA"],
    "epithelial":   ["EPCAM", "KRT18", "KRT19", "CLDN4", "MUC1"],
    "endothelial":  ["KDR", "ICAM1", "PECAM1", "CDH5", "ENG"],
    "myeloid":      ["CD14", "LYZ", "FCGR3A", "S100A8", "S100A9"],
}

# Filter to genes present in named set (need >= 2)
ct_markers_filtered = {}
for ct, markers in CELL_TYPE_MARKERS_EXPANDED.items():
    present = [g for g in markers if g in gene_to_named_idx]
    if len(present) >= 2:
        ct_markers_filtered[ct] = present

print("Cell-type marker genes (filtered):", flush=True)
all_ct_gene_idxs = set()
for ct, genes in ct_markers_filtered.items():
    idxs = [gene_to_named_idx[g] for g in genes]
    all_ct_gene_idxs.update(idxs)
    print(f"  {ct} (n={len(genes)}): {genes}", flush=True)

ct_gene_sets = {ct: [gene_to_named_idx[g] for g in genes]
                for ct, genes in ct_markers_filtered.items()}

# Build within-cell-type and cross-cell-type pairs
within_ct_pairs = []
cross_ct_pairs = []
for ct, gene_idxs in ct_gene_sets.items():
    for a in range(len(gene_idxs)):
        for b in range(a+1, len(gene_idxs)):
            within_ct_pairs.append((gene_idxs[a], gene_idxs[b], ct))
ct_list = list(ct_gene_sets.keys())
for ci in range(len(ct_list)):
    for cj in range(ci+1, len(ct_list)):
        for gi in ct_gene_sets[ct_list[ci]]:
            for gj in ct_gene_sets[ct_list[cj]]:
                cross_ct_pairs.append((gi, gj, ct_list[ci], ct_list[cj]))

print(f"  Within-cell-type pairs: {len(within_ct_pairs)}", flush=True)
print(f"  Cross-cell-type pairs: {len(cross_ct_pairs)}", flush=True)

within_idx1 = np.array([p[0] for p in within_ct_pairs])
within_idx2 = np.array([p[1] for p in within_ct_pairs])
cross_idx1 = np.array([p[0] for p in cross_ct_pairs])
cross_idx2 = np.array([p[1] for p in cross_ct_pairs])

# Contamination control: random same-sized non-marker named genes
non_marker_genes = [i for i in range(N_NAMED) if i not in all_ct_gene_idxs]
n_marker_genes = len(all_ct_gene_idxs)
print(f"  Marker genes: {n_marker_genes}, Non-marker genes: {len(non_marker_genes)}", flush=True)

h01_results = []
h01_pairwise = {}  # per cell-type-pair AUROC at layer 8

for layer in range(N_LAYERS):
    ue = get_unit_emb(layer)

    d_within = np.sqrt(np.sum((ue[within_idx1] - ue[within_idx2])**2, axis=1))
    d_cross = np.sqrt(np.sum((ue[cross_idx1] - ue[cross_idx2])**2, axis=1))

    mw_u, mw_p = mannwhitneyu(d_within, d_cross, alternative='less')
    effect = np.mean(d_within) - np.mean(d_cross)
    auroc = 1.0 - mw_u / (len(d_within) * len(d_cross))

    # Contamination control: random non-marker same-sized partition
    n_types = len(ct_gene_sets)
    ct_sizes = [len(v) for v in ct_gene_sets.values()]
    total_ct_size = sum(ct_sizes)
    ctrl_rng = np.random.default_rng(layer + 1000)

    # Sample same number of non-marker genes
    if len(non_marker_genes) >= total_ct_size:
        ctrl_sel = ctrl_rng.choice(non_marker_genes, size=total_ct_size, replace=False)
        ctrl_within = []
        ctrl_cross = []
        pos = 0
        ctrl_sets = []
        for sz in ct_sizes:
            ctrl_sets.append(list(ctrl_sel[pos:pos+sz]))
            pos += sz
        for ci2 in range(n_types):
            for a in range(len(ctrl_sets[ci2])):
                for b in range(a+1, len(ctrl_sets[ci2])):
                    ctrl_within.append((ctrl_sets[ci2][a], ctrl_sets[ci2][b]))
            for cj2 in range(ci2+1, n_types):
                for gi in ctrl_sets[ci2]:
                    for gj in ctrl_sets[cj2]:
                        ctrl_cross.append((gi, gj))

        if len(ctrl_within) > 0 and len(ctrl_cross) > 0:
            cw_i1 = np.array([p[0] for p in ctrl_within])
            cw_i2 = np.array([p[1] for p in ctrl_within])
            cc_i1 = np.array([p[0] for p in ctrl_cross])
            cc_i2 = np.array([p[1] for p in ctrl_cross])
            d_ctrl_within = np.sqrt(np.sum((ue[cw_i1] - ue[cw_i2])**2, axis=1))
            d_ctrl_cross = np.sqrt(np.sum((ue[cc_i1] - ue[cc_i2])**2, axis=1))
            ctrl_mw_u, ctrl_mw_p = mannwhitneyu(d_ctrl_within, d_ctrl_cross, alternative='less')
            ctrl_auroc = 1.0 - ctrl_mw_u / (len(d_ctrl_within) * len(d_ctrl_cross))
            ctrl_effect = np.mean(d_ctrl_within) - np.mean(d_ctrl_cross)
        else:
            ctrl_auroc = 0.5
            ctrl_effect = 0.0
            ctrl_mw_p = 1.0
    else:
        ctrl_auroc = None
        ctrl_effect = None
        ctrl_mw_p = None

    h01_results.append({
        "layer": layer,
        "within_mean_dist": float(np.mean(d_within)),
        "cross_mean_dist": float(np.mean(d_cross)),
        "effect": float(effect),
        "auroc": float(auroc),
        "mw_p": float(mw_p),
        "control_auroc": float(ctrl_auroc) if ctrl_auroc is not None else None,
        "control_effect": float(ctrl_effect) if ctrl_effect is not None else None,
        "control_mw_p": float(ctrl_mw_p) if ctrl_mw_p is not None else None,
        "specificity_auroc_delta": float(auroc - ctrl_auroc) if ctrl_auroc is not None else None,
    })
    ctrl_str = f"{ctrl_auroc:.3f}" if ctrl_auroc is not None else "N/A"
    print(f"  L{layer:2d}: AUROC={auroc:.3f} effect={effect:.4f} mw_p={mw_p:.3e}"
          f" | ctrl_AUROC={ctrl_str}",
          flush=True)

    # Per-cell-type pair AUROC at layer 8
    if layer == 8:
        for ci2 in range(len(ct_list)):
            for cj2 in range(ci2+1, len(ct_list)):
                ct_a, ct_b = ct_list[ci2], ct_list[cj2]
                pair_cross = [(gi, gj) for gi, gj, ca, cb in cross_ct_pairs
                              if (ca == ct_a and cb == ct_b) or (ca == ct_b and cb == ct_a)]
                if len(pair_cross) == 0:
                    continue
                pc_i1 = np.array([p[0] for p in pair_cross])
                pc_i2 = np.array([p[1] for p in pair_cross])
                d_pair = np.sqrt(np.sum((ue[pc_i1] - ue[pc_i2])**2, axis=1))

                # vs random same-sized pair
                rand_sel = rng.choice(len(all_named_pairs), size=len(pair_cross), replace=False)
                rp_i1 = np.array([all_named_pairs[k][0] for k in rand_sel])
                rp_i2 = np.array([all_named_pairs[k][1] for k in rand_sel])
                d_rand = np.sqrt(np.sum((ue[rp_i1] - ue[rp_i2])**2, axis=1))
                pw_u, pw_p = mannwhitneyu(d_pair, d_rand, alternative='less')
                pw_auroc = 1.0 - pw_u / (len(d_pair) * len(d_rand))
                h01_pairwise[f"{ct_a}_vs_{ct_b}"] = {
                    "n_pairs": len(pair_cross),
                    "mean_cross_dist": float(np.mean(d_pair)),
                    "mean_rand_dist": float(np.mean(d_rand)),
                    "auroc_vs_random": float(pw_auroc),
                    "mw_p": float(pw_p),
                }

aurocs1 = [r["auroc"] for r in h01_results]
ctrl_aurocs = [r["control_auroc"] for r in h01_results if r["control_auroc"] is not None]
delta_aurocs = [r["specificity_auroc_delta"] for r in h01_results if r["specificity_auroc_delta"] is not None]

# Layer curve functional form: linear regression on logit(AUROC) vs layer
logit_aurocs = [np.log(a / (1.0 - a + 1e-6) + 1e-6) for a in aurocs1]
slope, intercept, rval, pval, stderr = linregress(list(range(N_LAYERS)), logit_aurocs)

print(f"\n  Cell-type AUROC: mean={np.mean(aurocs1):.3f} range=[{np.min(aurocs1):.3f},{np.max(aurocs1):.3f}]")
print(f"  Control AUROC:   mean={np.mean(ctrl_aurocs):.3f} range=[{np.min(ctrl_aurocs):.3f},{np.max(ctrl_aurocs):.3f}]")
print(f"  Specificity delta AUROC: mean={np.mean(delta_aurocs):.3f}")
print(f"  Layer curve linear fit (logit AUROC vs layer): slope={slope:.4f} r={rval:.3f} p={pval:.4f}")
print(f"  Cell-type pair AUROCs at L8: {h01_pairwise}", flush=True)

n_ct_types = len(ct_markers_filtered)
h01_out = {
    "hypothesis": "H01_cell_type_marker_expansion_contamination_control",
    "n_cell_types": n_ct_types,
    "cell_types": {ct: genes for ct, genes in ct_markers_filtered.items()},
    "n_within_pairs": len(within_ct_pairs),
    "n_cross_pairs": len(cross_ct_pairs),
    "n_marker_genes": n_marker_genes,
    "n_non_marker_genes": len(non_marker_genes),
    "per_layer": h01_results,
    "pairwise_auroc_layer8": h01_pairwise,
    "layer_curve_fit": {
        "slope_logit_auroc_vs_layer": float(slope),
        "r_value": float(rval),
        "p_value": float(pval),
    },
    "summary": {
        "auroc_mean": float(np.mean(aurocs1)),
        "auroc_range": [float(np.min(aurocs1)), float(np.max(aurocs1))],
        "control_auroc_mean": float(np.mean(ctrl_aurocs)) if ctrl_aurocs else None,
        "specificity_delta_mean": float(np.mean(delta_aurocs)) if delta_aurocs else None,
        "n_layers_sig_mw": sum(1 for r in h01_results if r["mw_p"] < 0.05),
        "n_layers_control_not_sig": sum(1 for r in h01_results
                                        if r["control_mw_p"] is not None and r["control_mw_p"] >= 0.05),
    }
}
json.dump(h01_out, open(ITER_DIR / "h01_cell_type_expansion_contamination.json", "w"), indent=2)
print("  -> Saved h01_cell_type_expansion_contamination.json", flush=True)


# ═══════════════════════════════════════════════════════════════════════════════
# H02: GO Biological Process proximity (new method using mygene)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== H02: GO Biological Process proximity ===", flush=True)

try:
    import mygene
    mg = mygene.MyGeneInfo()

    # Fetch GO BP annotations for named genes (batch query)
    print(f"  Fetching GO BP annotations for {N_NAMED} genes ...", flush=True)

    # Query in batches of 200
    gene_go_bp = {}  # gene -> set of GO BP IDs
    batch_size = 200
    n_successful = 0

    for batch_start in range(0, N_NAMED, batch_size):
        batch_genes = named_genes[batch_start:batch_start + batch_size]
        try:
            results = mg.querymany(batch_genes, scopes='symbol', fields='go.BP', species='human',
                                   returnall=False, verbose=False)
            for r in results:
                gene_sym = r.get('query', '')
                if 'go' in r and 'BP' in r['go']:
                    bp_terms = r['go']['BP']
                    if isinstance(bp_terms, dict):
                        bp_terms = [bp_terms]
                    go_ids = {t['id'] for t in bp_terms if 'id' in t}
                    gene_go_bp[gene_sym] = go_ids
                    n_successful += 1
        except Exception as e:
            print(f"  Warning: batch {batch_start} failed: {e}", flush=True)

    print(f"  Genes with GO BP annotations: {len(gene_go_bp)}/{N_NAMED}", flush=True)

    # Compute GO BP Jaccard for named-gene pairs with annotations
    go_pairs = []  # (i, j, jaccard)
    genes_with_go = [g for g in named_genes if g in gene_go_bp and len(gene_go_bp[g]) > 0]
    genes_with_go_idx = {g: gene_to_named_idx[g] for g in genes_with_go}
    go_gene_list = sorted(genes_with_go)
    print(f"  Genes with non-empty GO BP: {len(go_gene_list)}", flush=True)

    # Compute Jaccard for all pairs (can be large, so cap)
    MAX_GO_PAIRS = 20000
    n_go = len(go_gene_list)
    pair_count = 0
    for gi in range(n_go):
        for gj in range(gi+1, n_go):
            ga, gb = go_gene_list[gi], go_gene_list[gj]
            ia, ib = genes_with_go_idx[ga], genes_with_go_idx[gb]
            set_a, set_b = gene_go_bp[ga], gene_go_bp[gb]
            inter = len(set_a & set_b)
            union = len(set_a | set_b)
            jacc = inter / union if union > 0 else 0.0
            go_pairs.append((ia, ib, jacc))
            pair_count += 1
            if pair_count >= MAX_GO_PAIRS:
                break
        if pair_count >= MAX_GO_PAIRS:
            break

    print(f"  GO pairs computed: {len(go_pairs)}", flush=True)

    # Filter to pairs with Jaccard > 0 (some shared GO BP terms)
    go_pairs_pos = [(i, j, jac) for i, j, jac in go_pairs if jac > 0]
    go_pairs_zero = [(i, j, jac) for i, j, jac in go_pairs if jac == 0]
    print(f"  Pairs with Jaccard > 0: {len(go_pairs_pos)}", flush=True)
    print(f"  Pairs with Jaccard = 0: {len(go_pairs_zero)}", flush=True)

    go_valid = len(go_pairs_pos) >= 10 and len(go_pairs_zero) >= 10

    if go_valid:
        go_i1 = np.array([p[0] for p in go_pairs])
        go_i2 = np.array([p[1] for p in go_pairs])
        go_jac = np.array([p[2] for p in go_pairs])

        go_pos_i1 = np.array([p[0] for p in go_pairs_pos])
        go_pos_i2 = np.array([p[1] for p in go_pairs_pos])
        go_pos_jac = np.array([p[2] for p in go_pairs_pos])

        # High vs low GO Jaccard split (median split on positive pairs)
        med_jac = float(np.median(go_pos_jac))
        go_high = [(i, j) for i, j, jac in go_pairs_pos if jac >= med_jac]
        go_low = [(i, j) for i, j, jac in go_pairs_pos if jac < med_jac]
        print(f"  GO high Jaccard (>={med_jac:.3f}): {len(go_high)}, low: {len(go_low)}", flush=True)

        h02_results = []
        for layer in range(N_LAYERS):
            ue = get_unit_emb(layer)

            # Spearman(GO_jaccard, distance) for all pairs
            d_all_go = np.sqrt(np.sum((ue[go_i1] - ue[go_i2])**2, axis=1))
            rho_all, p_all = spearmanr(go_jac, d_all_go)

            # For positive pairs only
            if len(go_pairs_pos) >= 10:
                d_pos_go = np.sqrt(np.sum((ue[go_pos_i1] - ue[go_pos_i2])**2, axis=1))
                rho_pos, p_pos = spearmanr(go_pos_jac, d_pos_go)
            else:
                rho_pos, p_pos = 0.0, 1.0

            # AUROC: high-Jaccard vs zero-Jaccard pairs
            if len(go_high) >= 5 and len(go_pairs_zero) >= 5:
                gh_i1 = np.array([p[0] for p in go_high])
                gh_i2 = np.array([p[1] for p in go_high])
                gz_i1 = np.array([p[0] for p in go_pairs_zero[:len(go_high)*3]])
                gz_i2 = np.array([p[1] for p in go_pairs_zero[:len(go_high)*3]])
                d_high = np.sqrt(np.sum((ue[gh_i1] - ue[gh_i2])**2, axis=1))
                d_zero = np.sqrt(np.sum((ue[gz_i1] - ue[gz_i2])**2, axis=1))
                mw_u, mw_p = mannwhitneyu(d_high, d_zero, alternative='less')
                auroc_high_vs_zero = 1.0 - mw_u / (len(d_high) * len(d_zero))
                effect_hvsz = float(np.mean(d_high) - np.mean(d_zero))
            else:
                auroc_high_vs_zero = 0.5
                mw_p = 1.0
                effect_hvsz = 0.0

            h02_results.append({
                "layer": layer,
                "spearman_rho_all": float(rho_all),
                "spearman_p_all": float(p_all),
                "spearman_rho_positive": float(rho_pos),
                "spearman_p_positive": float(p_pos),
                "auroc_high_vs_zero_jaccard": float(auroc_high_vs_zero),
                "mw_p_high_vs_zero": float(mw_p),
                "effect_high_vs_zero": effect_hvsz,
            })
            print(f"  L{layer:2d}: rho_all={rho_all:.3f}(p={p_all:.3e})"
                  f" rho_pos={rho_pos:.3f}(p={p_pos:.3e})"
                  f" AUROC_hi_vs_0={auroc_high_vs_zero:.3f}",
                  flush=True)

        rhos_all = [r["spearman_rho_all"] for r in h02_results]
        aurocs_h2 = [r["auroc_high_vs_zero_jaccard"] for r in h02_results]
        n_sig = sum(1 for r in h02_results if r["spearman_p_all"] < 0.05)
        n_auroc_sig = sum(1 for r in h02_results if r["mw_p_high_vs_zero"] < 0.05)
        print(f"\n  GO Spearman rho mean={np.mean(rhos_all):.3f}, sig layers={n_sig}/12")
        print(f"  GO AUROC mean={np.mean(aurocs_h2):.3f}, sig layers={n_auroc_sig}/12")

        h02_status = "tested"
        h02_direction = "positive" if np.mean(rhos_all) < -0.01 and n_sig >= 6 else "inconclusive"
        h02_result_value = f"Spearman rho mean={np.mean(rhos_all):.3f}, AUROC mean={np.mean(aurocs_h2):.3f}, {n_sig}/12 sig layers"
    else:
        print("  Insufficient GO pairs for analysis", flush=True)
        h02_results = []
        h02_status = "blocked"
        h02_direction = "inconclusive"
        h02_result_value = f"Insufficient pairs: {len(go_pairs_pos)} positive, {len(go_pairs_zero)} zero"
        n_sig = 0
        rhos_all = [0.0]*12
        aurocs_h2 = [0.5]*12
        med_jac = 0.0

    h02_out = {
        "hypothesis": "H02_GO_BP_proximity",
        "n_genes_with_go_bp": len(gene_go_bp),
        "n_genes_with_nonempty_go": len(genes_with_go),
        "n_go_pairs_total": len(go_pairs),
        "n_go_pairs_positive_jaccard": len(go_pairs_pos),
        "n_go_pairs_zero_jaccard": len(go_pairs_zero),
        "median_positive_jaccard": float(med_jac) if go_valid else 0.0,
        "per_layer": h02_results,
        "summary": {
            "spearman_rho_all_mean": float(np.mean(rhos_all)),
            "auroc_high_vs_zero_mean": float(np.mean(aurocs_h2)),
            "n_layers_spearman_sig": n_sig,
            "n_layers_auroc_sig": n_auroc_sig if go_valid else 0,
            "direction": h02_direction,
        }
    }
    json.dump(h02_out, open(ITER_DIR / "h02_go_bp_proximity.json", "w"), indent=2)
    print("  -> Saved h02_go_bp_proximity.json", flush=True)

except ImportError:
    print("  mygene not available — using cached GO Jaccard from iter_0015", flush=True)
    # Fallback: use iter_0015 h03_go_vs_ppi_driver.json data
    h03_15 = json.load(open(ITER15_DIR / "h03_go_vs_ppi_driver.json"))
    h02_results = []
    for ld in h03_15["layers"]:
        layer = ld["layer"]
        go_only_obs = ld["groups"]["GO-only"]["obs"]
        neither_obs = ld["groups"]["Neither"]["obs"]
        # Fraction in top-K neighbors is proxy for proximity
        # GO-only = pairs with GO Jaccard > median but not in STRING
        h02_results.append({
            "layer": layer,
            "go_only_topK_frac": go_only_obs,
            "neither_topK_frac": neither_obs,
            "effect_go_vs_neither": go_only_obs - neither_obs,
        })

    effects_15 = [r["effect_go_vs_neither"] for r in h02_results]
    h02_status = "tested"
    h02_direction = "positive" if np.mean(effects_15) > 0 else "negative"
    h02_result_value = f"GO-only vs Neither topK enrichment: mean_effect={np.mean(effects_15):.3f} (from iter_0015 fallback)"

    h02_out = {
        "hypothesis": "H02_GO_BP_proximity",
        "fallback": True,
        "source": "iter_0015/h03_go_vs_ppi_driver.json",
        "per_layer": h02_results,
        "summary": {
            "effect_mean": float(np.mean(effects_15)),
            "direction": h02_direction,
        }
    }
    json.dump(h02_out, open(ITER_DIR / "h02_go_bp_proximity.json", "w"), indent=2)
    print("  -> Saved h02_go_bp_proximity.json (fallback)", flush=True)


# ═══════════════════════════════════════════════════════════════════════════════
# H03: TRRUST directional split (Activation vs Repression proximity)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== H03: TRRUST directional split (Activation vs Repression) ===", flush=True)

trrust_pairs_by_dir = {"Activation": [], "Repression": [], "Unknown": []}
with open(TRRUST_PATH) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 3:
            continue
        tf, target, direction = parts[0], parts[1], parts[2]
        if tf in gene_to_named_idx and target in gene_to_named_idx:
            i, j = gene_to_named_idx[tf], gene_to_named_idx[target]
            if i != j:
                key = (min(i,j), max(i,j))
                # Exclude pairs that are also in STRING (same as iter_0022 H01 approach)
                dir_key = direction if direction in trrust_pairs_by_dir else "Unknown"
                trrust_pairs_by_dir[dir_key].append((i, j))

for d, ps in trrust_pairs_by_dir.items():
    print(f"  TRRUST {d}: {len(ps)} pairs (all, incl. STRING overlap)", flush=True)

# Filter to only exclusive (not in STRING) for cleaner test
trrust_excl_by_dir = {}
for direction, pairs in trrust_pairs_by_dir.items():
    excl = [(i,j) for i,j in pairs if (min(i,j), max(i,j)) not in string_set]
    trrust_excl_by_dir[direction] = excl
    print(f"  TRRUST-exclusive {direction}: {len(excl)} pairs", flush=True)

# Main test: activation vs repression embedding distance comparison
act_pairs = trrust_excl_by_dir.get("Activation", [])
rep_pairs = trrust_excl_by_dir.get("Repression", [])
print(f"  Testing Activation (N={len(act_pairs)}) vs Repression (N={len(rep_pairs)})", flush=True)

can_test_dir = len(act_pairs) >= 5 and len(rep_pairs) >= 5
print(f"  Can test directional split: {can_test_dir}", flush=True)

h03_results = []
for layer in range(N_LAYERS):
    ue = get_unit_emb(layer)

    result = {"layer": layer}

    # Overall TRRUST exclusive vs non-STRING (reproduce iter_0022 with all dirs)
    all_excl = act_pairs + rep_pairs + trrust_excl_by_dir.get("Unknown", [])
    if len(all_excl) > 0:
        te_i1 = np.array([p[0] for p in all_excl])
        te_i2 = np.array([p[1] for p in all_excl])
        d_te = np.sqrt(np.sum((ue[te_i1] - ue[te_i2])**2, axis=1))
        d_ns = np.sqrt(np.sum((ue[ns_idx1] - ue[ns_idx2])**2, axis=1))
        mw_all, p_all = mannwhitneyu(d_te, d_ns, alternative='less')
        result["trrust_excl_all_effect"] = float(np.mean(d_te) - np.mean(d_ns))
        result["trrust_excl_all_auroc"] = float(1.0 - mw_all / (len(d_te) * len(d_ns)))
        result["trrust_excl_all_p"] = float(p_all)

    # Directional split
    if can_test_dir:
        act_i1 = np.array([p[0] for p in act_pairs])
        act_i2 = np.array([p[1] for p in act_pairs])
        rep_i1 = np.array([p[0] for p in rep_pairs])
        rep_i2 = np.array([p[1] for p in rep_pairs])

        d_act = np.sqrt(np.sum((ue[act_i1] - ue[act_i2])**2, axis=1))
        d_rep = np.sqrt(np.sum((ue[rep_i1] - ue[rep_i2])**2, axis=1))

        # Activation vs non-STRING
        mw_act_ns, p_act_ns = mannwhitneyu(d_act, d_ns, alternative='less')
        auroc_act = 1.0 - mw_act_ns / (len(d_act) * len(d_ns))

        # Repression vs non-STRING
        mw_rep_ns, p_rep_ns = mannwhitneyu(d_rep, d_ns, alternative='less')
        auroc_rep = 1.0 - mw_rep_ns / (len(d_rep) * len(d_ns))

        # Activation vs Repression (two-sided: are they different?)
        mw_act_rep, p_act_rep = mannwhitneyu(d_act, d_rep, alternative='two-sided')

        result.update({
            "act_mean_dist": float(np.mean(d_act)),
            "rep_mean_dist": float(np.mean(d_rep)),
            "act_vs_nonstring_auroc": float(auroc_act),
            "rep_vs_nonstring_auroc": float(auroc_rep),
            "act_vs_nonstring_p": float(p_act_ns),
            "rep_vs_nonstring_p": float(p_rep_ns),
            "act_vs_rep_p": float(p_act_rep),
            "act_vs_rep_effect": float(np.mean(d_act) - np.mean(d_rep)),
        })
        print(f"  L{layer:2d}: act_AUROC={auroc_act:.3f}(p={p_act_ns:.3e})"
              f" rep_AUROC={auroc_rep:.3f}(p={p_rep_ns:.3e})"
              f" act_vs_rep_p={p_act_rep:.3e} effect={np.mean(d_act)-np.mean(d_rep):.4f}",
              flush=True)
    else:
        print(f"  L{layer:2d}: insufficient pairs for directional split", flush=True)

    h03_results.append(result)

# Summary
if can_test_dir:
    act_aurocs = [r["act_vs_nonstring_auroc"] for r in h03_results]
    rep_aurocs = [r["rep_vs_nonstring_auroc"] for r in h03_results]
    act_rep_ps = [r["act_vs_rep_p"] for r in h03_results]
    act_rep_effects = [r["act_vs_rep_effect"] for r in h03_results]
    print(f"\n  Activation AUROC (vs non-STRING): mean={np.mean(act_aurocs):.3f}")
    print(f"  Repression AUROC (vs non-STRING): mean={np.mean(rep_aurocs):.3f}")
    print(f"  Act vs Rep: effect mean={np.mean(act_rep_effects):.4f}")
    print(f"  Act vs Rep: {sum(1 for p in act_rep_ps if p < 0.05)}/12 layers sig")
    h03_result_value = f"Act AUROC={np.mean(act_aurocs):.3f}, Rep AUROC={np.mean(rep_aurocs):.3f}, act-rep effect={np.mean(act_rep_effects):.4f}"
    h03_direction = "positive" if abs(np.mean(act_rep_effects)) > 0.01 and sum(1 for p in act_rep_ps if p < 0.05) >= 4 else "inconclusive"
else:
    act_aurocs = rep_aurocs = [0.5]*12
    act_rep_effects = [0.0]*12
    act_rep_ps = [1.0]*12
    h03_result_value = f"Act N={len(act_pairs)}, Rep N={len(rep_pairs)}: insufficient for directional split"
    h03_direction = "inconclusive"

h03_out = {
    "hypothesis": "H03_TRRUST_directional_split",
    "n_trrust_excl_activation": len(act_pairs),
    "n_trrust_excl_repression": len(rep_pairs),
    "n_trrust_excl_unknown": len(trrust_excl_by_dir.get("Unknown", [])),
    "can_test_directional_split": can_test_dir,
    "per_layer": h03_results,
    "summary": {
        "activation_auroc_vs_nonstring_mean": float(np.mean(act_aurocs)),
        "repression_auroc_vs_nonstring_mean": float(np.mean(rep_aurocs)),
        "activation_vs_repression_effect_mean": float(np.mean(act_rep_effects)),
        "n_layers_act_sig_vs_nonstring": sum(1 for r in h03_results
                                              if r.get("act_vs_nonstring_p", 1.0) < 0.05),
        "n_layers_rep_sig_vs_nonstring": sum(1 for r in h03_results
                                              if r.get("rep_vs_nonstring_p", 1.0) < 0.05),
        "n_layers_act_vs_rep_sig": sum(1 for p in act_rep_ps if p < 0.05),
        "direction": h03_direction,
    }
}
json.dump(h03_out, open(ITER_DIR / "h03_trrust_directional_split.json", "w"), indent=2)
print("  -> Saved h03_trrust_directional_split.json", flush=True)


# ─── Final summary ─────────────────────────────────────────────────────────────
print("\n=== FINAL SUMMARY ===", flush=True)
print(f"H01: Cell-type expansion AUROC mean={np.mean(aurocs1):.3f}"
      f" (ctrl={np.mean(ctrl_aurocs):.3f if ctrl_aurocs else 0:.3f})"
      f" layer_slope={slope:.4f}(p={pval:.4f})", flush=True)
print(f"H02: GO BP proximity rho_mean={np.mean(rhos_all):.3f}", flush=True)
print(f"H03: TRRUST act={np.mean(act_aurocs):.3f} rep={np.mean(rep_aurocs):.3f}"
      f" act-rep_effect={np.mean(act_rep_effects):.4f}", flush=True)
print("Done.", flush=True)
