"""iter_0027 screen: H01 PC1 identity, H02 TRRUST polarity, H03 Dorothea threshold sweep"""
import numpy as np
import json
import csv
from pathlib import Path
from scipy.stats import spearmanr, mannwhitneyu
from scipy.stats import rankdata as rank_stats
from sklearn.decomposition import PCA
from collections import defaultdict
import re
import warnings
warnings.filterwarnings("ignore")

PROJECT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
               "/subproject_41_claude_topology_hypothesis_screening_autoloop")
ITER_DIR = PROJECT / "iterations" / "iter_0027"
CYCLE1 = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
              "/subproject_38_geometric_residual_stream_interpretability"
              "/implementation/outputs/cycle1_main")
TRRUST_PATH = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
                   "/single_cell_mechinterp/external/networks/trrust_human.tsv")
DOROTHEA_PATH = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
                     "/single_cell_mechinterp/external/networks/dorothea_human.tsv")
ITER_DIR.mkdir(parents=True, exist_ok=True)

# ─── Load embeddings ─────────────────────────────────────────────────────────
print("Loading embeddings...")
emb = np.load(CYCLE1 / "layer_gene_embeddings.npy")   # [12, 4803, 512]
n_layers, n_genes, n_dim = emb.shape
print(f"  Shape: {emb.shape}")

# Build vocab from cycle1_edge_dataset.tsv or gene list
gene_list_path = CYCLE1 / "gene_list.txt"
if gene_list_path.exists():
    all_genes = [l.strip() for l in open(gene_list_path) if l.strip()]
else:
    edge_path = CYCLE1 / "cycle1_edge_dataset.tsv"
    genes_seen = set()
    with open(edge_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if "source_gene" in row: genes_seen.add(row["source_gene"])
            if "target_gene" in row: genes_seen.add(row["target_gene"])
    all_genes = sorted(genes_seen)
print(f"  Vocab size: {len(all_genes)}")
gene2idx = {g: i for i, g in enumerate(all_genes)}

# Named genes: HGNC symbol pattern
named_genes_in_vocab = [g for g in all_genes
                        if re.match(r'^[A-Za-z][A-Za-z0-9\-]+$', g)
                        and not g.startswith("ENS")]
named_idxs = [gene2idx[g] for g in named_genes_in_vocab]
named_emb = emb[:, named_idxs, :]   # [12, n_named, 512]
print(f"  Named genes: {len(named_genes_in_vocab)}, emb: {named_emb.shape}")

# ─── Load TRRUST ──────────────────────────────────────────────────────────────
trrust_tfs = set()
trrust_pairs = []
with open(TRRUST_PATH) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 3:
            tf, target, direction = parts[0], parts[1], parts[2]
            trrust_tfs.add(tf)
            trrust_pairs.append((tf, target, direction))

# ─── Load Dorothea ────────────────────────────────────────────────────────────
dorothea_tfs = set()
dorothea_pairs = []   # (tf, target, conf_score)
conf_map = {"A": 1.0, "B": 0.75, "C": 0.5, "D": 0.25}
named_set = set(named_genes_in_vocab)
with open(DOROTHEA_PATH) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        tf = row.get("tf") or row.get("source", "")
        tgt = row.get("target", "")
        conf = row.get("confidence", "D")
        dorothea_tfs.add(tf)
        score = conf_map.get(str(conf).strip(), 0.25)
        if tf in gene2idx and tgt in gene2idx:
            dorothea_pairs.append((tf, tgt, score))
all_known_tfs = trrust_tfs | dorothea_tfs
print(f"  Known TFs (TRRUST+Dorothea): {len(all_known_tfs)}")

# ─────────────────────────────────────────────────────────────────────────────
# H01: PC1 Biological Identity at L11
# ─────────────────────────────────────────────────────────────────────────────
print("\n=== H01: PC1 Identity at L11 ===")
n_named = len(named_genes_in_vocab)

L11_emb = named_emb[11]
L11_c = L11_emb - L11_emb.mean(axis=0)
pca11 = PCA(n_components=10)
pca11.fit(L11_c)
pc1_scores = pca11.transform(L11_c)[:, 0]
ev_ratio = pca11.explained_variance_ratio_
print(f"  PC1 EV ratio: {ev_ratio[0]*100:.1f}%  (top5: {[f'{x*100:.1f}%' for x in ev_ratio[:5]]})")

gene_to_named_idx = {g: i for i, g in enumerate(named_genes_in_vocab)}

# Features
is_tf = np.array([1.0 if g in all_known_tfs else 0.0 for g in named_genes_in_vocab])
# Hub: TRRUST target count per TF
hub_degree = defaultdict(int)
for tf, target, direction in trrust_pairs:
    hub_degree[tf] += 1
hub_arr = np.array([float(hub_degree.get(g, 0)) for g in named_genes_in_vocab])

hla_genes = {"HLA-A","HLA-B","HLA-C","HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","HLA-E","HLA-F","HLA-G"}
ap1_genes = {"FOS","FOSB","JUN","JUNB","JUND","FOSL1","FOSL2","ATF3","ATF4"}
is_hla = np.array([1.0 if g in hla_genes else 0.0 for g in named_genes_in_vocab])
is_ap1 = np.array([1.0 if g in ap1_genes else 0.0 for g in named_genes_in_vocab])

n_tf, n_nontf = int(is_tf.sum()), n_named - int(is_tf.sum())
n_hla_found = int(is_hla.sum())
n_ap1_found = int(is_ap1.sum())
print(f"  TF={n_tf}, non-TF={n_nontf}, HLA={n_hla_found}, AP1={n_ap1_found}")

def auroc_mw(a, b):
    stat, pval = mannwhitneyu(a, b, alternative="two-sided")
    return stat / (len(a)*len(b)), pval

results_h01 = {
    "hypothesis": "H01_PC1_L11_identity",
    "pc1_ev_ratio": float(ev_ratio[0]),
    "pc1_top5_ev": [float(x) for x in ev_ratio[:5]],
}

# TF vs non-TF
auroc_tf, pval_tf = auroc_mw(pc1_scores[is_tf==1], pc1_scores[is_tf==0])
results_h01["tf_vs_nontf"] = {"n_tf": n_tf, "n_nontf": n_nontf,
                               "auroc": float(auroc_tf), "p_value": float(pval_tf),
                               "mean_tf": float(pc1_scores[is_tf==1].mean()),
                               "mean_nontf": float(pc1_scores[is_tf==0].mean())}
print(f"  TF vs non-TF: AUROC={auroc_tf:.3f}, p={pval_tf:.4f}")

# HLA
if n_hla_found >= 2:
    auroc_hla, pval_hla = auroc_mw(pc1_scores[is_hla==1], pc1_scores[is_hla==0])
    results_h01["hla_vs_nonhla"] = {"n_hla": n_hla_found, "auroc": float(auroc_hla),
                                     "p_value": float(pval_hla),
                                     "mean_hla": float(pc1_scores[is_hla==1].mean()),
                                     "mean_nonhla": float(pc1_scores[is_hla==0].mean())}
    print(f"  HLA vs non-HLA: AUROC={auroc_hla:.3f}, p={pval_hla:.4f}")

# AP1
if n_ap1_found >= 2:
    auroc_ap1, pval_ap1 = auroc_mw(pc1_scores[is_ap1==1], pc1_scores[is_ap1==0])
    results_h01["ap1_vs_nonap1"] = {"n_ap1": n_ap1_found, "auroc": float(auroc_ap1),
                                     "p_value": float(pval_ap1),
                                     "mean_ap1": float(pc1_scores[is_ap1==1].mean()),
                                     "mean_nonap1": float(pc1_scores[is_ap1==0].mean())}
    print(f"  AP1 vs non-AP1: AUROC={auroc_ap1:.3f}, p={pval_ap1:.4f}")

# Hub degree correlation
hub_nonzero = hub_arr > 0
if hub_nonzero.sum() >= 10:
    rho_hub, p_hub = spearmanr(hub_arr[hub_nonzero], pc1_scores[hub_nonzero])
    results_h01["hub_degree_spearman"] = {"rho": float(rho_hub), "p_value": float(p_hub),
                                           "n": int(hub_nonzero.sum())}
    print(f"  Hub degree (TRRUST n_targets) vs PC1: rho={rho_hub:.3f}, p={p_hub:.4f}, n={hub_nonzero.sum()}")

# Top/bottom genes on PC1
sorted_idxs = np.argsort(pc1_scores)
top20 = [(named_genes_in_vocab[i], float(pc1_scores[i])) for i in sorted_idxs[-20:][::-1]]
bot20 = [(named_genes_in_vocab[i], float(pc1_scores[i])) for i in sorted_idxs[:20]]
results_h01["top20_pc1_genes"] = top20
results_h01["bottom20_pc1_genes"] = bot20
print(f"  Top PC1: {[g for g,_ in top20[:10]]}")
print(f"  Bottom PC1: {[g for g,_ in bot20[:10]]}")

# TF AUROC on PC1 per layer
tf_auroc_per_layer = []
for layer in range(n_layers):
    Lc = named_emb[layer] - named_emb[layer].mean(axis=0)
    pca_l = PCA(n_components=1).fit(Lc)
    pc1_l = pca_l.transform(Lc)[:, 0]
    auroc_l, pval_l = auroc_mw(pc1_l[is_tf==1], pc1_l[is_tf==0])
    tf_auroc_per_layer.append({"layer": layer, "auroc": float(auroc_l), "pval": float(pval_l)})
results_h01["tf_auroc_per_layer"] = tf_auroc_per_layer
_tfl = [round(x["auroc"], 3) for x in tf_auroc_per_layer]
print(f"  TF AUROC on PC1 per layer: {_tfl}")

out = ITER_DIR / "h01_pc1_identity.json"
json.dump(results_h01, open(out, "w"), indent=2)
print(f"  → {out}")


# ─────────────────────────────────────────────────────────────────────────────
# H02: TRRUST Activation/Repression Polarity in Embedding Space
# ─────────────────────────────────────────────────────────────────────────────
print("\n=== H02: TRRUST Activator/Repressor Polarity ===")

# Build named-gene TRRUST pairs
tf_act = defaultdict(list)
tf_rep = defaultdict(list)
act_pairs = []
rep_pairs = []
for tf, target, direction in trrust_pairs:
    if tf in named_set and target in named_set:
        if direction == "Activation":
            tf_act[tf].append(target)
            act_pairs.append((tf, target))
        elif direction == "Repression":
            tf_rep[tf].append(target)
            rep_pairs.append((tf, target))

print(f"  Activation pairs (named): {len(act_pairs)}")
print(f"  Repression pairs (named): {len(rep_pairs)}")

act_targets_unique = list(set(t for _, t in act_pairs))
rep_targets_unique = list(set(t for _, t in rep_pairs))

results_h02 = {
    "hypothesis": "H02_TRRUST_polarity",
    "n_act_pairs": len(act_pairs),
    "n_rep_pairs": len(rep_pairs),
}

# Test 1: Mean PC1 loading of activation targets vs repression targets at L11
act_idxs = [gene_to_named_idx[g] for g in act_targets_unique if g in gene_to_named_idx]
rep_idxs = [gene_to_named_idx[g] for g in rep_targets_unique if g in gene_to_named_idx]

if len(act_idxs) >= 3 and len(rep_idxs) >= 3:
    act_pc1 = pc1_scores[act_idxs]
    rep_pc1 = pc1_scores[rep_idxs]
    auroc_pol, pval_pol = auroc_mw(act_pc1, rep_pc1)
    results_h02["target_polarity_L11_pc1"] = {
        "n_act": len(act_idxs), "n_rep": len(rep_idxs),
        "mean_act_pc1": float(act_pc1.mean()), "mean_rep_pc1": float(rep_pc1.mean()),
        "auroc": float(auroc_pol), "p_value": float(pval_pol)
    }
    print(f"  Act vs Rep targets PC1 (L11): AUROC={auroc_pol:.3f}, p={pval_pol:.4f}")
    print(f"    mean_act={act_pc1.mean():.3f}, mean_rep={rep_pc1.mean():.3f}")

# Test 2: Within-TF direction polarity (TFs with >=2 act and >=2 rep targets)
tfs_both = [tf for tf in tf_act if len(tf_act[tf]) >= 2 and len(tf_rep[tf]) >= 2
            and tf in gene_to_named_idx]
print(f"  TFs with >=2 act and >=2 rep named targets: {tfs_both}")
within_tf = []
for tf in tfs_both:
    act_ts = [t for t in tf_act[tf] if t in gene_to_named_idx]
    rep_ts = [t for t in tf_rep[tf] if t in gene_to_named_idx]
    if len(act_ts) < 2 or len(rep_ts) < 2:
        continue
    delta = pc1_scores[[gene_to_named_idx[t] for t in act_ts]].mean() - \
            pc1_scores[[gene_to_named_idx[t] for t in rep_ts]].mean()
    within_tf.append({"tf": tf, "n_act": len(act_ts), "n_rep": len(rep_ts),
                       "delta_act_minus_rep": float(delta),
                       "act_targets": act_ts, "rep_targets": rep_ts})
    print(f"    {tf}: n_act={len(act_ts)}, n_rep={len(rep_ts)}, delta={delta:.3f}")
results_h02["within_tf_polarity"] = within_tf

# Test 3: Polarity AUROC per layer (target act vs rep on PC1)
polarity_per_layer = []
for layer in range(n_layers):
    Lc = named_emb[layer] - named_emb[layer].mean(axis=0)
    pc1_l = PCA(n_components=1).fit(Lc).transform(Lc)[:, 0]
    if len(act_idxs) >= 3 and len(rep_idxs) >= 3:
        a_s = pc1_l[act_idxs]
        r_s = pc1_l[rep_idxs]
        auroc_l, pval_l = auroc_mw(a_s, r_s)
        polarity_per_layer.append({"layer": layer, "auroc": float(auroc_l),
                                    "pval": float(pval_l),
                                    "mean_act": float(a_s.mean()), "mean_rep": float(r_s.mean())})
    else:
        polarity_per_layer.append({"layer": layer, "auroc": 0.5, "pval": 1.0})
results_h02["polarity_per_layer"] = polarity_per_layer
_pol = [round(x["auroc"], 3) for x in polarity_per_layer]
print(f"  Polarity AUROC per layer: {_pol}")

# Test 4: L2 distance — are activation pairs closer than repression pairs? (supplementary)
if len(act_pairs) >= 5 and len(rep_pairs) >= 5:
    act_dists_l11 = []
    rep_dists_l11 = []
    L11_e = named_emb[11]
    for tf, tgt in act_pairs:
        if tf in gene_to_named_idx and tgt in gene_to_named_idx:
            d = np.linalg.norm(L11_e[gene_to_named_idx[tf]] - L11_e[gene_to_named_idx[tgt]])
            act_dists_l11.append(float(d))
    for tf, tgt in rep_pairs:
        if tf in gene_to_named_idx and tgt in gene_to_named_idx:
            d = np.linalg.norm(L11_e[gene_to_named_idx[tf]] - L11_e[gene_to_named_idx[tgt]])
            rep_dists_l11.append(float(d))
    if len(act_dists_l11) >= 3 and len(rep_dists_l11) >= 3:
        stat_d, pval_d = mannwhitneyu(act_dists_l11, rep_dists_l11, alternative="two-sided")
        auroc_d = stat_d / (len(act_dists_l11) * len(rep_dists_l11))
        results_h02["l2_dist_act_vs_rep_L11"] = {
            "n_act": len(act_dists_l11), "n_rep": len(rep_dists_l11),
            "mean_act_dist": float(np.mean(act_dists_l11)),
            "mean_rep_dist": float(np.mean(rep_dists_l11)),
            "auroc_closer_is_act": float(1 - auroc_d),  # AUROC of act being closer
            "p_value": float(pval_d)
        }
        print(f"  L2 dist act vs rep (L11): mean_act={np.mean(act_dists_l11):.2f}, mean_rep={np.mean(rep_dists_l11):.2f}, p={pval_d:.4f}")

out = ITER_DIR / "h02_trrust_polarity.json"
json.dump(results_h02, open(out, "w"), indent=2)
print(f"  → {out}")


# ─────────────────────────────────────────────────────────────────────────────
# H03: Dorothea Confidence Threshold Sensitivity Sweep + PR elbow alignment
# ─────────────────────────────────────────────────────────────────────────────
print("\n=== H03: Threshold Sensitivity Sweep (Dorothea confidence) ===")

# Dorothea pairs, named only
doro_named = [(tf, tgt, s) for tf, tgt, s in dorothea_pairs
              if tf in named_set and tgt in named_set and tf != tgt]
print(f"  Dorothea named pairs: {len(doro_named)}")

# Threshold sweep: 5 confidence thresholds
thresholds = [0.25, 0.50, 0.625, 0.75, 1.00]   # maps to D-and-above, C-and-above, B-and-above, B-only, A-only
threshold_labels = ["D+", "C+", "B+/above-midpoint", "B+", "A-only"]

rng = np.random.default_rng(42)
results_h03 = {
    "hypothesis": "H03_dorothea_threshold_sweep",
    "n_dorothea_named_total": len(doro_named),
    "thresholds": thresholds,
    "threshold_labels": threshold_labels,
}

# Precompute null pairs once
n_null_base = 500
null_i1 = rng.integers(0, n_named, n_null_base * 2)
null_i2 = rng.integers(0, n_named, n_null_base * 2)
valid_null = null_i1 != null_i2
null_i1 = null_i1[valid_null][:n_null_base]
null_i2 = null_i2[valid_null][:n_null_base]

threshold_results = []
for thresh, label in zip(thresholds, threshold_labels):
    pos_pairs = [(tf, tgt, s) for tf, tgt, s in doro_named if s >= thresh]
    if len(pos_pairs) < 3:
        threshold_results.append({"threshold": thresh, "label": label, "n_pos": len(pos_pairs), "best_auroc": None})
        continue
    n_pos = len(pos_pairs)
    p_i1 = np.array([gene_to_named_idx[tf] for tf, tgt, s in pos_pairs])
    p_i2 = np.array([gene_to_named_idx[tgt] for tf, tgt, s in pos_pairs])

    layer_aurocs = []
    for layer in range(n_layers):
        Le = named_emb[layer]
        d_pos = np.linalg.norm(Le[p_i1] - Le[p_i2], axis=1)
        d_null = np.linalg.norm(Le[null_i1] - Le[null_i2], axis=1)
        labels_bin = np.concatenate([np.ones(n_pos), np.zeros(len(null_i1))])
        scores_bin = -np.concatenate([d_pos, d_null])
        ranks = rank_stats(scores_bin)
        pos_rank_sum = ranks[labels_bin == 1].sum()
        auroc = (pos_rank_sum - n_pos*(n_pos+1)/2) / (n_pos * len(null_i1))
        layer_aurocs.append(float(auroc))

    best_auroc = max(layer_aurocs)
    best_layer = int(np.argmax(layer_aurocs))
    threshold_results.append({
        "threshold": thresh, "label": label, "n_pos": n_pos,
        "best_auroc": float(best_auroc), "best_layer": best_layer,
        "per_layer_auroc": layer_aurocs
    })
    print(f"  {label} (>=score {thresh}): n={n_pos}, best_AUROC={best_auroc:.3f} @ L{best_layer}")

results_h03["threshold_results"] = threshold_results

# PR collapse data from iter_0026
pr_by_layer = {0:21.07,1:14.12,2:11.19,3:10.11,4:7.73,5:6.88,
               6:6.12,7:5.35,8:4.61,9:3.22,10:2.06,11:1.68}

# Elbow alignment: compare best-AUROC layer vs PR at that layer
best_layers = [r["best_layer"] for r in threshold_results if r.get("best_layer") is not None]
best_aurocs = [r["best_auroc"] for r in threshold_results if r.get("best_auroc") is not None]
pr_at_best = [pr_by_layer[l] for l in best_layers]
thresh_vals = [r["threshold"] for r in threshold_results if r.get("best_layer") is not None]

if len(best_layers) >= 3:
    rho_th, p_th = spearmanr(thresh_vals, best_layers)
    rho_pr, p_pr = spearmanr(pr_at_best, best_aurocs)
    results_h03["spearman_threshold_vs_best_layer"] = {"rho": float(rho_th), "p": float(p_th)}
    results_h03["spearman_pr_vs_best_auroc"] = {"rho": float(rho_pr), "p": float(p_pr)}
    results_h03["pr_at_best_auroc_layers"] = pr_at_best
    results_h03["best_auroc_layers"] = best_layers
    results_h03["best_aurocs"] = best_aurocs
    print(f"  Spearman(threshold, best_layer): rho={rho_th:.3f}, p={p_th:.4f}")
    print(f"  Spearman(PR_at_best_layer, best_auroc): rho={rho_pr:.3f}, p={p_pr:.4f}")
    print(f"  Best layers: {best_layers}, PR there: {[f'{x:.1f}' for x in pr_at_best]}")

out = ITER_DIR / "h03_threshold_sweep.json"
json.dump(results_h03, open(out, "w"), indent=2)
print(f"  → {out}")

print("\n=== All 3 hypotheses complete ===")
