"""H-E only: TRRUST signed regulation in SV3/SV4"""
import numpy as np
import json
import pickle
from pathlib import Path
from collections import defaultdict
from scipy.stats import mannwhitneyu
import csv

ITER_DIR = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0017"
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

emb = np.load(EMB_PATH)
N_LAYERS, N_GENES, N_DIM = emb.shape
gene_to_emb_idx = {}
with open(EDGE_PATH) as f:
    rdr = csv.DictReader(f, delimiter="\t")
    for row in rdr:
        src, si = row["source"], row.get("source_idx", "")
        tgt, ti = row["target"], row.get("target_idx", "")
        if src and si and src not in gene_to_emb_idx:
            gene_to_emb_idx[src] = int(si)
        if tgt and ti and tgt not in gene_to_emb_idx:
            gene_to_emb_idx[tgt] = int(ti)

named_genes = sorted(gene_to_emb_idx.keys())
gene_indices = np.array([gene_to_emb_idx[g] for g in named_genes])
gene_to_local = {g: i for i, g in enumerate(named_genes)}
emb_named = emb[:, gene_indices, :]
emb_mc = emb_named - emb_named.mean(axis=1, keepdims=True)

def get_sv_projections(layer, sv_idx, K=52):
    U, S, Vt = np.linalg.svd(emb_mc[layer], full_matrices=False)
    sv_proj = U[:, sv_idx] * S[sv_idx]
    sorted_idx = np.argsort(sv_proj)
    top_genes = set(named_genes[i] for i in sorted_idx[-K:])
    bot_genes = set(named_genes[i] for i in sorted_idx[:K])
    return sv_proj, top_genes, bot_genes

# Load TRRUST
act_pairs, rep_pairs = [], []
with open(TRRUST_PATH) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        tf, target, reg_type = parts[0], parts[1], parts[2]
        if tf in gene_to_local and target in gene_to_local:
            if "Activation" in reg_type:
                act_pairs.append((tf, target))
            if "Repression" in reg_type:
                rep_pairs.append((tf, target))
act_pairs = list(set(act_pairs))
rep_pairs = list(set(rep_pairs))
print(f"TRRUST act={len(act_pairs)}, rep={len(rep_pairs)}")

he_results = {}
for ax_name, sv_idx in [("SV3", 2), ("SV4", 3)]:
    layer_results = {}
    for layer in range(12):
        sv_proj, top_genes, bot_genes = get_sv_projections(layer, sv_idx, K=52)
        gene_proj = {g: float(sv_proj[gene_to_local[g]]) for g in named_genes}

        act_abs = [abs(gene_proj[tf] - gene_proj[tgt]) for tf, tgt in act_pairs]
        rep_abs = [abs(gene_proj[tf] - gene_proj[tgt]) for tf, tgt in rep_pairs]

        act_copole = [1 if np.sign(gene_proj[tf]) == np.sign(gene_proj[tgt]) else 0 for tf, tgt in act_pairs]
        rep_copole = [1 if np.sign(gene_proj[tf]) == np.sign(gene_proj[tgt]) else 0 for tf, tgt in rep_pairs]
        act_cr = float(np.mean(act_copole))
        rep_cr = float(np.mean(rep_copole))

        stat, p_mwu = mannwhitneyu(act_abs, rep_abs, alternative="less")
        layer_results[layer] = {
            "act_copole_rate": round(act_cr, 4),
            "rep_copole_rate": round(rep_cr, 4),
            "delta": round(act_cr - rep_cr, 4),
            "mwu_p": round(float(p_mwu), 6),
            "significant": p_mwu < 0.05
        }

    n_sig = sum(1 for r in layer_results.values() if r["significant"])
    mean_act = float(np.mean([r["act_copole_rate"] for r in layer_results.values()]))
    mean_rep = float(np.mean([r["rep_copole_rate"] for r in layer_results.values()]))
    print(f"{ax_name}: n_sig={n_sig}/12, mean_act={mean_act:.3f}, mean_rep={mean_rep:.3f}")
    he_results[ax_name] = {
        "sv_idx": sv_idx, "n_act": len(act_pairs), "n_rep": len(rep_pairs),
        "n_layers_significant": n_sig,
        "mean_act_copole": round(mean_act, 4),
        "mean_rep_copole": round(mean_rep, 4),
        "mean_delta": round(mean_act - mean_rep, 4),
        "layer_results": layer_results
    }

with open(ITER_DIR / "h_e_trrust_sv3_sv4_signed.json", "w") as f:
    json.dump(he_results, f, indent=2)
print("Saved h_e_trrust_sv3_sv4_signed.json")

# Update iter0017_results.json
with open(ITER_DIR / "iter0017_results.json") as f:
    r = json.load(f)
r["H_E_trrust_sv3_sv4_signed"] = {
    ax: {"n_layers_sig": he_results[ax]["n_layers_significant"],
         "mean_act_copole": he_results[ax]["mean_act_copole"],
         "mean_rep_copole": he_results[ax]["mean_rep_copole"],
         "mean_delta": he_results[ax]["mean_delta"]}
    for ax in ["SV3","SV4"]
}
with open(ITER_DIR / "iter0017_results.json", "w") as f:
    json.dump(r, f, indent=2)
print("Updated iter0017_results.json")
print(json.dumps(r["H_E_trrust_sv3_sv4_signed"], indent=2))
