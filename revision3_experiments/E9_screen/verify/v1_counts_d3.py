"""Independent verifier (A6 / E9): record counts, family counts, and the D3 gene-row misalignment.
Reads only raw inputs (iteration JSONs, gene_list.txt, h5ad var index, edge dataset, embedding). No e9_lib.
"""
import csv, json, re
from pathlib import Path
import numpy as np, h5py

P41 = Path(__file__).resolve().parents[3]
REPO = P41.parent
ITER = P41 / "iterations"
CYC1 = REPO / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main"
H5K = REPO / "single_cell_mechinterp/outputs/tabula_sapiens_processed.h5ad"
out = {}

# ---- record counts ----
dirs = sorted(d for d in ITER.iterdir() if d.is_dir() and d.name.startswith("iter_"))
files = [d / "executor_hypothesis_screen.json" for d in dirs if (d / "executor_hypothesis_screen.json").exists()]
recs = []
for f in files:
    for h in json.load(open(f)).get("hypotheses", []):
        recs.append((f.parent.name, h.get("id"), h.get("family", "unknown")))
num = lambda it: int(it[5:9])
r183 = [r for r in recs if num(r[0]) <= 63]
out["n_iteration_dirs"] = len(dirs)
out["dirs_without_json"] = [d.name for d in dirs if not (d / "executor_hypothesis_screen.json").exists()]
out["n_json_files"] = len(files)
out["n_records"] = len(recs)
out["n_records_excl_incomplete"] = sum("incomplete" not in r[0] for r in recs)
out["n_183"] = len(r183)
out["n_json_files_1_63"] = len({r[0] for r in r183})
out["families_183"] = sorted({r[2] for r in r183})
out["families_all"] = sorted({r[2] for r in recs})

# ---- D3: gene_list.txt vs h5ad var; misaligned vs true mapping ----
lines = [l.rstrip("\n").strip() for l in open(CYC1 / "gene_list.txt")]
with h5py.File(H5K, "r") as f:
    var = [g.decode() if isinstance(g, bytes) else str(g) for g in f["var"]["_index"][:]]
nonblank = [(i, g) for i, g in enumerate(lines) if g]
out["gene_list_lines"] = len(lines); out["gene_list_nonblank"] = len(nonblank)
out["gene_list_line_i_equals_var_i"] = all(var[i] == g for i, g in nonblank)
ed = list(csv.DictReader(open(CYC1 / "cycle1_edge_dataset.tsv"), delimiter="\t"))
true = {}
for r in ed:
    true.setdefault(r["source"], int(r["source_idx"])); true.setdefault(r["target"], int(r["target_idx"]))
out["edge_idx_consistent_with_var"] = all(var[i] == g for g, i in true.items())
out["n_named"] = len(true)
mis = {g: k for k, g in enumerate([g for _, g in nonblank])}   # as in iters 27, 29-40
emb = np.load(CYC1 / "layer_gene_embeddings.npy", mmap_mode="r")
nz = np.linalg.norm(np.asarray(emb[0], dtype=np.float64), axis=1) > 0
out["kidney_A_n_observed"] = int(nz.sum()); out["kidney_A_G"] = int(emb.shape[1])
out["n_named_on_own_row_under_misaligned"] = int(sum(mis[g] == true[g] for g in true if g in mis))
out["misaligned_rows_range"] = [min(mis.values()), max(mis.values())]
zero_mis = sorted(g for g in true if not nz[mis[g]]); zero_true = sorted(g for g in true if not nz[true[g]])
out["n_zero_under_misaligned"] = len(zero_mis); out["zero_under_misaligned"] = zero_mis
out["n_zero_under_true"] = len(zero_true); out["zero_under_true"] = zero_true
out["n_observed_under_true"] = len(true) - len(zero_true)
# what iter_0030 actually called OOV
txt = (ITER / "iter_0030/executor_iteration_report.md").read_text()
out["iter30_report_mentions_195"] = "195" in txt

# ---- iter_0015 exact Spearman p, n=5, rho=1 ----
from itertools import permutations
from scipy.stats import spearmanr
x = np.arange(5); cnt = 0
for p in permutations(range(5)):
    if abs(spearmanr(x, p)[0]) >= 1 - 1e-12: cnt += 1
out["iter15_exact_two_sided_p_n5_rho1"] = cnt / 120
out["iter15_scipy_p"] = float(spearmanr([1, 2, 3, 4, 5], [1.00, 1.48, 2.09, 3.18, 5.02])[1])

json.dump(out, open(Path(__file__).with_suffix(".json"), "w"), indent=1)
print(json.dumps(out, indent=1))
