"""Independent verifier (A6 / E9): zero-row genes in the SV1 top pole (iter_0006), immune TF/target label counts,
and per-layer z of the iter_0015 SV4/SV5 co-pole test against the brainstormer's success criterion."""
import csv, json
from pathlib import Path
import numpy as np

P41 = Path(__file__).resolve().parents[3]
REPO = P41.parent
S38 = REPO / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs"
out = {}
def named(run):
    g2i = {}; src = set()
    for r in csv.DictReader(open(run / "cycle1_edge_dataset.tsv"), delimiter="\t"):
        g2i.setdefault(r["source"], int(r["source_idx"])); g2i.setdefault(r["target"], int(r["target_idx"])); src.add(r["source"])
    return g2i, src
E = np.load(S38 / "cycle1_main/layer_gene_embeddings.npy", mmap_mode="r")
nz = np.linalg.norm(np.asarray(E[0], dtype=np.float64), axis=1) > 0
g2i, _ = named(S38 / "cycle1_main")
X = np.asarray(E[11], dtype=np.float64); X = X - X.mean(0)
U = np.linalg.svd(X, full_matrices=False)[0]
order = sorted(g2i, key=lambda g: -U[g2i[g], 0]); top = set(order[:52])
zero = [g for g in g2i if not nz[g2i[g]]]
out["kidA_zero_named_in_SV1_top52"] = [sum(g in top for g in zero), len(zero)]
EI = np.load(S38 / "cycle4_immune_main/layer_gene_embeddings.npy", mmap_mode="r")
nzI = np.linalg.norm(np.asarray(EI[0], dtype=np.float64), axis=1) > 0
gI, srcI = named(S38 / "cycle4_immune_main")
obsn = [g for g in gI if nzI[gI[g]]]
out["immune_named_observed_tf_targetonly"] = [len(gI), len(obsn), sum(g in srcI for g in obsn), sum(g not in srcI for g in obsn)]
d = json.load(open(P41 / "iterations/iter_0015/h01_sv45_ppi_copole.json"))
for k in ("SV4", "SV5"):
    z = [l["svd_axes"][k]["z_score"] for l in d["layers"]]; p = [l["svd_axes"][k]["emp_p"] for l in d["layers"]]
    out[f"iter15_{k}"] = dict(n_z_gt_3=sum(x > 3 for x in z), n_p_lt_05=sum(x < 0.05 for x in p), min_z_above_3=min(x for x in z if x > 3))
json.dump(out, open(Path(__file__).with_suffix(".json"), "w"), indent=1)
print(json.dumps(out, indent=1))
