"""Independent re-computation of the iter_0006 'SV1 extracellular OR = 6.37' ledger number and of its zero-row
dependence, written from the iter_0006 source (run_h02_svd_biology) without importing E11 code.
Output: verify/verify_ledger.json"""
import csv, json, pickle
from pathlib import Path
import numpy as np
from scipy.stats import fisher_exact

HERE = Path(__file__).resolve().parent; REPO = HERE.parents[3]
S38 = REPO / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main"
E = np.load(S38 / "layer_gene_embeddings.npy", mmap_mode="r")
G2GO = pickle.load(open(REPO / "single_cell_mechinterp/data/perturb/gene2go_all.pkl", "rb"))
g2i = {}
for r in csv.DictReader(open(S38 / "cycle1_edge_dataset.tsv"), delimiter="\t"):
    g2i[r["source"]] = int(r["source_idx"]); g2i[r["target"]] = int(r["target_idx"])
names = list(g2i)
X = np.asarray(E[11]); Xc = X - X.mean(0, keepdims=True)
U, s, _ = np.linalg.svd(Xc, full_matrices=False)
obs = np.linalg.norm(np.asarray(E[0], np.float64), axis=1) > 0


def quart(gs, proj):
    o = sorted(zip(gs, proj), key=lambda t: -t[1]); n25 = max(5, len(o) // 4)
    top = {g for g, _ in o[:n25]}; bot = {g for g, _ in o[-n25:]}
    tg = {g for g in gs if "GO:0005615" in G2GO.get(g, ())}
    a, c = len(tg & top), len(tg & bot)
    OR, p = fisher_exact([[a, n25 - a], [c, n25 - c]], alternative="greater")
    return dict(n=len(gs), n25=n25, a=a, c=c, OR=float(OR), p=float(p),
                zero_in_top=int(sum(not obs[g2i[g]] for g in top)), zero_in_bot=int(sum(not obs[g2i[g]] for g in bot)))


res = dict(n_named=len(names), n_named_zero=int(sum(not obs[g2i[g]] for g in names)),
           r_SV1_vs_observed=float(np.corrcoef(U[:, 0], obs.astype(float))[0, 1]),
           own=quart(names, [U[g2i[g], 0] for g in names]))
keep = [g for g in names if obs[g2i[g]]]
res["own_drop_zero_rows"] = quart(keep, [U[g2i[g], 0] for g in keep])
res["own_drop_zero_rows_flipped"] = quart(keep, [-U[g2i[g], 0] for g in keep])
zr = [g for g in names if not obs[g2i[g]]]
res["ext_rate_zero_named"] = float(np.mean(["GO:0005615" in G2GO.get(g, ()) for g in zr]))
res["ext_rate_obs_named"] = float(np.mean(["GO:0005615" in G2GO.get(g, ()) for g in keep]))
(HERE / "verify_ledger.json").write_text(json.dumps(res, indent=1))
print(json.dumps(res, indent=1))
