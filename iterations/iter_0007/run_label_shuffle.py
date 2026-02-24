import numpy as np, json, csv, pickle, sys
from pathlib import Path
from scipy.stats import fisher_exact
from collections import defaultdict

ITER_DIR = Path('/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0007')
EMB_PATH = Path('/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main/layer_gene_embeddings.npy')
EDGE_PATH = Path('/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main/cycle1_edge_dataset.tsv')
GENE2GO_PKL = Path('/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/single_cell_mechinterp/data/perturb/gene2go_all.pkl')

RNG = np.random.default_rng(99)
N_SHUFFLE = 500
MIN_TERM = 5; MAX_TERM = 200

print("Loading...", flush=True)
emb = np.load(EMB_PATH)
gene2idx = {}
with open(EDGE_PATH) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        gene2idx[row['source'].strip().upper()] = int(row['source_idx'])
        gene2idx[row['target'].strip().upper()] = int(row['target_idx'])
gene_list = list(gene2idx.keys())
print(f"N genes: {len(gene_list)}", flush=True)

with open(GENE2GO_PKL, "rb") as f:
    gene2go_raw = pickle.load(f)
gene_set = set(gene_list)
term2genes = defaultdict(set)
for gene, terms in gene2go_raw.items():
    g = gene.strip().upper()
    if g not in gene_set: continue
    for t in terms:
        if isinstance(t, dict):
            go_id = t.get("GO_ID", t.get("go_id", ""))
            aspect = t.get("aspect", t.get("Category", "P"))
        else:
            go_id = str(t); aspect = "P"
        if go_id.startswith("GO:") and "P" in str(aspect):
            term2genes[go_id].add(g)
term2genes = {t: gs for t, gs in term2genes.items() if MIN_TERM <= len(gs & gene_set) <= MAX_TERM}
print(f"GO terms: {len(term2genes)}", flush=True)

def go_enrich(top, bot, t2g):
    top_s, bot_s = set(top), set(bot)
    results = []
    for term, ann in t2g.items():
        a = len(ann & top_s); b = len(ann & bot_s)
        c = len(top_s)-a; d = len(bot_s)-b
        if a+b == 0: continue
        odds, p = fisher_exact([[a,b],[c,d]], alternative="greater")
        results.append({"term":term,"a":a,"b":b,"or":odds,"p":p})
    results.sort(key=lambda x: x["p"])
    return results

# SV1 label-shuffle null
print("\n=== SV1 label-shuffle null ===", flush=True)
L11 = emb[11]
L11_c = L11 - L11.mean(axis=0)
U, S, Vt = np.linalg.svd(L11_c, full_matrices=False)
gene_arr = np.array(gene_list)
sv1_vals = np.array([float(U[gene2idx[g], 0] * S[0]) for g in gene_list])
sorted_pos = np.argsort(sv1_vals)
n = len(gene_list); q = n//4
obs_top = gene_arr[sorted_pos[-q:]].tolist()
obs_bot = gene_arr[sorted_pos[:q]].tolist()
obs_res = go_enrich(obs_top, obs_bot, term2genes)
obs_sv1_p = obs_res[0]["p"] if obs_res else 1.0
obs_sv1_term = obs_res[0]["term"] if obs_res else "none"
print(f"  Observed: p={obs_sv1_p:.5f} ({obs_sv1_term})", flush=True)

null_sv1_ps = []
for rep in range(N_SHUFFLE):
    perm = RNG.permutation(n)
    shuf = gene_arr[perm]
    top_g = shuf[sorted_pos[-q:]].tolist()
    bot_g = shuf[sorted_pos[:q]].tolist()
    res = go_enrich(top_g, bot_g, term2genes)
    null_sv1_ps.append(res[0]["p"] if res else 1.0)
    if rep % 100 == 0:
        print(f"  rep {rep}: null_p={null_sv1_ps[-1]:.5f}", flush=True)

null_sv1_ps = np.array(null_sv1_ps)
sv1_null_p5 = float(np.percentile(null_sv1_ps, 5))
sv1_emp_p = float(np.mean(null_sv1_ps <= obs_sv1_p))
sv1_pass = bool(obs_sv1_p < sv1_null_p5)
print(f"  obs_p={obs_sv1_p:.5f}, null_p5={sv1_null_p5:.5f}, emp_p={sv1_emp_p:.4f}, PASS={sv1_pass}", flush=True)
np.save(ITER_DIR / "h01_sv1_label_shuffle_null_ps.npy", null_sv1_ps)

# Drift label-shuffle null
print("\n=== Drift label-shuffle null ===", flush=True)
L0, L11 = emb[0], emb[11]
drift_all = np.linalg.norm(L11 - L0, axis=1)
named_drift_vals = np.array([float(drift_all[gene2idx[g]]) for g in gene_list])
sorted_drift_pos = np.argsort(named_drift_vals)
top50 = gene_arr[sorted_drift_pos[-n//2:]].tolist()
bot50 = gene_arr[sorted_drift_pos[:n//2]].tolist()
obs_drift_res = go_enrich(top50, bot50, term2genes)
obs_drift_p = obs_drift_res[0]["p"] if obs_drift_res else 1.0
obs_drift_term = obs_drift_res[0]["term"] if obs_drift_res else "none"
obs_drift_or = obs_drift_res[0]["or"] if obs_drift_res else 1.0
obs_drift_nsig = sum(1 for r in obs_drift_res if r["p"] < 0.05)
print(f"  Observed: p={obs_drift_p:.5f} ({obs_drift_term}, OR={obs_drift_or:.2f}), n_sig={obs_drift_nsig}", flush=True)

null_drift_ps = []
for rep in range(N_SHUFFLE):
    perm = RNG.permutation(n)
    shuf = gene_arr[perm]
    top_g = shuf[sorted_drift_pos[-n//2:]].tolist()
    bot_g = shuf[sorted_drift_pos[:n//2]].tolist()
    res = go_enrich(top_g, bot_g, term2genes)
    null_drift_ps.append(res[0]["p"] if res else 1.0)
    if rep % 100 == 0:
        print(f"  rep {rep}: null_p={null_drift_ps[-1]:.5f}", flush=True)

null_drift_ps = np.array(null_drift_ps)
drift_null_p5 = float(np.percentile(null_drift_ps, 5))
drift_emp_p = float(np.mean(null_drift_ps <= obs_drift_p))
drift_pass = bool(obs_drift_p < drift_null_p5)
print(f"  obs_p={obs_drift_p:.5f}, null_p5={drift_null_p5:.5f}, emp_p={drift_emp_p:.4f}, PASS={drift_pass}", flush=True)
np.save(ITER_DIR / "h03_drift_label_shuffle_null_ps.npy", null_drift_ps)

summary = {
    "sv1_obs_p": float(obs_sv1_p), "sv1_obs_term": obs_sv1_term,
    "sv1_null_p5": sv1_null_p5, "sv1_null_mean": float(null_sv1_ps.mean()),
    "sv1_empirical_p": sv1_emp_p, "sv1_pass": sv1_pass,
    "drift_obs_p": float(obs_drift_p), "drift_obs_term": obs_drift_term,
    "drift_obs_or": float(obs_drift_or), "drift_obs_nsig_p05": obs_drift_nsig,
    "drift_null_p5": drift_null_p5, "drift_null_mean": float(null_drift_ps.mean()),
    "drift_empirical_p": drift_emp_p, "drift_pass": drift_pass,
    "n_shuffle": N_SHUFFLE
}
with open(ITER_DIR / "label_shuffle_null_summary.json", "w") as f:
    json.dump(summary, f, indent=2)

print("\n=== DONE ===", flush=True)
print(json.dumps(summary, indent=2), flush=True)
