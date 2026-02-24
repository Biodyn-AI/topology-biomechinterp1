"""
iter_0025 H01 FIX: Multi-predictor joint model with correct data loading.
STRING from string_ppi_score04_cache.json, Dorothea with correct column names.
"""

import numpy as np
import json
import csv
from pathlib import Path
from scipy.stats import spearmanr, rankdata
import warnings
warnings.filterwarnings("ignore")

PROJECT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
               "/subproject_41_claude_topology_hypothesis_screening_autoloop")
ITER_DIR = PROJECT / "iterations" / "iter_0025"
ITER15_DIR = PROJECT / "iterations" / "iter_0015"
CYCLE1 = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
              "/subproject_38_geometric_residual_stream_interpretability"
              "/implementation/outputs/cycle1_main")
DOROTHEA_PATH = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
                     "/single_cell_mechinterp/external/networks/dorothea_human.tsv")

# Load embeddings
print("Loading embeddings ...", flush=True)
emb = np.load(CYCLE1 / "layer_gene_embeddings.npy")  # [12, 4803, 512]
N_LAYERS = emb.shape[0]

with open(CYCLE1 / "gene_list.txt") as f:
    vocab_genes = [line.strip() for line in f]
gene_to_emb_idx = {g: i for i, g in enumerate(vocab_genes) if g}

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
print(f"  Named genes: {N_NAMED}", flush=True)

# Pre-compute L2 distances
named_emb = emb[:, named_idx, :]
def pairwise_l2(X):
    diff = X[:, None, :] - X[None, :, :]
    dist_matrix = np.sqrt((diff**2).sum(axis=-1))
    triu_i, triu_j = np.triu_indices(len(X), k=1)
    return dist_matrix[triu_i, triu_j], triu_i, triu_j

print("Computing pairwise L2 ...", flush=True)
layer_dists = []
triu_i_all = triu_j_all = None
for L in range(N_LAYERS):
    d, ti, tj = pairwise_l2(named_emb[L])
    layer_dists.append(d)
    triu_i_all = ti
    triu_j_all = tj
layer_dists = np.array(layer_dists)  # [12, N_pairs]
N_PAIRS = layer_dists.shape[1]
print(f"  Pairs: {N_PAIRS}", flush=True)

# Load STRING pairs (named gene pairs only)
print("Loading STRING data ...", flush=True)
string_cache = json.load(open(ITER15_DIR / "string_ppi_score04_cache.json"))
string_score_map = {}
for p in string_cache['pairs']:
    g1, g2 = p['g1'], p['g2']
    sc = float(p['score'])
    if g1 in gene_to_named_idx and g2 in gene_to_named_idx:
        key = (min(g1, g2), max(g1, g2))
        string_score_map[key] = sc
print(f"  STRING pairs (named genes): {len(string_score_map)}", flush=True)

# Load Dorothea (correct column names: source, target, confidence)
print("Loading Dorothea data ...", flush=True)
confidence_levels = {'A': 5, 'B': 4, 'C': 3, 'D': 2, 'E': 1}
dorothea_conf_map = {}
with open(DOROTHEA_PATH) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        tf = row.get('source', '')
        target = row.get('target', '')
        conf = row.get('confidence', 'E')
        conf_val = confidence_levels.get(conf.upper(), 1)
        if tf in gene_to_named_idx and target in gene_to_named_idx:
            key = (min(tf, target), max(tf, target))
            dorothea_conf_map[key] = max(dorothea_conf_map.get(key, 0), conf_val)
print(f"  Dorothea pairs (named genes): {len(dorothea_conf_map)}", flush=True)

# Load GO annotations (from mygene - reuse logic from iter_0024)
print("Fetching GO annotations ...", flush=True)
try:
    import mygene
    mg = mygene.MyGeneInfo()
    results = mg.querymany(named_genes, scopes='symbol', fields='go.BP,go.CC',
                           species='human', returnall=False, verbose=False)
    gene_go_bp = {}
    gene_go_cc = {}
    for r in results:
        if 'notfound' in r and r['notfound']:
            continue
        sym = r.get('query', '')
        if sym not in gene_to_named_idx:
            continue
        go_bp = r.get('go', {}).get('BP', [])
        go_cc = r.get('go', {}).get('CC', [])
        if isinstance(go_bp, dict): go_bp = [go_bp]
        if isinstance(go_cc, dict): go_cc = [go_cc]
        gene_go_bp[sym] = set(t['id'] for t in go_bp if 'id' in t)
        gene_go_cc[sym] = set(t['id'] for t in go_cc if 'id' in t)
    print(f"  GO BP: {len(gene_go_bp)} genes, GO CC: {len(gene_go_cc)} genes", flush=True)
    go_available = True
except Exception as e:
    print(f"  GO fetch failed: {e}", flush=True)
    go_available = False
    gene_go_bp = {}
    gene_go_cc = {}

def jaccard(s1, s2):
    if not s1 or not s2:
        return 0.0
    return len(s1 & s2) / len(s1 | s2)

# Build feature matrix
print("Building feature matrix ...", flush=True)
feature_names = ['STRING_score', 'Dorothea_conf', 'GO_BP_jaccard', 'GO_CC_jaccard']
pair_features = np.zeros((N_PAIRS, 4))

for k in range(N_PAIRS):
    gi = named_genes[triu_i_all[k]]
    gj = named_genes[triu_j_all[k]]
    key = (min(gi, gj), max(gi, gj))

    pair_features[k, 0] = string_score_map.get(key, 0.0)
    pair_features[k, 1] = dorothea_conf_map.get(key, 0.0)

    if go_available:
        pair_features[k, 2] = jaccard(gene_go_bp.get(gi, set()), gene_go_bp.get(gj, set()))
        pair_features[k, 3] = jaccard(gene_go_cc.get(gi, set()), gene_go_cc.get(gj, set()))

print(f"  Feature matrix shape: {pair_features.shape}", flush=True)
for fi, fn in enumerate(feature_names):
    print(f"  {fn}: {(pair_features[:,fi]>0).sum()} non-zero", flush=True)

# Partial R^2 via rank regression residuals
def partial_spearman_r2(y, X, feature_idx):
    n_feat = X.shape[1]
    other_idx = [i for i in range(n_feat) if i != feature_idx]

    y_rank = rankdata(y).astype(float)
    feat_rank = rankdata(X[:, feature_idx]).astype(float)

    if len(other_idx) == 0:
        rho, p = spearmanr(y_rank, feat_rank)
        return rho**2, rho, p

    X_other = np.column_stack([rankdata(X[:, j]).astype(float) for j in other_idx])
    X_other_aug = np.column_stack([np.ones(len(X_other)), X_other])

    try:
        beta_y, _, _, _ = np.linalg.lstsq(X_other_aug, y_rank, rcond=None)
        y_resid = y_rank - X_other_aug @ beta_y
        beta_f, _, _, _ = np.linalg.lstsq(X_other_aug, feat_rank, rcond=None)
        feat_resid = feat_rank - X_other_aug @ beta_f
    except Exception:
        y_resid = y_rank
        feat_resid = feat_rank

    rho, p = spearmanr(y_resid, feat_resid)
    return rho**2, rho, p

print("Computing partial R^2 per layer ...", flush=True)
h01_per_layer = []
for L in range(N_LAYERS):
    y = -layer_dists[L]  # negative distance: higher = more similar

    partial_r2 = {}
    partial_rho = {}
    partial_p = {}
    uni_rho = {}
    uni_p = {}

    for fi, fname in enumerate(feature_names):
        r2, rho, p = partial_spearman_r2(y, pair_features, fi)
        partial_r2[fname] = float(r2)
        partial_rho[fname] = float(rho)
        partial_p[fname] = float(p)
        rho_uni, p_uni = spearmanr(y, pair_features[:, fi])
        uni_rho[fname] = float(rho_uni)
        uni_p[fname] = float(p_uni)

    h01_per_layer.append({
        'layer': L,
        'partial_r2': partial_r2,
        'partial_rho': partial_rho,
        'partial_p': partial_p,
        'univariate_rho': uni_rho,
        'univariate_p': uni_p
    })

    if L in (0, 4, 8, 11):
        pr2 = partial_r2
        print(f"  L{L}: partial R^2: STRING={pr2['STRING_score']:.5f}, "
              f"Dor={pr2['Dorothea_conf']:.5f}, "
              f"BP={pr2['GO_BP_jaccard']:.5f}, "
              f"CC={pr2['GO_CC_jaccard']:.5f} | "
              f"total={sum(pr2.values()):.5f}", flush=True)

# VIF computation
print("Computing VIF at layer 8 ...", flush=True)
try:
    X8 = np.column_stack([rankdata(pair_features[:, fi]).astype(float) for fi in range(4)])
    X8_norm = (X8 - X8.mean(0)) / (X8.std(0) + 1e-10)
    corr_mat = np.corrcoef(X8_norm.T)
    vif_vals = np.diag(np.linalg.inv(corr_mat))
    vif_dict = {feature_names[i]: float(vif_vals[i]) for i in range(4)}
    print(f"  VIF: {vif_dict}", flush=True)
except Exception as e:
    print(f"  VIF failed: {e}", flush=True)
    vif_dict = {}

# Summary
layer8_pr2 = h01_per_layer[8]['partial_r2']
best_predictor = max(layer8_pr2, key=layer8_pr2.get)
total_pr2_L8 = sum(layer8_pr2.values())
uni_rho_L8 = h01_per_layer[8]['univariate_rho']

# Determine direction
sig_predictors = sum(1 for fname in feature_names
                     if abs(h01_per_layer[8]['partial_rho'].get(fname, 0)) > 0.01
                     and h01_per_layer[8]['partial_p'].get(fname, 1) < 0.01)

print(f"\nSummary:", flush=True)
print(f"  Best predictor at L8: {best_predictor} (partial R^2={layer8_pr2[best_predictor]:.5f})", flush=True)
print(f"  Total partial R^2 at L8: {total_pr2_L8:.5f}", flush=True)
print(f"  Sig predictors at L8: {sig_predictors}", flush=True)
print(f"  Univariate Spearman at L8: STRING={uni_rho_L8['STRING_score']:.4f}, "
      f"Dor={uni_rho_L8['Dorothea_conf']:.4f}, "
      f"BP={uni_rho_L8['GO_BP_jaccard']:.4f}, "
      f"CC={uni_rho_L8['GO_CC_jaccard']:.4f}", flush=True)

result = {
    'hypothesis': 'H01_multi_predictor_joint_model',
    'n_pairs_total': N_PAIRS,
    'feature_names': feature_names,
    'n_nonzero': {fn: int((pair_features[:, fi] > 0).sum())
                  for fi, fn in enumerate(feature_names)},
    'vif_layer8': vif_dict,
    'layer8_partial_r2': layer8_pr2,
    'layer8_partial_rho': h01_per_layer[8]['partial_rho'],
    'layer8_partial_p': h01_per_layer[8]['partial_p'],
    'layer8_univariate_rho': uni_rho_L8,
    'layer8_univariate_p': h01_per_layer[8]['univariate_p'],
    'best_predictor_layer8': best_predictor,
    'total_partial_r2_layer8': float(total_pr2_L8),
    'n_sig_predictors_layer8': sig_predictors,
    'per_layer': h01_per_layer,
    'summary_direction': 'positive' if sig_predictors >= 2 else ('inconclusive' if sig_predictors == 1 else 'negative')
}

with open(ITER_DIR / 'h01_multi_predictor_joint_model.json', 'w') as f:
    json.dump(result, f, indent=2)
print("Saved h01_multi_predictor_joint_model.json", flush=True)
