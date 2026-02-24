"""
iter_0025 Multi-Hypothesis Screen

H01 (manifold_distance, new_method: Multi-predictor joint model):
    For all named-gene pairs, combine STRING score, Dorothea confidence,
    GO CC Jaccard, GO BP Jaccard into a joint regression predicting -L2_distance.
    Compute partial R^2 for each predictor at each layer to quantify independent
    contributions. Also test VIF for multicollinearity.
    This is the keystone experiment consolidating all 5 biological anchors.

H02 (manifold_distance, new_method: Layer-resolved encoding timeline):
    Compile existing per-layer Spearman curves from prior iteration artifacts
    (STRING, TRRUST, Dorothea, GO CC, GO BP). For each anchor, identify its
    peak layer and test if different biological dimensions peak at different depths.
    Uses only existing artifacts — cost = low.

H03 (null_sensitivity, new_method: Chromosomal proximity negative control):
    Fetch chromosomal locations for the 209 named genes via mygene API.
    Test AUROC for same-chromosome vs different-chromosome pairs in L2 distance.
    Expected null: chromosome should NOT predict embedding proximity more than
    chance if the signal is functional, not genomic. This validates that prior
    positive results (STRING, GO, Dorothea) are not genomic artifacts.
"""

import numpy as np
import json
import csv
from pathlib import Path
from scipy.stats import spearmanr, mannwhitneyu, rankdata
import warnings
warnings.filterwarnings("ignore")

# ─── Paths ────────────────────────────────────────────────────────────────────
PROJECT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
               "/subproject_41_claude_topology_hypothesis_screening_autoloop")
ITER_DIR = PROJECT / "iterations" / "iter_0025"
ITER22_DIR = PROJECT / "iterations" / "iter_0022"
ITER23_DIR = PROJECT / "iterations" / "iter_0023"
ITER24_DIR = PROJECT / "iterations" / "iter_0024"
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

# ─── Pre-compute L2 distances ────────────────────────────────────────────────
print("Pre-computing pairwise L2 distances ...", flush=True)
named_emb = emb[:, named_idx, :]  # [12, N_NAMED, 512]

def pairwise_l2(X):
    """X: [N, D] → upper-triangle L2 distances as flat array"""
    diff = X[:, None, :] - X[None, :, :]
    dist_matrix = np.sqrt((diff**2).sum(axis=-1))
    triu_i, triu_j = np.triu_indices(len(X), k=1)
    return dist_matrix[triu_i, triu_j], triu_i, triu_j

layer_dists = []
triu_i_all = triu_j_all = None
for L in range(N_LAYERS):
    d, ti, tj = pairwise_l2(named_emb[L])
    layer_dists.append(d)
    triu_i_all = ti
    triu_j_all = tj

layer_dists = np.array(layer_dists)  # [12, N_pairs]
N_PAIRS = layer_dists.shape[1]
print(f"  Pairwise L2 shape: {layer_dists.shape}", flush=True)

# Create pair-key lookup for fast assignment
pair_key_to_idx = {}
for k in range(N_PAIRS):
    i, j = triu_i_all[k], triu_j_all[k]
    gi, gj = named_genes[i], named_genes[j]
    pair_key_to_idx[(gi, gj)] = k
    pair_key_to_idx[(gj, gi)] = k


# =============================================================================
# H01: Multi-predictor joint model (partial R^2 + VIF)
# =============================================================================
print("\n=== H01: Multi-predictor joint model ===", flush=True)

# ── Load STRING scores for named gene pairs ──────────────────────────────────
# Reconstruct from edge_dataset.tsv which has string scores
string_scores = {}
with open(EDGE_PATH) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        g1, g2 = row['source'], row['target']
        # Score from STRING (combined_score field)
        try:
            sc = float(row.get('combined_score', row.get('string_score', 0)))
        except (ValueError, TypeError):
            sc = 0.0
        if g1 in gene_to_named_idx and g2 in gene_to_named_idx:
            key = (min(g1, g2), max(g1, g2))
            string_scores[key] = sc

print(f"  STRING pairs loaded: {len(string_scores)}", flush=True)

# ── Load Dorothea confidence ──────────────────────────────────────────────────
dorothea_conf_map = {}  # (g1, g2) -> numeric confidence score
confidence_levels = {'A': 5, 'B': 4, 'C': 3, 'D': 2, 'E': 1}

if DOROTHEA_PATH.exists():
    with open(DOROTHEA_PATH) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            tf = row.get('tf', row.get('TF', ''))
            target = row.get('target', row.get('Target', ''))
            conf = row.get('confidence', row.get('Confidence', 'E'))
            conf_val = confidence_levels.get(conf.upper(), 1)
            if tf in gene_to_named_idx and target in gene_to_named_idx:
                key = (min(tf, target), max(tf, target))
                # Take max confidence if multiple entries
                dorothea_conf_map[key] = max(
                    dorothea_conf_map.get(key, 0), conf_val)
    print(f"  Dorothea pairs: {len(dorothea_conf_map)}", flush=True)
else:
    print("  WARNING: Dorothea file not found, skipping Dorothea predictor", flush=True)

# ── Fetch GO CC and BP Jaccard similarities ───────────────────────────────────
print("  Fetching GO annotations via mygene ...", flush=True)
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
        if isinstance(go_bp, dict):
            go_bp = [go_bp]
        if isinstance(go_cc, dict):
            go_cc = [go_cc]
        gene_go_bp[sym] = set(t['id'] for t in go_bp if 'id' in t)
        gene_go_cc[sym] = set(t['id'] for t in go_cc if 'id' in t)

    print(f"  GO BP genes: {len(gene_go_bp)}, GO CC genes: {len(gene_go_cc)}", flush=True)
    go_available = True
except Exception as e:
    print(f"  WARNING: GO fetch failed: {e}", flush=True)
    go_available = False
    gene_go_bp = {}
    gene_go_cc = {}

def jaccard(s1, s2):
    if not s1 or not s2:
        return np.nan
    return len(s1 & s2) / len(s1 | s2)

# ── Build feature matrix for all pairs ───────────────────────────────────────
print("  Building feature matrix ...", flush=True)

# For each pair: string_score, dorothea_conf, go_cc_jaccard, go_bp_jaccard
# Use 0 for missing STRING/Dorothea, NaN for missing GO
pair_features = []
pair_dist_layer8 = []
valid_pair_idx = []

for k in range(N_PAIRS):
    gi = named_genes[triu_i_all[k]]
    gj = named_genes[triu_j_all[k]]
    key = (min(gi, gj), max(gi, gj))

    str_sc = string_scores.get(key, 0.0)
    dor_conf = dorothea_conf_map.get(key, 0.0)

    # GO Jaccard
    if go_available:
        bp_j = jaccard(gene_go_bp.get(gi, set()), gene_go_bp.get(gj, set()))
        cc_j = jaccard(gene_go_cc.get(gi, set()), gene_go_cc.get(gj, set()))
    else:
        bp_j = np.nan
        cc_j = np.nan

    pair_features.append([str_sc, dor_conf, bp_j if not np.isnan(bp_j) else 0.0,
                          cc_j if not np.isnan(cc_j) else 0.0])
    valid_pair_idx.append(k)

pair_features = np.array(pair_features)  # [N_PAIRS, 4]
feature_names = ['STRING_score', 'Dorothea_conf', 'GO_BP_jaccard', 'GO_CC_jaccard']
print(f"  Feature matrix shape: {pair_features.shape}", flush=True)
print(f"  Non-zero STRING: {(pair_features[:,0]>0).sum()}", flush=True)
print(f"  Non-zero Dorothea: {(pair_features[:,1]>0).sum()}", flush=True)
print(f"  Non-zero GO BP: {(pair_features[:,2]>0).sum()}", flush=True)
print(f"  Non-zero GO CC: {(pair_features[:,3]>0).sum()}", flush=True)

# ── Compute partial R^2 at each layer ─────────────────────────────────────────
def partial_spearman_r2(y, X, feature_idx, n_perms=0):
    """
    Partial R^2 of feature feature_idx after partialling out other features.
    Uses Spearman rank correlation, partialling via residuals of rank regression.
    """
    n_feat = X.shape[1]
    other_idx = [i for i in range(n_feat) if i != feature_idx]

    # Rank y
    y_rank = rankdata(y)

    # Rank the target feature
    feat_rank = rankdata(X[:, feature_idx])

    # If no other features to partial out, just compute Spearman
    if len(other_idx) == 0:
        rho, p = spearmanr(y_rank, feat_rank)
        return rho**2, rho, p

    # Partial out other features using OLS on ranks
    # Residualize y_rank on other features
    X_other = np.column_stack([rankdata(X[:, j]) for j in other_idx])
    X_other_aug = np.column_stack([np.ones(len(X_other)), X_other])

    # OLS: y_rank_resid = y_rank - X_other @ beta
    try:
        beta_y, _, _, _ = np.linalg.lstsq(X_other_aug, y_rank, rcond=None)
        y_resid = y_rank - X_other_aug @ beta_y
    except Exception:
        y_resid = y_rank

    # Residualize feat_rank on other features
    try:
        beta_f, _, _, _ = np.linalg.lstsq(X_other_aug, feat_rank, rcond=None)
        feat_resid = feat_rank - X_other_aug @ beta_f
    except Exception:
        feat_resid = feat_rank

    rho, p = spearmanr(y_resid, feat_resid)
    return rho**2, rho, p

# Negative distance as target (higher = closer = more related)
h01_per_layer = []
print("  Computing partial R^2 at each layer ...", flush=True)

for L in range(N_LAYERS):
    y = -layer_dists[L]  # negative distance: higher = closer

    # Pairwise Spearman for each feature alone
    results_L = {'layer': L}

    # Joint model partial R^2
    partial_r2 = {}
    partial_rho = {}
    partial_p = {}

    for fi, fname in enumerate(feature_names):
        r2, rho, p = partial_spearman_r2(y, pair_features, fi)
        partial_r2[fname] = r2
        partial_rho[fname] = rho
        partial_p[fname] = p

    # Also compute univariate Spearman for comparison
    uni_rho = {}
    uni_p = {}
    for fi, fname in enumerate(feature_names):
        rho, p = spearmanr(y, pair_features[:, fi])
        uni_rho[fname] = rho
        uni_p[fname] = p

    results_L['partial_r2'] = partial_r2
    results_L['partial_rho'] = partial_rho
    results_L['partial_p'] = partial_p
    results_L['univariate_rho'] = uni_rho
    results_L['univariate_p'] = uni_p

    h01_per_layer.append(results_L)

    if L in (0, 8, 11):
        print(f"  L{L} partial R^2: STRING={partial_r2['STRING_score']:.4f}, "
              f"Dor={partial_r2['Dorothea_conf']:.4f}, "
              f"BP={partial_r2['GO_BP_jaccard']:.4f}, "
              f"CC={partial_r2['GO_CC_jaccard']:.4f}", flush=True)

# ── VIF at layer 8 ────────────────────────────────────────────────────────────
print("  Computing VIF at layer 8 ...", flush=True)
try:
    from numpy.linalg import inv
    X8 = np.column_stack([rankdata(pair_features[:, fi]) for fi in range(4)])
    X8_norm = (X8 - X8.mean(0)) / (X8.std(0) + 1e-8)
    corr_mat = np.corrcoef(X8_norm.T)
    # VIF = diagonal of inv(corr_mat)
    vif_vals = np.diag(inv(corr_mat))
    vif_dict = {feature_names[i]: float(vif_vals[i]) for i in range(4)}
    print(f"  VIF: {vif_dict}", flush=True)
except Exception as e:
    print(f"  VIF computation failed: {e}", flush=True)
    vif_dict = {}

# Find which feature has highest partial R^2 at each layer
layer8_partial_r2 = h01_per_layer[8]['partial_r2']
best_predictor = max(layer8_partial_r2, key=layer8_partial_r2.get)
total_partial_r2_L8 = sum(layer8_partial_r2.values())

h01_result = {
    'hypothesis': 'H01_multi_predictor_joint_model',
    'n_pairs': N_PAIRS,
    'feature_names': feature_names,
    'n_nonzero_string': int((pair_features[:,0]>0).sum()),
    'n_nonzero_dorothea': int((pair_features[:,1]>0).sum()),
    'n_nonzero_go_bp': int((pair_features[:,2]>0).sum()),
    'n_nonzero_go_cc': int((pair_features[:,3]>0).sum()),
    'vif_layer8': vif_dict,
    'layer8_partial_r2': layer8_partial_r2,
    'layer8_partial_rho': h01_per_layer[8]['partial_rho'],
    'layer8_partial_p': h01_per_layer[8]['partial_p'],
    'layer8_univariate_rho': h01_per_layer[8]['univariate_rho'],
    'layer8_univariate_p': h01_per_layer[8]['univariate_p'],
    'best_predictor_layer8': best_predictor,
    'total_partial_r2_layer8': float(total_partial_r2_L8),
    'per_layer': h01_per_layer,
    'summary_direction': 'positive' if total_partial_r2_L8 > 0.01 else 'inconclusive'
}

with open(ITER_DIR / 'h01_multi_predictor_joint_model.json', 'w') as f:
    json.dump(h01_result, f, indent=2)
print(f"  Saved h01_multi_predictor_joint_model.json", flush=True)
print(f"  Best predictor at L8: {best_predictor}", flush=True)
print(f"  Total partial R^2 at L8: {total_partial_r2_L8:.4f}", flush=True)


# =============================================================================
# H02: Layer-resolved encoding timeline (from existing artifacts)
# =============================================================================
print("\n=== H02: Layer-resolved encoding timeline ===", flush=True)

# Load per-layer Spearman from iter22 (STRING AUROC / Spearman)
iter22_data = {}
with open(ITER22_DIR / 'h01_string_auroc_trrust_exclusive.json') as f:
    d22 = json.load(f)
for entry in d22.get('per_layer', []):
    L = entry['layer']
    iter22_data[L] = {
        'string_spearman': entry.get('string_spearman_rho', np.nan),
        'string_auroc': entry.get('string_auroc', np.nan)
    }

# Load per-layer from iter23 (TRRUST directional split + GO BP)
with open(ITER23_DIR / 'h02_go_bp_proximity.json') as f:
    d23_go = json.load(f)
iter23_go_data = {}
for entry in d23_go.get('per_layer', []):
    L = entry['layer']
    iter23_go_data[L] = {
        'go_bp_spearman': entry.get('spearman_rho', np.nan)
    }

with open(ITER23_DIR / 'h03_trrust_directional_split.json') as f:
    d23_trrust = json.load(f)
iter23_trrust_data = {}
for entry in d23_trrust.get('per_layer', []):
    L = entry['layer']
    iter23_trrust_data[L] = {
        'trrust_activation_auroc': entry.get('activation_auroc', np.nan)
    }

# Load per-layer from iter24 (Dorothea confidence + GO CC/BP)
with open(ITER24_DIR / 'h01_dorothea_confidence.json') as f:
    d24_dor = json.load(f)
iter24_dor_data = {}
for entry in d24_dor.get('per_layer', []):
    L = entry['layer']
    iter24_dor_data[L] = {
        'dorothea_auroc': entry.get('auroc_high_vs_bg', np.nan)
    }

with open(ITER24_DIR / 'h03_go_ontology_comparison.json') as f:
    d24_go = json.load(f)
iter24_go_data = {}
for entry in d24_go.get('per_layer', []):
    L = entry['layer']
    iter24_go_data[L] = {
        'go_cc_spearman': entry.get('spearman_cc', np.nan),
        'go_bp_spearman_24': entry.get('spearman_bp', np.nan)
    }

# Also use H01 current iteration for multi-predictor partial R^2
# Compile a timeline table
timeline = []
for L in range(N_LAYERS):
    row = {'layer': L}

    # STRING
    row['string_spearman'] = iter22_data.get(L, {}).get('string_spearman', np.nan)
    row['string_auroc'] = iter22_data.get(L, {}).get('string_auroc', np.nan)

    # TRRUST activation
    row['trrust_activation_auroc'] = iter23_trrust_data.get(L, {}).get('trrust_activation_auroc', np.nan)

    # Dorothea
    row['dorothea_auroc'] = iter24_dor_data.get(L, {}).get('dorothea_auroc', np.nan)

    # GO CC
    row['go_cc_spearman'] = iter24_go_data.get(L, {}).get('go_cc_spearman', np.nan)

    # GO BP (use iter24 for consistency)
    row['go_bp_spearman'] = iter24_go_data.get(L, {}).get('go_bp_spearman_24', np.nan)

    # H01 partial R^2 total
    if L < len(h01_per_layer):
        row['joint_partial_r2_total'] = sum(h01_per_layer[L]['partial_r2'].values())
        row['partial_r2_string'] = h01_per_layer[L]['partial_r2']['STRING_score']
        row['partial_r2_dorothea'] = h01_per_layer[L]['partial_r2']['Dorothea_conf']
        row['partial_r2_go_bp'] = h01_per_layer[L]['partial_r2']['GO_BP_jaccard']
        row['partial_r2_go_cc'] = h01_per_layer[L]['partial_r2']['GO_CC_jaccard']

    timeline.append(row)

# Find peak layer for each anchor
anchors = {
    'STRING (Spearman)': [row['string_spearman'] for row in timeline],
    'TRRUST activation (AUROC)': [row['trrust_activation_auroc'] for row in timeline],
    'Dorothea A/B (AUROC)': [row['dorothea_auroc'] for row in timeline],
    'GO CC (Spearman)': [row['go_cc_spearman'] for row in timeline],
    'GO BP (Spearman)': [row['go_bp_spearman'] for row in timeline],
}

peak_layers = {}
for anchor, vals in anchors.items():
    valid = [(i, v) for i, v in enumerate(vals) if not np.isnan(v)]
    if valid:
        peak_layer = max(valid, key=lambda x: abs(x[1]))[0]
        peak_val = max(valid, key=lambda x: abs(x[1]))[1]
        peak_layers[anchor] = {'peak_layer': peak_layer, 'peak_value': float(peak_val)}

print(f"  Peak layers per anchor:", flush=True)
for anchor, info in peak_layers.items():
    print(f"    {anchor}: L{info['peak_layer']} (val={info['peak_value']:.4f})", flush=True)

# Test diversity: do different anchors peak at different layers?
peak_layer_list = [v['peak_layer'] for v in peak_layers.values()]
peak_diversity = len(set(peak_layer_list))
peak_span = max(peak_layer_list) - min(peak_layer_list) if peak_layer_list else 0

# Compute correlation between anchor curves to test independence
anchor_arr = np.array([[v if not np.isnan(v) else 0.0 for v in vals]
                       for _, vals in anchors.items()])  # [n_anchors, 12]
anchor_corr = np.corrcoef(anchor_arr)
mean_inter_anchor_corr = float(np.mean(anchor_corr[np.triu_indices(len(anchors), k=1)]))

h02_result = {
    'hypothesis': 'H02_layer_resolved_encoding_timeline',
    'timeline': timeline,
    'peak_layers': peak_layers,
    'peak_diversity_n_unique': peak_diversity,
    'peak_span_layers': peak_span,
    'mean_inter_anchor_correlation': mean_inter_anchor_corr,
    'n_anchors': len(anchors),
    'summary': f"Peak span {peak_span} layers, {peak_diversity}/{len(anchors)} unique peaks",
    'summary_direction': 'positive' if peak_span >= 3 else 'inconclusive'
}

with open(ITER_DIR / 'h02_layer_encoding_timeline.json', 'w') as f:
    json.dump(h02_result, f, indent=2)
print(f"  Saved h02_layer_encoding_timeline.json", flush=True)
print(f"  Peak span: {peak_span} layers across {peak_diversity} unique peak layers", flush=True)
print(f"  Inter-anchor correlation: {mean_inter_anchor_corr:.4f}", flush=True)


# =============================================================================
# H03: Chromosomal proximity negative control
# =============================================================================
print("\n=== H03: Chromosomal proximity negative control ===", flush=True)

# Fetch chromosomal locations via mygene
print("  Fetching chromosomal data via mygene ...", flush=True)
gene_chrom = {}

try:
    import mygene
    mg = mygene.MyGeneInfo()
    results = mg.querymany(named_genes, scopes='symbol',
                           fields='genomic_pos,chrom,symbol',
                           species='human', returnall=False, verbose=False)

    for r in results:
        if 'notfound' in r and r['notfound']:
            continue
        sym = r.get('query', '')
        if sym not in gene_to_named_idx:
            continue

        # Try genomic_pos first
        gpos = r.get('genomic_pos', None)
        if gpos:
            if isinstance(gpos, list):
                gpos = gpos[0]  # Take first hit
            chrom = str(gpos.get('chr', '')).replace('chr', '').upper()
            if chrom:
                gene_chrom[sym] = chrom
        elif r.get('chrom'):
            gene_chrom[sym] = str(r['chrom']).replace('chr', '').upper()

    print(f"  Chromosomal data retrieved for {len(gene_chrom)}/{N_NAMED} genes", flush=True)
    chr_available = len(gene_chrom) >= 50

except Exception as e:
    print(f"  WARNING: Chromosome fetch failed: {e}", flush=True)
    chr_available = False

if chr_available:
    # Classify pairs as same-chromosome vs different-chromosome
    same_chr_dists = {L: [] for L in range(N_LAYERS)}
    diff_chr_dists = {L: [] for L in range(N_LAYERS)}
    n_same = 0
    n_diff = 0

    for k in range(N_PAIRS):
        gi = named_genes[triu_i_all[k]]
        gj = named_genes[triu_j_all[k]]

        chr_i = gene_chrom.get(gi, None)
        chr_j = gene_chrom.get(gj, None)

        if chr_i is None or chr_j is None:
            continue

        if chr_i == chr_j:
            n_same += 1
            for L in range(N_LAYERS):
                same_chr_dists[L].append(layer_dists[L, k])
        else:
            n_diff += 1
            for L in range(N_LAYERS):
                diff_chr_dists[L].append(layer_dists[L, k])

    print(f"  Same-chromosome pairs: {n_same}, Different-chromosome pairs: {n_diff}", flush=True)

    h03_per_layer = []
    for L in range(N_LAYERS):
        s = np.array(same_chr_dists[L])
        d = np.array(diff_chr_dists[L])

        if len(s) < 5 or len(d) < 5:
            h03_per_layer.append({'layer': L, 'auroc': 0.5, 'mw_p': 1.0,
                                  'mean_same': np.nan, 'mean_diff': np.nan})
            continue

        # AUROC: does same-chr have LOWER distance (more similar)?
        # Positive effect would mean genomic proximity = embedding proximity
        stat, p = mannwhitneyu(d, s, alternative='greater')  # test if diff > same
        n1, n2 = len(d), len(s)
        auroc = stat / (n1 * n2)

        h03_per_layer.append({
            'layer': L,
            'auroc_diff_gt_same': float(auroc),
            'mw_p': float(p),
            'mean_same_chr': float(s.mean()),
            'mean_diff_chr': float(d.mean()),
            'n_same': int(len(s)),
            'n_diff': int(len(d))
        })

        if L in (0, 8, 11):
            print(f"  L{L}: AUROC(diff>same)={auroc:.4f}, p={p:.4e}, "
                  f"mean_same={s.mean():.3f}, mean_diff={d.mean():.3f}", flush=True)

    layer8_auroc_chr = h03_per_layer[8]['auroc_diff_gt_same']
    n_sig_chr = sum(1 for e in h03_per_layer if e.get('mw_p', 1.0) < 0.05)

    # Determine direction: if AUROC ~ 0.5, chromosomal proximity is null (good control)
    # If AUROC >> 0.5, chromosomal proximity IS predictive (potential confound)
    direction = 'negative' if abs(layer8_auroc_chr - 0.5) < 0.05 else 'positive'

    h03_result = {
        'hypothesis': 'H03_chromosomal_proximity_negative_control',
        'n_genes_with_chrom': len(gene_chrom),
        'n_same_chr_pairs': n_same,
        'n_diff_chr_pairs': n_diff,
        'layer8_auroc_diff_gt_same': float(layer8_auroc_chr),
        'n_layers_significant': n_sig_chr,
        'per_layer': h03_per_layer,
        'interpretation': (
            'NULL: chromosomal co-location does NOT predict embedding proximity (positive control for specificity)'
            if abs(layer8_auroc_chr - 0.5) < 0.05
            else f'CONFOUND WARNING: chromosomal proximity predicts embedding proximity (AUROC={layer8_auroc_chr:.3f})'
        ),
        'summary_direction': direction
    }
else:
    print("  FALLBACK: using known chromosomal assignments from literature", flush=True)
    # Minimal fallback: assign random chromosomes to test framework
    known_chr = {
        'TP53': '17', 'BRCA1': '17', 'BRCA2': '13', 'MYC': '8',
        'EGFR': '7', 'PTEN': '10', 'RB1': '13', 'APC': '5',
        'KRAS': '12', 'NRAS': '1', 'BRAF': '7', 'PIK3CA': '3',
        'CDK4': '12', 'CDK6': '7', 'CCND1': '11', 'CDKN2A': '9',
        'MDM2': '12', 'MDM4': '1', 'ATM': '11', 'CHEK2': '22',
        'ERBB2': '17', 'ERBB3': '12', 'MET': '7', 'ALK': '2',
        'ROS1': '6', 'RET': '10', 'FGFR1': '8', 'FGFR2': '10',
        'FGFR3': '4', 'NTRK1': '1', 'NTRK2': '9', 'NTRK3': '15',
    }
    gene_chrom = {g: known_chr[g] for g in named_genes if g in known_chr}
    print(f"  Known chr for {len(gene_chrom)} genes (fallback)", flush=True)

    h03_result = {
        'hypothesis': 'H03_chromosomal_proximity_negative_control',
        'status': 'blocked_insufficient_data',
        'fallback_genes': len(gene_chrom),
        'summary_direction': 'inconclusive'
    }

with open(ITER_DIR / 'h03_chromosomal_proximity_control.json', 'w') as f:
    json.dump(h03_result, f, indent=2)
print(f"  Saved h03_chromosomal_proximity_control.json", flush=True)


# =============================================================================
# Summary
# =============================================================================
print("\n=== SUMMARY ===", flush=True)
print(f"H01 (multi-predictor): total partial R^2 L8={total_partial_r2_L8:.4f}, "
      f"best={best_predictor}", flush=True)
print(f"H02 (timeline): peak span={peak_span}L, {peak_diversity} unique peaks, "
      f"inter-anchor corr={mean_inter_anchor_corr:.3f}", flush=True)
h03_layer8 = h03_result.get('layer8_auroc_diff_gt_same', 'N/A')
print(f"H03 (chromosomal control): L8 AUROC={h03_layer8}", flush=True)

print("\nDone.", flush=True)
