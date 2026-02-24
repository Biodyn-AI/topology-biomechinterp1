"""
iter_0005 - Multi-Hypothesis Screen

H01: GO term embedding clustering (full gene vocabulary approach using 209-gene index map)
     - For each GO Biological Process term with >= 5 genes in our vocabulary,
       compute mean pairwise cosine distance within the GO group vs a
       null distribution (random same-size groups from 209-gene pool).
     - Test across all 12 layers; report z-score and fraction of GO terms
       showing significant clustering.

H02: Per-gene residual drift analysis + GO enrichment
     - Compute L2 norm of (layer_11 - layer_0) per gene = "residual drift".
     - Test whether top-drift genes are enriched in specific GO categories
       vs bottom-drift genes (Fisher's exact test on GO terms).

H03: Effective rank (spectral complexity) per layer
     - Compute effective rank = exp(entropy of squared singular values)
       for the embedding submatrix of the 209 known genes at each layer.
     - Compare to TwoNN ID from iter_0004.

Data: scGPT lung embeddings, cycle1, 12 layers x 4803 genes x 512 dims.
      Known gene->index mapping: 209 genes from cycle1_edge_dataset.tsv
"""

import numpy as np
import json
import csv
import os
import sys
import urllib.request
import urllib.parse
import time
from pathlib import Path
from scipy.spatial.distance import cdist
from scipy.stats import fisher_exact, pearsonr, spearmanr
from collections import defaultdict

ITER_DIR = Path(__file__).parent
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

RNG = np.random.default_rng(42)
N_NULL_REPS = 20
MIN_TERM_SIZE = 5
MAX_TERM_SIZE = 200


def cosine_dist_matrix(X):
    norms = np.linalg.norm(X, axis=1, keepdims=True)
    norms = np.where(norms < 1e-10, 1.0, norms)
    X_n = X / norms
    sim = X_n @ X_n.T
    return 1.0 - np.clip(sim, -1, 1)


def mean_pairwise_cosine_dist(X_norm, idxs):
    sub = X_norm[idxs]
    D = cosine_dist_matrix(sub)
    n = len(idxs)
    mask = np.triu(np.ones((n, n), dtype=bool), k=1)
    return float(D[mask].mean()) if mask.sum() > 0 else np.nan


def effective_rank(X):
    """exp(entropy of normalized squared singular values)."""
    _, S, _ = np.linalg.svd(X, full_matrices=False)
    s2 = S ** 2
    s2_norm = s2 / s2.sum()
    s2_norm = s2_norm[s2_norm > 1e-12]
    H = -np.sum(s2_norm * np.log(s2_norm))
    return float(np.exp(H))


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
print("Loading embeddings...", flush=True)
emb = np.load(EMB_PATH)  # (12, 4803, 512)
N_LAYERS, N_TOTAL_GENES, DIM = emb.shape
print(f"  embedding shape: {N_LAYERS} x {N_TOTAL_GENES} x {DIM}")

# Build gene -> embedding index from edge dataset
gene_to_emb_idx = {}
with open(EDGE_PATH) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        gene_to_emb_idx[row['source']] = int(row['source_idx'])
        gene_to_emb_idx[row['target']] = int(row['target_idx'])

gene_list = sorted(gene_to_emb_idx.keys())
N_GENES = len(gene_list)
gene_emb_idxs = np.array([gene_to_emb_idx[g] for g in gene_list])
print(f"  Known genes: {N_GENES}, embedding idx range: {gene_emb_idxs.min()}-{gene_emb_idxs.max()}")


# ---------------------------------------------------------------------------
# H03: Effective rank per layer (fast)
# ---------------------------------------------------------------------------
print("\n=== H03: Effective rank per layer ===", flush=True)

twonn_path = ITER_DIR.parent / "iter_0004" / "h01_intrinsic_dim_per_layer.csv"
twonn_data = {}
with open(twonn_path) as f:
    reader = csv.DictReader(f)
    for row in reader:
        twonn_data[int(row['layer'])] = float(row['intrinsic_dim_twonn'])

h03_rows = []
for layer in range(N_LAYERS):
    X_sub = emb[layer][gene_emb_idxs]  # (209, 512)
    er = effective_rank(X_sub)
    id_twonn = twonn_data.get(layer, np.nan)
    print(f"  Layer {layer:2d}: effective_rank={er:.3f}, TwoNN_ID={id_twonn:.3f}")
    h03_rows.append({'layer': layer, 'effective_rank': er, 'twonn_id': id_twonn})

h03_path = ITER_DIR / "h03_effective_rank_per_layer.csv"
with open(h03_path, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['layer', 'effective_rank', 'twonn_id'])
    w.writeheader()
    w.writerows(h03_rows)

er_vals = np.array([r['effective_rank'] for r in h03_rows])
id_vals = np.array([r['twonn_id'] for r in h03_rows])
pr, pp = pearsonr(er_vals, id_vals)
sr, sp = spearmanr(er_vals, id_vals)
print(f"  Pearson r={pr:.3f} (p={pp:.4f}), Spearman r={sr:.3f} (p={sp:.4f})")
h03_summary = {
    'er_range': [float(er_vals.min()), float(er_vals.max())],
    'er_layer0': float(er_vals[0]),
    'er_layer11': float(er_vals[-1]),
    'er_monotone_decrease': bool(all(er_vals[i] >= er_vals[i+1] for i in range(len(er_vals)-1))),
    'pearson_r_vs_twonn': float(pr),
    'pearson_p': float(pp),
    'spearman_r_vs_twonn': float(sr),
    'spearman_p': float(sp),
}


# ---------------------------------------------------------------------------
# H02: Per-gene residual drift
# ---------------------------------------------------------------------------
print("\n=== H02: Per-gene residual drift ===", flush=True)

layer0_sub = emb[0][gene_emb_idxs]    # (209, 512)
layer11_sub = emb[11][gene_emb_idxs]   # (209, 512)
drift = np.linalg.norm(layer11_sub - layer0_sub, axis=1)  # (209,)
drift_sorted_idx = np.argsort(drift)[::-1]

print(f"  Drift range: {drift.min():.3f}-{drift.max():.3f}, mean={drift.mean():.3f}, std={drift.std():.3f}")

h02_gene_drift_path = ITER_DIR / "h02_gene_residual_drift.csv"
with open(h02_gene_drift_path, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['gene', 'emb_idx', 'drift_l2', 'drift_rank'])
    for rank, local_idx in enumerate(drift_sorted_idx):
        w.writerow([gene_list[local_idx], gene_emb_idxs[local_idx], float(drift[local_idx]), rank])

top_n = min(50, N_GENES // 4)
bot_n = top_n
top_genes = set(gene_list[i] for i in drift_sorted_idx[:top_n])
bot_genes = set(gene_list[i] for i in drift_sorted_idx[-bot_n:])
print(f"  Top-{top_n} drift genes: {sorted(top_genes)[:10]}...")
print(f"  Bot-{bot_n} drift genes: {sorted(bot_genes)[:10]}...")


# ---------------------------------------------------------------------------
# Fetch GO annotations
# ---------------------------------------------------------------------------
def fetch_go_annotations(gene_symbols, species='human'):
    """Fetch GO BP annotations using the mygene Python client."""
    import mygene
    mg = mygene.MyGeneInfo()
    go_dict = defaultdict(set)
    all_syms = list(gene_symbols)
    results = mg.querymany(all_syms, scopes='symbol', fields='go.BP', species=species, verbose=False)
    for hit in results:
        sym = hit.get('query', '')
        go_bp = hit.get('go', {}).get('BP', [])
        if isinstance(go_bp, dict):
            go_bp = [go_bp]
        for term in go_bp:
            if isinstance(term, dict) and 'id' in term:
                go_dict[sym].add(term['id'])
    return go_dict

print("\nFetching GO annotations...", flush=True)
try:
    go_annotations = fetch_go_annotations(gene_list, species='human')
    n_annotated = sum(1 for g in gene_list if go_annotations.get(g))
    print(f"  Annotated: {n_annotated}/{N_GENES}")
    go_fetch_ok = True
except Exception as e:
    print(f"  GO fetch failed: {e}")
    go_fetch_ok = False
    go_annotations = {}

# Build GO term -> gene indices (local indices into our 209-gene list)
go_to_local_idxs = defaultdict(list)
for local_idx, gene in enumerate(gene_list):
    for term in go_annotations.get(gene, set()):
        go_to_local_idxs[term].append(local_idx)

valid_terms = {t: idxs for t, idxs in go_to_local_idxs.items()
               if MIN_TERM_SIZE <= len(idxs) <= MAX_TERM_SIZE}
print(f"  Valid GO BP terms ({MIN_TERM_SIZE}-{MAX_TERM_SIZE} genes): {len(valid_terms)}")


# ---------------------------------------------------------------------------
# H02 GO enrichment: top vs bottom drift genes
# ---------------------------------------------------------------------------
h02_enrichment_rows = []
if go_fetch_ok and len(valid_terms) > 0:
    top_local_idxs = set(drift_sorted_idx[:top_n].tolist())
    bot_local_idxs = set(drift_sorted_idx[-bot_n:].tolist())

    for term, t_idxs in valid_terms.items():
        t_set = set(t_idxs)
        a = len(t_set & top_local_idxs)
        b = top_n - a
        c = len(t_set & bot_local_idxs)
        d = bot_n - c
        if a + c == 0:
            continue
        _, p = fisher_exact([[a, b], [c, d]], alternative='greater')
        odds = (a * d) / max((b * c), 1e-9)
        h02_enrichment_rows.append({
            'go_term': term,
            'n_term_genes': len(t_idxs),
            'n_in_top': a,
            'n_in_bot': c,
            'odds_ratio': float(odds),
            'p_fisher': float(p),
        })

    h02_enrichment_rows.sort(key=lambda x: x['p_fisher'])
    h02_enrich_path = ITER_DIR / "h02_go_enrichment_drift.csv"
    with open(h02_enrich_path, 'w', newline='') as f:
        fieldnames = ['go_term', 'n_term_genes', 'n_in_top', 'n_in_bot', 'odds_ratio', 'p_fisher']
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(h02_enrichment_rows)

    n_sig = sum(1 for r in h02_enrichment_rows if r['p_fisher'] < 0.05)
    n_total = len(h02_enrichment_rows)
    n_fdr = sum(1 for i, r in enumerate(h02_enrichment_rows)
                if r['p_fisher'] < 0.05 * (i + 1) / n_total) if n_total > 0 else 0
    print(f"  GO enrichment tests: {n_total}, sig(p<0.05): {n_sig}, FDR-adj ~{n_fdr}")
    if h02_enrichment_rows:
        t = h02_enrichment_rows[0]
        print(f"  Top hit: {t['go_term']} OR={t['odds_ratio']:.2f} p={t['p_fisher']:.4f} (top={t['n_in_top']}, bot={t['n_in_bot']})")
    h02_summary = {
        'n_go_terms_tested': n_total,
        'n_sig_p05': n_sig,
        'n_sig_fdr': n_fdr,
        'drift_mean': float(drift.mean()),
        'drift_std': float(drift.std()),
        'top_hit_term': h02_enrichment_rows[0]['go_term'] if h02_enrichment_rows else None,
        'top_hit_or': float(h02_enrichment_rows[0]['odds_ratio']) if h02_enrichment_rows else None,
        'top_hit_p': float(h02_enrichment_rows[0]['p_fisher']) if h02_enrichment_rows else None,
    }
else:
    h02_summary = {'error': 'no GO data or no valid terms', 'n_go_terms_tested': 0}


# ---------------------------------------------------------------------------
# H01: GO term clustering in embedding space across layers
# ---------------------------------------------------------------------------
print("\n=== H01: GO term embedding clustering ===", flush=True)

if go_fetch_ok and len(valid_terms) > 0:
    # Sample up to 80 terms for speed
    term_keys = list(valid_terms.keys())
    if len(term_keys) > 80:
        rng2 = np.random.default_rng(42)
        term_keys = list(rng2.choice(term_keys, size=80, replace=False))
    print(f"  Testing {len(term_keys)} GO terms across {N_LAYERS} layers, {N_NULL_REPS} null reps")

    h01_rows = []
    for layer in range(N_LAYERS):
        X_sub = emb[layer][gene_emb_idxs]  # (209, 512)
        # Normalize
        norms = np.linalg.norm(X_sub, axis=1, keepdims=True)
        norms = np.where(norms < 1e-10, 1.0, norms)
        X_norm = X_sub / norms

        # Observed: mean pairwise cosine dist per GO term
        obs_dists = []
        term_sizes = []
        for term in term_keys:
            t_idxs = valid_terms[term]
            if len(t_idxs) < MIN_TERM_SIZE:
                continue
            sub = X_norm[t_idxs]
            D = cosine_dist_matrix(sub)
            n = len(t_idxs)
            mask = np.triu(np.ones((n, n), dtype=bool), k=1)
            obs_dists.append(float(D[mask].mean()) if mask.sum() > 0 else np.nan)
            term_sizes.append(n)

        obs_mean = np.nanmean(obs_dists)

        # Null: random groups of same sizes
        rng3 = np.random.default_rng(999 + layer)
        null_means = []
        for _ in range(N_NULL_REPS):
            rep = []
            for sz in term_sizes:
                rand_idxs = rng3.choice(N_GENES, size=sz, replace=False)
                sub = X_norm[rand_idxs]
                D = cosine_dist_matrix(sub)
                mask = np.triu(np.ones((sz, sz), dtype=bool), k=1)
                rep.append(float(D[mask].mean()) if mask.sum() > 0 else np.nan)
            null_means.append(np.nanmean(rep))

        null_mean = np.mean(null_means)
        null_std = np.std(null_means)
        z = (obs_mean - null_mean) / (null_std + 1e-12)
        frac_below = sum(1 for d in obs_dists if not np.isnan(d) and d < null_mean) / max(1, len(obs_dists))

        print(f"  Layer {layer:2d}: obs={obs_mean:.4f}, null={null_mean:.4f}, z={z:.3f}, frac_below={frac_below:.3f} (n_terms={len(obs_dists)})")
        h01_rows.append({
            'layer': layer,
            'n_terms': len(obs_dists),
            'obs_mean_cosine_dist': float(obs_mean),
            'null_mean': float(null_mean),
            'null_std': float(null_std),
            'z_score': float(z),
            'frac_terms_below_null': float(frac_below),
        })

    h01_path = ITER_DIR / "h01_go_clustering_by_layer.csv"
    with open(h01_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['layer', 'n_terms', 'obs_mean_cosine_dist', 'null_mean', 'null_std', 'z_score', 'frac_terms_below_null'])
        w.writeheader()
        w.writerows(h01_rows)

    z_vals = [r['z_score'] for r in h01_rows]
    mean_z = float(np.mean(z_vals))
    n_sig = sum(1 for z in z_vals if z < -1.96)
    print(f"\n  Summary: mean_z={mean_z:.3f}, layers with z<-1.96: {n_sig}/12")
    h01_summary = {
        'n_terms': len(term_keys),
        'mean_z': mean_z,
        'n_sig_layers_z_lt_neg1p96': n_sig,
        'z_range': [float(min(z_vals)), float(max(z_vals))],
        'z_by_layer': z_vals,
    }
else:
    h01_summary = {'error': 'no GO data', 'mean_z': None}
    print("  Skipping (no GO data)")


# ---------------------------------------------------------------------------
# Save combined results JSON
# ---------------------------------------------------------------------------
results = {
    'iteration': 'iter_0005',
    'H01_go_clustering': h01_summary,
    'H02_residual_drift_enrichment': h02_summary,
    'H03_effective_rank': h03_summary,
}
results_path = ITER_DIR / "iter0005_results.json"
with open(results_path, 'w') as f:
    json.dump(results, f, indent=2)

print(f"\n=== Complete. Results in {ITER_DIR} ===")
print(json.dumps(results, indent=2))
