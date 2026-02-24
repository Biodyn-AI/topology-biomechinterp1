"""
iter_0040 Multi-Hypothesis Screen

H01 (manifold_distance / new_method): Full 12-layer GC attractor onset scan + IRF4 plasma-bridge
    For each layer L0..L11 in cycle4_immune:
    - rank of each GC-TF (BATF, BACH2, BCL6, PAX5) among genes nearest to B-cell centroid
    - rank of IRF4 (plasma-bridge candidate)
    - detect "attractor onset layer" = first layer where GC rank < 20 (top-20 proximity)

H02 (cross_model_alignment / new_family): Geneformer token embedding B-cell precision@10
    Use Geneformer input token embeddings (cycle12 data) to test if B-cell markers cluster
    in Geneformer's own embedding space. Compare to scGPT L0 result.

H03 (manifold_distance / refinement): GC-plasma subspace principal angle across 12 layers
    SVD on 4-gene GC cluster vs plasma cluster at each layer.
    Track principal angle (scipy.linalg.subspace_angles) across L0..L11.
    Test: does GC-plasma angle increase (diverge) with layer depth?
    CD19+BLK minimal anchor variant for precision@10 across all layers.
"""

import numpy as np
import json
import csv
from pathlib import Path
from scipy.stats import spearmanr
from scipy.linalg import subspace_angles
import anndata as ad

# ─── Paths ────────────────────────────────────────────────────────────────────
ROOT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work")
CYCLE1 = ROOT / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main"
CYCLE4 = ROOT / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle4_immune_main"
IMMUNE_H5AD = ROOT / "single_cell_mechinterp/outputs/tabula_sapiens_immune_subset_hpn_processed.h5ad"
GF_DIR = ROOT / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle12_geneformer_immune_bootstrap"
ITER_DIR = ROOT / "subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0040"
ITER_DIR.mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng(42)

# ─── Load cycle1 embeddings ───────────────────────────────────────────────────
print("Loading cycle1 embeddings...", flush=True)
emb1 = np.load(CYCLE1 / "layer_gene_embeddings.npy")   # [12, 4803, 512]
with open(CYCLE1 / "gene_list.txt") as f:
    vocab1 = [line.strip() for line in f if line.strip()]
gene_to_idx1 = {g: i for i, g in enumerate(vocab1)}
print(f"  cycle1 shape: {emb1.shape}, genes: {len(vocab1)}", flush=True)

# ─── Load cycle4_immune embeddings ───────────────────────────────────────────
print("Loading cycle4_immune embeddings...", flush=True)
emb4 = np.load(CYCLE4 / "layer_gene_embeddings.npy")   # [12, 4941, 512]
adata4 = ad.read_h5ad(IMMUNE_H5AD)
gene_list4 = list(adata4.var_names)
gene_to_idx4 = {g: i for i, g in enumerate(gene_list4)}
print(f"  cycle4 shape: {emb4.shape}, genes: {len(gene_list4)}", flush=True)

# Filter out zero-norm genes from cycle4
norms0_4 = np.linalg.norm(emb4[0], axis=1)  # [N_GENES]
valid_mask4 = norms0_4 > 1e-8
valid_genes4 = [g for g, v in zip(gene_list4, valid_mask4) if v]
valid_idx4 = np.array([gene_to_idx4[g] for g in valid_genes4])
gene_to_valid4 = {g: i for i, g in enumerate(valid_genes4)}
emb4_valid = emb4[:, valid_idx4, :]  # [12, N_valid, 512]
print(f"  cycle4 valid (non-zero): {len(valid_genes4)}", flush=True)

# ─── Gene panels ────────────────────────────────────────────────────────────
GC_TFS = ['BATF', 'BACH2', 'BCL6', 'PAX5']
PLASMA_BRIDGE = ['IRF4']
PLASMA_MARKERS = ['JCHAIN', 'SDC1', 'PRDM1']
BCELL_ANCHOR = ['MS4A1', 'CD79A', 'BLK']        # plasma-exclusive B-cell markers
BCELL_FULL = ['MS4A1', 'CD19', 'CD79A', 'BLK', 'PRDM1']
BCELL_MINIMAL = ['CD19', 'BLK']                  # minimal anchor from iter_0039

# ─── Utility functions ────────────────────────────────────────────────────────
def get_centroid(emb_layer, genes, gene_map):
    """Compute mean embedding vector for a gene set."""
    found = [g for g in genes if g in gene_map]
    if not found:
        return None, []
    idxs = [gene_map[g] for g in found]
    return emb_layer[idxs].mean(axis=0), found

def rank_in_neighbors(emb_layer, query_vec, target_genes, gene_map):
    """Return rank (1-based) of each target gene among all genes by L2 dist to query_vec."""
    dists = np.linalg.norm(emb_layer - query_vec, axis=1)  # [N_genes]
    sorted_idx = np.argsort(dists)
    rank_of = {}
    for gene in target_genes:
        if gene not in gene_map:
            rank_of[gene] = None
        else:
            gidx = gene_map[gene]
            rank_of[gene] = int(np.where(sorted_idx == gidx)[0][0]) + 1  # 1-based
    return rank_of

def precision_at_k(emb_layer, anchor_genes, target_genes, gene_map, k=10, exclude_anchor=True):
    """Precision@k: fraction of k-NN of centroid that are target genes."""
    centroid_vec, found_anchor = get_centroid(emb_layer, anchor_genes, gene_map)
    if centroid_vec is None:
        return None, 0
    dists = np.linalg.norm(emb_layer - centroid_vec, axis=1)
    sorted_genes = [list(gene_map.keys())[i] for i in np.argsort(dists)]
    if exclude_anchor:
        sorted_genes = [g for g in sorted_genes if g not in set(anchor_genes)]
    target_set = set(target_genes) & set(gene_map.keys())
    hits = sum(1 for g in sorted_genes[:k] if g in target_set)
    return hits / k, len(target_set)


# ═══════════════════════════════════════════════════════════════════════════════
# H01: Full 12-layer GC attractor onset scan + IRF4 plasma-bridge (cycle4)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== H01: Full 12-layer GC attractor onset scan ===", flush=True)

h01_layers = []
TARGET_GENES_H01 = GC_TFS + PLASMA_BRIDGE + PLASMA_MARKERS

# Check availability
for g in TARGET_GENES_H01 + BCELL_ANCHOR:
    print(f"  {g}: {'in-vocab' if g in gene_to_valid4 else 'OOV'}", flush=True)

bcell_anchor_available = [g for g in BCELL_ANCHOR if g in gene_to_valid4]
print(f"\nUsing B-cell anchor (cycle4): {bcell_anchor_available}", flush=True)

for layer in range(12):
    E = emb4_valid[layer]  # [N_valid, 512]

    # Compute B-cell centroid
    centroid_vec, found_anchor = get_centroid(E, BCELL_ANCHOR, gene_to_valid4)
    if centroid_vec is None:
        print(f"  L{layer}: anchor not found", flush=True)
        continue

    # Rank all target genes
    ranks = rank_in_neighbors(E, centroid_vec, TARGET_GENES_H01, gene_to_valid4)

    # Compute centroid distance from L0
    if layer == 0:
        centroid_L0 = centroid_vec.copy()
    drift_from_L0 = float(np.linalg.norm(centroid_vec - centroid_L0))

    layer_result = {
        'layer': layer,
        'anchor_genes': bcell_anchor_available,
        'ranks': ranks,
        'centroid_drift_from_L0': drift_from_L0,
    }
    h01_layers.append(layer_result)

    gc_ranks_str = ', '.join([f"{g}:{ranks.get(g,'OOV')}" for g in GC_TFS])
    irf4_rank = ranks.get('IRF4', 'OOV')
    print(f"  L{layer:2d}: GC ranks [{gc_ranks_str}], IRF4={irf4_rank}, drift={drift_from_L0:.2f}", flush=True)

# Find attractor onset: first layer where any GC-TF rank <= 20
attractor_onset = None
for entry in h01_layers:
    gc_ranks = [entry['ranks'].get(g) for g in GC_TFS if entry['ranks'].get(g) is not None]
    if gc_ranks and min(gc_ranks) <= 20:
        attractor_onset = entry['layer']
        break

print(f"\nH01 summary: GC attractor onset layer = {attractor_onset}", flush=True)

# Spearman rho for GC-TF mean rank across layers (should decrease)
layers_list = [e['layer'] for e in h01_layers]
mean_gc_ranks = []
for entry in h01_layers:
    gc_r = [entry['ranks'].get(g) for g in GC_TFS if entry['ranks'].get(g) is not None]
    mean_gc_ranks.append(np.mean(gc_r) if gc_r else np.nan)

irf4_ranks_by_layer = [entry['ranks'].get('IRF4') for entry in h01_layers]
irf4_valid = [(l, r) for l, r in zip(layers_list, irf4_ranks_by_layer) if r is not None]

if len([r for r in mean_gc_ranks if not np.isnan(r)]) >= 3:
    rho_gc, p_gc = spearmanr(layers_list, mean_gc_ranks)
    print(f"H01: GC mean rank Spearman rho vs layer = {rho_gc:.3f} (p={p_gc:.4f})", flush=True)
else:
    rho_gc, p_gc = None, None

if len(irf4_valid) >= 3:
    l_irf4, r_irf4 = zip(*irf4_valid)
    rho_irf4, p_irf4 = spearmanr(l_irf4, r_irf4)
    print(f"H01: IRF4 rank Spearman rho vs layer = {rho_irf4:.3f} (p={p_irf4:.4f})", flush=True)
else:
    rho_irf4, p_irf4 = None, None

# Save H01
h01_result = {
    'layers': h01_layers,
    'attractor_onset_layer': attractor_onset,
    'gc_mean_rank_spearman_rho': rho_gc,
    'gc_mean_rank_spearman_p': p_gc,
    'irf4_rank_spearman_rho': rho_irf4,
    'irf4_rank_spearman_p': p_irf4,
    'gc_mean_ranks_by_layer': mean_gc_ranks,
    'irf4_ranks_by_layer': irf4_ranks_by_layer,
}
with open(ITER_DIR / "h01_gc_attractor_scan.json", "w") as f:
    json.dump(h01_result, f, indent=2, default=lambda x: float(x) if isinstance(x, (np.floating, np.integer)) else x)
print("H01 saved.", flush=True)


# ═══════════════════════════════════════════════════════════════════════════════
# H02: Geneformer token embedding B-cell precision@10
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== H02: Geneformer cross-model B-cell precision@10 ===", flush=True)

try:
    # Load Geneformer gene token map
    gf_gene_token = {}
    with open(GF_DIR / "geneformer_gene_token_map.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            gf_gene_token[row['gene']] = int(row['token_id'])
    print(f"  Geneformer token map: {len(gf_gene_token)} genes", flush=True)

    # Load the Geneformer model token embeddings
    # Try to load from Hugging Face cache or model path
    gf_model_paths = [
        ROOT / "single_cell_mechinterp/external/Geneformer",
        Path.home() / ".cache/huggingface/hub/models--ctheodoris--Geneformer",
        Path("/tmp/geneformer_model"),
    ]

    gf_emb_matrix = None
    for mpath in gf_model_paths:
        if mpath.exists():
            print(f"  Checking Geneformer path: {mpath}", flush=True)
            # Look for pytorch_model.bin or model.safetensors
            for fname in ['pytorch_model.bin', 'model.safetensors']:
                candidates = list(mpath.rglob(fname))
                if candidates:
                    print(f"  Found model file: {candidates[0]}", flush=True)
                    import torch
                    state = torch.load(candidates[0], map_location='cpu', weights_only=False)
                    # Geneformer uses BertForMaskedLM; token embeddings at bert.embeddings.word_embeddings.weight
                    for key in state:
                        if 'word_embeddings.weight' in key or 'token_type_embeddings' in key:
                            print(f"  Key: {key} shape: {state[key].shape}", flush=True)
                    if 'bert.embeddings.word_embeddings.weight' in state:
                        gf_emb_matrix = state['bert.embeddings.word_embeddings.weight'].numpy()
                        print(f"  Loaded GF word embeddings: {gf_emb_matrix.shape}", flush=True)
                    break
            if gf_emb_matrix is not None:
                break

    if gf_emb_matrix is None:
        print("  Geneformer model weights not found locally. Attempting huggingface download...", flush=True)
        try:
            from transformers import AutoModel
            import torch
            model = AutoModel.from_pretrained("ctheodoris/Geneformer", trust_remote_code=True)
            gf_emb_matrix = model.embeddings.word_embeddings.weight.detach().numpy()
            print(f"  Downloaded GF word embeddings: {gf_emb_matrix.shape}", flush=True)
        except Exception as e2:
            print(f"  Could not load Geneformer: {e2}", flush=True)
            gf_emb_matrix = None

    if gf_emb_matrix is not None:
        # Build gene -> embedding mapping
        gf_gene_list = sorted(gf_gene_token.keys())
        gf_valid_genes = [g for g in gf_gene_list if gf_gene_token.get(g, -1) < gf_emb_matrix.shape[0]]
        gf_emb_by_gene = {g: gf_emb_matrix[gf_gene_token[g]] for g in gf_valid_genes}

        # Create embedding matrix for precision@k computation
        gf_emb_arr = np.stack([gf_emb_by_gene[g] for g in gf_valid_genes], axis=0)  # [N, D]
        gf_gene_to_local = {g: i for i, g in enumerate(gf_valid_genes)}

        print(f"  GF embedding matrix: {gf_emb_arr.shape}", flush=True)

        # Check availability
        bcell_gf = [g for g in BCELL_FULL if g in gf_gene_to_local]
        gc_gf = [g for g in GC_TFS if g in gf_gene_to_local]
        print(f"  B-cell markers in GF: {bcell_gf}", flush=True)
        print(f"  GC-TFs in GF: {gc_gf}", flush=True)

        # Precision@10 for B-cell centroid → GC-TFs
        p10_gf, n_target = precision_at_k(gf_emb_arr, BCELL_FULL, GC_TFS, gf_gene_to_local, k=10)

        # Null: 200 random anchor sets of same size
        n_anchor = len(bcell_gf)
        null_p10 = []
        for _ in range(200):
            rand_anchor = rng.choice(gf_valid_genes, size=n_anchor, replace=False).tolist()
            p, _ = precision_at_k(gf_emb_arr, rand_anchor, GC_TFS, gf_gene_to_local, k=10)
            if p is not None:
                null_p10.append(p)
        null_p10 = np.array(null_p10)
        z_p10 = (p10_gf - null_p10.mean()) / (null_p10.std() + 1e-10) if p10_gf is not None else None
        pctile = float(np.mean(null_p10 < p10_gf)) if p10_gf is not None else None

        print(f"  GF B-cell precision@10 = {p10_gf:.3f}, null mean={null_p10.mean():.3f}, z={z_p10:.2f}, pctile={pctile:.3f}", flush=True)

        h02_result = {
            'status': 'tested',
            'gf_emb_shape': list(gf_emb_matrix.shape),
            'bcell_in_gf': bcell_gf,
            'gc_tfs_in_gf': gc_gf,
            'precision_at_10': float(p10_gf) if p10_gf is not None else None,
            'null_mean': float(null_p10.mean()),
            'null_std': float(null_p10.std()),
            'z_score': float(z_p10) if z_p10 is not None else None,
            'null_percentile': float(pctile) if pctile is not None else None,
        }
    else:
        h02_result = {'status': 'blocked', 'reason': 'Geneformer model weights not accessible'}

except Exception as e:
    print(f"  H02 error: {e}", flush=True)
    h02_result = {'status': 'blocked', 'reason': str(e)}

with open(ITER_DIR / "h02_geneformer_bcell.json", "w") as f:
    json.dump(h02_result, f, indent=2)
print("H02 saved.", flush=True)


# ═══════════════════════════════════════════════════════════════════════════════
# H03: GC-plasma subspace principal angle + CD19+BLK minimal anchor precision@10
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== H03: GC-plasma subspace principal angles + minimal anchor scan ===", flush=True)

# Part A: Subspace principal angles per layer (cycle4)
# GC subspace: BATF, BACH2, BCL6, PAX5
# Plasma subspace: JCHAIN, SDC1 + whatever plasma markers are available
# B-cell subspace: MS4A1, CD79A, BLK (anchor)

GC_SET = ['BATF', 'BACH2', 'BCL6', 'PAX5']
PLASMA_SET = ['JCHAIN', 'SDC1', 'PRDM1', 'IRF4']
BCELL_SET = ['MS4A1', 'CD79A', 'BLK']

gc_available = [g for g in GC_SET if g in gene_to_valid4]
plasma_available = [g for g in PLASMA_SET if g in gene_to_valid4]
bcell_available = [g for g in BCELL_SET if g in gene_to_valid4]

print(f"  GC genes available: {gc_available}", flush=True)
print(f"  Plasma genes available: {plasma_available}", flush=True)
print(f"  B-cell genes available: {bcell_available}", flush=True)

h03_layers = []
h03_precision = []

for layer in range(12):
    E = emb4_valid[layer]  # [N_valid, 512]

    # --- Part A: Subspace principal angles ---
    def get_gene_matrix(genes, gene_map, emb):
        found = [g for g in genes if g in gene_map]
        if len(found) < 2:
            return None, found
        idxs = [gene_map[g] for g in found]
        return emb[idxs], found  # [n_genes, D]

    gc_mat, gc_found = get_gene_matrix(gc_available, gene_to_valid4, E)
    plasma_mat, plasma_found = get_gene_matrix(plasma_available, gene_to_valid4, E)
    bcell_mat, bcell_found = get_gene_matrix(bcell_available, gene_to_valid4, E)

    angles_gc_plasma = None
    angles_gc_bcell = None

    if gc_mat is not None and plasma_mat is not None:
        # Use SVD-based subspace representation
        # Subspace spanned by gene embeddings: take top min(n,d) singular vectors
        try:
            angles_rad = subspace_angles(gc_mat.T, plasma_mat.T)
            angles_gc_plasma = float(np.degrees(angles_rad[0]))  # smallest principal angle
        except Exception as e:
            print(f"  L{layer} subspace_angles error: {e}", flush=True)

    if gc_mat is not None and bcell_mat is not None:
        try:
            angles_rad2 = subspace_angles(gc_mat.T, bcell_mat.T)
            angles_gc_bcell = float(np.degrees(angles_rad2[0]))
        except Exception as e:
            print(f"  L{layer} subspace_angles2 error: {e}", flush=True)

    # --- Part B: Precision@10 for CD19+BLK minimal anchor -> GC-TFs ---
    bcell_minimal_avail = [g for g in BCELL_MINIMAL if g in gene_to_valid4]
    p10_minimal = None
    if bcell_minimal_avail:
        gene_list_local = list(gene_to_valid4.keys())
        gene_local_map = {g: i for i, g in enumerate(gene_list_local)}
        E_full_check = emb4_valid[layer]
        p10_minimal, _ = precision_at_k(E_full_check, bcell_minimal_avail, GC_TFS, gene_to_valid4, k=10)

    # --- Part C: Full B-cell anchor precision@10 for GC-TFs ---
    p10_full = None
    if bcell_available:
        p10_full, _ = precision_at_k(E, bcell_available, GC_TFS, gene_to_valid4, k=10)

    layer_result = {
        'layer': layer,
        'gc_plasma_angle_deg': angles_gc_plasma,
        'gc_bcell_angle_deg': angles_gc_bcell,
        'precision_at_10_minimal_anchor': p10_minimal,
        'precision_at_10_full_anchor': p10_full,
        'gc_found': gc_found,
        'plasma_found': plasma_found,
    }
    h03_layers.append(layer_result)
    print(f"  L{layer:2d}: GC-plasma angle={angles_gc_plasma:.1f}°, GC-bcell angle={angles_gc_bcell:.1f}°, p@10(minimal)={p10_minimal}, p@10(full)={p10_full}", flush=True)

# Spearman trend for GC-plasma angle across layers
angles_gcp = [e['gc_plasma_angle_deg'] for e in h03_layers]
angles_gcb = [e['gc_bcell_angle_deg'] for e in h03_layers]
p10_minimal_all = [e['precision_at_10_minimal_anchor'] for e in h03_layers]
p10_full_all = [e['precision_at_10_full_anchor'] for e in h03_layers]
layers_h03 = [e['layer'] for e in h03_layers]

rho_angle_gcp, p_angle_gcp = spearmanr(layers_h03, [a for a in angles_gcp if a is not None])
rho_p10_minimal, p_p10_minimal = spearmanr(layers_h03, [p for p in p10_minimal_all if p is not None])
rho_p10_full, p_p10_full = spearmanr(layers_h03, [p for p in p10_full_all if p is not None])

print(f"\nH03 summary:", flush=True)
print(f"  GC-plasma angle Spearman rho vs layer: {rho_angle_gcp:.3f} (p={p_angle_gcp:.4f})", flush=True)
print(f"  p@10 minimal anchor Spearman rho vs layer: {rho_p10_minimal:.3f} (p={p_p10_minimal:.4f})", flush=True)
print(f"  p@10 full anchor Spearman rho vs layer: {rho_p10_full:.3f} (p={p_p10_full:.4f})", flush=True)

# Null comparison for precision@10 at L11 with minimal anchor
E_L11 = emb4_valid[11]
n_minimal = len([g for g in BCELL_MINIMAL if g in gene_to_valid4])
null_p10_minimal = []
if n_minimal > 0:
    for _ in range(200):
        rand_anchor = rng.choice(list(gene_to_valid4.keys()), size=n_minimal, replace=False).tolist()
        p, _ = precision_at_k(E_L11, rand_anchor, GC_TFS, gene_to_valid4, k=10)
        if p is not None:
            null_p10_minimal.append(p)
    null_p10_minimal = np.array(null_p10_minimal)
    p10_minimal_L11 = p10_minimal_all[11] if p10_minimal_all[11] is not None else 0.0
    z_minimal_L11 = (p10_minimal_L11 - null_p10_minimal.mean()) / (null_p10_minimal.std() + 1e-10)
    pctile_minimal_L11 = float(np.mean(null_p10_minimal < p10_minimal_L11))
    print(f"  L11 minimal anchor p@10={p10_minimal_L11:.3f}, null={null_p10_minimal.mean():.3f}±{null_p10_minimal.std():.3f}, z={z_minimal_L11:.2f}, pctile={pctile_minimal_L11:.3f}", flush=True)
else:
    null_p10_minimal = np.array([])
    z_minimal_L11 = None
    pctile_minimal_L11 = None

h03_result = {
    'layers': h03_layers,
    'gc_plasma_angle_spearman_rho': float(rho_angle_gcp),
    'gc_plasma_angle_spearman_p': float(p_angle_gcp),
    'p10_minimal_anchor_spearman_rho': float(rho_p10_minimal),
    'p10_minimal_anchor_spearman_p': float(p_p10_minimal),
    'p10_full_anchor_spearman_rho': float(rho_p10_full),
    'p10_full_anchor_spearman_p': float(p_p10_full),
    'null_p10_minimal_L11_mean': float(null_p10_minimal.mean()) if len(null_p10_minimal) > 0 else None,
    'null_p10_minimal_L11_std': float(null_p10_minimal.std()) if len(null_p10_minimal) > 0 else None,
    'p10_minimal_L11_z': float(z_minimal_L11) if z_minimal_L11 is not None else None,
    'p10_minimal_L11_pctile': float(pctile_minimal_L11) if pctile_minimal_L11 is not None else None,
}
with open(ITER_DIR / "h03_subspace_angles.json", "w") as f:
    json.dump(h03_result, f, indent=2, default=lambda x: float(x) if isinstance(x, (np.floating, np.integer)) else x)
print("H03 saved.", flush=True)


# ═══════════════════════════════════════════════════════════════════════════════
# Compile summary CSV
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== Writing summary CSV ===", flush=True)

summary_rows = []
# H01 summary
for entry in h01_layers:
    gc_r = {g: entry['ranks'].get(g) for g in GC_TFS}
    irf4_r = entry['ranks'].get('IRF4')
    summary_rows.append({
        'hypothesis': 'H01',
        'layer': entry['layer'],
        'metric': 'gc_mean_rank',
        'value': np.mean([v for v in gc_r.values() if v is not None]),
        'irf4_rank': irf4_r,
        'drift_from_L0': entry['centroid_drift_from_L0'],
        **{f'rank_{g}': v for g, v in gc_r.items()},
    })

# H03 summary
for entry in h03_layers:
    summary_rows.append({
        'hypothesis': 'H03',
        'layer': entry['layer'],
        'metric': 'gc_plasma_angle_deg',
        'value': entry['gc_plasma_angle_deg'],
        'gc_bcell_angle_deg': entry['gc_bcell_angle_deg'],
        'p10_minimal': entry['precision_at_10_minimal_anchor'],
        'p10_full': entry['precision_at_10_full_anchor'],
    })

import pandas as pd
df = pd.DataFrame(summary_rows)
df.to_csv(ITER_DIR / "iter_0040_summary.csv", index=False)
print(f"Summary CSV written: {len(df)} rows", flush=True)

print("\n=== All iter_0040 experiments complete ===", flush=True)
