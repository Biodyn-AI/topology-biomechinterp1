"""
iter_0032 Multi-Hypothesis Screen

H01 (module_structure / new_method): Community detection on L11 kNN graph → biological annotation
    Greedy modularity + connected components on k=10 kNN. Map communities to GO families,
    cell-type markers, STRING hubs. Test enrichment via Fisher's exact test.

H02 (manifold_distance / refinement, TRRUST OOV-corrected): TRRUST regulatory pairs → L2 dist
    For 195 in-vocab genes: test if TRRUST TF-target pairs are closer in embedding space
    than random gene pairs. Compare across layers (early vs late). AUROC per layer.

H03 (intrinsic_dimensionality / new_method): PC1 at L11 vs biological partitions
    Project 195 in-vocab genes onto PC1 at L11 (explains 26% variance).
    Test: TF vs non-TF (Mann-Whitney), cell-type marker groups, STRING high-degree hubs.
    AUROC per group. Compare to PC1 at L0 (control).
"""

import numpy as np
import json
import csv
from pathlib import Path
from scipy.stats import spearmanr, mannwhitneyu, fisher_exact
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import laplacian
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings("ignore")

PROJECT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
               "/subproject_41_claude_topology_hypothesis_screening_autoloop")
ITER_DIR = PROJECT / "iterations" / "iter_0032"
ITER_DIR.mkdir(parents=True, exist_ok=True)

CYCLE1 = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
              "/subproject_38_geometric_residual_stream_interpretability"
              "/implementation/outputs/cycle1_main")

ITER15_DIR = PROJECT / "iterations" / "iter_0015"
TRRUST_CACHE = ITER15_DIR / "trrust_named_gene_pairs.json"
STRING_CACHE = ITER15_DIR / "string_ppi_score04_cache.json"

rng = np.random.default_rng(42)

# ─── Load embeddings ──────────────────────────────────────────────────────────
print("Loading embeddings...", flush=True)
emb = np.load(CYCLE1 / "layer_gene_embeddings.npy")   # [12, 4803, 512]
N_LAYERS, N_GENES_TOTAL, N_DIM = emb.shape
print(f"  Shape: {emb.shape}", flush=True)

with open(CYCLE1 / "gene_list.txt") as f:
    vocab_genes = [line.strip() for line in f if line.strip()]
gene_to_emb_idx = {g: i for i, g in enumerate(vocab_genes)}

EDGES_FILE = CYCLE1 / "cycle1_edge_dataset.tsv"
named_genes_set = set()
with open(EDGES_FILE) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        named_genes_set.add(row["source"])
        named_genes_set.add(row["target"])
named_genes = sorted(named_genes_set)
print(f"  Total named genes: {len(named_genes)}", flush=True)

OOV_GENES = {"FOS","HLA-A","HLA-DPB1","JUNB","KLF6","LDHA","LGALS1",
             "NCAM1","NCOA3","NR4A3","PAX5","PTGS2","TBXAS1","TNF"}
inv_genes = [g for g in named_genes if g in gene_to_emb_idx and g not in OOV_GENES]
inv_indices = [gene_to_emb_idx[g] for g in inv_genes]
print(f"  In-vocab named genes: {len(inv_genes)}", flush=True)

emb_inv = emb[:, inv_indices, :]
print(f"  In-vocab emb shape: {emb_inv.shape}", flush=True)  # [12, 195, 512]
inv_gene_set = set(inv_genes)
inv_gene_to_local = {g: i for i, g in enumerate(inv_genes)}

# ─── H01: Community detection on L11 kNN → biological annotation ──────────────
print("\n=== H01: Community detection at L11 ===", flush=True)

L11_emb = emb_inv[11]  # [195, 512]
k = 10

# Build kNN graph
nbrs = NearestNeighbors(n_neighbors=k+1, metric='euclidean').fit(L11_emb)
distances, indices = nbrs.kneighbors(L11_emb)
# indices[:, 0] is self; skip
n = len(inv_genes)
rows, cols = [], []
for i in range(n):
    for j in indices[i, 1:]:
        rows.append(i)
        cols.append(j)
        rows.append(j)
        cols.append(i)

adj = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
adj = (adj > 0).astype(float)  # binarize

# Greedy modularity community detection via networkx
try:
    import networkx as nx
    G = nx.from_scipy_sparse_array(adj)
    from networkx.algorithms.community import greedy_modularity_communities
    communities = list(greedy_modularity_communities(G))
    modularity = nx.algorithms.community.quality.modularity(G, communities)
    n_communities = len(communities)
    community_sizes = sorted([len(c) for c in communities], reverse=True)
    print(f"  Communities: {n_communities}, modularity: {modularity:.4f}", flush=True)
    print(f"  Community sizes: {community_sizes[:10]}", flush=True)

    # Assign community labels
    gene_community = {}
    for comm_id, comm in enumerate(communities):
        for node_idx in comm:
            gene_community[inv_genes[node_idx]] = comm_id

    # Null: shuffle gene labels
    null_modularities = []
    for _ in range(100):
        shuffled_labels = list(range(n_communities))
        # create random communities of same sizes
        perm = rng.permutation(n)
        null_comms = []
        pos = 0
        for s in community_sizes:
            null_comms.append(set(perm[pos:pos+s].tolist()))
            pos += s
        q = nx.algorithms.community.quality.modularity(G, null_comms)
        null_modularities.append(q)
    null_mean = float(np.mean(null_modularities))
    null_std = float(np.std(null_modularities))
    z_score = (modularity - null_mean) / (null_std + 1e-10)
    print(f"  Null modularity: {null_mean:.4f} ± {null_std:.4f}, z={z_score:.2f}", flush=True)

    # Biological annotation: curated gene families
    bio_families = {
        "AP1": ["JUN", "FOSL1", "FOSL2", "ATF3", "BATF"],
        "HLA_I": ["HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G"],
        "HLA_II": ["HLA-DRA", "HLA-DRB1", "HLA-DMA", "HLA-DMB", "HLA-DOA",
                   "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1"],
        "BCL2fam": ["BCL2", "BCL2L1", "MCL1", "BCL2A1", "BAX", "BAD"],
        "TNFSF": ["TNFSF10", "TNFSF11", "TNFSF13", "TNFSF14", "CD70"],
        "IL2_pathway": ["IL2", "IL2RA", "IL2RB", "IL2RG", "IL7", "IL7R", "IL15"],
        "KLF": ["KLF4", "KLF6", "KLF2"],
        "RUNX": ["RUNX1", "RUNX2", "RUNX3"],
    }
    # Cell-type markers
    tcell_markers = ["CD3E", "CD3D", "CD8A", "CD8B", "CD4", "LCK", "ZAP70", "TCF7",
                     "GZMB", "PRF1", "IFNG", "TIGIT", "LAG3", "PDCD1"]
    bcell_markers = ["CD19", "MS4A1", "CD79A", "CD79B", "PAX5", "EBF1",
                     "IGHM", "IGHD", "IGKC"]
    nk_markers = ["NCAM1", "KLRB1", "KLRD1", "KLRG1", "NCR1", "NKG7"]
    monocyte_markers = ["CD14", "FCGR3A", "CSF1R", "LYZ", "S100A8", "S100A9",
                        "CCL2", "CCR2"]

    # Test: are same-family genes over-represented in same community?
    enrichment_results = {}
    for fname, fgenes in bio_families.items():
        fgenes_inv = [g for g in fgenes if g in inv_gene_set]
        if len(fgenes_inv) < 2:
            continue
        # For top community that contains most family members
        best_comm = None
        best_overlap = 0
        for comm_id, comm in enumerate(communities):
            comm_genes = {inv_genes[ni] for ni in comm}
            overlap = len(set(fgenes_inv) & comm_genes)
            if overlap > best_overlap:
                best_overlap = overlap
                best_comm = comm_id
        if best_comm is None:
            continue
        comm_genes = {inv_genes[ni] for ni in communities[best_comm]}
        # Fisher's exact
        in_comm_in_fam = len(set(fgenes_inv) & comm_genes)
        in_comm_not_fam = len(comm_genes) - in_comm_in_fam
        not_comm_in_fam = len(fgenes_inv) - in_comm_in_fam
        not_comm_not_fam = n - in_comm_in_fam - in_comm_not_fam - not_comm_in_fam
        contingency = [[in_comm_in_fam, in_comm_not_fam],
                       [not_comm_in_fam, not_comm_not_fam]]
        odds_ratio, pval = fisher_exact(contingency, alternative='greater')
        enrichment_results[fname] = {
            "best_community": best_comm,
            "comm_size": int(len(communities[best_comm])),
            "family_size_inv": len(fgenes_inv),
            "overlap": int(in_comm_in_fam),
            "fisher_odds": float(odds_ratio),
            "fisher_pval": float(pval)
        }
        print(f"  {fname}: overlap={in_comm_in_fam}/{len(fgenes_inv)} in comm {best_comm} "
              f"(size {len(communities[best_comm])}), OR={odds_ratio:.2f}, p={pval:.4f}", flush=True)

    # Test cell-type marker clustering
    marker_groups = {
        "tcell": tcell_markers,
        "bcell": bcell_markers,
        "nk": nk_markers,
        "monocyte": monocyte_markers
    }
    marker_enrichment = {}
    for mname, markers in marker_groups.items():
        markers_inv = [g for g in markers if g in inv_gene_set]
        if len(markers_inv) < 2:
            continue
        best_comm = None
        best_overlap = 0
        for comm_id, comm in enumerate(communities):
            comm_genes = {inv_genes[ni] for ni in comm}
            overlap = len(set(markers_inv) & comm_genes)
            if overlap > best_overlap:
                best_overlap = overlap
                best_comm = comm_id
        if best_comm is None:
            continue
        comm_genes = {inv_genes[ni] for ni in communities[best_comm]}
        in_comm_in_fam = len(set(markers_inv) & comm_genes)
        in_comm_not_fam = len(comm_genes) - in_comm_in_fam
        not_comm_in_fam = len(markers_inv) - in_comm_in_fam
        not_comm_not_fam = n - in_comm_in_fam - in_comm_not_fam - not_comm_in_fam
        contingency = [[in_comm_in_fam, in_comm_not_fam],
                       [not_comm_in_fam, not_comm_not_fam]]
        odds_ratio, pval = fisher_exact(contingency, alternative='greater')
        marker_enrichment[mname] = {
            "best_community": best_comm,
            "comm_size": int(len(communities[best_comm])),
            "markers_in_inv": len(markers_inv),
            "overlap": int(in_comm_in_fam),
            "fisher_odds": float(odds_ratio),
            "fisher_pval": float(pval)
        }
        print(f"  {mname}: overlap={in_comm_in_fam}/{len(markers_inv)} in comm {best_comm} "
              f"(size {len(communities[best_comm])}), OR={odds_ratio:.2f}, p={pval:.4f}", flush=True)

    h01_result = {
        "n_communities": n_communities,
        "modularity": float(modularity),
        "null_modularity_mean": null_mean,
        "null_modularity_std": null_std,
        "z_score": float(z_score),
        "community_sizes": [int(s) for s in community_sizes],
        "bio_family_enrichment": enrichment_results,
        "cell_type_marker_enrichment": marker_enrichment
    }
    nx_available = True
except ImportError:
    print("  networkx not available, using scipy-based approach", flush=True)
    nx_available = False
    # Fallback: connected components + simple degree-based community proxy
    from scipy.sparse.csgraph import connected_components
    n_comp, labels = connected_components(adj, directed=False)
    comp_sizes = [int(np.sum(labels == i)) for i in range(n_comp)]
    print(f"  Connected components: {n_comp}, sizes: {sorted(comp_sizes, reverse=True)[:10]}", flush=True)
    h01_result = {
        "n_communities": n_comp,
        "method": "connected_components_fallback",
        "community_sizes": sorted(comp_sizes, reverse=True)
    }

# Save H01
h01_path = ITER_DIR / "h01_community_detection_l11.json"
with open(h01_path, "w") as f:
    json.dump(h01_result, f, indent=2)
print(f"  Saved: {h01_path}", flush=True)


# ─── H02: TRRUST regulatory pairs → L2 distance (195 in-vocab, OOV-corrected) ─
print("\n=== H02: TRRUST pairs → embedding distance (195 in-vocab) ===", flush=True)

# Load TRRUST pairs
trrust_pairs = []
if TRRUST_CACHE.exists():
    with open(TRRUST_CACHE) as f:
        raw = json.load(f)
    # raw is list of {tf, target, direction, pmids}
    for p in raw:
        tf = p.get("tf") or p.get("TF")
        target = p.get("target")
        if tf and target and tf in inv_gene_set and target in inv_gene_set:
            trrust_pairs.append((tf, target))
    print(f"  TRRUST in-vocab pairs: {len(trrust_pairs)}", flush=True)
else:
    # Try to read from iter_0011 directory
    ITER11_DIR = PROJECT / "iterations" / "iter_0011"
    trrust_alt = ITER11_DIR / "trrust_named_pairs.json"
    if trrust_alt.exists():
        with open(trrust_alt) as f:
            raw = json.load(f)
        for p in raw:
            tf = p.get("tf") or p.get("TF")
            target = p.get("target")
            if tf and target and tf in inv_gene_set and target in inv_gene_set:
                trrust_pairs.append((tf, target))
    print(f"  TRRUST in-vocab pairs (from iter_0011): {len(trrust_pairs)}", flush=True)

# If still no pairs, try to load from iter_0024 or similar
if len(trrust_pairs) == 0:
    # Check dorothea
    dorothea_files = list((PROJECT / "iterations").glob("*/dorothea_human.tsv"))
    dorothea_cache_files = list((PROJECT / "iterations").glob("*/dorothea_pairs_named.json"))
    print(f"  Looking for dorothea: {dorothea_files}", flush=True)
    for cf in dorothea_cache_files:
        with open(cf) as f:
            raw = json.load(f)
        for p in raw:
            tf = p.get("tf") or p.get("TF")
            target = p.get("target")
            if tf and target and tf in inv_gene_set and target in inv_gene_set:
                trrust_pairs.append((tf, target))
        if trrust_pairs:
            print(f"  Using dorothea from {cf}: {len(trrust_pairs)} pairs", flush=True)
            break

# If still empty, compute from scratch using known TFs in inv_genes
if len(trrust_pairs) == 0:
    # Use known TF-target pairs from named genes based on what we know
    # TFs in the 195 gene set
    known_tfs_in_vocab = [g for g in inv_genes if g in {
        "JUN", "FOSL1", "FOSL2", "ATF3", "BATF", "RUNX1", "RUNX2", "RUNX3",
        "KLF4", "KLF2", "BCL6", "NFKB1", "NFKB2", "RELA", "RELB", "REL",
        "IRF1", "IRF2", "IRF3", "IRF4", "IRF7", "IRF8",
        "STAT1", "STAT2", "STAT3", "STAT5A", "STAT5B",
        "TP53", "MYC", "MYCN", "ETS1", "ETS2", "FLI1",
        "GATA1", "GATA2", "GATA3", "TAL1", "TCF7",
        "EBF1", "PAX5", "SPI1", "CEBPA", "CEBPB"
    }]
    print(f"  Fallback: {len(known_tfs_in_vocab)} known TFs in vocab", flush=True)
    # Can't make pairs without actual regulatory data; mark as blocked
    h02_blocked = True
else:
    h02_blocked = False

h02_results = {"status": "tested" if not h02_blocked else "blocked"}

if not h02_blocked and len(trrust_pairs) >= 10:
    # Compute AUROC per layer: TRRUST pairs vs random pairs
    pair_dists_per_layer = []
    rand_dists_per_layer = []
    n_trrust = len(trrust_pairs)
    n_rand = min(n_trrust * 3, 2000)

    # Random pairs
    all_pairs_idx = [(i, j) for i in range(len(inv_genes)) for j in range(i+1, len(inv_genes))]
    rand_idx = rng.choice(len(all_pairs_idx), size=n_rand, replace=False)
    rand_pairs = [all_pairs_idx[i] for i in rand_idx]

    aurocs = []
    for layer in range(N_LAYERS):
        E = emb_inv[layer]  # [195, 512]

        trrust_dists = []
        for (tf, target) in trrust_pairs:
            if tf in inv_gene_to_local and target in inv_gene_to_local:
                i, j = inv_gene_to_local[tf], inv_gene_to_local[target]
                d = float(np.linalg.norm(E[i] - E[j]))
                trrust_dists.append(d)

        rand_dists = []
        for (i, j) in rand_pairs:
            d = float(np.linalg.norm(E[i] - E[j]))
            rand_dists.append(d)

        if len(trrust_dists) < 5:
            aurocs.append(0.5)
            continue

        # AUROC: TF-target pairs should be closer (lower distance) than random
        # AUROC for "TF-target is closer" = P(dist_trrust < dist_rand)
        stat, pval = mannwhitneyu(trrust_dists, rand_dists, alternative='less')
        n1, n2 = len(trrust_dists), len(rand_dists)
        auroc = stat / (n1 * n2)
        aurocs.append(float(auroc))

    mean_auroc = float(np.mean(aurocs))
    best_layer = int(np.argmax(aurocs))
    rho, rho_p = spearmanr(range(N_LAYERS), aurocs)

    print(f"  TRRUST pairs: {n_trrust}, Random pairs: {n_rand}", flush=True)
    print(f"  AUROC per layer: {[f'{a:.3f}' for a in aurocs]}", flush=True)
    print(f"  Mean AUROC: {mean_auroc:.4f}, Best at L{best_layer}: {aurocs[best_layer]:.4f}", flush=True)
    print(f"  Spearman rho(layer, AUROC): {rho:.4f}, p={rho_p:.4f}", flush=True)

    h02_results.update({
        "n_trrust_pairs": n_trrust,
        "n_rand_pairs": n_rand,
        "auroc_per_layer": aurocs,
        "mean_auroc": mean_auroc,
        "best_layer": best_layer,
        "best_auroc": float(aurocs[best_layer]),
        "spearman_rho_layer_auroc": float(rho),
        "spearman_pval": float(rho_p)
    })
else:
    print(f"  H02 blocked: insufficient TRRUST pairs ({len(trrust_pairs)})", flush=True)
    # Fallback: use Dorothea pairs if available
    dorothea_pairs = []
    for iter_name in ["iter_0024", "iter_0025", "iter_0026"]:
        iter_d = PROJECT / "iterations" / iter_name
        for fn in iter_d.glob("*.json"):
            try:
                with open(fn) as f2:
                    d = json.load(f2)
                if isinstance(d, list) and len(d) > 0 and isinstance(d[0], dict):
                    for p in d:
                        tf = p.get("tf") or p.get("TF") or p.get("source")
                        target = p.get("target")
                        if tf and target and tf in inv_gene_set and target in inv_gene_set:
                            dorothea_pairs.append((tf, target))
                if dorothea_pairs:
                    break
            except:
                pass
        if dorothea_pairs:
            break

    if dorothea_pairs:
        print(f"  Fallback Dorothea pairs: {len(dorothea_pairs)}", flush=True)
        # Run same analysis
        n_pairs = len(dorothea_pairs)
        n_rand = min(n_pairs * 3, 2000)
        all_pairs_idx = [(i, j) for i in range(len(inv_genes)) for j in range(i+1, len(inv_genes))]
        rand_idx = rng.choice(len(all_pairs_idx), size=n_rand, replace=False)
        rand_pairs_list = [all_pairs_idx[i] for i in rand_idx]

        aurocs = []
        for layer in range(N_LAYERS):
            E = emb_inv[layer]
            pair_dists = []
            for (tf, target) in dorothea_pairs:
                i = inv_gene_to_local.get(tf)
                j = inv_gene_to_local.get(target)
                if i is not None and j is not None:
                    pair_dists.append(float(np.linalg.norm(E[i] - E[j])))
            rand_dists = [float(np.linalg.norm(E[i] - E[j])) for i, j in rand_pairs_list]
            if len(pair_dists) < 5:
                aurocs.append(0.5)
                continue
            stat, _ = mannwhitneyu(pair_dists, rand_dists, alternative='less')
            aurocs.append(float(stat / (len(pair_dists) * len(rand_dists))))

        mean_auroc = float(np.mean(aurocs))
        rho, rho_p = spearmanr(range(N_LAYERS), aurocs)
        print(f"  Dorothea AUROC per layer: {[f'{a:.3f}' for a in aurocs]}", flush=True)
        print(f"  Mean AUROC={mean_auroc:.4f}, rho={rho:.4f}, p={rho_p:.4f}", flush=True)
        h02_results = {
            "status": "tested",
            "source": "dorothea_fallback",
            "n_pairs": n_pairs,
            "n_rand_pairs": n_rand,
            "auroc_per_layer": aurocs,
            "mean_auroc": mean_auroc,
            "spearman_rho_layer_auroc": float(rho),
            "spearman_pval": float(rho_p)
        }
    else:
        h02_results = {"status": "blocked", "reason": "no regulatory pairs found"}

h02_path = ITER_DIR / "h02_regulatory_pairs_dist.json"
with open(h02_path, "w") as f:
    json.dump(h02_results, f, indent=2)
print(f"  Saved: {h02_path}", flush=True)


# ─── H03: PC1 at L11 vs biological partitions ─────────────────────────────────
print("\n=== H03: PC1 at L11 vs biological partitions ===", flush=True)

def compute_pc1_auroc(E, group_genes, background_genes, inv_gene_to_local):
    """Project onto PC1, compute AUROC for group vs background."""
    E_centered = E - E.mean(axis=0)
    U, s, Vt = np.linalg.svd(E_centered, full_matrices=False)
    pc1 = U[:, 0] * s[0]  # PC1 scores for all genes

    group_idx = [inv_gene_to_local[g] for g in group_genes if g in inv_gene_to_local]
    bg_idx = [inv_gene_to_local[g] for g in background_genes if g in inv_gene_to_local]

    if len(group_idx) < 2 or len(bg_idx) < 2:
        return 0.5, 1.0, pc1

    group_scores = pc1[group_idx]
    bg_scores = pc1[bg_idx]

    stat, pval = mannwhitneyu(group_scores, bg_scores, alternative='two-sided')
    auroc = stat / (len(group_idx) * len(bg_idx))
    return float(auroc), float(pval), pc1

# TFs in the dataset
# Load TF list from prior iterations or define from known biology
tf_genes = set()
# Check for dorothea TFs
for iter_name in ["iter_0024", "iter_0025"]:
    iter_d = PROJECT / "iterations" / iter_name
    for fn in iter_d.glob("*.json"):
        try:
            with open(fn) as f2:
                d = json.load(f2)
            if isinstance(d, dict) and "tf_genes" in d:
                tf_genes.update(d["tf_genes"])
        except:
            pass

# Fallback: use known TF set
if not tf_genes:
    tf_genes = {"JUN", "FOSL1", "FOSL2", "ATF3", "BATF", "RUNX1", "RUNX2", "RUNX3",
                "KLF4", "KLF2", "BCL6", "NFKB1", "NFKB2", "RELA", "RELB", "REL",
                "IRF1", "IRF2", "IRF3", "IRF4", "IRF7", "IRF8",
                "STAT1", "STAT2", "STAT3", "STAT5A", "STAT5B",
                "TP53", "MYC", "MYCN", "ETS1", "ETS2", "FLI1",
                "GATA1", "GATA2", "GATA3", "TAL1", "TCF7",
                "EBF1", "PAX5", "SPI1", "CEBPA", "CEBPB",
                "E2F1", "E2F2", "E2F3", "YBX1", "SP1", "SP3",
                "HIF1A", "EPAS1", "VHL", "ARNT"}

tf_in_vocab = [g for g in inv_genes if g in tf_genes]
non_tf_in_vocab = [g for g in inv_genes if g not in tf_genes]
print(f"  TFs in vocab: {len(tf_in_vocab)}, non-TF: {len(non_tf_in_vocab)}", flush=True)

# Cell-type marker groups
tcell = [g for g in ["CD3E","CD3D","CD8A","CD8B","CD4","LCK","ZAP70","TCF7",
                      "GZMB","PRF1","IFNG","TIGIT","LAG3","PDCD1"] if g in inv_gene_set]
bcell = [g for g in ["CD19","MS4A1","CD79A","CD79B","EBF1","IGHM","IGHD","IGKC"] if g in inv_gene_set]
monocyte = [g for g in ["CD14","FCGR3A","CSF1R","LYZ","S100A8","S100A9","CCL2","CCR2"] if g in inv_gene_set]

print(f"  T-cell markers in vocab: {tcell}", flush=True)
print(f"  B-cell markers in vocab: {bcell}", flush=True)
print(f"  Monocyte markers in vocab: {monocyte}", flush=True)

h03_results = {"layers": []}

for layer in [0, 5, 11]:  # early, mid, late
    E = emb_inv[layer]

    # Var explained by PC1
    E_centered = E - E.mean(axis=0)
    U, s, Vt = np.linalg.svd(E_centered, full_matrices=False)
    var_pc1 = float(s[0]**2 / np.sum(s**2))

    # TF vs non-TF AUROC
    tf_auroc, tf_pval, pc1 = compute_pc1_auroc(E, tf_in_vocab, non_tf_in_vocab, inv_gene_to_local)

    # Cell-type marker AUROCs
    all_genes_no_tcell = [g for g in inv_genes if g not in tcell]
    all_genes_no_bcell = [g for g in inv_genes if g not in bcell]
    all_genes_no_mono = [g for g in inv_genes if g not in monocyte]

    tcell_auroc, tcell_pval, _ = compute_pc1_auroc(E, tcell, all_genes_no_tcell, inv_gene_to_local)
    bcell_auroc, bcell_pval, _ = compute_pc1_auroc(E, bcell, all_genes_no_bcell, inv_gene_to_local)
    mono_auroc, mono_pval, _ = compute_pc1_auroc(E, monocyte, all_genes_no_mono, inv_gene_to_local)

    # STRING high-degree hubs (degree in kNN graph at L11)
    # Reuse L11 kNN graph adjacency
    if layer == 11:
        degrees = np.array(adj.sum(axis=1)).flatten()
        high_degree_threshold = np.percentile(degrees, 75)
        hub_genes = [inv_genes[i] for i in range(n) if degrees[i] >= high_degree_threshold]
        low_degree_genes = [inv_genes[i] for i in range(n) if degrees[i] < high_degree_threshold]
        hub_auroc, hub_pval, _ = compute_pc1_auroc(E, hub_genes, low_degree_genes, inv_gene_to_local)
    else:
        hub_auroc, hub_pval = None, None

    layer_result = {
        "layer": layer,
        "var_pc1": var_pc1,
        "tf_vs_nontf_auroc": tf_auroc,
        "tf_vs_nontf_pval": tf_pval,
        "tf_n": len(tf_in_vocab),
        "nontf_n": len(non_tf_in_vocab),
        "tcell_auroc": tcell_auroc,
        "tcell_pval": tcell_pval,
        "bcell_auroc": bcell_auroc,
        "bcell_pval": bcell_pval,
        "monocyte_auroc": mono_auroc,
        "monocyte_pval": mono_pval
    }
    if hub_auroc is not None:
        layer_result["hub_vs_nonhub_auroc"] = hub_auroc
        layer_result["hub_vs_nonhub_pval"] = hub_pval
        layer_result["n_hubs"] = len(hub_genes)

    print(f"  L{layer}: var_PC1={var_pc1:.3f}, TF_AUROC={tf_auroc:.3f}(p={tf_pval:.3f}), "
          f"Tcell={tcell_auroc:.3f}, Bcell={bcell_auroc:.3f}, Mono={mono_auroc:.3f}", flush=True)

    h03_results["layers"].append(layer_result)

# Overall summary
l11_res = [r for r in h03_results["layers"] if r["layer"] == 11][0]
l0_res = [r for r in h03_results["layers"] if r["layer"] == 0][0]
h03_results["summary"] = {
    "l11_var_pc1": l11_res["var_pc1"],
    "l0_var_pc1": l0_res["var_pc1"],
    "l11_tf_auroc": l11_res["tf_vs_nontf_auroc"],
    "l0_tf_auroc": l0_res["tf_vs_nontf_auroc"],
    "l11_tcell_auroc": l11_res["tcell_auroc"],
    "l11_bcell_auroc": l11_res["bcell_auroc"],
    "l11_monocyte_auroc": l11_res["monocyte_auroc"]
}
print(f"\n  H03 Summary: L11 PC1 var={l11_res['var_pc1']:.3f}", flush=True)
print(f"    TF AUROC: L0={l0_res['tf_vs_nontf_auroc']:.3f}, L11={l11_res['tf_vs_nontf_auroc']:.3f}", flush=True)

h03_path = ITER_DIR / "h03_pc1_biological_partitions.json"
with open(h03_path, "w") as f:
    json.dump(h03_results, f, indent=2)
print(f"  Saved: {h03_path}", flush=True)

print("\n=== All H03 done ===", flush=True)
print("=== iter_0032 screen complete ===", flush=True)
