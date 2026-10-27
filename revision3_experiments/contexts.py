"""Round-3 data layer: every scGPT residual-stream context, correctly labelled.

The round-2 primary matrix `cycle1_main` (12 x 4803 x 512) was described in the
manuscript as Tabula Sapiens immune cells. Its gene axis comes from
`tabula_sapiens_processed.h5ad`, which is the Tabula Sapiens KIDNEY dataset
(11,376 cells, tissue_in_publication == "Kidney" for every cell; 82% kidney
epithelium), and in the 256 cells behind the embedding the proximal-tubule
markers LRP2/CUBN are expressed in 137/114 cells against 43 for PTPRC. It is
therefore labelled `kidney_A` here. The genuinely immune runs share the 4,941-gene
axis of `tabula_sapiens_immune_subset_hpn_processed.h5ad`.

Each context = one scGPT forward-pass extraction: 256 sampled cells, per-layer
mean-pooled gene embeddings (12 x G x 512), gene axis = the h5ad var index.
"""
from __future__ import annotations
import gzip, json
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[2]
S38 = REPO / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs"
S41 = REPO / "subproject_41_claude_topology_hypothesis_screening_autoloop/iterations"
SC = REPO / "single_cell_mechinterp/outputs"
ICE = SC / "invariant_causal_edges"

H5_IMMUNE = SC / "tabula_sapiens_immune_subset_hpn_processed.h5ad"     # 20,000 x 4,941
H5_KIDNEY_A = SC / "tabula_sapiens_processed.h5ad"                      # 11,376 x 4,803 (TS kidney)
H5_KIDNEY_B = ICE / "kidney/processed.h5ad"                             #  8,000 x 8,143
H5_LUNG = ICE / "lung/processed.h5ad"                                   # 15,000 x 8,181
H5_EXT_LUNG = ICE / "external_lung/processed.h5ad"                      #  9,000 x 8,229 (Krasnow)

def _c(emb, h5, tissue, group, note=""):
    return dict(emb=emb, h5ad=h5, tissue=tissue, group=group, note=note)

CONTEXTS = {
    # --- immune: primary -----------------------------------------------------------
    "immune":          _c(S38/"cycle4_immune_main/layer_gene_embeddings.npy",  H5_IMMUNE, "immune", "primary"),
    "immune_seed43":   _c(S38/"cycle4_immune_seed43/layer_gene_embeddings.npy", H5_IMMUNE, "immune", "resample", "independent 256-cell draw"),
    "immune_seed44":   _c(S38/"cycle4_immune_seed44/layer_gene_embeddings.npy", H5_IMMUNE, "immune", "resample", "independent 256-cell draw"),
    "immune_TSP14":    _c(S38/"cycle10_immune_donor_TSP14_seed42/layer_gene_embeddings.npy", H5_IMMUNE, "immune", "donor"),
    "immune_TSP21":    _c(S38/"cycle10_immune_donor_TSP21_seed42/layer_gene_embeddings.npy", H5_IMMUNE, "immune", "donor"),
    "immune_TSP25":    _c(S38/"cycle10_immune_donor_TSP25_seed42/layer_gene_embeddings.npy", H5_IMMUNE, "immune", "donor"),
    "immune_Bcell":    _c(S38/"cycle11_immune_celltype_b_cell/layer_gene_embeddings.npy", H5_IMMUNE, "immune", "celltype"),
    "immune_CD4T":     _c(S38/"cycle11_immune_celltype_cd4_t/layer_gene_embeddings.npy", H5_IMMUNE, "immune", "celltype"),
    "immune_CD8T":     _c(S38/"cycle11_immune_celltype_cd8_t/layer_gene_embeddings.npy", H5_IMMUNE, "immune", "celltype"),
    # --- kidney ----------------------------------------------------------------------
    "kidney_A":        _c(S38/"cycle1_main/layer_gene_embeddings.npy", H5_KIDNEY_A, "kidney", "tissue",
                          "4,803-gene HVG axis; the matrix analysed in the submitted version"),
    "kidney_B":        _c(S41/"iter_0067/kidney_scgpt_run/layer_gene_embeddings.npy", H5_KIDNEY_B, "kidney", "tissue",
                          "8,143-gene HVG axis"),
    "kidney_ckpt_whole_human": _c(S41/"iter_0087/kidney_scgpt_whole_human_run/layer_gene_embeddings.npy", H5_KIDNEY_B, "kidney", "checkpoint"),
    "kidney_ckpt_kidney":      _c(S41/"iter_0087/kidney_scgpt_kidney_run/layer_gene_embeddings.npy", H5_KIDNEY_B, "kidney", "checkpoint"),
    "kidney_ckpt_brain":       _c(S41/"iter_0087/kidney_scgpt_brain_run/layer_gene_embeddings.npy", H5_KIDNEY_B, "kidney", "checkpoint"),
    # --- lung --------------------------------------------------------------------------
    "lung":            _c(S38/"cycle6_lung_main/layer_gene_embeddings.npy", H5_LUNG, "lung", "tissue", "Tabula Sapiens lung"),
    "lung_external":   _c(S38/"cycle7_external_lung_main/layer_gene_embeddings.npy", H5_EXT_LUNG, "lung", "tissue",
                          "Krasnow lung atlas, independent study"),
}

def h5ad_genes(path):
    import h5py
    with h5py.File(path, "r") as f:
        raw = f["var"]["_index"][:]
    return [g.decode() if isinstance(g, bytes) else str(g) for g in raw]

def load_context(name, mmap=True):
    c = CONTEXTS[name]
    emb = np.load(c["emb"], mmap_mode="r" if mmap else None)
    genes = h5ad_genes(c["h5ad"])
    assert emb.shape[1] == len(genes), f"{name}: emb {emb.shape} vs {len(genes)} genes"
    return emb, genes

# --------------------------------------------------------------------------------------
# reference networks
# --------------------------------------------------------------------------------------
NET = REPO / "single_cell_mechinterp/external/networks"
STRING_LINKS = REPO / "subproject_61_two_axes_revision/data/string_links_v12.txt.gz"
STRING_INFO = REPO / "subproject_61_two_axes_revision/data/string_protein_info_v12.txt.gz"

def load_string_genomewide(min_score=700):
    """Genome-wide STRING v12 combined-score links mapped to gene symbols (undirected, deduplicated)."""
    ens2sym = {}
    with gzip.open(STRING_INFO, "rt") as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            ens2sym[p[0]] = p[1]
    pairs = {}
    with gzip.open(STRING_LINKS, "rt") as f:
        next(f)
        for line in f:
            a, b, s = line.split()
            s = int(s)
            if s < min_score: continue
            ga, gb = ens2sym.get(a), ens2sym.get(b)
            if not ga or not gb or ga == gb: continue
            k = (ga, gb) if ga < gb else (gb, ga)
            pairs[k] = max(pairs.get(k, 0), s)
    return pairs            # {(g1,g2): combined_score}

def load_trrust():
    import pandas as pd
    t = pd.read_csv(NET/"trrust_human.tsv", sep="\t", header=None, names=["TF","target","mode","pmid"])
    t["TF"] = t.TF.astype(str); t["target"] = t.target.astype(str)
    pairs = sorted({(a, b) for a, b in zip(t.TF, t.target) if a != b})
    return t, pairs

def load_dorothea(chipseq=False, levels=("A","B","C")):
    import pandas as pd
    d = pd.read_csv(NET/("dorothea_chipseq_human.tsv" if chipseq else "dorothea_human.tsv"), sep="\t")
    keep = d.confidence.astype(str).apply(lambda s: any(l in s.split(";") for l in levels))
    d = d[keep]
    return sorted({(a, b) for a, b in zip(d.source.astype(str), d.target.astype(str)) if a != b})
