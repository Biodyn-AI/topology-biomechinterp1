#!/usr/bin/env python3
"""E9 — Reproducible derivation of the 4,803 / 195 gene vocabulary.

Reviewer 2 question: "Why are 4,803 genes considered? What is the filtering criterion?"

Trace:
  1. Tabula Sapiens immune subset is preprocessed via Scanpy (HVG selection,
     log1p, normalisation). The output `tabula_sapiens_processed.h5ad` has
     n_vars = 4803 — all 4,803 entries are flagged `highly_variable=True`.
  2. The 4,803 `_index` entries are HGNC symbols and are 100% covered by
     scGPT's whole-human vocabulary (60,694 gene tokens).
  3. The paper's "209 named genes" are the union of TFs and targets in the
     immune-curated regulatory edge list (`cycle1_edge_dataset.tsv`, 29,948 edges
     over 209 unique gene names). 14 of these have zero embedding norm at every
     layer (out-of-vocabulary at the input layer of scGPT in the immune-cell
     forward pass) → 195 in-vocabulary named genes used for biological tests.

This script runs the audit and emits a CSV showing the gene counts at each step.
"""
from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]  # biodyn-work
H5AD = REPO_ROOT / "single_cell_mechinterp/outputs/tabula_sapiens_processed.h5ad"
SCGPT_VOCAB = REPO_ROOT / "single_cell_mechinterp/external/scGPT_checkpoints/whole-human/vocab.json"
TRRUST = REPO_ROOT / "single_cell_mechinterp/external/networks/trrust_human.tsv"
STRING_NAMED = REPO_ROOT / "subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0012/string_ppi_named_genes.json"
CYCLE_TSV = REPO_ROOT / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main/cycle1_edge_dataset.tsv"
EMBEDDINGS = REPO_ROOT / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main/layer_gene_embeddings.npy"

OUT = Path(__file__).resolve().parent


def _h5ad_gene_index() -> list[str]:
    with h5py.File(H5AD, "r") as f:
        raw = f["var"]["_index"][:]
        is_hvg = f["var"]["highly_variable"][:]
    assert is_hvg.sum() == len(raw), "Expected all 4,803 entries flagged HVG"
    return [g.decode() if isinstance(g, bytes) else str(g) for g in raw]


def _vocab() -> set[str]:
    with open(SCGPT_VOCAB) as f:
        return set(json.load(f).keys())


def _trrust_gene_set() -> set[str]:
    t = pd.read_csv(TRRUST, sep="\t", header=None, names=["TF", "target", "mode", "pmid"])
    return set(t.TF.astype(str)) | set(t.target.astype(str))


def _string_named() -> set[str]:
    if not STRING_NAMED.exists():
        return set()
    obj = json.loads(STRING_NAMED.read_text())
    if isinstance(obj, list):
        if obj and isinstance(obj[0], list):  # list of [g1, g2, score]
            genes: set[str] = set()
            for row in obj:
                if isinstance(row, list) and len(row) >= 2:
                    genes.add(str(row[0]))
                    genes.add(str(row[1]))
            return genes
        return set(str(x) for x in obj)
    return set()


def _cycle_named() -> set[str]:
    if not CYCLE_TSV.exists():
        return set()
    edge = pd.read_csv(CYCLE_TSV, sep="\t")
    if "source" in edge.columns and "target" in edge.columns:
        return set(edge.source.astype(str)) | set(edge.target.astype(str))
    if "tf" in edge.columns and "target" in edge.columns:
        return set(edge.tf.astype(str)) | set(edge.target.astype(str))
    return set()


def _zero_norm_genes(gene_index: list[str]) -> set[str]:
    if not EMBEDDINGS.exists():
        return set()
    emb = np.load(EMBEDDINGS, mmap_mode="r")  # (12, 4803, 512)
    norms = np.linalg.norm(emb, axis=2)  # (12, 4803)
    # zero across ALL layers => OOV
    is_zero_all = (norms == 0).all(axis=0)
    return {g for g, z in zip(gene_index, is_zero_all) if z}


def main() -> None:
    rows = []
    gene_index = _h5ad_gene_index()
    gene_set = set(gene_index)
    vocab = _vocab()
    trrust = _trrust_gene_set()
    string_named = _string_named()
    cycle_named = _cycle_named()

    rows.append(("Tabula Sapiens preprocessed h5ad (HVG)", len(gene_index)))
    rows.append(("scGPT whole-human vocab tokens", len(vocab)))
    rows.append(("h5ad ∩ scGPT vocab", len(gene_set & vocab)))
    rows.append(("TRRUST unique genes (TFs + targets)", len(trrust)))
    rows.append(("TRRUST ∩ h5ad", len(trrust & gene_set)))
    rows.append(("STRING named-genes file (immune subset)", len(string_named)))
    rows.append(("STRING-named ∩ h5ad", len(string_named & gene_set)))
    rows.append(("cycle1_edge_dataset gene set (paper '209 named')", len(cycle_named)))
    rows.append(("cycle1 ∩ h5ad", len(cycle_named & gene_set)))

    if EMBEDDINGS.exists():
        zero_norm = _zero_norm_genes(gene_index)
        rows.append(("h5ad genes with zero embedding norm at every layer (OOV in cell forward pass)", len(zero_norm)))
        cycle_in_vocab = (cycle_named & gene_set) - zero_norm
        rows.append(("cycle1 ∩ h5ad ∩ embedded (paper's 195 in-vocab)", len(cycle_in_vocab)))

    out_df = pd.DataFrame(rows, columns=["step", "count"])
    out_df.to_csv(OUT / "vocab_derivation.csv", index=False)
    print(out_df.to_string(index=False))


if __name__ == "__main__":
    main()
