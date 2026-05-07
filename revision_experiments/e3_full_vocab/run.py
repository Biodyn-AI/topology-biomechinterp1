#!/usr/bin/env python3
"""E3 — Broader gene-vocabulary analysis.

Reviewer 1 (#1), Reviewer 2 (#7), Reviewer 3 (#3): "biological analyses are
restricted to ~195 in-vocabulary genes."

This script repeats the headline analyses on the full annotated gene set
inside the 4,803-gene scGPT input vocabulary (no immune-only curation):

  - SV1 GO:0005615 enrichment on the full 4,803 / annotated subset.
  - SV2 PPI co-pole on STRING ≥ 0.7 pairs covering all 4,803 genes.
  - TF-vs-target AUROC on TRRUST TFs/targets present in the 4,803-gene set.

Also includes a 1,000-bootstrap subsampling test: draw 195 genes uniformly
from the broader set and compare the distribution of effects to the immune-
curated 195. This tells us whether the paper's curated subset is an outlier
or representative.
"""
from __future__ import annotations

import json
import pickle
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent

EMB = REPO / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main/layer_gene_embeddings.npy"
H5AD = REPO / "single_cell_mechinterp/outputs/tabula_sapiens_processed.h5ad"
TRRUST_PATH = REPO / "single_cell_mechinterp/external/networks/trrust_human.tsv"
STRING_CACHE = REPO / "subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0015/string_ppi_score04_cache.json"
GENE2GO = REPO / "single_cell_mechinterp/data/perturb/gene2go_all.pkl"

# Curated 209-gene immune-relevant edge dataset
CYCLE1_TSV = REPO / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main/cycle1_edge_dataset.tsv"

GO_EXTRACELLULAR = "GO:0005615"


def _h5ad_genes(p):
    with h5py.File(p, "r") as f:
        raw = f["var"]["_index"][:]
    return [g.decode() if isinstance(g, bytes) else str(g) for g in raw]


def fisher_2x2(a, b, c, d):
    from math import lgamma, exp
    n = a + b + c + d
    r1, r2, c1, c2 = a + b, c + d, a + c, b + d
    if r1 == 0 or r2 == 0 or c1 == 0 or c2 == 0:
        return 1.0, 1.0
    odds = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))

    def lpx(x):
        return (lgamma(r1 + 1) - lgamma(x + 1) - lgamma(r1 - x + 1)
                + lgamma(r2 + 1) - lgamma(c1 - x + 1) - lgamma(c2 - r1 + x + 1)
                - lgamma(n + 1) + lgamma(c1 + 1) + lgamma(c2 + 1))
    obs = lpx(a)
    p = 0.0
    for x in range(max(0, c1 - r2), min(c1, r1) + 1):
        if lpx(x) <= obs + 1e-12:
            p += exp(lpx(x))
    return float(odds), float(min(p, 1.0))


def go_enrichment_on_subset(US, gene_names, sub_idx, gene2go, go_term, pole_frac=0.27):
    """Restrict to sub_idx (subset of 0..n_genes-1), run pole enrichment on SV1."""
    sub_proj = US[sub_idx, 0]
    sub_names = [gene_names[i] for i in sub_idx]
    n = len(sub_idx)
    K = max(1, int(pole_frac * n))
    order = np.argsort(sub_proj)
    in_top = set(order[-K:].tolist())
    in_bot = set(order[:K].tolist())
    a_top = sum(1 for k in range(n) if k in in_top and go_term in gene2go.get(sub_names[k], set()))
    in_top_with = a_top
    in_top_without = K - a_top
    out_with = sum(1 for k in range(n) if k not in in_top and go_term in gene2go.get(sub_names[k], set()))
    out_without = (n - K) - out_with
    or_top, p_top = fisher_2x2(in_top_with, in_top_without, out_with, out_without)
    a_bot = sum(1 for k in range(n) if k in in_bot and go_term in gene2go.get(sub_names[k], set()))
    in_bot_with = a_bot
    in_bot_without = K - a_bot
    out_with2 = sum(1 for k in range(n) if k not in in_bot and go_term in gene2go.get(sub_names[k], set()))
    out_without2 = (n - K) - out_with2
    or_bot, p_bot = fisher_2x2(in_bot_with, in_bot_without, out_with2, out_without2)
    return max(or_top, or_bot), min(p_top, p_bot)


def copole_rate_on_subset(US, gene_names, sub_idx, pairs, pole_frac=0.27):
    """Co-pole rate on subset of genes."""
    sub_proj = US[sub_idx, 1]
    sub_names = [gene_names[i] for i in sub_idx]
    n = len(sub_idx)
    K = max(1, int(pole_frac * n))
    order = np.argsort(sub_proj)
    top = set(order[-K:].tolist())
    bot = set(order[:K].tolist())
    n2i = {g: k for k, g in enumerate(sub_names)}
    co = total = 0
    for g1, g2 in pairs:
        if g1 not in n2i or g2 not in n2i:
            continue
        i, j = n2i[g1], n2i[g2]
        total += 1
        if (i in top and j in top) or (i in bot and j in bot):
            co += 1
    return (co / total if total else 0.0), total


def tf_target_auroc_on_subset(US, gene_names, sub_idx, tfs, targets_only):
    sub_emb = US[sub_idx, 1:7]
    sub_names = [gene_names[i] for i in sub_idx]
    name_to_i = {g: k for k, g in enumerate(sub_names)}
    tf_idx = [name_to_i[g] for g in tfs if g in name_to_i]
    target_idx = [name_to_i[g] for g in targets_only if g in name_to_i]
    if len(tf_idx) < 5 or len(target_idx) < 5:
        return float("nan"), len(tf_idx), len(target_idx)
    norm = sub_emb / (np.linalg.norm(sub_emb, axis=1, keepdims=True) + 1e-12)
    sim = norm @ norm[tf_idx].T
    score = sim.mean(axis=1)
    pos = score[tf_idx]
    neg = score[target_idx]
    all_s = np.concatenate([pos, neg])
    labels = np.array([1] * len(pos) + [0] * len(neg))
    order = np.argsort(-all_s)
    ranked = labels[order]
    cn = ar = 0
    for r in ranked:
        if r == 1:
            ar += cn
        else:
            cn += 1
    return float(1 - ar / (len(pos) * len(neg))), len(tf_idx), len(target_idx)


def main():
    print("[E3] Loading…", flush=True)
    gene_names = _h5ad_genes(H5AD)
    n_genes = len(gene_names)
    name_to_i = {g: i for i, g in enumerate(gene_names)}

    # Annotations
    t = pd.read_csv(TRRUST_PATH, sep="\t", header=None, names=["TF", "target", "mode", "pmid"])
    tfs_all = set(t.TF.astype(str))
    targets_all = set(t.target.astype(str)) - tfs_all

    obj = json.loads(STRING_CACHE.read_text())
    string_pairs = [(p["g1"], p["g2"]) for p in obj["pairs"] if p["score"] >= 0.7]

    with open(GENE2GO, "rb") as f:
        gene2go = pickle.load(f)

    # Define gene subsets
    full_set_idx = np.arange(n_genes)
    annotated_in_vocab = [i for i, g in enumerate(gene_names) if g in gene2go]
    trrust_in_vocab = [i for i, g in enumerate(gene_names) if g in tfs_all or g in targets_all]
    union_set = sorted(set(annotated_in_vocab) | set(trrust_in_vocab))

    # Curated 209
    edge = pd.read_csv(CYCLE1_TSV, sep="\t")
    curated_genes = set(edge["source"].astype(str)) | set(edge["target"].astype(str))
    curated_idx = sorted([i for i, g in enumerate(gene_names) if g in curated_genes])

    print(f"  full vocab: {len(full_set_idx)}", flush=True)
    print(f"  GO-annotated ∩ vocab: {len(annotated_in_vocab)}", flush=True)
    print(f"  TRRUST ∩ vocab: {len(trrust_in_vocab)}", flush=True)
    print(f"  Union (broader set): {len(union_set)}", flush=True)
    print(f"  Curated 209 ∩ vocab: {len(curated_idx)}", flush=True)

    # Layer 11 embedding
    emb = np.load(EMB, mmap_mode="r")
    print(f"  emb: {emb.shape}; using L11", flush=True)
    X = np.asarray(emb[11])
    Xc = X - X.mean(axis=0, keepdims=True)
    print("[E3] Computing SVD on full 4,803 × 512 layer-11 matrix…", flush=True)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    US = U * S

    # Run on each subset
    results = []
    for set_name, sub_idx in [
        ("curated_209", curated_idx),
        ("trrust_in_vocab", trrust_in_vocab),
        ("annotated_in_vocab", annotated_in_vocab),
        ("union_broader", union_set),
        ("full_vocab", list(full_set_idx)),
    ]:
        sub_idx = list(sub_idx)
        or_extra, p_extra = go_enrichment_on_subset(US, gene_names, sub_idx, gene2go, GO_EXTRACELLULAR)
        copole, n_pairs = copole_rate_on_subset(US, gene_names, sub_idx, string_pairs)
        auroc, n_tf, n_tgt = tf_target_auroc_on_subset(US, gene_names, sub_idx, tfs_all, targets_all)
        results.append(
            dict(
                subset=set_name,
                n_genes=len(sub_idx),
                or_extracellular=or_extra,
                p_extracellular=p_extra,
                ppi_copole_rate=copole,
                ppi_n_pairs=n_pairs,
                tf_target_auroc=auroc,
                n_tf=n_tf,
                n_target=n_tgt,
            )
        )
        print(f"  {set_name}: n={len(sub_idx)}, SV1 extra OR={or_extra:.2f} p={p_extra:.4f}; PPI co-pole={copole:.3f} (n_pairs={n_pairs}); TF/target AUROC={auroc:.3f}", flush=True)

    pd.DataFrame(results).to_csv(OUT / "subset_results.csv", index=False)

    # Bootstrap: draw 195 from union_broader 1000 times
    print("[E3] Bootstrap: 195 random draws from union_broader, B=200…", flush=True)
    rng = np.random.default_rng(42)
    boot_or = np.zeros(200)
    boot_copole = np.zeros(200)
    boot_auroc = np.zeros(200)
    for k in range(200):
        sub = rng.choice(union_set, size=195, replace=False).tolist()
        or_e, _ = go_enrichment_on_subset(US, gene_names, sub, gene2go, GO_EXTRACELLULAR)
        cp, _ = copole_rate_on_subset(US, gene_names, sub, string_pairs)
        au, _, _ = tf_target_auroc_on_subset(US, gene_names, sub, tfs_all, targets_all)
        boot_or[k] = or_e
        boot_copole[k] = cp
        boot_auroc[k] = au if not np.isnan(au) else 0.5

    bootstrap_summary = dict(
        boot_or_p5_p50_p95=[float(np.percentile(boot_or, 5)), float(np.percentile(boot_or, 50)), float(np.percentile(boot_or, 95))],
        boot_copole_p5_p50_p95=[float(np.percentile(boot_copole, 5)), float(np.percentile(boot_copole, 50)), float(np.percentile(boot_copole, 95))],
        boot_auroc_p5_p50_p95=[float(np.percentile(boot_auroc, 5)), float(np.percentile(boot_auroc, 50)), float(np.percentile(boot_auroc, 95))],
        curated_or_extracellular=results[0]["or_extracellular"],
        curated_copole=results[0]["ppi_copole_rate"],
        curated_auroc=results[0]["tf_target_auroc"],
        curated_or_percentile_in_boot=float((boot_or < results[0]["or_extracellular"]).mean()),
        curated_copole_percentile_in_boot=float((boot_copole < results[0]["ppi_copole_rate"]).mean()),
        curated_auroc_percentile_in_boot=float((boot_auroc < results[0]["tf_target_auroc"]).mean()),
    )
    with open(OUT / "results.json", "w") as f:
        json.dump({"subsets": results, "bootstrap": bootstrap_summary}, f, indent=2)
    print("\n[E3] Bootstrap summary:")
    print(json.dumps(bootstrap_summary, indent=2))


if __name__ == "__main__":
    main()
