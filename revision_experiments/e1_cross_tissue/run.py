#!/usr/bin/env python3
"""E1 — Cross-tissue replication of headline findings.

Replicates the paper's headline analyses on kidney scGPT embeddings:
  (1) SV1 GO:0005615 (extracellular space) co-pole enrichment
  (2) SV2 STRING ≥ 0.7 PPI co-pole rate (12-layer mean and per-layer)
  (3) TF-vs-target AUROC on the joint SV2-SV7 subspace

For each finding, the script also reports the immune (paper-original) value as
a comparison anchor. Outputs a single results.json plus per-layer CSVs for
manuscript Table S7.
"""
from __future__ import annotations

import json
import pickle
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent

# Embedding sources -----------------------------------------------------------
EMB_SOURCES = {
    "immune": REPO / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main/layer_gene_embeddings.npy",
    "kidney": REPO / "subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0067/kidney_scgpt_run/layer_gene_embeddings.npy",
}
H5AD_SOURCES = {
    "immune": REPO / "single_cell_mechinterp/outputs/tabula_sapiens_processed.h5ad",
    "kidney": REPO / "single_cell_mechinterp/outputs/invariant_causal_edges/kidney/processed.h5ad",
}

# Annotation sources ---------------------------------------------------------
TRRUST_PATH = REPO / "single_cell_mechinterp/external/networks/trrust_human.tsv"
STRING_CACHE = REPO / "subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0015/string_ppi_score04_cache.json"
GENE2GO = REPO / "single_cell_mechinterp/data/perturb/gene2go_all.pkl"

GO_EXTRACELLULAR = "GO:0005615"  # extracellular space
GO_CYTOSOL = "GO:0005829"        # cytosol


def _h5ad_genes(p: Path) -> list[str]:
    with h5py.File(p, "r") as f:
        raw = f["var"]["_index"][:]
    return [g.decode() if isinstance(g, bytes) else str(g) for g in raw]


def _load_trrust() -> tuple[set[str], set[str], set[tuple[str, str]]]:
    t = pd.read_csv(TRRUST_PATH, sep="\t", header=None, names=["TF", "target", "mode", "pmid"])
    tfs = set(t.TF.astype(str))
    targets = set(t.target.astype(str)) - tfs  # target-only
    pairs = set(zip(t.TF.astype(str), t.target.astype(str)))
    return tfs, targets, pairs


def _load_string_pairs(min_score: float = 0.7) -> list[tuple[str, str]]:
    obj = json.loads(STRING_CACHE.read_text())
    return [(p["g1"], p["g2"]) for p in obj["pairs"] if p["score"] >= min_score]


def _load_gene2go() -> dict[str, set[str]]:
    with open(GENE2GO, "rb") as f:
        d = pickle.load(f)
    return {k: set(v) if not isinstance(v, set) else v for k, v in d.items()}


# ---------- Analysis primitives ---------------------------------------------

def svd_gene_projections(emb_layer: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Mean-centre rows, run thin SVD; return U @ S (n_genes x k)."""
    X = emb_layer - emb_layer.mean(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    return U * S, S


def fisher_2x2(a: int, b: int, c: int, d: int) -> tuple[float, float]:
    """Fisher exact two-sided. Returns (OR, p)."""
    from math import comb, lgamma, exp

    n = a + b + c + d
    r1 = a + b
    r2 = c + d
    c1 = a + c
    c2 = b + d
    if r1 == 0 or r2 == 0 or c1 == 0 or c2 == 0:
        return 1.0, 1.0
    odds = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    # log p_x using lgamma
    def lpx(x: int) -> float:
        return (lgamma(r1 + 1) - lgamma(x + 1) - lgamma(r1 - x + 1)
                + lgamma(r2 + 1) - lgamma(c1 - x + 1) - lgamma(c2 - r1 + x + 1)
                - lgamma(n + 1) + lgamma(c1 + 1) + lgamma(c2 + 1))
    obs_lp = lpx(a)
    p = 0.0
    for x in range(max(0, c1 - r2), min(c1, r1) + 1):
        if lpx(x) <= obs_lp + 1e-12:
            p += exp(lpx(x))
    return float(odds), float(min(p, 1.0))


def go_pole_enrichment(
    proj: np.ndarray,           # (n_genes,)
    gene_names: list[str],
    gene2go: dict[str, set[str]],
    go_term: str,
    pole_frac: float = 0.27,
    pole: str = "top",
) -> tuple[float, float, int, int]:
    """Test whether genes in the top/bottom `pole_frac` of `proj` are enriched
    for `go_term`."""
    n = len(gene_names)
    K = max(1, int(pole_frac * n))
    order = np.argsort(proj)
    if pole == "top":
        pole_idx = set(order[-K:].tolist())
    else:
        pole_idx = set(order[:K].tolist())
    in_pole_with = sum(1 for i in pole_idx if go_term in gene2go.get(gene_names[i], set()))
    in_pole_without = K - in_pole_with
    out_with = sum(1 for i in range(n) if i not in pole_idx and go_term in gene2go.get(gene_names[i], set()))
    out_without = (n - K) - out_with
    odds, p = fisher_2x2(in_pole_with, in_pole_without, out_with, out_without)
    return odds, p, in_pole_with, in_pole_with + out_with


def copole_rate(proj: np.ndarray, gene_names: list[str], pairs: Iterable[tuple[str, str]], pole_frac: float = 0.27) -> tuple[float, int]:
    """Fraction of pairs where both genes fall in same pole (top OR bottom)."""
    name_to_i = {g: i for i, g in enumerate(gene_names)}
    n = len(gene_names)
    K = max(1, int(pole_frac * n))
    order = np.argsort(proj)
    top = set(order[-K:].tolist())
    bot = set(order[:K].tolist())

    co = total = 0
    for g1, g2 in pairs:
        if g1 not in name_to_i or g2 not in name_to_i:
            continue
        i, j = name_to_i[g1], name_to_i[g2]
        total += 1
        if (i in top and j in top) or (i in bot and j in bot):
            co += 1
    return (co / total if total else 0.0), total


def copole_z(proj: np.ndarray, gene_names: list[str], pairs: list[tuple[str, str]], n_perm: int = 200, pole_frac: float = 0.27, rng: np.random.Generator = None) -> dict:
    rng = rng or np.random.default_rng(42)
    obs_rate, n_pairs = copole_rate(proj, gene_names, pairs, pole_frac)
    null_rates = np.zeros(n_perm)
    perm_idx = np.arange(len(gene_names))
    name_to_i = {g: i for i, g in enumerate(gene_names)}
    n = len(gene_names)
    K = max(1, int(pole_frac * n))
    for k in range(n_perm):
        # permute gene name -> index assignment
        rng.shuffle(perm_idx)
        permuted_names = [gene_names[perm_idx[i]] for i in range(n)]
        # use fixed proj (so the poles are same set), but gene labels shuffled
        n2i = {g: i for i, g in enumerate(permuted_names)}
        order = np.argsort(proj)
        top = set(order[-K:].tolist())
        bot = set(order[:K].tolist())
        co = total = 0
        for g1, g2 in pairs:
            if g1 not in n2i or g2 not in n2i:
                continue
            i, j = n2i[g1], n2i[g2]
            total += 1
            if (i in top and j in top) or (i in bot and j in bot):
                co += 1
        null_rates[k] = co / total if total else 0.0
    null_mean = null_rates.mean()
    null_std = null_rates.std(ddof=1) if null_rates.std(ddof=1) > 0 else 1e-9
    z = (obs_rate - null_mean) / null_std
    emp_p = float((null_rates >= obs_rate).mean())
    return dict(obs=obs_rate, null_mean=float(null_mean), z=float(z), emp_p=emp_p, n_pairs=n_pairs, n_perm=n_perm)


def tf_target_auroc(proj_subspace: np.ndarray, gene_names: list[str], tfs: set[str], targets_only: set[str], n_perm: int = 100, rng: np.random.Generator = None) -> dict:
    """Mean cosine similarity to known TFs as a TF-score; AUROC against a
    target-only label set, with permutation null."""
    rng = rng or np.random.default_rng(42)
    name_to_i = {g: i for i, g in enumerate(gene_names)}
    tf_idx = [name_to_i[g] for g in tfs if g in name_to_i]
    target_idx = [name_to_i[g] for g in targets_only if g in name_to_i]
    if len(tf_idx) < 5 or len(target_idx) < 5:
        return dict(auroc=float("nan"), perm_p=float("nan"), n_tf=len(tf_idx), n_target=len(target_idx))

    # Cosine similarity of every gene to the TF set's mean cosine
    norm_proj = proj_subspace / (np.linalg.norm(proj_subspace, axis=1, keepdims=True) + 1e-12)
    sim = norm_proj @ norm_proj[tf_idx].T  # (n, n_tf)
    score = sim.mean(axis=1)  # higher score => more TF-like

    # AUROC
    pos = score[tf_idx]
    neg = score[target_idx]
    all_scores = np.concatenate([pos, neg])
    labels = np.array([1] * len(pos) + [0] * len(neg))
    order = np.argsort(-all_scores)
    ranked = labels[order]
    n_pos, n_neg = len(pos), len(neg)
    cumneg = 0
    auroc_num = 0
    for r in ranked:
        if r == 1:
            auroc_num += cumneg
        else:
            cumneg += 1
    auroc = 1.0 - auroc_num / (n_pos * n_neg)

    # Permutation null
    perm_aurocs = np.zeros(n_perm)
    combined = np.array(tf_idx + target_idx)
    n_tf = len(tf_idx)
    for k in range(n_perm):
        perm = rng.permutation(combined)
        perm_tf = perm[:n_tf]
        perm_targ = perm[n_tf:]
        sim_p = norm_proj @ norm_proj[perm_tf].T
        score_p = sim_p.mean(axis=1)
        pos_p = score_p[perm_tf]
        neg_p = score_p[perm_targ]
        ap = np.concatenate([pos_p, neg_p])
        al = np.array([1] * n_tf + [0] * len(perm_targ))
        order_p = np.argsort(-ap)
        ranked_p = al[order_p]
        cn = ar = 0
        for r in ranked_p:
            if r == 1:
                ar += cn
            else:
                cn += 1
        perm_aurocs[k] = 1.0 - ar / (n_tf * len(perm_targ))
    perm_p = float((perm_aurocs >= auroc).mean())
    return dict(auroc=float(auroc), perm_p=perm_p, perm_mean=float(perm_aurocs.mean()), n_tf=n_tf, n_target=len(target_idx))


def main() -> None:
    print("[E1] Loading annotations…")
    tfs, targets_only, trrust_pairs = _load_trrust()
    string_pairs = _load_string_pairs(0.7)
    gene2go = _load_gene2go()
    print(f"  TRRUST TFs: {len(tfs)}, target-only: {len(targets_only)}, pairs: {len(trrust_pairs)}")
    print(f"  STRING ≥ 0.7 pairs: {len(string_pairs)}")
    print(f"  gene2go entries: {len(gene2go)}")

    results: dict[str, dict] = {}

    # Skip immune (already done) — focus on kidney
    for tissue in ("kidney",):
        print(f"\n[E1] === Tissue: {tissue} ===")
        emb_path = EMB_SOURCES[tissue]
        h5_path = H5AD_SOURCES[tissue]
        if not emb_path.exists():
            print(f"  Missing embeddings at {emb_path}")
            continue
        if not h5_path.exists():
            print(f"  Missing h5ad at {h5_path}")
            continue
        emb = np.load(emb_path, mmap_mode="r")  # (12, n_genes, 512)
        gene_names = _h5ad_genes(h5_path)
        n_layers, n_genes, d = emb.shape
        print(f"  embedding shape: {emb.shape}; gene names: {len(gene_names)}")

        # SVD layer-by-layer; record SV1 enrichment, SV2 PPI co-pole z, TF/target AUROC on SV2-7
        per_layer: list[dict] = []
        rng = np.random.default_rng(42)
        for li in range(n_layers):
            X = emb[li]
            US, S = svd_gene_projections(X)
            sv1, sv2 = US[:, 0], US[:, 1]
            sv2_7 = US[:, 1:7]

            # SV1 GO:0005615 (top pole) and GO:0005829 (bottom pole)
            or_extra_top, p_extra_top, n_hits_extra, n_term_extra = go_pole_enrichment(sv1, gene_names, gene2go, GO_EXTRACELLULAR, pole="top")
            or_extra_bot, p_extra_bot, _, _ = go_pole_enrichment(sv1, gene_names, gene2go, GO_EXTRACELLULAR, pole="bottom")
            or_cyto_bot, p_cyto_bot, n_hits_cyto, n_term_cyto = go_pole_enrichment(sv1, gene_names, gene2go, GO_CYTOSOL, pole="bottom")

            # Allow either polarity (different tissues might flip SV1 sign)
            or_extra = max(or_extra_top, or_extra_bot)
            p_extra = min(p_extra_top, p_extra_bot)

            # SV2 PPI co-pole
            ppi_z = copole_z(sv2, gene_names, string_pairs, n_perm=200, rng=rng)

            # TF-vs-target AUROC on SV2-7
            tf_auroc = tf_target_auroc(sv2_7, gene_names, tfs, targets_only, n_perm=100, rng=rng)

            per_layer.append(
                dict(
                    layer=li,
                    or_GO0005615=or_extra,
                    p_GO0005615=p_extra,
                    n_hits_extracellular=n_hits_extra,
                    n_term_extracellular_total=n_term_extra,
                    or_GO0005829_bottom=or_cyto_bot,
                    p_GO0005829_bottom=p_cyto_bot,
                    ppi_obs_rate=ppi_z["obs"],
                    ppi_null_rate=ppi_z["null_mean"],
                    ppi_z=ppi_z["z"],
                    ppi_emp_p=ppi_z["emp_p"],
                    ppi_n_pairs=ppi_z["n_pairs"],
                    tf_target_auroc=tf_auroc["auroc"],
                    tf_target_perm_p=tf_auroc["perm_p"],
                    tf_target_n_tf=tf_auroc["n_tf"],
                    tf_target_n_target=tf_auroc["n_target"],
                )
            )

            print(
                f"  L{li:02d}: SV1 extra OR={or_extra:.2f} p={p_extra:.4f}; "
                f"SV2 PPI z={ppi_z['z']:.2f} (obs {ppi_z['obs']:.3f} vs null {ppi_z['null_mean']:.3f}); "
                f"TF-target AUROC={tf_auroc['auroc']:.3f} (n_tf={tf_auroc['n_tf']}, n_target={tf_auroc['n_target']})"
            )

        df = pd.DataFrame(per_layer)
        df.to_csv(OUT / f"per_layer_{tissue}.csv", index=False)

        # Aggregate
        results[tissue] = dict(
            n_genes=int(n_genes),
            n_layers=int(n_layers),
            sv1_extracellular_OR_layer11=float(df.iloc[-1].or_GO0005615),
            sv1_extracellular_p_layer11=float(df.iloc[-1].p_GO0005615),
            sv1_extracellular_OR_max=float(df.or_GO0005615.max()),
            ppi_z_mean=float(df.ppi_z.mean()),
            ppi_z_layers_significant=int((df.ppi_emp_p <= 0.05).sum()),
            tf_target_auroc_mean=float(df.tf_target_auroc.mean()),
            tf_target_auroc_max=float(df.tf_target_auroc.max()),
            tf_target_layers_significant=int((df.tf_target_perm_p <= 0.05).sum()),
            tf_target_n_tf=int(df.tf_target_n_tf.iloc[0]) if len(df) else 0,
            tf_target_n_target=int(df.tf_target_n_target.iloc[0]) if len(df) else 0,
        )
        print(f"  AGG {tissue}: SV1 OR(L11)={results[tissue]['sv1_extracellular_OR_layer11']:.2f}, "
              f"PPI mean z={results[tissue]['ppi_z_mean']:.2f} ({results[tissue]['ppi_z_layers_significant']}/12 sig), "
              f"TF/target AUROC mean={results[tissue]['tf_target_auroc_mean']:.3f}")

    # Save aggregate
    with open(OUT / "results.json", "w") as f:
        json.dump(results, f, indent=2)

    print("\n[E1] Done. Results in results.json and per_layer_*.csv")


if __name__ == "__main__":
    main()
