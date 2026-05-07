#!/usr/bin/env python3
"""E2 — Cross-model replication.

Reviewer 1 (#1): "single model checkpoint."

Compares headline analyses across:
  - scGPT layer 11 (the paper's primary)
  - Geneformer V1-10M
  - Geneformer V2-104M
  - ESM2 protein-LM (cross-modality reference; covered separately by E10)

For each model, we compute on the shared gene set:
  - SV1 GO:0005615 enrichment
  - SV2 STRING ≥ 0.7 PPI co-pole z
  - TF-vs-target AUROC on SV2-7
  - CKA agreement between scGPT-L11 and the model
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

MODELS = {
    "scgpt_L11": REPO / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main/layer_gene_embeddings.npy",
    "geneformer_v1_10m": REPO / "subproject_55_geneformer_v1_10m_gpl_replication/embeddings/geneformer_v1_10m_gene_embeddings.npz",
    "geneformer_v2_104m": REPO / "subproject_56_geneformer_v2_104m_gpl_replication/embeddings/geneformer_v2_104m_gene_embeddings.npz",
}

H5AD = REPO / "single_cell_mechinterp/outputs/tabula_sapiens_processed.h5ad"
TRRUST_PATH = REPO / "single_cell_mechinterp/external/networks/trrust_human.tsv"
STRING_CACHE = REPO / "subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0015/string_ppi_score04_cache.json"
GENE2GO = REPO / "single_cell_mechinterp/data/perturb/gene2go_all.pkl"

GO_EXTRACELLULAR = "GO:0005615"


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


def go_pole_enrichment(proj, gene_names, gene2go, go_term, pole_frac=0.27):
    n = len(gene_names)
    K = max(1, int(pole_frac * n))
    order = np.argsort(proj)
    in_top = set(order[-K:].tolist())
    in_bot = set(order[:K].tolist())
    a_top = sum(1 for i in in_top if go_term in gene2go.get(gene_names[i], set()))
    in_top_with = a_top
    in_top_without = K - a_top
    out_with = sum(1 for i in range(n) if i not in in_top and go_term in gene2go.get(gene_names[i], set()))
    out_without = (n - K) - out_with
    or_top, p_top = fisher_2x2(in_top_with, in_top_without, out_with, out_without)
    a_bot = sum(1 for i in in_bot if go_term in gene2go.get(gene_names[i], set()))
    in_bot_with = a_bot
    in_bot_without = K - a_bot
    out_with2 = sum(1 for i in range(n) if i not in in_bot and go_term in gene2go.get(gene_names[i], set()))
    out_without2 = (n - K) - out_with2
    or_bot, p_bot = fisher_2x2(in_bot_with, in_bot_without, out_with2, out_without2)
    return max(or_top, or_bot), min(p_top, p_bot)


def copole_z(proj, gene_names, pairs, n_perm=100, pole_frac=0.27, rng=None):
    rng = rng or np.random.default_rng(42)
    n2i = {g: i for i, g in enumerate(gene_names)}
    n = len(gene_names)
    K = max(1, int(pole_frac * n))
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
    obs = co / total if total else 0.0
    null_rates = np.zeros(n_perm)
    perm_idx = np.arange(n)
    for k in range(n_perm):
        rng.shuffle(perm_idx)
        n2i_p = {gene_names[perm_idx[i]]: i for i in range(n)}
        co_n = total_n = 0
        for g1, g2 in pairs:
            if g1 not in n2i_p or g2 not in n2i_p:
                continue
            i, j = n2i_p[g1], n2i_p[g2]
            total_n += 1
            if (i in top and j in top) or (i in bot and j in bot):
                co_n += 1
        null_rates[k] = co_n / total_n if total_n else 0.0
    null_mean = null_rates.mean()
    null_std = null_rates.std(ddof=1) if null_rates.std(ddof=1) > 0 else 1e-9
    return dict(obs=obs, null_mean=float(null_mean), z=float((obs - null_mean) / null_std), n_pairs=total)


def tf_target_auroc(proj_subspace, gene_names, tfs, targets_only):
    name_to_i = {g: i for i, g in enumerate(gene_names)}
    tf_idx = [name_to_i[g] for g in tfs if g in name_to_i]
    target_idx = [name_to_i[g] for g in targets_only if g in name_to_i]
    if len(tf_idx) < 5 or len(target_idx) < 5:
        return float("nan"), len(tf_idx), len(target_idx)
    norm_proj = proj_subspace / (np.linalg.norm(proj_subspace, axis=1, keepdims=True) + 1e-12)
    sim = norm_proj @ norm_proj[tf_idx].T
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


def cka(X, Y):
    Xc = X - X.mean(axis=0, keepdims=True)
    Yc = Y - Y.mean(axis=0, keepdims=True)
    K = Xc @ Xc.T
    L = Yc @ Yc.T
    n = K.shape[0]
    H = np.eye(n) - np.ones((n, n)) / n
    Kc = H @ K @ H
    Lc = H @ L @ H
    return float(np.sum(Kc * Lc) / (np.sqrt(np.sum(Kc * Kc) * np.sum(Lc * Lc)) + 1e-12))


def load_model(name, path):
    if name == "scgpt_L11":
        emb = np.load(path, mmap_mode="r")
        scg_l11 = np.asarray(emb[11])
        with h5py.File(H5AD, "r") as f:
            raw = f["var"]["_index"][:]
        gene_names = [g.decode() if isinstance(g, bytes) else str(g) for g in raw]
        return scg_l11, gene_names
    else:
        n = np.load(path, allow_pickle=True)
        return n["embeddings"], list(n["gene_names"])


def main():
    print("[E2] Loading models…", flush=True)
    model_data = {}
    for name, path in MODELS.items():
        emb, names = load_model(name, path)
        print(f"  {name}: shape={emb.shape}, n_genes={len(names)}", flush=True)
        model_data[name] = (emb, names)

    # Annotations
    t = pd.read_csv(TRRUST_PATH, sep="\t", header=None, names=["TF", "target", "mode", "pmid"])
    tfs = set(t.TF.astype(str))
    targets_only = set(t.target.astype(str)) - tfs
    obj = json.loads(STRING_CACHE.read_text())
    string_pairs = [(p["g1"], p["g2"]) for p in obj["pairs"] if p["score"] >= 0.7]
    with open(GENE2GO, "rb") as f:
        g2go = pickle.load(f)

    # Common shared-gene set: intersection of all model gene_names
    common = None
    for _, (emb, names) in model_data.items():
        s = set(names)
        common = s if common is None else common & s
    common = sorted(common)
    print(f"  Shared gene set across all models: {len(common)}", flush=True)

    rng = np.random.default_rng(42)
    rows = []
    sub_embs = {}
    for name, (emb, names) in model_data.items():
        idx = {g: i for i, g in enumerate(names)}
        sub = emb[[idx[g] for g in common]]
        sub_embs[name] = sub
        Xc = sub - sub.mean(axis=0, keepdims=True)
        print(f"[E2] SVD on {name} (sub shape {sub.shape})…", flush=True)
        U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
        US = U * S
        sv1 = US[:, 0]
        sv2 = US[:, 1] if US.shape[1] > 1 else US[:, 0]
        sv27 = US[:, 1:7] if US.shape[1] >= 7 else US[:, 1:]

        or_extra, p_extra = go_pole_enrichment(sv1, common, g2go, GO_EXTRACELLULAR)
        ppi = copole_z(sv2, common, string_pairs, n_perm=100, rng=rng)
        au, n_tf, n_tgt = tf_target_auroc(sv27, common, tfs, targets_only)

        rows.append(
            dict(
                model=name,
                shape=str(sub.shape),
                or_GO0005615=or_extra,
                p_GO0005615=p_extra,
                ppi_obs_rate=ppi["obs"],
                ppi_null_rate=ppi["null_mean"],
                ppi_z=ppi["z"],
                ppi_n_pairs=ppi["n_pairs"],
                tf_target_auroc=au,
                n_tf=n_tf,
                n_target=n_tgt,
            )
        )
        print(f"  {name}: SV1 extra OR={or_extra:.2f} p={p_extra:.4f}; PPI z={ppi['z']:.2f}; TF/target AUROC={au:.3f}", flush=True)

    # CKA: scgpt_L11 vs each Geneformer
    print("[E2] CKA comparisons…", flush=True)
    sg = sub_embs["scgpt_L11"]
    if sg.shape[0] > 3000:
        sub_idx = rng.choice(sg.shape[0], size=3000, replace=False)
    else:
        sub_idx = np.arange(sg.shape[0])
    cka_rows = []
    for name in ("geneformer_v1_10m", "geneformer_v2_104m"):
        cka_v = cka(sg[sub_idx].astype(np.float32), sub_embs[name][sub_idx].astype(np.float32))
        cka_rows.append(dict(comparison=f"scgpt_L11_vs_{name}", cka=cka_v))
        print(f"  {name}: CKA={cka_v:.3f}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "results.csv", index=False)
    cka_df = pd.DataFrame(cka_rows)
    cka_df.to_csv(OUT / "cka.csv", index=False)
    summary = dict(per_model=rows, cka=cka_rows, n_shared_genes=len(common))
    with open(OUT / "results.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("\n[E2] Summary:")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
