"""
iter_0007 - Multi-Hypothesis Screen

H01: SV1 feature-shuffle null (100 permutations)
H02: Layer-wise SVD trajectory (12 layers)
H03: Drift feature-shuffle null + SV1 negative-pole GO enrichment

Data: layer_gene_embeddings.npy [12, 4803, 512], cycle1_edge_dataset.tsv (209 named genes).
GO annotations from local pkl cache.
"""

import numpy as np
import json
import csv
import pickle
import sys
from pathlib import Path
from scipy.stats import fisher_exact
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
GENE2GO_PKL = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/single_cell_mechinterp/data/perturb/gene2go_all.pkl"
)

RNG = np.random.default_rng(42)
N_SHUFFLE = 100
MIN_TERM_SIZE = 5
MAX_TERM_SIZE = 200


def load_embeddings():
    emb = np.load(EMB_PATH)
    print(f"Embeddings shape: {emb.shape}", flush=True)
    return emb


def load_gene_index_map():
    gene2idx = {}
    with open(EDGE_PATH) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene2idx[row['source'].strip().upper()] = int(row['source_idx'])
            gene2idx[row['target'].strip().upper()] = int(row['target_idx'])
    print(f"Loaded {len(gene2idx)} named genes with indices", flush=True)
    return gene2idx


def load_go_annotations(gene_list):
    """Load GO BP annotations from pkl cache and build term->genes mapping."""
    print(f"Loading GO annotations from pkl cache...", flush=True)
    with open(GENE2GO_PKL, "rb") as f:
        gene2go = pickle.load(f)

    gene_set = set(gene_list)
    term2genes = defaultdict(set)

    for gene, terms in gene2go.items():
        g = gene.strip().upper()
        if g not in gene_set:
            continue
        # terms may be a set/list of GO IDs or dicts
        for t in terms:
            if isinstance(t, dict):
                go_id = t.get("GO_ID", t.get("go_id", ""))
                aspect = t.get("aspect", t.get("Category", "P"))
            else:
                go_id = str(t)
                aspect = "P"  # assume BP
            if go_id.startswith("GO:") and "P" in str(aspect):
                term2genes[go_id].add(g)

    # Filter to terms with enough coverage
    valid = {t: gs for t, gs in term2genes.items()
             if MIN_TERM_SIZE <= len(gs & gene_set) <= MAX_TERM_SIZE}
    print(f"  Total GO terms: {len(term2genes)}, valid (size {MIN_TERM_SIZE}-{MAX_TERM_SIZE}): {len(valid)}", flush=True)
    return valid


def go_enrichment_fisher(top_genes, bottom_genes, term2genes):
    top_set = set(top_genes)
    bot_set = set(bottom_genes)
    results = []
    for term, annotated in term2genes.items():
        a = len(annotated & top_set)
        b = len(annotated & bot_set)
        c = len(top_set) - a
        d = len(bot_set) - b
        if a + b == 0:
            continue
        odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        results.append({"term": term, "a": a, "b": b, "c": c, "d": d,
                         "odds_ratio": odds, "p_value": p})
    results.sort(key=lambda x: x["p_value"])
    return results


# =========================================================
# H01: SV1 feature-shuffle null
# =========================================================

def run_h01_sv1_shuffle_null(emb, gene_list, gene2idx, term2genes):
    print("\n=== H01: SV1 feature-shuffle null (100 perms) ===", flush=True)

    L11 = emb[11]
    L11_c = L11 - L11.mean(axis=0)
    U, S, Vt = np.linalg.svd(L11_c, full_matrices=False)

    named_proj = {g: float(U[gene2idx[g], 0] * S[0]) for g in gene_list}
    sorted_genes = sorted(named_proj, key=named_proj.get)
    n = len(sorted_genes)
    q = n // 4
    top_g = sorted_genes[-q:]
    bot_g = sorted_genes[:q]

    obs_res = go_enrichment_fisher(top_g, bot_g, term2genes)
    obs_top_p = obs_res[0]["p_value"] if obs_res else 1.0
    obs_term = obs_res[0]["term"] if obs_res else "none"
    obs_or = obs_res[0]["odds_ratio"] if obs_res else 1.0
    n_sig_p05 = sum(1 for r in obs_res if r["p_value"] < 0.05)
    print(f"  Observed: p={obs_top_p:.5f} ({obs_term}, OR={obs_or:.2f}), n_sig_p05={n_sig_p05}", flush=True)

    null_ps = []
    for rep in range(N_SHUFFLE):
        perm = RNG.permutation(L11_c.shape[1])
        U_p, S_p, _ = np.linalg.svd(L11_c[:, perm], full_matrices=False)
        pp = {g: float(U_p[gene2idx[g], 0] * S_p[0]) for g in gene_list}
        s_p = sorted(pp, key=pp.get)
        res = go_enrichment_fisher(s_p[-q:], s_p[:q], term2genes)
        null_ps.append(res[0]["p_value"] if res else 1.0)
        if rep % 20 == 0:
            print(f"  Rep {rep}/100: null_p={null_ps[-1]:.4f}", flush=True)

    null_ps = np.array(null_ps)
    null_p5 = float(np.percentile(null_ps, 5))
    pct_rank = float(np.mean(null_ps >= obs_top_p))
    print(f"  Null p5={null_p5:.5f}, obs pct_rank={pct_rank:.3f}, pass={obs_top_p < null_p5}", flush=True)

    np.save(ITER_DIR / "h01_sv1_null_top_ps.npy", null_ps)
    with open(ITER_DIR / "h01_sv1_obs_enrichment.csv", "w", newline="") as f:
        if obs_res:
            w = csv.DictWriter(f, fieldnames=obs_res[0].keys())
            w.writeheader(); w.writerows(obs_res[:50])

    return {
        "obs_top_p": obs_top_p, "obs_term": obs_term, "obs_or": obs_or,
        "obs_n_sig_p05": n_sig_p05,
        "n_mapped": n, "n_shuffle": N_SHUFFLE, "n_go_terms": len(term2genes),
        "null_mean": float(null_ps.mean()), "null_p5": null_p5,
        "null_p25": float(np.percentile(null_ps, 25)),
        "null_p50": float(null_ps.median()) if hasattr(null_ps, 'median') else float(np.median(null_ps)),
        "pct_rank": pct_rank, "pass": bool(obs_top_p < null_p5),
        "obs_top5": obs_res[:5] if obs_res else [],
        "null_ps": null_ps.tolist()
    }, top_g, bot_g


# =========================================================
# H02: Layer-wise SVD trajectory
# =========================================================

def run_h02_layerwise_svd(emb, gene_list, gene2idx, term2genes):
    print("\n=== H02: Layer-wise SVD trajectory ===", flush=True)
    rows = []
    gene_sv1_by_layer = {}

    for layer in range(12):
        L = emb[layer]
        L_c = L - L.mean(axis=0)
        U, S, Vt = np.linalg.svd(L_c, full_matrices=False)

        total_var = float((S**2).sum())
        sv1_var = float(S[0]**2) / total_var
        sv1_sv2_ratio = float(S[0] / S[1]) if S[1] > 0 else float("inf")
        cumvar5 = float(np.cumsum(S**2)[4] / total_var) if len(S) > 4 else float(np.cumsum(S**2)[-1] / total_var)

        p = (S**2) / total_var
        p = p[p > 0]
        er = float(np.exp(-np.sum(p * np.log(p))))

        named_proj = {g: float(U[gene2idx[g], 0] * S[0]) for g in gene_list}
        gene_sv1_by_layer[layer] = named_proj

        sorted_g = sorted(named_proj, key=named_proj.get)
        n = len(sorted_g)
        q = max(5, n // 4)
        res = go_enrichment_fisher(sorted_g[-q:], sorted_g[:q], term2genes)
        top_p = res[0]["p_value"] if res else 1.0
        top_term = res[0]["term"] if res else "none"
        top_or = res[0]["odds_ratio"] if res else 1.0

        row = {"layer": layer, "sv1_var": sv1_var, "sv1_sv2_ratio": sv1_sv2_ratio,
               "cumvar5": cumvar5, "eff_rank": er,
               "top_go_p": top_p, "top_go_term": top_term, "top_go_or": top_or,
               "s1": float(S[0]), "s2": float(S[1])}
        rows.append(row)
        print(f"  L{layer:2d}: SV1/SV2={sv1_sv2_ratio:.2f}, sv1_var={sv1_var:.3f}, "
              f"ER={er:.2f}, GO_p={top_p:.4f} ({top_term})", flush=True)

    with open(ITER_DIR / "h02_layerwise_svd_trajectory.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader(); w.writerows(rows)

    ratios = [r["sv1_sv2_ratio"] for r in rows]
    return {
        "rows": rows,
        "layer_2x": next((i for i, r in enumerate(ratios) if r >= 2.0), None),
        "layer_3x": next((i for i, r in enumerate(ratios) if r >= 3.0), None),
        "layer_5x": next((i for i, r in enumerate(ratios) if r >= 5.0), None),
        "sv1_sv2_l0": ratios[0], "sv1_sv2_l11": ratios[11],
        "sv1_var_l0": rows[0]["sv1_var"], "sv1_var_l11": rows[11]["sv1_var"],
    }, gene_sv1_by_layer


# =========================================================
# H03: Drift shuffle null + SV1 negative-pole
# =========================================================

def run_h03_drift_null_sv1_bottom(emb, gene_list, gene2idx, term2genes, gene_sv1_by_layer):
    print("\n=== H03: Drift shuffle null + SV1 bottom-pole enrichment ===", flush=True)

    L0, L11 = emb[0], emb[11]
    drift = np.linalg.norm(L11 - L0, axis=1)
    named_drift = {g: float(drift[gene2idx[g]]) for g in gene_list}

    sorted_d = sorted(named_drift, key=named_drift.get)
    n = len(sorted_d)
    q50 = n // 2
    top_drift = sorted_d[-q50:]
    bot_drift = sorted_d[:q50]

    obs_res = go_enrichment_fisher(top_drift, bot_drift, term2genes)
    obs_p = obs_res[0]["p_value"] if obs_res else 1.0
    obs_term = obs_res[0]["term"] if obs_res else "none"
    obs_or = obs_res[0]["odds_ratio"] if obs_res else 1.0
    obs_n_sig = sum(1 for r in obs_res if r["p_value"] < 0.05)
    print(f"  Observed drift: p={obs_p:.5f} ({obs_term}, OR={obs_or:.2f}), n_sig={obs_n_sig}", flush=True)

    # Shuffle null
    L0_c = L0 - L0.mean(axis=0)
    L11_c = L11 - L11.mean(axis=0)
    null_ps = []
    for rep in range(N_SHUFFLE):
        perm = RNG.permutation(L0.shape[1])
        drift_p = np.linalg.norm(L11_c[:, perm] - L0_c[:, perm], axis=1)
        named_dp = {g: float(drift_p[gene2idx[g]]) for g in gene_list}
        sorted_dp = sorted(named_dp, key=named_dp.get)
        nq = len(sorted_dp) // 2
        res = go_enrichment_fisher(sorted_dp[-nq:], sorted_dp[:nq], term2genes)
        null_ps.append(res[0]["p_value"] if res else 1.0)
        if rep % 20 == 0:
            print(f"  Drift rep {rep}/100: null_p={null_ps[-1]:.4f}", flush=True)

    null_ps = np.array(null_ps)
    null_p5 = float(np.percentile(null_ps, 5))
    pct_rank = float(np.mean(null_ps >= obs_p))
    print(f"  Drift null p5={null_p5:.5f}, pct_rank={pct_rank:.3f}, pass={obs_p < null_p5}", flush=True)
    np.save(ITER_DIR / "h03_drift_null_ps.npy", null_ps)
    with open(ITER_DIR / "h03_drift_obs_enrichment.csv", "w", newline="") as f:
        if obs_res:
            w = csv.DictWriter(f, fieldnames=obs_res[0].keys())
            w.writeheader(); w.writerows(obs_res[:50])

    # SV1 bottom-pole enrichment (layer 11)
    sv1_l11 = gene_sv1_by_layer.get(11, {})
    sorted_sv1 = sorted(sv1_l11, key=sv1_l11.get)
    n_sv1 = len(sorted_sv1); q = n_sv1 // 4
    sv1_bot = sorted_sv1[:q]
    sv1_top = sorted_sv1[-q:]
    bot_res = go_enrichment_fisher(sv1_bot, sv1_top, term2genes)
    bot_p = bot_res[0]["p_value"] if bot_res else 1.0
    bot_term = bot_res[0]["term"] if bot_res else "none"
    bot_or = bot_res[0]["odds_ratio"] if bot_res else 1.0
    bot_n_sig = sum(1 for r in bot_res if r["p_value"] < 0.05)
    print(f"  SV1 bottom-pole: p={bot_p:.5f} ({bot_term}, OR={bot_or:.2f}), n_sig={bot_n_sig}", flush=True)
    with open(ITER_DIR / "h03_sv1_bottom_pole_enrichment.csv", "w", newline="") as f:
        if bot_res:
            w = csv.DictWriter(f, fieldnames=bot_res[0].keys())
            w.writeheader(); w.writerows(bot_res[:50])

    return {
        "obs_drift_p": obs_p, "obs_drift_term": obs_term, "obs_drift_or": obs_or,
        "obs_drift_n_sig_p05": obs_n_sig,
        "null_drift_mean": float(null_ps.mean()), "null_drift_p5": null_p5,
        "drift_pct_rank": pct_rank, "drift_pass": bool(obs_p < null_p5),
        "sv1_bot_p": bot_p, "sv1_bot_term": bot_term, "sv1_bot_or": bot_or,
        "sv1_bot_n_sig": bot_n_sig,
        "sv1_bot_genes": sv1_bot[:10],
        "obs_drift_top5": obs_res[:5] if obs_res else [],
        "sv1_bot_top5": bot_res[:5] if bot_res else [],
        "null_drift_ps": null_ps.tolist()
    }


# =========================================================
# Main
# =========================================================

def main():
    emb = load_embeddings()
    gene2idx = load_gene_index_map()
    gene_list = list(gene2idx.keys())

    term2genes = load_go_annotations(gene_list)

    if len(term2genes) < 5:
        print("CRITICAL: GO annotations unavailable. Aborting.", flush=True)
        sys.exit(1)

    print(f"Using {len(term2genes)} GO terms for enrichment.", flush=True)

    h01, sv1_top, sv1_bot = run_h01_sv1_shuffle_null(emb, gene_list, gene2idx, term2genes)
    h02, gene_sv1_by_layer = run_h02_layerwise_svd(emb, gene_list, gene2idx, term2genes)
    h03 = run_h03_drift_null_sv1_bottom(emb, gene_list, gene2idx, term2genes, gene_sv1_by_layer)

    results = {
        "iteration": "iter_0007",
        "H01_sv1_shuffle_null": {k: v for k, v in h01.items() if k != "null_ps"},
        "H02_layerwise_svd": {k: v for k, v in h02.items() if k != "rows"},
        "H02_layerwise_svd_rows": h02["rows"],
        "H03_drift_null_sv1_bottom": {k: v for k, v in h03.items() if k not in ("null_drift_ps",)},
    }

    with open(ITER_DIR / "iter0007_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    print("\n=== FINAL SUMMARY ===", flush=True)
    print(f"H01 SV1 shuffle null: obs_p={h01['obs_top_p']:.5f}, null_p5={h01['null_p5']:.5f}, pass={h01['pass']}", flush=True)
    print(f"H02 SVD trajectory: SV1/SV2 L0={h02['sv1_sv2_l0']:.2f}->L11={h02['sv1_sv2_l11']:.2f}, "
          f"2x@{h02['layer_2x']}, 5x@{h02['layer_5x']}", flush=True)
    print(f"H03 drift null: obs_p={h03['obs_drift_p']:.5f}, null_p5={h03['null_drift_p5']:.5f}, pass={h03['drift_pass']}", flush=True)
    print(f"H03 SV1 bottom-pole: p={h03['sv1_bot_p']:.5f} ({h03['sv1_bot_term']}, OR={h03['sv1_bot_or']:.2f})", flush=True)
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
