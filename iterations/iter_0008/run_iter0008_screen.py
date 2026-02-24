"""
iter_0008 - Multi-Hypothesis Screen

H01: SV2/SV3 axes GO enrichment + label-shuffle null at layer-11
H02: Layer-specific compartment transients — layer-3 mitochondrion, layer-7/8 ER lumen vs gene-label null
H03: Signal peptide (UniProt keyword) vs SV1 top-pole enrichment (Fisher test + null)

Data: layer_gene_embeddings.npy [12, 4803, 512], cycle1_edge_dataset.tsv (209 named genes)
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

ITER_DIR = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0008"
)
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
N_SHUFFLE = 500
MIN_TERM = 5
MAX_TERM = 200

# ---------------------------------------------------------------------------
# Load shared data
# ---------------------------------------------------------------------------
print("Loading embeddings...", flush=True)
emb = np.load(EMB_PATH)  # [12, 4803, 512]
print(f"  Shape: {emb.shape}", flush=True)

gene2idx = {}
with open(EDGE_PATH) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        gene2idx[row['source'].strip().upper()] = int(row['source_idx'])
        gene2idx[row['target'].strip().upper()] = int(row['target_idx'])
gene_list = sorted(gene2idx.keys())
gene_arr = np.array(gene_list)
n_genes = len(gene_list)
print(f"  Named genes: {n_genes}", flush=True)

print("Loading GO annotations...", flush=True)
with open(GENE2GO_PKL, "rb") as f:
    gene2go_raw = pickle.load(f)
gene_set = set(gene_list)
term2genes = defaultdict(set)
for gene, terms in gene2go_raw.items():
    g = gene.strip().upper()
    if g not in gene_set:
        continue
    for t in terms:
        if isinstance(t, dict):
            go_id = t.get("GO_ID", t.get("go_id", ""))
            aspect = t.get("aspect", t.get("Category", "P"))
        else:
            go_id = str(t)
            aspect = "P"
        if go_id.startswith("GO:") and "P" in str(aspect):
            term2genes[go_id].add(g)
term2genes_bp = {t: gs for t, gs in term2genes.items() if MIN_TERM <= len(gs & gene_set) <= MAX_TERM}
print(f"  GO BP terms (size {MIN_TERM}-{MAX_TERM}): {len(term2genes_bp)}", flush=True)

# Also build CC (cellular component) annotations for subcellular compartment tests
term2genes_cc = defaultdict(set)
for gene, terms in gene2go_raw.items():
    g = gene.strip().upper()
    if g not in gene_set:
        continue
    for t in terms:
        if isinstance(t, dict):
            go_id = t.get("GO_ID", t.get("go_id", ""))
            aspect = t.get("aspect", t.get("Category", "C"))
        else:
            go_id = str(t)
            aspect = "C"
        if go_id.startswith("GO:") and "C" in str(aspect):
            term2genes_cc[go_id].add(g)
term2genes_cc = {t: gs for t, gs in term2genes_cc.items() if MIN_TERM <= len(gs & gene_set) <= MAX_TERM}
print(f"  GO CC terms (size {MIN_TERM}-{MAX_TERM}): {len(term2genes_cc)}", flush=True)

# Combine BP+CC for broader search
term2genes_all = {**term2genes_bp, **term2genes_cc}


def go_enrich_top(top_genes, bg_genes, t2g):
    """Enrichment of top_genes vs bg_genes (rest)."""
    top_s = set(top_genes)
    bg_s = set(bg_genes) - top_s
    results = []
    for term, ann in t2g.items():
        a = len(ann & top_s)
        b = len(ann & bg_s)
        c = len(top_s) - a
        d = len(bg_s) - b
        if a + b == 0:
            continue
        odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        results.append({"term": term, "a": a, "b": b, "c": c, "d": d, "odds_ratio": odds, "p_value": p})
    results.sort(key=lambda x: x["p_value"])
    return results


def label_shuffle_null(sv_vals, gene_arr_local, t2g, n_shuffle, rng, q_frac=0.25):
    """Gene-label shuffle null: shuffle gene labels, keep SV positions fixed."""
    n = len(gene_arr_local)
    q = int(n * q_frac)
    sorted_pos = np.argsort(sv_vals)
    null_ps = []
    for _ in range(n_shuffle):
        perm = rng.permutation(n)
        shuf = gene_arr_local[perm]
        top_g = shuf[sorted_pos[-q:]].tolist()
        bg_g = shuf.tolist()
        res = go_enrich_top(top_g, bg_g, t2g)
        null_ps.append(res[0]["p_value"] if res else 1.0)
    return np.array(null_ps)


# ===========================================================================
# H01: SV2 and SV3 axes at layer-11, GO enrichment + gene-label shuffle null
# ===========================================================================
print("\n=== H01: SV2/SV3 axes at layer-11 ===", flush=True)
L11 = emb[11]
L11_c = L11 - L11.mean(axis=0)
U, S, Vt = np.linalg.svd(L11_c, full_matrices=False)

# Extract SV1/SV2/SV3 gene projections
sv_projs = {}
for sv_idx, sv_name in [(0, "sv1"), (1, "sv2"), (2, "sv3")]:
    vals = np.array([float(U[gene2idx[g], sv_idx] * S[sv_idx]) for g in gene_list])
    sv_projs[sv_name] = vals

q = n_genes // 4

h01_results = {}
for sv_name in ["sv1", "sv2", "sv3"]:
    vals = sv_projs[sv_name]
    sorted_pos = np.argsort(vals)
    top_genes = gene_arr[sorted_pos[-q:]].tolist()
    bg_genes = gene_arr.tolist()

    obs_res = go_enrich_top(top_genes, bg_genes, term2genes_all)
    obs_p = obs_res[0]["p_value"] if obs_res else 1.0
    obs_term = obs_res[0]["term"] if obs_res else "none"
    obs_or = obs_res[0]["odds_ratio"] if obs_res else float("nan")
    print(f"  {sv_name} top: p={obs_p:.5f} term={obs_term} OR={obs_or:.2f}", flush=True)

    # Also check bottom pole
    bot_genes = gene_arr[sorted_pos[:q]].tolist()
    bot_res = go_enrich_top(bot_genes, bg_genes, term2genes_all)
    bot_p = bot_res[0]["p_value"] if bot_res else 1.0
    bot_term = bot_res[0]["term"] if bot_res else "none"
    bot_or = bot_res[0]["odds_ratio"] if bot_res else float("nan")
    print(f"  {sv_name} bot: p={bot_p:.5f} term={bot_term} OR={bot_or:.2f}", flush=True)

    h01_results[sv_name] = {
        "top_p": obs_p, "top_term": obs_term, "top_or": float(obs_or),
        "top_top5": obs_res[:5] if obs_res else [],
        "bot_p": bot_p, "bot_term": bot_term, "bot_or": float(bot_or),
        "bot_top5": bot_res[:5] if bot_res else [],
    }

# Gene-label shuffle null for SV2 and SV3 (top pole)
print("  Running label-shuffle null for SV2 and SV3 top poles...", flush=True)
h01_null = {}
for sv_name in ["sv2", "sv3"]:
    vals = sv_projs[sv_name]
    null_ps = label_shuffle_null(vals, gene_arr, term2genes_all, N_SHUFFLE, RNG)
    emp_p = float((null_ps <= h01_results[sv_name]["top_p"]).mean())
    h01_null[sv_name] = {
        "null_mean": float(null_ps.mean()),
        "null_p5": float(np.percentile(null_ps, 5)),
        "empirical_p": emp_p,
        "obs_p": h01_results[sv_name]["top_p"],
    }
    np.save(ITER_DIR / f"h01_{sv_name}_null_ps.npy", null_ps)
    print(f"  {sv_name}: null_mean={h01_null[sv_name]['null_mean']:.5f}, empirical_p={emp_p:.3f}", flush=True)

# Save H01 enrichment tables
for sv_name in ["sv1", "sv2", "sv3"]:
    rows = h01_results[sv_name]["top_top5"] + h01_results[sv_name]["bot_top5"]
    fpath = ITER_DIR / f"h01_{sv_name}_enrichment.csv"
    if rows:
        with open(fpath, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            for r in rows:
                # make OR serializable
                r2 = dict(r)
                if r2.get("odds_ratio") == float("inf"):
                    r2["odds_ratio"] = "inf"
                writer.writerow(r2)


# ===========================================================================
# H02: Layer-specific compartment transients — mitochondrion (layer-3), ER lumen (layer-7/8)
# Target terms: GO:0005739 (mitochondrion), GO:0005788 (ER lumen)
# ===========================================================================
print("\n=== H02: Layer-specific compartment transients ===", flush=True)

TARGET_TERMS = {
    "mito_l3": {"layer": 3, "term": "GO:0005739", "label": "mitochondrion"},
    "er_lumen_l7": {"layer": 7, "term": "GO:0005788", "label": "ER lumen"},
    "er_lumen_l8": {"layer": 8, "term": "GO:0005788", "label": "ER lumen"},
}

# Combine BP+CC (already done above in term2genes_all)
# Also need raw CC-only for these specific terms
term2genes_cc_raw = defaultdict(set)
for gene, terms in gene2go_raw.items():
    g = gene.strip().upper()
    if g not in gene_set:
        continue
    for t in terms:
        if isinstance(t, dict):
            go_id = t.get("GO_ID", t.get("go_id", ""))
            # any aspect OK — we want CC terms
        else:
            go_id = str(t)
        if go_id.startswith("GO:"):
            term2genes_cc_raw[go_id].add(g)

h02_results = {}

for key, cfg in TARGET_TERMS.items():
    layer_idx = cfg["layer"]
    target_term = cfg["term"]
    layer_label = cfg["label"]

    L = emb[layer_idx]
    Lc = L - L.mean(axis=0)
    U_l, S_l, Vt_l = np.linalg.svd(Lc, full_matrices=False)
    sv1_vals = np.array([float(U_l[gene2idx[g], 0] * S_l[0]) for g in gene_list])

    sorted_pos = np.argsort(sv1_vals)
    q = n_genes // 4
    top_genes = gene_arr[sorted_pos[-q:]].tolist()
    bg_genes = gene_arr.tolist()

    # Check if target term has annotated genes in our set
    ann_genes = term2genes_cc_raw.get(target_term, set()) & gene_set
    print(f"  {key} layer-{layer_idx} {target_term}: {len(ann_genes)} annotated genes in set", flush=True)

    if len(ann_genes) < 3:
        print(f"  Skipping {key}: too few annotated genes", flush=True)
        h02_results[key] = {"status": "skipped", "reason": "too_few_annotated", "n_annotated": len(ann_genes)}
        continue

    # Fisher test: target term enriched in top-25% SV1?
    a = len(set(top_genes) & ann_genes)
    b = len(ann_genes) - a
    c = len(top_genes) - a
    d = n_genes - len(top_genes) - b
    odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    print(f"  {key}: a={a} b={b} c={c} d={d} OR={odds:.2f} p={p:.5f}", flush=True)

    # Gene-label shuffle null for this specific test
    null_ps = []
    for rep in range(N_SHUFFLE):
        perm = RNG.permutation(n_genes)
        shuf = gene_arr[perm]
        top_shuf = set(shuf[sorted_pos[-q:]].tolist())
        aa = len(top_shuf & ann_genes)
        bb = len(ann_genes) - aa
        cc = q - aa
        dd = n_genes - q - bb
        _, pp = fisher_exact([[aa, bb], [cc, dd]], alternative="greater")
        null_ps.append(pp)
    null_arr = np.array(null_ps)
    emp_p = float((null_arr <= p).mean())
    np.save(ITER_DIR / f"h02_{key}_null_ps.npy", null_arr)

    print(f"  {key}: null_mean={null_arr.mean():.4f} emp_p={emp_p:.3f}", flush=True)
    h02_results[key] = {
        "layer": layer_idx, "term": target_term, "label": layer_label,
        "n_annotated": len(ann_genes), "a": a, "b": b, "c": c, "d": d,
        "odds_ratio": float(odds), "p_value": float(p),
        "null_mean": float(null_arr.mean()), "null_p5": float(np.percentile(null_arr, 5)),
        "empirical_p": emp_p,
        "pass": emp_p < 0.05,
    }

# Save H02 results table
h02_rows = []
for key, res in h02_results.items():
    if "odds_ratio" in res:
        h02_rows.append({"key": key, **res})
if h02_rows:
    with open(ITER_DIR / "h02_layer_compartment_transients.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=h02_rows[0].keys())
        writer.writeheader()
        for r in h02_rows:
            writer.writerow(r)


# ===========================================================================
# H03: Signal peptide genes vs SV1 top-pole (layer-11)
# Use curated signal-peptide gene list from literature + UniProt keywords overlap
# Proxy: genes annotated to GO:0005615 (extracellular space) + GO:0005576 (extracellular region)
#        as signal-peptide proxies (proteins secreted need signal peptide)
# ===========================================================================
print("\n=== H03: Signal peptide proxy vs SV1 top-pole ===", flush=True)

# Signal peptide proxy: genes in extracellular space (GO:0005615) OR extracellular region (GO:0005576)
SIGNAL_PEPTIDE_TERMS = {"GO:0005615", "GO:0005576", "GO:0005788"}  # extracell space, region, ER lumen
signal_genes = set()
for t in SIGNAL_PEPTIDE_TERMS:
    signal_genes |= (term2genes_cc_raw.get(t, set()) & gene_set)
# Also check term2genes_all
for t in SIGNAL_PEPTIDE_TERMS:
    signal_genes |= (term2genes_all.get(t, set()) & gene_set)
print(f"  Signal-peptide proxy genes: {len(signal_genes)}", flush=True)

# SV1 at layer-11
sv1_vals = sv_projs["sv1"]
sorted_pos = np.argsort(sv1_vals)
q = n_genes // 4
top_genes_sv1 = set(gene_arr[sorted_pos[-q:]].tolist())
bg_rest = gene_set - top_genes_sv1

a = len(signal_genes & top_genes_sv1)
b = len(signal_genes - top_genes_sv1)
c = len(top_genes_sv1) - a
d = len(bg_rest) - b
odds_sp, p_sp = fisher_exact([[a, b], [c, d]], alternative="greater")
print(f"  Signal peptide: a={a} b={b} OR={odds_sp:.2f} p={p_sp:.5f}", flush=True)

# Null for H03
null_ps_sp = []
for rep in range(N_SHUFFLE):
    perm = RNG.permutation(n_genes)
    shuf = gene_arr[perm]
    top_shuf = set(shuf[sorted_pos[-q:]].tolist())
    aa = len(signal_genes & top_shuf)
    bb = len(signal_genes - top_shuf)
    cc = q - aa
    dd = n_genes - q - bb
    _, pp = fisher_exact([[aa, bb], [cc, dd]], alternative="greater")
    null_ps_sp.append(pp)
null_sp_arr = np.array(null_ps_sp)
emp_p_sp = float((null_sp_arr <= p_sp).mean())
np.save(ITER_DIR / "h03_signal_peptide_null_ps.npy", null_sp_arr)

print(f"  Signal peptide null: mean={null_sp_arr.mean():.5f} emp_p={emp_p_sp:.3f}", flush=True)

h03_result = {
    "n_signal_genes": len(signal_genes),
    "proxy_terms": sorted(SIGNAL_PEPTIDE_TERMS),
    "a": a, "b": b, "c": c, "d": d,
    "odds_ratio": float(odds_sp),
    "p_value": float(p_sp),
    "null_mean": float(null_sp_arr.mean()),
    "null_p5": float(np.percentile(null_sp_arr, 5)),
    "empirical_p": emp_p_sp,
    "pass": emp_p_sp < 0.05,
    "signal_genes_in_top": sorted(signal_genes & top_genes_sv1),
}

# Also compute SV1 score for signal vs non-signal genes
signal_sv1_scores = [sv1_vals[i] for i, g in enumerate(gene_list) if g in signal_genes]
nonsig_sv1_scores = [sv1_vals[i] for i, g in enumerate(gene_list) if g not in signal_genes]
from scipy.stats import mannwhitneyu
stat, p_mw = mannwhitneyu(signal_sv1_scores, nonsig_sv1_scores, alternative="greater")
h03_result["signal_median_sv1"] = float(np.median(signal_sv1_scores))
h03_result["nonsig_median_sv1"] = float(np.median(nonsig_sv1_scores))
h03_result["mw_stat"] = float(stat)
h03_result["mw_p"] = float(p_mw)
print(f"  MW test signal vs non-signal SV1: stat={stat:.1f} p={p_mw:.5f}", flush=True)
print(f"  Median SV1: signal={np.median(signal_sv1_scores):.3f} non-signal={np.median(nonsig_sv1_scores):.3f}", flush=True)

# ===========================================================================
# Save all results
# ===========================================================================
results = {
    "iteration": "iter_0008",
    "H01_sv2_sv3_axes": h01_results,
    "H01_null": h01_null,
    "H02_layer_transients": h02_results,
    "H03_signal_peptide": h03_result,
}

with open(ITER_DIR / "iter0008_results.json", "w") as f:
    import math
    def default_serial(obj):
        if isinstance(obj, (np.integer,)): return int(obj)
        if isinstance(obj, (np.floating,)):
            if math.isinf(obj): return "inf"
            if math.isnan(obj): return None
            return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        if isinstance(obj, set): return sorted(obj)
        raise TypeError(f"Not serializable: {type(obj)}")
    json.dump(results, f, indent=2, default=default_serial)

print("\n=== Done. Results saved to iter0008_results.json ===", flush=True)
