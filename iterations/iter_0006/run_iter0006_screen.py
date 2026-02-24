"""
iter_0006 - Multi-Hypothesis Screen

H01: Full-vocabulary residual drift + GO enrichment
     - Compute L2(layer_11 - layer_0) for ALL 4803 scGPT gene positions.
     - Map 209 named genes to their vocab indices.
     - For the named genes, compute their drift PERCENTILE relative to the
       full 4803-gene distribution (context-normalized drift).
     - Run GO enrichment (Fisher exact) comparing top-50 vs bottom-50 drift-ranked
       named genes, as well as top-25% vs bottom-25%.
     - Hypothesis: high-drift named genes enriched in regulatory/immune terms (replicate
       iter_0005 H02 finding, now with full-vocab context normalization).

H02: SVD biology of the layer-11 embedding subspace
     - Compute SVD of the full 4803×512 layer-11 embedding matrix.
     - Project 209 named genes onto top-3 singular vectors.
     - Run GO enrichment for top vs bottom quartile genes on each SV.
     - Also: compute explained variance by top-k components and compare to
       effective rank from iter_0005 (ER=1.28 at layer 11).

H03: Effective-rank curvature and breakpoint analysis
     - Using ER data from iter_0005 (layers 0..11).
     - Compute first and second finite differences.
     - Find the layer of maximum compression rate (min of first diff).
     - Cross-correlate ER curve with TwoNN curve.
     - Novel: compute "compression half-life" = layer at which ER drops to 50% of layer-0.

Data: layer_gene_embeddings.npy [12, 4803, 512], cycle1_edge_dataset.tsv (209 named genes).
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
ER_PATH = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
    "/subproject_41_claude_topology_hypothesis_screening_autoloop"
    "/iterations/iter_0005/h03_effective_rank_per_layer.csv"
)

RNG = np.random.default_rng(42)
N_NULL_REPS = 50  # More reps for better power
MIN_TERM_SIZE = 5
MAX_TERM_SIZE = 200
GO_API_URL = "https://quickgo.ebi.ac.uk/services/annotation/downloadSearch"

# =========================================================
# Data loading
# =========================================================

def load_embeddings():
    emb = np.load(EMB_PATH)  # [12, 4803, 512]
    print(f"Embeddings shape: {emb.shape}", flush=True)
    return emb


def load_gene_index_map():
    """Load the 209 named genes and their embedding indices."""
    gene2idx = {}
    with open(EDGE_PATH) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene2idx[row['source']] = int(row['source_idx'])
            gene2idx[row['target']] = int(row['target_idx'])
    print(f"Loaded {len(gene2idx)} named genes with indices", flush=True)
    return gene2idx


def load_er_data():
    """Load effective rank data from iter_0005."""
    rows = []
    with open(ER_PATH) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({
                'layer': int(row['layer']),
                'er': float(row['effective_rank']),
                'twonn': float(row['twonn_id'])
            })
    return rows


# =========================================================
# GO annotation fetching (QuickGO API)
# =========================================================

def fetch_go_bp_annotations(gene_list, max_per_term=500, retries=3):
    """
    Fetch GO Biological Process annotations for gene list via QuickGO.
    Returns dict: go_term_id -> list of gene symbols in our list.
    """
    gene_set = set(gene_list)
    term2genes = defaultdict(set)

    # Use EBI QuickGO gene annotations endpoint
    genes_str = ",".join(list(gene_set)[:500])  # API limit
    params = {
        "geneProductId": genes_str,
        "aspect": "biological_process",
        "taxonId": "9606",
        "limit": "500",
        "downloadLimit": "10000",
    }
    url = f"https://www.ebi.ac.uk/QuickGO/services/annotation/search?{urllib.parse.urlencode(params)}"
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read())
            if "results" not in data:
                break
            for entry in data["results"]:
                gene_sym = entry.get("geneProductId", "")
                go_id = entry.get("goId", "")
                # geneProductId can be UniProt; skip non-symbol
                if gene_sym in gene_set and go_id.startswith("GO:"):
                    term2genes[go_id].add(gene_sym)
            print(f"  QuickGO: fetched {len(term2genes)} GO terms, attempt {attempt+1}", flush=True)
            break
        except Exception as e:
            print(f"  QuickGO attempt {attempt+1} failed: {e}", flush=True)
            time.sleep(2)

    return term2genes


def fetch_go_via_uniprot_idmapping(gene_list, retries=2):
    """
    Fetch GO Biological Process terms via UniProt ID mapping API (gene symbol -> GO terms).
    Returns gene2go: dict gene_symbol -> set of GO BP term ids.
    """
    gene2go = defaultdict(set)
    # Use UniProt ID mapping: gene names -> UniProt IDs, then fetch GO annotations
    # Step 1: Map gene symbols to UniProt
    url_map = "https://rest.uniprot.org/idmapping/run"
    genes_str = " ".join(gene_list)
    data_bytes = urllib.parse.urlencode({
        "from": "Gene_Name",
        "to": "UniProtKB",
        "ids": genes_str,
        "taxId": "9606"
    }).encode()
    job_id = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url_map, data=data_bytes, method="POST")
            with urllib.request.urlopen(req, timeout=30) as resp:
                result = json.loads(resp.read())
            job_id = result.get("jobId")
            break
        except Exception as e:
            print(f"  UniProt map attempt {attempt+1} failed: {e}", flush=True)
            time.sleep(2)
    if not job_id:
        return gene2go

    # Poll for completion
    for _ in range(15):
        time.sleep(2)
        status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
        try:
            with urllib.request.urlopen(status_url, timeout=15) as resp:
                status = json.loads(resp.read())
            if status.get("jobStatus") in ("RUNNING", None):
                continue
            break
        except:
            continue

    # Fetch results
    results_url = f"https://rest.uniprot.org/idmapping/stream/{job_id}?format=json&fields=go"
    try:
        with urllib.request.urlopen(results_url, timeout=30) as resp:
            results = json.loads(resp.read())
        for item in results.get("results", []):
            gene_sym = item.get("from", "")
            for go_entry in item.get("to", {}).get("uniProtkbId", []):
                pass  # complex parsing; skip
            # Parse Gene Ontology from the to field
            to_field = item.get("to", {})
            if isinstance(to_field, dict):
                for go_ref in to_field.get("geneOntology", []):
                    aspect = go_ref.get("aspect", "")
                    go_id = go_ref.get("id", "")
                    if "biological" in aspect.lower() and go_id:
                        gene2go[gene_sym].add(go_id)
    except Exception as e:
        print(f"  UniProt results fetch failed: {e}", flush=True)

    return gene2go


def fetch_go_from_gene2go_pkl():
    """Try loading pre-cached GO annotations from the single_cell_mechinterp repo."""
    pkl_path = Path(
        "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
        "/single_cell_mechinterp/data/perturb/gene2go_all.pkl"
    )
    if not pkl_path.exists():
        return {}
    import pickle
    with open(pkl_path, 'rb') as f:
        gene2go = pickle.load(f)
    print(f"  Loaded gene2go cache: {len(gene2go)} entries", flush=True)
    return gene2go


# =========================================================
# H01: Full-vocabulary residual drift
# =========================================================

def run_h01_full_vocab_drift(emb, gene2idx):
    """
    Compute drift for all 4803 genes.
    Rank named genes in full-vocab context.
    Run GO enrichment on high vs low drift named genes.
    """
    print("\n=== H01: Full-vocab residual drift ===", flush=True)

    layer0 = emb[0]   # [4803, 512]
    layer11 = emb[11]  # [4803, 512]
    drift_all = np.linalg.norm(layer11 - layer0, axis=1)  # [4803]

    # Full-vocab stats
    drift_mean = float(drift_all.mean())
    drift_std = float(drift_all.std())
    drift_median = float(np.median(drift_all))

    # Named genes only
    named_genes = list(gene2idx.keys())
    named_indices = [gene2idx[g] for g in named_genes]
    named_drifts = drift_all[named_indices]

    # Percentile ranks in full-vocab context
    all_sorted = np.sort(drift_all)
    percentiles = np.searchsorted(all_sorted, named_drifts) / len(all_sorted)

    gene_drift_data = [
        {"gene": g, "emb_idx": gene2idx[g],
         "drift_l2": float(named_drifts[i]),
         "drift_pct_fullvocab": float(percentiles[i])}
        for i, g in enumerate(named_genes)
    ]
    gene_drift_data.sort(key=lambda x: -x["drift_l2"])

    # Save drift data
    drift_out = ITER_DIR / "h01_fullvocab_drift.csv"
    with open(drift_out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['gene','emb_idx','drift_l2','drift_pct_fullvocab'])
        writer.writeheader()
        writer.writerows(gene_drift_data)
    print(f"  Saved {len(gene_drift_data)} named gene drifts to {drift_out}", flush=True)

    # Also save full-vocab drift distribution summary
    pct_vals = [1, 5, 25, 50, 75, 95, 99]
    full_vocab_summary = {
        "n_genes": int(len(drift_all)),
        "mean": drift_mean,
        "std": drift_std,
        "median": drift_median,
        "percentiles": {str(p): float(np.percentile(drift_all, p)) for p in pct_vals},
        "named_gene_mean_drift": float(named_drifts.mean()),
        "named_gene_median_pct": float(np.median(percentiles)),
    }
    print(f"  Full vocab drift: mean={drift_mean:.2f}, std={drift_std:.2f}", flush=True)
    print(f"  Named gene median percentile in full vocab: {np.median(percentiles)*100:.1f}%", flush=True)

    # GO enrichment
    print("  Loading GO annotations...", flush=True)
    gene2go = fetch_go_from_gene2go_pkl()

    # Filter to our named genes
    gene2go_named = {g: gene2go.get(g, set()) for g in named_genes if g in gene2go}
    print(f"  Named genes with GO annotations: {len(gene2go_named)}", flush=True)

    # Build term -> genes mapping
    term2genes_named = defaultdict(set)
    for gene, terms in gene2go_named.items():
        for t in terms:
            term2genes_named[t].add(gene)

    # Filter terms by size
    valid_terms = [t for t, gs in term2genes_named.items()
                   if MIN_TERM_SIZE <= len(gs) <= MAX_TERM_SIZE]
    print(f"  Valid GO BP terms (5-200 named genes): {len(valid_terms)}", flush=True)

    # High drift = top 50 named genes, Low drift = bottom 50
    top50 = set(g["gene"] for g in gene_drift_data[:50])
    bot50 = set(g["gene"] for g in gene_drift_data[-50:])
    # Alternative: top/bottom 25% of named genes
    n25 = max(5, len(named_genes) // 4)
    top25pct = set(g["gene"] for g in gene_drift_data[:n25])
    bot25pct = set(g["gene"] for g in gene_drift_data[-n25:])

    enrichment_rows = []
    for term in valid_terms:
        term_genes = term2genes_named[term]
        # Top50 vs Bot50
        a = len(term_genes & top50)
        b = len(top50) - a
        c = len(term_genes & bot50)
        d = len(bot50) - c
        if a + c == 0:
            continue
        OR50, p50 = fisher_exact([[a, b], [c, d]], alternative='greater')
        enrichment_rows.append({
            "go_term": term,
            "n_term_in_top50": a,
            "n_term_in_bot50": c,
            "OR_top50_vs_bot50": float(OR50),
            "p_top50_vs_bot50": float(p50),
        })

    # Sort by p-value
    enrichment_rows.sort(key=lambda x: x["p_top50_vs_bot50"])

    n_sig_p05 = sum(1 for r in enrichment_rows if r["p_top50_vs_bot50"] < 0.05)
    n_sig_p01 = sum(1 for r in enrichment_rows if r["p_top50_vs_bot50"] < 0.01)
    n_total = len(enrichment_rows)

    # Benjamini-Hochberg FDR
    if n_total > 0:
        pvals = [r["p_top50_vs_bot50"] for r in enrichment_rows]
        bh_thresh = [0.05 * (i+1) / n_total for i in range(n_total)]
        n_sig_fdr = sum(1 for p, t in zip(pvals, bh_thresh) if p <= t)
        for i, r in enumerate(enrichment_rows):
            r["bh_fdr_sig"] = pvals[i] <= bh_thresh[i]
    else:
        n_sig_fdr = 0

    print(f"  Enrichment: {n_total} terms tested, {n_sig_p05} p<0.05, {n_sig_p01} p<0.01, {n_sig_fdr} FDR<0.05", flush=True)
    if enrichment_rows:
        top_hit = enrichment_rows[0]
        print(f"  Top hit: {top_hit['go_term']} OR={top_hit['OR_top50_vs_bot50']:.2f} p={top_hit['p_top50_vs_bot50']:.4f}", flush=True)

    # Save enrichment
    enrich_out = ITER_DIR / "h01_go_enrichment_fullvocab_drift.csv"
    if enrichment_rows:
        with open(enrich_out, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(enrichment_rows[0].keys()))
            writer.writeheader()
            writer.writerows(enrichment_rows)
        print(f"  Saved enrichment to {enrich_out}", flush=True)

    return {
        "full_vocab_summary": full_vocab_summary,
        "n_terms_tested": n_total,
        "n_sig_p05": n_sig_p05,
        "n_sig_p01": n_sig_p01,
        "n_sig_fdr": n_sig_fdr,
        "top_hit_term": enrichment_rows[0]["go_term"] if enrichment_rows else None,
        "top_hit_or": enrichment_rows[0]["OR_top50_vs_bot50"] if enrichment_rows else None,
        "top_hit_p": enrichment_rows[0]["p_top50_vs_bot50"] if enrichment_rows else None,
    }


# =========================================================
# H02: SVD biology of the layer-11 embedding
# =========================================================

def run_h02_svd_biology(emb, gene2idx):
    """
    SVD of the full 4803×512 layer-11 embedding matrix.
    Project named genes onto top singular vectors.
    GO enrichment of top vs bottom quartile on SV1.
    """
    print("\n=== H02: SVD biology of layer-11 embeddings ===", flush=True)

    layer11 = emb[11]  # [4803, 512]
    # Center
    layer11_c = layer11 - layer11.mean(axis=0, keepdims=True)

    print("  Computing SVD...", flush=True)
    U, S, Vt = np.linalg.svd(layer11_c, full_matrices=False)
    # U: [4803, 512], S: [512], Vt: [512, 512]

    explained_var_ratio = (S**2) / (S**2).sum()
    cumvar_top5 = float(explained_var_ratio[:5].sum())
    cumvar_top20 = float(explained_var_ratio[:20].sum())
    sv_top5 = [float(s) for s in S[:5]]

    # Effective rank from SVD (should match iter_0005 ER)
    ent = -np.sum(explained_var_ratio * np.log(explained_var_ratio + 1e-15))
    er_from_svd = float(np.exp(ent))

    print(f"  Top-5 singular values: {[f'{s:.2f}' for s in sv_top5]}", flush=True)
    print(f"  Cumulative var top-5: {cumvar_top5:.3f}, top-20: {cumvar_top20:.3f}", flush=True)
    print(f"  Effective rank (SVD): {er_from_svd:.3f} (iter_0005 ER = 1.279)", flush=True)

    # Named gene projections onto top SVs
    named_genes = list(gene2idx.keys())
    named_indices = [gene2idx[g] for g in named_genes]
    named_U = U[named_indices, :]  # [209, 512]

    sv_projections = []
    for i, g in enumerate(named_genes):
        sv_projections.append({
            "gene": g,
            "sv1_loading": float(named_U[i, 0]),
            "sv2_loading": float(named_U[i, 1]),
            "sv3_loading": float(named_U[i, 2]),
            "sv1_abs": float(abs(named_U[i, 0])),
        })

    sv_projections.sort(key=lambda x: -x["sv1_abs"])

    sv_out = ITER_DIR / "h02_svd_gene_projections.csv"
    with open(sv_out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(sv_projections[0].keys()))
        writer.writeheader()
        writer.writerows(sv_projections)
    print(f"  Saved {len(sv_projections)} gene SV projections to {sv_out}", flush=True)

    # SV1 top/bottom gene summary
    top_sv1 = [r["gene"] for r in sv_projections[:20]]
    bot_sv1 = [r["gene"] for r in sv_projections[-20:]]
    print(f"  Top-20 genes by |SV1| loading: {top_sv1[:10]}...", flush=True)
    print(f"  Bottom-20 genes by |SV1| loading: {bot_sv1[:10]}...", flush=True)

    # GO enrichment: high vs low SV1 magnitude
    print("  Loading GO annotations for SV1 enrichment...", flush=True)
    gene2go = fetch_go_from_gene2go_pkl()
    gene2go_named = {g: gene2go.get(g, set()) for g in named_genes if g in gene2go}
    term2genes_named = defaultdict(set)
    for gene, terms in gene2go_named.items():
        for t in terms:
            term2genes_named[t].add(gene)
    valid_terms = [t for t, gs in term2genes_named.items()
                   if MIN_TERM_SIZE <= len(gs) <= MAX_TERM_SIZE]
    print(f"  Valid GO BP terms: {len(valid_terms)}", flush=True)

    # Sort by SV1 loading (signed, not absolute) for direction test
    sv_projections_signed = sorted(sv_projections, key=lambda x: -x["sv1_loading"])
    n25 = max(5, len(sv_projections_signed) // 4)
    top_q = set(r["gene"] for r in sv_projections_signed[:n25])
    bot_q = set(r["gene"] for r in sv_projections_signed[-n25:])

    sv1_enrich_rows = []
    for term in valid_terms:
        term_genes = term2genes_named[term]
        a = len(term_genes & top_q)
        b = len(top_q) - a
        c = len(term_genes & bot_q)
        d = len(bot_q) - c
        if a + c == 0:
            continue
        OR, p = fisher_exact([[a, b], [c, d]], alternative='greater')
        sv1_enrich_rows.append({
            "go_term": term,
            "n_in_top_q": a,
            "n_in_bot_q": c,
            "OR": float(OR),
            "p": float(p),
        })

    sv1_enrich_rows.sort(key=lambda x: x["p"])
    if sv1_enrich_rows:
        pvals = [r["p"] for r in sv1_enrich_rows]
        n_total = len(sv1_enrich_rows)
        bh_thresh = [0.05 * (i+1) / n_total for i in range(n_total)]
        n_sig_fdr = sum(1 for p, t in zip(pvals, bh_thresh) if p <= t)
        for i, r in enumerate(sv1_enrich_rows):
            r["bh_fdr_sig"] = pvals[i] <= bh_thresh[i]
    else:
        n_sig_fdr = 0

    n_sig_p05_sv1 = sum(1 for r in sv1_enrich_rows if r["p"] < 0.05)
    print(f"  SV1 enrichment: {len(sv1_enrich_rows)} terms, {n_sig_p05_sv1} p<0.05, {n_sig_fdr} FDR<0.05", flush=True)
    if sv1_enrich_rows:
        top_hit = sv1_enrich_rows[0]
        print(f"  SV1 top enrichment: {top_hit['go_term']} OR={top_hit['OR']:.2f} p={top_hit['p']:.4f}", flush=True)

    sv1_enrich_out = ITER_DIR / "h02_svd_sv1_go_enrichment.csv"
    if sv1_enrich_rows:
        with open(sv1_enrich_out, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(sv1_enrich_rows[0].keys()))
            writer.writeheader()
            writer.writerows(sv1_enrich_rows)

    # Also save layer-0 vs layer-11 SVD comparison
    layer0 = emb[0]
    layer0_c = layer0 - layer0.mean(axis=0, keepdims=True)
    _, S0, _ = np.linalg.svd(layer0_c, full_matrices=False)
    ev0 = (S0**2) / (S0**2).sum()
    ent0 = -np.sum(ev0 * np.log(ev0 + 1e-15))
    er0 = float(np.exp(ent0))
    cumvar0_5 = float(ev0[:5].sum())
    cumvar0_20 = float(ev0[:20].sum())

    svd_comparison = [
        {"layer": 0, "er": er0, "cumvar_top5": cumvar0_5, "cumvar_top20": cumvar0_20,
         "sv1": float(S0[0]), "sv2": float(S0[1]), "sv3": float(S0[2])},
        {"layer": 11, "er": er_from_svd, "cumvar_top5": cumvar_top5, "cumvar_top20": cumvar_top20,
         "sv1": float(S[0]), "sv2": float(S[1]), "sv3": float(S[2])},
    ]
    comp_out = ITER_DIR / "h02_svd_layer_comparison.csv"
    with open(comp_out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(svd_comparison[0].keys()))
        writer.writeheader()
        writer.writerows(svd_comparison)

    print(f"  Layer-0 ER={er0:.3f} cumvar5={cumvar0_5:.3f}, Layer-11 ER={er_from_svd:.3f} cumvar5={cumvar_top5:.3f}", flush=True)

    return {
        "layer11_er_from_svd": er_from_svd,
        "layer0_er_from_svd": er0,
        "cumvar_top5_layer11": cumvar_top5,
        "cumvar_top20_layer11": cumvar_top20,
        "sv1_top5_genes": top_sv1[:5],
        "sv1_n_sig_p05": n_sig_p05_sv1,
        "sv1_n_sig_fdr": n_sig_fdr,
        "sv1_top_term": sv1_enrich_rows[0]["go_term"] if sv1_enrich_rows else None,
        "sv1_top_term_or": sv1_enrich_rows[0]["OR"] if sv1_enrich_rows else None,
        "sv1_top_term_p": sv1_enrich_rows[0]["p"] if sv1_enrich_rows else None,
    }


# =========================================================
# H03: Effective rank curvature breakpoint
# =========================================================

def run_h03_er_curvature(er_data):
    """
    Compute ER curvature (first/second finite differences) from iter_0005 data.
    Find the inflection point and compression half-life.
    Cross-correlate ER with TwoNN.
    """
    print("\n=== H03: ER curvature breakpoint ===", flush=True)

    layers = np.array([r['layer'] for r in er_data])
    er = np.array([r['er'] for r in er_data])
    twonn = np.array([r['twonn'] for r in er_data])

    # First differences (rate of ER compression)
    d_er = np.diff(er)  # 11 values
    d_twonn = np.diff(twonn)

    # Second differences (acceleration of compression)
    d2_er = np.diff(d_er)  # 10 values

    # Max compression rate layer
    max_compress_layer = int(layers[1:][np.argmin(d_er)])
    max_compress_rate = float(d_er.min())

    # Inflection point: where d2_er changes sign (or is closest to zero)
    # If all same sign, use argmax of |d2_er|
    sign_changes = np.where(np.diff(np.sign(d2_er)))[0]
    if len(sign_changes) > 0:
        inflection_layer = int(layers[2:][sign_changes[0]])
    else:
        inflection_layer = int(layers[2:][np.argmax(np.abs(d2_er))])

    # Compression half-life: layer at which ER drops to 50% of layer-0
    er0 = er[0]
    half_er = er0 / 2.0
    half_life_layer = None
    for i, e in enumerate(er):
        if e <= half_er:
            half_life_layer = int(layers[i])
            break

    # Cross-correlation ER vs TwoNN
    pearson_r, pearson_p = pearsonr(er, twonn)
    spearman_r, spearman_p = spearmanr(er, twonn)

    # Layer-by-layer ratios
    curvature_data = []
    for i in range(len(layers)):
        row = {
            "layer": int(layers[i]),
            "er": float(er[i]),
            "twonn": float(twonn[i]),
            "er_first_diff": float(d_er[i-1]) if i > 0 else None,
            "er_second_diff": float(d2_er[i-2]) if i > 1 else None,
        }
        curvature_data.append(row)

    curv_out = ITER_DIR / "h03_er_curvature.csv"
    with open(curv_out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['layer','er','twonn','er_first_diff','er_second_diff'])
        writer.writeheader()
        for row in curvature_data:
            writer.writerow({k: ('' if v is None else v) for k, v in row.items()})
    print(f"  Saved curvature data to {curv_out}", flush=True)

    print(f"  ER layer-0: {er0:.3f}, ER layer-11: {er[-1]:.3f}", flush=True)
    print(f"  Max compression rate at layer {max_compress_layer}: {max_compress_rate:.3f}", flush=True)
    print(f"  ER inflection point at layer {inflection_layer}", flush=True)
    print(f"  ER half-life layer: {half_life_layer}", flush=True)
    print(f"  Pearson r(ER,TwoNN) = {pearson_r:.3f}, p = {pearson_p:.4f}", flush=True)
    print(f"  Spearman r(ER,TwoNN) = {spearman_r:.3f}, p = {spearman_p:.4f}", flush=True)

    return {
        "er_layer0": float(er0),
        "er_layer11": float(er[-1]),
        "er_fold_change": float(er0 / er[-1]),
        "max_compress_layer": max_compress_layer,
        "max_compress_rate": max_compress_rate,
        "inflection_layer": inflection_layer,
        "half_life_layer": half_life_layer,
        "pearson_r_er_twonn": float(pearson_r),
        "pearson_p_er_twonn": float(pearson_p),
        "spearman_r_er_twonn": float(spearman_r),
        "spearman_p_er_twonn": float(spearman_p),
    }


# =========================================================
# Main
# =========================================================

def main():
    print("=== iter_0006 Multi-Hypothesis Screen ===", flush=True)
    emb = load_embeddings()
    gene2idx = load_gene_index_map()
    er_data = load_er_data()

    results = {"iteration": "iter_0006"}

    # H01
    try:
        h01 = run_h01_full_vocab_drift(emb, gene2idx)
        results["H01_fullvocab_drift"] = h01
    except Exception as e:
        print(f"H01 failed: {e}", flush=True)
        import traceback
        traceback.print_exc()
        results["H01_fullvocab_drift"] = {"error": str(e)}

    # H02
    try:
        h02 = run_h02_svd_biology(emb, gene2idx)
        results["H02_svd_biology"] = h02
    except Exception as e:
        print(f"H02 failed: {e}", flush=True)
        import traceback
        traceback.print_exc()
        results["H02_svd_biology"] = {"error": str(e)}

    # H03
    try:
        h03 = run_h03_er_curvature(er_data)
        results["H03_er_curvature"] = h03
    except Exception as e:
        print(f"H03 failed: {e}", flush=True)
        import traceback
        traceback.print_exc()
        results["H03_er_curvature"] = {"error": str(e)}

    # Save results
    results_out = ITER_DIR / "iter0006_results.json"
    with open(results_out, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {results_out}", flush=True)
    print(json.dumps(results, indent=2), flush=True)


if __name__ == "__main__":
    main()
