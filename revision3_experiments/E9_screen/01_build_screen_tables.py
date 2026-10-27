"""E9 / A6 steps 1-3 -- per-record provenance table (S2), per-family table (S1), context usage.

Sources (read-only): every iterations/iter_*/executor_hypothesis_screen.json; the .py scripts in each
iteration folder; executor_iteration_report.md (iterations 45-63 left no script, so the data paths written in
their reports and the record's own method text are used); round-2 audit (revision2_experiments/
master_p_round2.csv); this folder's audit_409_records.csv (from 02_audit_409.py).

Every cell that could not be determined programmatically or from the files is 'unknown'. Where a value is a
hand-assigned reading of the iteration's own files (iterations 1-63, which underlie the preprint), the column
`provenance_basis` says so.
"""
from __future__ import annotations
import json, re
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
P41 = HERE.parents[1]
ITER = P41 / "iterations"
R2 = P41 / "revision2_experiments"

# ------------------------------------------------------------------ context vocabulary
CTX_RE = [  # (regex on script text, context name)
    (r"cycle1_main_v2env_seed42", "kidney_A_v2env_seed42"), (r"cycle2_maxgenes1024", "kidney_A_maxgenes1024"),
    (r"cycle1_main(?!_v2env)", "kidney_A"), (r"cycle1_seed43", "kidney_A_seed43"), (r"cycle1_seed44", "kidney_A_seed44"),
    (r"cycle4_immune_main", "immune"), (r"cycle4_immune_seed43", "immune_seed43"), (r"cycle4_immune_seed44", "immune_seed44"),
    (r"cycle6_lung_main", "lung"), (r"cycle6_lung_seed43", "lung_seed43"), (r"cycle6_lung_seed44", "lung_seed44"),
    (r"cycle7_external_lung_main", "lung_external"), (r"cycle7_external_lung_seed43", "lung_external_seed43"),
    (r"cycle7_external_lung_seed44", "lung_external_seed44"),
    (r"cycle10_immune_donor_TSP14", "immune_TSP14"), (r"cycle10_immune_donor_TSP21", "immune_TSP21"),
    (r"cycle10_immune_donor_TSP25", "immune_TSP25"),
    (r"cycle11_immune_celltype_b_cell", "immune_Bcell"), (r"cycle11_immune_celltype_cd4_t", "immune_CD4T"),
    (r"cycle11_immune_celltype_cd8_t", "immune_CD8T"),
    (r"kidney_scgpt_whole_human_run", "kidney_ckpt_whole_human"), (r"kidney_scgpt_kidney_run", "kidney_ckpt_kidney"),
    (r"kidney_scgpt_brain_run", "kidney_ckpt_brain"), (r"kidney_scgpt_run(?!_)", "kidney_B"),
    (r"cycle12_geneformer_immune_bootstrap", "geneformer_immune"), (r"cycle12_geneformer_lung_bootstrap", "geneformer_lung"),
    (r"cycle12_geneformer_external_lung_bootstrap", "geneformer_lung_external"),
    (r"cycle12_geneformer_kidney_bootstrap", "geneformer_kidney"), (r"cycle12_geneformer_liver_bootstrap", "geneformer_liver"),
    (r"cycle12_geneformer_heart_bootstrap", "geneformer_heart"),
]
TISSUE = {"kidney_A": "kidney (TS kidney; labelled immune in the preprint)", "kidney_A_seed43": "kidney", "kidney_A_seed44": "kidney",
          "kidney_A_v2env_seed42": "kidney", "kidney_A_maxgenes1024": "kidney", "kidney_B": "kidney",
          "kidney_ckpt_whole_human": "kidney", "kidney_ckpt_kidney": "kidney", "kidney_ckpt_brain": "kidney",
          "immune": "immune", "immune_seed43": "immune", "immune_seed44": "immune", "immune_TSP14": "immune",
          "immune_TSP21": "immune", "immune_TSP25": "immune", "immune_Bcell": "immune (B cells)",
          "immune_CD4T": "immune (CD4 T)", "immune_CD8T": "immune (CD8 T)", "lung": "lung", "lung_seed43": "lung",
          "lung_seed44": "lung", "lung_external": "lung (Krasnow atlas)", "lung_external_seed43": "lung (Krasnow atlas)",
          "lung_external_seed44": "lung (Krasnow atlas)", "lung_attention": "lung (scGPT attention matrix)",
          "random_gaussian": "none (synthetic null)", "geneformer_static": "none (Geneformer token embeddings)"}
AXIS = {"kidney_A": 4803, "kidney_A_seed43": 4803, "kidney_A_seed44": 4803, "kidney_A_v2env_seed42": 4803,
        "kidney_A_maxgenes1024": 4803, "kidney_B": 8143, "kidney_ckpt_whole_human": 8143, "kidney_ckpt_kidney": 8143,
        "kidney_ckpt_brain": 8143, "lung": 8181, "lung_seed43": 8181, "lung_seed44": 8181, "lung_external": 8229,
        "lung_external_seed43": 8229, "lung_external_seed44": 8229, "lung_attention": 8181, "random_gaussian": 4803,
        "geneformer_static": 20275}
for k in ("immune", "immune_seed43", "immune_seed44", "immune_TSP14", "immune_TSP21", "immune_TSP25", "immune_Bcell",
          "immune_CD4T", "immune_CD8T"):
    AXIS[k] = 4941
for t in ("immune", "lung", "lung_external", "kidney", "liver", "heart"):
    TISSUE[f"geneformer_{t}"] = f"{t} (Geneformer, not scGPT)"; AXIS[f"geneformer_{t}"] = "unknown"


def script_contexts(it_dir):
    txt = "\n".join(p.read_text(errors="ignore") for p in it_dir.glob("*.py") if not p.name.startswith("._"))
    found = [c for rx, c in CTX_RE if re.search(rx, txt)]
    if "attention_scores" in txt and "invariant_causal_edges" in txt: found.append("lung_attention")
    return sorted(set(found)), txt


# ------------------------------------------------------------------ manual reading, iterations 1-63
SEEDS_I = "immune;immune_seed43;immune_seed44"
def ctx_1_63(it, hid, text):
    """Matrix actually loaded, iterations 1-63 (from each iteration's script; reports for 45-63)."""
    if it == 2: return "kidney_A;kidney_A_seed43;kidney_A_seed44", "script"
    if 4 <= it <= 38:
        if (it, hid) in {(18, "H03"), (19, "H03"), (20, "H02")}: return "kidney_A;lung_attention", "script"
        if (it, hid) == (19, "H01"): return "kidney_A;geneformer_static", "script"
        if (it, hid) == (18, "H01"): return "random_gaussian", "script (Gaussian replacement of the kidney_A matrix)"
        return "kidney_A", "script"
    if it == 39: return ("immune" if hid == "H02" else "kidney_A"), "script + record text"
    if it == 40: return ("geneformer_static" if hid == "H02" else "immune"), "script + record text"
    if 41 <= it <= 44: return "immune", "script"
    if 45 <= it <= 58:
        if (it, hid) in {(46, "H01"), (55, "H01"), (56, "H02"), (57, "H01"), (57, "H02"), (58, "H01"), (58, "H03")}:
            return SEEDS_I, "report + record text (no script)"
        return "immune", "report + record text (no script)"
    if it == 59:
        if hid == "H01": return "kidney_A;kidney_A_seed43;kidney_A_seed44", "report (no script); summed embeddings"
        if hid == "H03": return "immune;immune_seed43;immune_seed44", "re-analysis of iter_0058 output (no script)"
        return "kidney_A", "report (no script); summed embeddings"
    if it == 60: return "kidney_A", "report (no script); summed embeddings"
    if it == 61: return SEEDS_I, "report + record text (no script)"
    if it == 62:
        if hid == "H01": return SEEDS_I, "report + record text (no script)"
        if hid == "H03": return "kidney_A;immune", "report + record text (no script)"
        return "immune", "report + record text (no script)"
    if it == 63: return "immune", "report + record text (no script)"
    return "unknown", "unknown"


def mapping_1_63(it, hid, ctx):
    """Gene->row mapping used for kidney_A labels."""
    if "kidney_A" not in ctx: return "n/a (immune/other: h5ad var_names or edge-dataset indices)" if ctx != "unknown" else "unknown"
    if it in (2,) or (it == 4 and hid in ("H01", "H03")): return "n/a (unlabelled random gene sample)"
    if it == 27 or 29 <= it <= 40: return "MISALIGNED (gene_list.txt blank lines dropped: labels on rows 0-208)"
    if 4 <= it <= 26 or it == 28: return "true (edge-dataset indices or blank-preserving gene_list.txt)"
    if it in (59, 60, 62): return "true (edge-dataset source_idx/target_idx, per report)"
    return "unknown"


def geneset_1_63(it, hid, ctx, text, mapping):
    t = text
    if it == 2: return "300 random genes per seed (4,803 axis; PCA-20)"
    if it == 4 and hid in ("H01", "H03"): return "400 random genes (4,803 axis; PCA-20)"
    if it == 6 and hid == "H01": return "full 4,803 (drift) + 209 named genes for GO"
    if it in (6, 7) and ("4803" in t or "full" in t.lower()): return "full 4,803 matrix (2,831 all-zero rows); 209 named genes projected"
    if ctx.startswith("random_gaussian"): return "random Gaussian [12,4803,512]; 209 named-gene positions"
    if "MISALIGNED" in mapping:
        return ("209 labels on rows 0-208 (misaligned)" if it <= 30 else
                "'195 in-vocab' = rows 0-208 minus the 14 zero rows among them (misaligned labels)")
    if it == 61 and hid == "H01": return "1,496 genes non-zero in all 3 immune seeds"
    if it == 62 and hid == "H03": return "80 genes shared by the kidney_A and immune edge sets"
    if "kidney_A" in ctx and it in (59, 60): return "1,972 observed rows (summed embeddings); labels on the 209-gene edge set (288 pos / 864 neg rows)"
    if "kidney_A" in ctx and "true" in mapping:
        return "209 named genes, true rows (22 of 209 are never-observed all-zero rows)"
    if ctx.startswith("immune") or ctx.startswith("geneformer"):
        if it in (45, 46) and ("2000" in t or "n=2000" in t): return "2,000 random of 4,941 (zero rows included)"
        if "2039" in t or "nonzero" in t.lower() or it >= 47: return "2,039 observed; labels on the 295 observed of 361 edge-dataset genes"
        return "4,941 axis incl. 2,902 zero rows (ranks/centroids over all rows); marker panels"
    return "unknown"


SPECTRAL_RE = re.compile(r"\bSV\s?\d|\bPC\s?\d|\bSVD\b|singular|spectral|effective rank|eff_rank|participation ratio|\bPR\b|PCA", re.I)
def svdconv_1_63(it, hid, ctx, text, mapping):
    t = text
    if not SPECTRAL_RE.search(t): return "not spectral"
    if it in (2, 4) and "PCA" in t: return "not spectral (PCA-20 preprocessing of sampled genes)"
    if "PCA" in t and not re.search(r"SV\d|SVD|singular|PC\d", t): return "not spectral (PCA preprocessing)"
    if ctx.startswith("random_gaussian"): return "subset-fitted (209 rows of a Gaussian matrix)"
    if it in (6, 7) and ("4803" in t or "full" in t.lower()): return "full matrix with zero rows (4,803 rows, 2,831 all-zero)"
    if "MISALIGNED" in mapping: return "subset-fitted on misaligned rows (rows 0-208; 195 non-zero after iter_0030)"
    if "kidney_A" in ctx and "true" in mapping and it <= 28: return "subset-fitted (209 true rows, incl. 22 all-zero rows)"
    if "kidney_A" in ctx and it in (59, 60): return "observed non-zero rows only (1,972; summed, not mean, embeddings)"
    if (it, hid) in {(46, "H02")}: return "random 2,000-row subsample incl. zero rows (centred over the subsample)"
    if (it, hid) in {(46, "H01"), (45, "H02")}: return "not spectral (TwoNN on 2,000-row subsample incl. zero rows)"
    if ctx.startswith("immune") and ("nonzero" in t.lower() or "2039" in t or it >= 47): return "observed non-zero rows only (2,039)"
    if ctx.startswith("immune"): return "unknown (spectral wording; rows not stated)"
    return "unknown"


CLAIMS = {  # (iteration, id) -> ledger claim ids (E9_claims_ledger.csv)
    (2, "H01_graph_topology_cc"): "C57", (6, "H02"): "C08", (6, "H03"): "C61", (7, "H01"): "C10", (7, "H02"): "C03(related: kidney_A full-matrix 57.2->93.2%)",
    (7, "H03"): "C09", (8, "H01"): "C12", (8, "H02"): "C11", (9, "H01"): "C12", (9, "H02"): "C11", (9, "H03"): "C04;C16",
    (10, "H01"): "C12", (10, "H02"): "C13;C60", (10, "H03"): "C04;C16", (11, "H01"): "C13;C14", (11, "H03"): "C53",
    (12, "H03"): "C21", (13, "H01"): "C21", (14, "H01"): "C15", (14, "H03"): "C25", (15, "H01"): "C22", (15, "H02"): "C23",
    (15, "H03"): "C24", (16, "H02"): "C26", (18, "H03"): "C56", (19, "H03"): "C56(related)", (19, "H01"): "C50",
    (22, "H01"): "C27", (22, "H02"): "C18(initial 0.853)", (23, "H01"): "C18;C19", (26, "H01"): "C20", (30, "H04"): "C59",
    (31, "H03"): "C06", (34, "H03"): "C40", (35, "H02"): "C39(precursor)", (36, "H01"): "C39(precursor)", (36, "H02"): "C39(precursor)",
    (37, "H01"): "C39", (38, "H02"): "C44", (40, "H01"): "C42", (40, "H02"): "C51", (41, "H01"): "C47", (41, "H02"): "C42;C45",
    (41, "H03"): "C46", (42, "H02"): "C45", (42, "H03"): "C41", (43, "H01"): "C43", (45, "H02"): "C05", (46, "H01"): "C05",
    (46, "H02"): "C02;C03;C07", (47, "H01"): "contradicts C02 (not reported)", (47, "H02"): "contradicts C03 (not reported)",
    (51, "H01"): "C31", (51, "H02"): "C17", (52, "H03"): "C32", (53, "H01"): "C32", (55, "H02"): "C30",
    (56, "H01"): "C28;C29;C30", (57, "H01"): "C29", (58, "H03"): "C52", (59, "H03"): "C52", (61, "H01"): "C55",
    (61, "H03"): "C54", (62, "H01"): "C33;C34;C35", (62, "H02"): "C34;C36", (63, "H63-A"): "C38", (63, "H63-G"): "C37;C58",
}

PREPRINT_FAMILIES = ["Persistent homology", "Graph topology", "Intrinsic dimensionality", "Cross-model alignment",
                     "SVD biological axes", "PPI network encoding", "Cell-type/family clustering",
                     "Attention-SVD dissociation", "Regulatory geometry (SV2-SV4)", "Regulatory geometry (SV5-SV7)",
                     "Edge-level geometry", "Signed regulation", "B-cell attractor dynamics"]
def _family_rules(fam, s):
    tl = s.lower()
    R = [("Attention-SVD dissociation", r"attention"),
         ("Cross-model alignment", r"geneformer|cross-model|\bcka\b|procrustes"),
         ("Persistent homology", r"persistent homology|ripser|betti|\bh1\b|\bh0\b|0-dim ph|persistence"),
         ("Cell-type/family clustering", r"cell-type marker|cell type marker|gene famil|marker gene cluster"),
         ("B-cell attractor dynamics", r"b-cell|bcell|gc-tf|\bgc\b|batf|bach2|bcl6|pax5|prdm1|plasma|lineage|\btcr\b|t-cell|myeloid|attractor|\bnk\b|cell-type axis|cell-type axes|apc"),
         ("Signed regulation", r"activation|repression|signed regulation|activator"),
         ("Edge-level geometry", r"\bffl\b|feed-forward|edge-level|edge auroc|edge probability|signed displacement|edge cosine|auroc mechanism|sparsity vs depth"),
         ("Regulatory geometry (SV5-SV7)", r"sv5-7|sv5-10|sv\(k|sv2-7|joint 6d|joint sv2|\b6d\b|tf boundary|tf/target|tf vs target|dual-role|tf family|laplacian|directionality|bfs hierarchy|tf->target|sv energy|effective rank as|effective rank vs auroc"),
         ("Regulatory geometry (SV2-SV4)", r"sv2-4|sv2-sv4"),
         ("PPI network encoding", r"string|\bppi\b|dorothea|go co-annotation|go biological process proximity|go ontology|chromosom|hub|confidence|trrust pairs distance|multi-predictor|encoding timeline|co-regulation|precision@k"),
         ("SVD biological axes", r"\bsv\d|\bpc\d|\bgo\b|go:|compartment|enrichment|co-pole|signal peptide|drift|housekeeping|svd|pole"),
         ("Intrinsic dimensionality", r"twonn|effective rank|participation ratio|intrinsic dim|spectral|singular value|wasserstein|\bswd\b|zero-norm|\boov\b|spectrum|compression|\bid\b"),
         ("Graph topology", r"knn|graph|spectral gap|communit|clustering|transitivity|modularity|connected component|trajectory|bifurcation|slope|centrality")]
    for nm, rx in R:
        if re.search(rx, tl): return nm
    return "unmapped"


def preprint_family(fam, text, name=""):
    """Heuristic assignment of a record to the preprint's 13-family taxonomy (Table 'hypothesis_summary'):
    rules applied to the record name first, then to name+method+result."""
    for s in (name, text):
        r = _family_rules(fam, s or "")
        if r != "unmapped": return r
    return "unmapped"


SV_RE = re.compile(r"\b(?:SV|PC)\s?\d+(?:\s?[-–/+]\s?(?:SV|PC)?\s?\d+)*|\bSV\(k,\s?k\+2\)", re.I)
def sv_indices(text):
    seen = []
    for m in SV_RE.findall(text):
        m = re.sub(r"\s+", "", m).upper()
        if m not in seen: seen.append(m)
    return ";".join(seen) if seen else ""


REF_RX = [("TRRUST", r"trrust|tf-target|tf→target|tf->target|589 (?:valid |positive )|735 pos|288 pos|tf genes|target-only|tf/target|tf vs target|tf-vs-target|tf-target|tf boundary|tf \(n="),
          ("STRING", r"string"), ("DoRothEA", r"dorothea"), ("GO", r"\bgo\b|go:|gene ontology|compartment"),
          ("Geneformer", r"geneformer"), ("scGPT attention", r"attention"), ("housekeeping (Eisenberg)", r"housekeeping"),
          ("chromosome position", r"chromosom"), ("marker panels", r"marker|lineage|b-cell|t-cell|myeloid|gc-tf"),
          ("co-expression", r"coexpr|co-expression"), ("GENIE3/GRNBoost2", r"genie3|grnboost"), ("perturbation data", r"perturb")]
def references(text):
    tl = text.lower(); out = [n for n, rx in REF_RX if re.search(rx, tl)]
    m = re.search(r"score\s*(?:>=|≥)\s*0?\.?(\d+)", tl)
    if "STRING" in out and m: out[out.index("STRING")] = f"STRING (score>={m.group(0).split('=')[-1].strip() if '=' in m.group(0) else m.group(1)})"
    return ";".join(out) if out else "none stated"


NP_PATTERNS = [
    r"(?:N|n|N_null|n_null|n_perm|N_perm|N_shuffle)\s*=\s*([\d,]+)\s*(?:[\w\- ]{0,30}?)(?:perm|shuffl|null|label|rep|random|draw)",
    r"(?:perm\w*|shuffl\w*|null)[^.;]{0,25}?\(\s*N\s*=\s*([\d,]+)",
    r"(?:perm\w*|shuffl\w*|null)[^.;]{0,8}?N\s*=\s*([\d,]+)",
    r"([\d,]+)\s*(?:[\w\-]+\s){0,3}(?:permutations?|perms?|shuffles?|label shuffles|reps|replicates|nulls?|random gene sets|random draws)\b",
    r"\((?:N\s*=\s*)?([\d,]+)\s*reps?\)",
]
def nperm_broad(method):
    vals = []
    for rx in NP_PATTERNS:
        for m in re.finditer(rx, method or "", re.I):
            ctx = (method[max(0, m.start() - 15): m.end() + 5]).lower()
            if "bootstrap" in ctx: continue
            try: vals.append(int(m.group(1).replace(",", "")))
            except ValueError: pass
    vals = [v for v in vals if v >= 10]
    return ";".join(str(v) for v in sorted(set(vals))) if vals else "not stated"


PHASE2_TOK = [(r"external[_ -]lung|krasnow", "lung_external"), (r"(?<!external_)(?<!external )\blung\b|lung_main", "lung"),
              (r"b_cell|b-cell context|bcell context", "immune_Bcell"), (r"\bcd4", "immune_CD4T"), (r"\bcd8", "immune_CD8T"),
              (r"tsp14", "immune_TSP14"), (r"tsp21", "immune_TSP21"), (r"tsp25", "immune_TSP25"),
              (r"\bimmune\b|immune_main", "immune"), (r"\bkidney\b", "kidney*"), (r"liver", "geneformer_liver"),
              (r"heart", "geneformer_heart"), (r"seed43", "seed43"), (r"seed44", "seed44")]
def ctx_phase2(script_ctx, text):
    if not script_ctx: return "none (re-analysis of earlier iteration outputs / synthesis)", "script (no matrix loaded)"
    tl = text.lower(); hits = set()
    for rx, c in PHASE2_TOK:
        if re.search(rx, tl): hits.add(c)
    scg = [c for c in script_ctx if not c.startswith("geneformer") and c != "lung_attention"]
    chosen = set()
    for c in scg:
        base = c.replace("_seed43", "").replace("_seed44", "")
        fam = ("kidney*" if base.startswith("kidney") else base)
        if fam in hits or base in hits:
            if c.endswith("seed43") and "seed43" not in hits: continue
            if c.endswith("seed44") and "seed44" not in hits: continue
            chosen.add(c)
    if "geneformer" in tl or "cross-model" in tl or "cross_model" in tl:
        chosen |= {c for c in script_ctx if c.startswith("geneformer") and c.split("_", 1)[1].replace("lung_external", "lung_external") in
                   {h.replace("immune_Bcell", "immune") for h in hits} | {"immune", "lung"}}
    if chosen: return ";".join(sorted(chosen)), "script contexts filtered by record text (heuristic)"
    return ";".join(script_ctx), "script-level (record does not name its domain; uncertain)"


def svdconv_phase2(text, script_txt):
    if not SPECTRAL_RE.search(text): return "not spectral"
    if "build_svd_projections" in script_txt or re.search(r"_cnt\[.*\]\s*>\s*0", script_txt):
        return ("subset-fitted on the observed edge-dataset genes (count>0), mean-centred, randomized SVD "
                "(iter_0082 helper build_svd_projections; script-level, unverified per record)")
    if "nonzero" in text.lower() or "non-zero" in text.lower(): return "observed non-zero rows only (record text)"
    if re.search(r"norm\([^)]*\)\s*>\s*0|nonzero|non_zero|valid_mask|cnt\[[^\]]*\]\s*>\s*0|count[^\n]{0,40}>\s*0", script_txt):
        return "observed non-zero rows only (script-level; unverified per record)"
    return "unknown"


def main():
    aud = pd.read_csv(HERE / "audit_409_records.csv")
    r2 = pd.read_csv(R2 / "master_p_round2.csv")
    r2["k"] = r2.groupby("iter").cumcount()
    recs = []
    sc_cache = {}
    for jp in sorted(ITER.glob("iter_*/executor_hypothesis_screen.json")):
        it_dir = jp.parent; it = it_dir.name; n = int(it[5:9])
        if it_dir not in sc_cache: sc_cache[it_dir] = script_contexts(it_dir)
        sctx, stxt = sc_cache[it_dir]
        rep = it_dir / "executor_iteration_report.md"
        for k, h in enumerate(json.loads(jp.read_text()).get("hypotheses", [])):
            text = " ".join(str(h.get(x, "")) for x in ("name", "method", "result_value"))
            hid = h.get("id", "")
            if n <= 63:
                ctx, basis = ctx_1_63(n, hid, text)
                mp = mapping_1_63(n, hid, ctx)
                gs = geneset_1_63(n, hid, ctx, text, mp)
                sv = svdconv_1_63(n, hid, ctx, text, mp)
                basis = "hand reading of iteration files: " + basis
            else:
                ctx, basis = ctx_phase2(sctx, text)
                mp = "true (edge-dataset indices / var_names; script-level)" if "kidney_A" in ctx or "immune" in ctx else (
                    "n/a" if ctx.startswith("none") else "unknown")
                gs = "unknown (script-level; see the iteration script)"
                if "matched" in text.lower() or "common edge" in text.lower() or "shared" in text.lower():
                    gs = "edge set matched across models/contexts (record text); size unknown"
                sv = svdconv_phase2(text, stxt)
            ctxs = [c for c in ctx.split(";") if c]
            tissue = ";".join(sorted({TISSUE.get(c, "unknown") for c in ctxs})) if not ctx.startswith("none") else "n/a"
            axis = ";".join(sorted({str(AXIS.get(c, "unknown")) for c in ctxs})) if not ctx.startswith("none") else "n/a"
            claim = CLAIMS.get((n, hid), "")
            spectral = bool(SPECTRAL_RE.search(text))
            if not claim: claim = "none (spectral record, not cited)" if spectral else "non-spectral"
            if n > 63: claim = "n/a (after the preprint's 63 iterations)"
            a = aud[(aud.iteration == it) & (aud.record_index == k)].iloc[0]
            row = dict(iteration=it, record_index=k, hypothesis_id=hid, name=h.get("name", ""), json_family=h.get("family", ""),
                       in_183=n <= 63, aborted_duplicate_run="incomplete" in it,
                       matrix_loaded=ctx, true_tissue=tissue, gene_axis_size=axis, gene_row_mapping=mp, gene_set=gs,
                       svd_convention=sv, spectral=spectral, sv_indices_recorded=sv_indices(text) or "none",
                       references=references(text), n_perm_stated=nperm_broad(h.get("method", "")),
                       n_perm_round2_regex=(int(a.n_perm_method) if pd.notna(a.n_perm_method) else "not stated"),
                       decision=h.get("decision", ""), result_direction=h.get("result_direction", ""),
                       preprint_claim=claim, preprint_family_heuristic=(preprint_family(h.get("family", ""), text, h.get("name", "")) if n <= 63 else "n/a"),
                       provenance_basis=basis,
                       script_contexts=";".join(sctx) if sctx else "none",
                       has_script=any(not p.name.startswith("._") for p in it_dir.glob("*.py")))
            if n <= 63:
                rr = r2[(r2.iter == it) & (r2.k == k)].iloc[0]
                assert rr.hypothesis_id == hid
                row.update(p_round1=rr.p_round1, primary_p_round2=rr.p_round2, n_perm_round2=rr.n_perm,
                           resolution_limited_round2=rr.resolution_limited, bh_q_round2=rr.bh_q,
                           bh_pass_round2=bool(rr.bh_pass), bonferroni_pass_round2=bool(rr.bonferroni_pass))
            else:
                row.update(p_round1="n/a (not in 183)", primary_p_round2="n/a (not in 183)", n_perm_round2="n/a",
                           resolution_limited_round2="n/a", bh_q_round2="n/a", bh_pass_round2="n/a", bonferroni_pass_round2="n/a")
            row.update(p_this_audit=a.p_floored, bh_pass_all_records=bool(a.bh_pass) if pd.notna(a.p_floored) else False,
                       bonferroni_pass_all_records=bool(a.bonferroni_pass) if pd.notna(a.p_floored) else False)
            recs.append(row)
    S2 = pd.DataFrame(recs)
    S2.to_csv(HERE / "S2_hypotheses.csv", index=False)

    # ---------------- S1: per JSON family ----------------
    fam_rows = []
    for fam, g in S2.groupby("json_family"):
        g183 = g[g.in_183]
        svs = Counter(s for v in g[g.spectral].sv_indices_recorded for s in v.split(";") if s and s != "none")
        refs = Counter(s for v in g.references for s in v.split(";") if s != "none stated")
        kid183 = g183[g183.matrix_loaded.str.contains("kidney_A")]
        fam_rows.append(dict(
            json_family=fam, n_183=len(g183), n_all=len(g), n_all_excl_aborted_run=int((~g.aborted_duplicate_run).sum()),
            n_extractable_p_183=int(pd.to_numeric(g183.primary_p_round2, errors="coerce").notna().sum()),
            n_extractable_p_all=int(g.p_this_audit.notna().sum()),
            n_bh_pass_round2_183=int((g183.bh_pass_round2 == True).sum()),
            n_bonferroni_pass_round2_183=int((g183.bonferroni_pass_round2 == True).sum()),
            n_bh_pass_all_records=int(g.bh_pass_all_records.sum()), n_bonferroni_pass_all_records=int(g.bonferroni_pass_all_records.sum()),
            n_spectral_183=int(g183.spectral.sum()), n_spectral_all=int(g.spectral.sum()),
            sv_directions_bands=";".join(f"{k}({v})" for k, v in svs.most_common(12)),
            references=";".join(f"{k}({v})" for k, v in refs.most_common()),
            n_183_on_kidney_A=len(kid183),
            n_183_on_kidney_A_labelled_immune_by_preprint=len(kid183),
            n_183_on_kidney_A_whose_own_record_name_says_immune_or_lung=int(kid183.name.fillna("").str.lower().str.contains("immune|lung").sum()),
            n_183_misaligned_labels=int(g183.gene_row_mapping.str.startswith("MISALIGNED").sum()),
            n_183_on_immune=int(g183.matrix_loaded.str.contains("immune").sum()),
            preprint_families_heuristic=";".join(f"{k}({v})" for k, v in Counter(g183.preprint_family_heuristic).most_common())))
    S1 = pd.DataFrame(fam_rows).sort_values("n_all", ascending=False)
    tot = {c: (S1[c].sum() if S1[c].dtype != object else "") for c in S1.columns}; tot["json_family"] = "TOTAL"
    S1 = pd.concat([S1, pd.DataFrame([tot])], ignore_index=True)
    S1.to_csv(HERE / "S1_family_table.csv", index=False)
    # JSON family x preprint family mapping (183 records)
    s183 = S2[S2.in_183]
    ct = pd.crosstab(s183.json_family, s183.preprint_family_heuristic)
    ct = ct.reindex(columns=[c for c in PREPRINT_FAMILIES + ["unmapped"] if c in ct.columns])
    ct.to_csv(HERE / "S1_family_mapping_json_x_preprint.csv")

    # ---------------- context usage ----------------
    use = defaultdict(lambda: dict(iters=set(), records=0, records_183=0, svs=Counter(), claims=Counter()))
    for r in S2.itertuples():
        for c in str(r.matrix_loaded).split(";"):
            if not c or c.startswith("none"): continue
            u = use[c]; u["iters"].add(int(r.iteration[5:9])); u["records"] += 1; u["records_183"] += int(r.in_183)
            if r.spectral:
                for s in r.sv_indices_recorded.split(";"):
                    if s and s != "none": u["svs"][s] += 1
            if r.preprint_claim.startswith("C") or r.preprint_claim.startswith("contradicts"):
                for cl in r.preprint_claim.split(";"): u["claims"][cl] += 1
    # script-level loading (a context can be loaded by an iteration's script without a record naming it)
    script_iters = defaultdict(set)
    for it_dir, (sctx, _) in sc_cache.items():
        for c in sctx: script_iters[c].add(int(it_dir.name[5:9]))
    cu = []
    for c in sorted(set(use) | set(script_iters)):
        u = use.get(c, dict(iters=set(), records=0, records_183=0, svs=Counter(), claims=Counter()))
        allit = u["iters"] | script_iters.get(c, set())
        cu.append(dict(context=c, true_tissue=TISSUE.get(c, "unknown"), gene_axis=AXIS.get(c, "unknown"),
                       n_iterations_loading_script_level=len(script_iters.get(c, set())),
                       n_iterations_with_records_on_it=len(u["iters"]), n_iterations_any=len(allit),
                       first_iteration=min(allit) if allit else "", last_iteration=max(allit) if allit else "",
                       n_records=u["records"], n_records_183=u["records_183"],
                       sv_bands_recorded=";".join(f"{k}({v})" for k, v in u["svs"].most_common(10)),
                       preprint_claims_derived=";".join(sorted(u["claims"]))))
    pd.DataFrame(cu).sort_values("n_records", ascending=False).to_csv(HERE / "context_usage.csv", index=False)
    summ = dict(n_records=len(S2), n_183=int(S2.in_183.sum()),
                n_183_kidney_A=int(s183.matrix_loaded.str.contains("kidney_A").sum()),
                n_183_immune=int(s183.matrix_loaded.str.contains("immune").sum()),
                n_183_misaligned=int(s183.gene_row_mapping.str.startswith("MISALIGNED").sum()),
                n_183_spectral=int(s183.spectral.sum()),
                n_183_by_svd_convention=s183.svd_convention.value_counts().to_dict(),
                n_all_by_svd_convention=S2.svd_convention.value_counts().to_dict(),
                n_183_preprint_family=s183.preprint_family_heuristic.value_counts().to_dict(),
                json_families_183=sorted(s183.json_family.unique()), json_families_all=sorted(S2.json_family.unique()),
                n_records_with_claim=int(S2.preprint_claim.str.match(r"^C\d").sum()),
                n_iterations_without_script=int(sum(not v for v in S2.drop_duplicates("iteration").has_script)))
    (HERE / "screen_tables_summary.json").write_text(json.dumps(summ, indent=1, default=str))
    print(json.dumps(summ, indent=1, default=str))


if __name__ == "__main__":
    main()
