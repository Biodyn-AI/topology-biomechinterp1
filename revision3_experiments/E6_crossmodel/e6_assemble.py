"""E6 step 3: assemble the deliverable CSVs and results.json from work/."""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import kendalltau
sys.path.insert(0, str(Path(__file__).resolve().parent))
import e6lib as L
S = L.S
H = L.HERE; W = L.WORK

COLS = ["set", "embedding", "model", "label", "kind", "ctx", "layer", "d", "n_genes", "task", "estimator", "metric",
        "n_pos", "n_neg", "n_tf", "value", "ci_lo", "ci_hi", "ci_method", "ci_lo_pct", "ci_hi_pct", "boot_mean", "boot_sd",
        "upper95_one_sided", "p", "p_b", "p_N", "p_method", "p_raw", "null_mean", "null_sd", "null_q975",
        "best_band", "best_band_z", "best_z_band", "best_z", "p_best_z", "n_bands_padj05",
        "value_vocab_centred", "value_detection_stratified"]


def load_rows():
    rows = []
    for d in sorted(W.iterdir()):
        if d.name.startswith("._"): continue
        if not d.is_dir(): continue
        for f in sorted(x for x in d.glob("*.json") if not x.name.startswith("._")):
            rows += json.loads(f.read_text())["rows"]
    df = pd.DataFrame(rows)
    df = df[[c for c in COLS if c in df.columns] + [c for c in df.columns if c not in COLS]]
    # --- resume-agent corrections (see SUMMARY.txt) -------------------------------------------
    # (1) label: the UCE protein file is ESM2-15B per-protein (5,120-d), not ESM2-3B (2,560-d)
    df["label"] = df["label"].astype(str).str.replace(ESM_OLD, ESM_LABEL, regex=False)
    # (2) estimator (iii): replace max-T rows by the exact-sign 1-D-band recomputation (e6_patch1d.py)
    df = apply_p1d(df)
    # (3) CONVENTION3_ERRATA item 3: delete-d scaling of every vertex-subsampling interval
    vs = df.ci_method.astype(str).str.startswith("vertex subsampling")
    n = df.loc[vs, "n_genes"].astype(float); m = np.floor(0.5 * n); f = np.sqrt(n / (n - m))
    df["ci_lo_deleted"] = np.nan; df["ci_hi_deleted"] = np.nan; df["deleted_scale_factor"] = np.nan
    df.loc[vs, "ci_lo_deleted"] = df.loc[vs, "value"] - f * (df.loc[vs, "value"] - df.loc[vs, "ci_lo"])
    df.loc[vs, "ci_hi_deleted"] = df.loc[vs, "value"] + f * (df.loc[vs, "ci_hi"] - df.loc[vs, "value"])
    df.loc[vs, "deleted_scale_factor"] = f
    df.loc[vs, "ci_method"] = df.loc[vs, "ci_method"] + "; ci_lo/ci_hi = pre-registered sqrt(m/n); ci_*_deleted = delete-d sqrt(m/(n-m))"
    # (4) permutation resolution at alpha = 0.05 (stats3.assert_resolution) for every row with a p-value
    def res_ok(N):
        if not np.isfinite(N): return np.nan
        try: S.assert_resolution(int(N), 0.05); return True
        except ValueError: return False
    df["resolution_ok_alpha05"] = df["p_N"].astype(float).map(res_ok)
    return df


ESM_OLD = "ESM2 (UCE protein file; 'ESM2-3B' in submitted version)"
ESM_LABEL = "ESM2-15B (UCE per-protein)"
P1D_NOTE = "1-D bands rescored with exact sign normalisation (e6_patch1d.py)"


def load_p1d(setname, emb):
    f = W / setname / f"{emb}.p1d.npz"
    return np.load(f) if f.exists() else None


def apply_p1d(df):
    df = df.copy()
    for c in ("value_pre1dfix", "best_band_pre1dfix", "p_pre1dfix", "p1d_fixed"):
        df[c] = np.nan if c != "best_band_pre1dfix" else None
    df["p1d_fixed"] = False
    idx = df.index[df.estimator == "bestband"]
    for i in idx:
        r = df.loc[i]; z = load_p1d(r.set, r.embedding)
        if z is None or f"tobs_{r.task}" not in z.files:
            raise SystemExit(f"missing 1-D patch for {r.set} {r.embedding} {r.task}: run e6_patch1d.py first")
        to = z[f"tobs_{r.task}"].astype(float); tn = z[f"tnull_{r.task}"].astype(float)
        zz, pa, pr = S.max_t_adjust(to, tn); jb = int(np.argmax(to)); jz = int(np.argmax(zz))
        df.loc[i, ["value_pre1dfix", "p_pre1dfix"]] = [r.value, r.p]; df.loc[i, "best_band_pre1dfix"] = r.best_band
        df.loc[i, ["value", "best_band_z", "p", "p_raw", "best_z", "p_best_z", "n_bands_padj05"]] = [
            float(to[jb]), float(zz[jb]), float(pa[jb]), float(pr[jb]), float(zz[jz]), float(pa[jz]), int((pa < 0.05).sum())]
        df.loc[i, "best_band"] = L.BAND_NAMES[jb]; df.loc[i, "best_z_band"] = L.BAND_NAMES[jz]
        df.loc[i, "p1d_fixed"] = True; df.loc[i, "p_method"] = str(r.p_method) + "; " + P1D_NOTE
    return df


def boots(setname, emb, task, est="full"):
    z = np.load(W / setname / f"{emb}.npz"); return z[f"boot_{task}_{est}"].astype(float)


def crosslayer_maxt(setname, ctx):
    out = []
    for task, metric in (("node_lambert", "AUROC (LOO mean cosine to TF set)"), ("node_trrust", "AUROC (LOO mean cosine to TF set)"),
                         ("ppi", "pair-cosine AUROC"), ("edge", "within-TF stratified AUROC")):
        to, tn = [], []
        for l in range(12):
            f = W / setname / f"{ctx}@L{l}.npz"
            if not f.exists(): return out
            z = load_p1d(setname, f"{ctx}@L{l}")
            if z is None: raise SystemExit(f"missing 1-D patch for {setname} {ctx}@L{l}")
            if f"tobs_{task}" not in z.files: break
            to.append(z[f"tobs_{task}"]); tn.append(z[f"tnull_{task}"].astype(float))
        if len(to) < 12: continue
        to = np.concatenate(to); tn = np.concatenate(tn, 1)
        zz, pa, pr = S.max_t_adjust(to, tn)
        j = int(np.argmax(to)); jz = int(np.argmax(zz)); nb = len(L.BANDS)
        out.append(dict(set=setname, embedding=f"{ctx}@all_layers", model="scgpt_contextual", label=f"scGPT contextual ({ctx}), 12 layers",
                        kind="contextual residual stream", ctx=ctx, layer=-1, d=512, task=task, estimator="bestband_alllayers",
                        metric=metric + " (best of 78 bands x 12 layers)", value=float(to[j]),
                        best_band=f"L{j // nb} {L.BAND_NAMES[j % nb]}", best_band_z=float(zz[j]), p=float(pa[j]), p_raw=float(pr[j]),
                        p_N=int(tn.shape[0]), p_method="Westfall-Young max-T over 936 (layer, band) tests; " + P1D_NOTE,
                        best_z_band=f"L{jz // nb} {L.BAND_NAMES[jz % nb]}", best_z=float(zz[jz]), p_best_z=float(pa[jz]),
                        n_bands_padj05=int((pa < 0.05).sum())))
    return out


def delta_row(setname, task, a, b, va, vb, ba, bb, la, lb):
    d = va - vb; db = ba - bb; shift = d - float(np.mean(db)); lo, hi = L.ci(db)
    lo, hi = lo + shift, hi + shift
    return dict(set=setname, task=task, estimator="full", model_a=a, label_a=la, model_b=b, label_b=lb,
                value_a=va, value_b=vb, delta=d, ci_lo=lo, ci_hi=hi, B=int(db.size),
                frac_boot_delta_le0=float(np.mean(db + shift <= 0)), ci_excludes_0=bool(lo > 0 or hi < 0),
                ci_method="paired gene bootstrap (identical gene draws for both embeddings), percentile re-centred")


def main():
    df = load_rows()
    if "set" not in df: raise SystemExit("no rows")
    tierA = df[(df.set == "tierA") & ~df.embedding.str.contains("@rp")]
    tierA.to_csv(H / "tierA_results.csv", index=False)
    tb = [df[df.set.str.startswith("tierB_")]]
    for ctx in L.CORE_CTX: tb.append(pd.DataFrame(crosslayer_maxt(f"tierB_{ctx}", ctx)))
    tierB = pd.concat(tb, ignore_index=True); tierB.to_csv(H / "tierB_results.csv", index=False)

    # ---------------- legacy reproduction
    sub = pd.read_csv(L.R3.parent / "revision2_experiments/crossmodel_round2.csv")
    rep = pd.read_csv(H / "legacy_submitted_repro_raw.csv")
    br = pd.read_csv(W / "legacy_bridge.csv")
    br["label"] = br["label"].astype(str).str.replace(ESM_OLD, ESM_LABEL, regex=False)
    fdd = np.sqrt(186 / (186 - 93))
    br["ppi_functional_ci_lo_deleted"] = br.ppi_functional_auroc_full - fdd * (br.ppi_functional_auroc_full - br.ppi_functional_ci_lo)
    br["ppi_functional_ci_hi_deleted"] = br.ppi_functional_auroc_full + fdd * (br.ppi_functional_ci_hi - br.ppi_functional_auroc_full)
    leg = df[df.set.isin(["legacy186", "legacy178"])]
    name_map = {"scGPT L11": "kidney_A@L11", "Geneformer V1-10M": "gf_v1_10m", "Geneformer V2-104M": "gf_v2_104m", "ESM2-3B": "esm2_uce"}
    LR = []
    for _, r in sub.iterrows():
        rr = rep[rep.model == r.model].iloc[0]; e = name_map[r.model]
        bb = br[br.embedding == e].iloc[0]
        def g(task, est="full", s="legacy186"):
            x = leg[(leg.set == s) & (leg.embedding == e) & (leg.task == task) & (leg.estimator == est)]
            return x.iloc[0] if len(x) else None
        nt, nl, pp, ed = g("node_trrust"), g("node_lambert"), g("ppi"), g("edge")
        sv = g("node_trrust", "SV2-SV7")
        LR.append(dict(model=r.model, label_corrected=(ESM_LABEL if e == "esm2_uce" else r.model), embedding=e, d=int(r.dim),
                       submitted_tf_target_auroc=r.tf_target_auroc, reproduced_tf_target_auroc=rr.tf_target_auroc,
                       submitted_copole_z=r.ppi_z, reproduced_copole_z=rr.ppi_z,
                       reproduced_exactly=bool(abs(r.tf_target_auroc - rr.tf_target_auroc) < 1e-12 and abs(r.ppi_z - rr.ppi_z) < 1e-12),
                       bridge_subsetSVD_SV2_SV7_auroc=bb.submitted_style_tf_target_auroc_subsetSVD_SV2_SV7,
                       est_i_tf_target_auroc=nt.value, est_i_tf_target_ci_lo=nt.ci_lo, est_i_tf_target_ci_hi=nt.ci_hi, est_i_tf_target_p=nt.p,
                       vocabSVD_SV2_SV7_tf_target_auroc=sv.value if sv is not None else np.nan,
                       est_i_lambert_auroc=nl.value, est_i_lambert_ci_lo=nl.ci_lo, est_i_lambert_ci_hi=nl.ci_hi,
                       est_i_ppi_physical_auroc=pp.value, est_i_ppi_physical_ci_lo=pp.ci_lo, est_i_ppi_physical_ci_hi=pp.ci_hi,
                       est_i_ppi_physical_ci_lo_deleted=pp.ci_lo_deleted, est_i_ppi_physical_ci_hi_deleted=pp.ci_hi_deleted,
                       n_physical_pairs=int(pp.n_pos),
                       est_i_ppi_functional_auroc=bb.ppi_functional_auroc_full, est_i_ppi_functional_ci_lo=bb.ppi_functional_ci_lo,
                       est_i_ppi_functional_ci_hi=bb.ppi_functional_ci_hi, n_functional_pairs=int(bb.functional_string_pairs),
                       est_i_ppi_functional_ci_lo_deleted=bb.ppi_functional_ci_lo_deleted, est_i_ppi_functional_ci_hi_deleted=bb.ppi_functional_ci_hi_deleted,
                       est_i_edge_auroc=ed.value, est_i_edge_ci_lo=ed.ci_lo, est_i_edge_ci_hi=ed.ci_hi, n_edge_pos=int(ed.n_pos)))
    LRdf = pd.DataFrame(LR)
    # every other legacy-set embedding (estimator (i))
    extra = []
    for e in sorted(set(leg[leg.set == "legacy186"].embedding)):
        x = leg[(leg.set == "legacy186") & (leg.embedding == e) & (leg.estimator == "full")]
        d = dict(model=x.label.iloc[0], embedding=e, d=int(x.d.iloc[0]))
        tn = {"node_trrust": "tf_target", "node_lambert": "lambert", "ppi": "ppi_physical", "edge": "edge"}
        for _, r in x.iterrows():
            t_ = tn[r.task]
            d[f"est_i_{t_}_auroc"] = r.value; d[f"est_i_{t_}_ci_lo"] = r.ci_lo; d[f"est_i_{t_}_ci_hi"] = r.ci_hi
            if t_ == "ppi_physical": d["est_i_ppi_physical_ci_lo_deleted"] = r.ci_lo_deleted; d["est_i_ppi_physical_ci_hi_deleted"] = r.ci_hi_deleted
        bb = br[br.embedding == e]
        if len(bb): d["bridge_subsetSVD_SV2_SV7_auroc"] = bb.submitted_style_tf_target_auroc_subsetSVD_SV2_SV7.iloc[0]
        extra.append(d)
    LRall = pd.concat([LRdf.assign(section="submitted Table 1 models"),
                       pd.DataFrame(extra).assign(section="all embeddings on legacy186, estimator (i)")], ignore_index=True)
    LRall.to_csv(H / "legacy_reproduction.csv", index=False)

    # ---------------- projection control
    P = []
    for fam, base in (("esm2_uce@rp1024", "esm2_uce"), ("scgpt_static@rp102", "scgpt_static")):
        for task in ("node_lambert", "node_trrust", "ppi", "edge"):
            x = df[(df.set == "tierA") & df.embedding.str.startswith(fam) & (df.task == task) & (df.estimator == "full")]
            f = tierA[(tierA.embedding == base) & (tierA.task == task) & (tierA.estimator == "full")]
            if not len(x) or not len(f): continue
            P.append(dict(family=fam, base_model=base, task=task, d_full=int(f.d.iloc[0]), d_projected=int(x.d.iloc[0]),
                          full_value=float(f.value.iloc[0]), full_ci_lo=float(f.ci_lo.iloc[0]), full_ci_hi=float(f.ci_hi.iloc[0]),
                          n_seeds=len(x), proj_mean=float(x.value.mean()), proj_sd=float(x.value.std(ddof=1)),
                          proj_min=float(x.value.min()), proj_max=float(x.value.max()),
                          proj_minus_full=float(x.value.mean() - f.value.iloc[0]),
                          seeds_values=";".join(f"{v:.4f}" for v in x.sort_values("embedding").value),
                          proj_ci_lo_mean=float(x.ci_lo.mean()), proj_ci_hi_mean=float(x.ci_hi.mean())))
    Pdf = pd.DataFrame(P)
    if len(Pdf):
        sc = Pdf[Pdf.base_model == "scgpt_static"].set_index("task")
        Pdf["scgpt_static_full_value"] = Pdf.task.map(sc.full_value)
    Pdf.to_csv(H / "projection_control.csv", index=False)

    # ---------------- pairwise deltas (node tasks, estimator (i), paired gene bootstrap)
    D = []
    def val(s, e, t):
        x = df[(df.set == s) & (df.embedding == e) & (df.task == t) & (df.estimator == "full")]
        return (float(x.value.iloc[0]), x.label.iloc[0]) if len(x) else (None, None)
    for task in ("node_lambert", "node_trrust"):
        ks = L.STATIC_KEYS
        for i in range(len(ks)):
            for j in range(len(ks)):
                if i == j: continue
                va, la = val("tierA", ks[i], task); vb, lb = val("tierA", ks[j], task)
                if va is None or vb is None or i > j: continue
                D.append(delta_row("tierA", task, ks[i], ks[j], va, vb, boots("tierA", ks[i], task), boots("tierA", ks[j], task), la, lb))
        for ctx in L.CORE_CTX:
            s = f"tierB_{ctx}"
            for l in range(12):
                c = f"{ctx}@L{l}"; vb, lb = val(s, c, task)
                if vb is None: continue
                for k in ks:
                    va, la = val(s, k, task)
                    if va is None: continue
                    D.append(delta_row(s, task, k, c, va, vb, boots(s, k, task), boots(s, c, task), la, lb))
        for s in ("legacy186", "legacy178"):
            c = "kidney_A@L11"; vb, lb = val(s, c, task)
            for k in ks:
                va, la = val(s, k, task)
                if va is None or vb is None: continue
                D.append(delta_row(s, task, k, c, va, vb, boots(s, k, task), boots(s, c, task), la, lb))
    Ddf = pd.DataFrame(D); Ddf.to_csv(H / "pairwise_deltas.csv", index=False)

    # ---------------- selection bias + Kendall tau
    sb = pd.read_csv(W / "selection_bias_core.csv"); smeta = json.loads((W / "selection_meta.json").read_text())
    ta = tierA[(tierA.task == "node_trrust") & (tierA.estimator == "full")].set_index("embedding")
    sb["label"] = sb["label"].astype(str).str.replace(ESM_OLD, ESM_LABEL, regex=False)
    sb["row_type"] = "subset_null"
    sb["tierA_value"] = sb.model.map(ta.value); sb["tierA_ci_lo"] = sb.model.map(ta.ci_lo); sb["tierA_ci_hi"] = sb.model.map(ta.ci_hi)
    # cross-check legacy values against e6_run
    for _, r in sb.iterrows():
        v, _ = val(r.variant, r.model, "node_trrust")
        if v is not None: assert abs(v - r.legacy_value) < 1e-6, (r.variant, r.model, v, r.legacy_value)
    T = []
    rng = np.random.default_rng(L.SEED + 909)
    for v in ("legacy186", "legacy178"):
        x = sb[sb.variant == v]
        ms = [m for m in x.model if (W / v / f"{m}.npz").exists() and (W / "tierA" / f"{m}.npz").exists()]
        if len(ms) < 3: continue
        x = x.set_index("model").loc[ms]
        tau, p = kendalltau(x.legacy_value, x.tierA_value)
        bl = np.array([boots(v, m, "node_trrust") for m in ms]); bt = np.array([boots("tierA", m, "node_trrust") for m in ms])
        taus = np.array([kendalltau(bl[:, r], bt[:, r])[0] for r in range(bl.shape[1])])
        lo, hi = L.ci(taus)
        tau_nl, _ = kendalltau(x.null_mean, x.legacy_value); tau_nt, _ = kendalltau(x.null_mean, x.tierA_value)
        T.append(dict(row_type="kendall_tau", variant=v, model="ALL", n_models=len(ms), kendall_tau_legacy_vs_tierA=tau,
                      kendall_p_asymptotic=p, tau_ci_lo=lo, tau_ci_hi=hi, tau_boot_B=int(taus.size),
                      tau_ci_method="paired gene bootstrap draws (legacy and Tier A resampled independently), percentile",
                      kendall_tau_subsetmean_vs_legacy=tau_nl, kendall_tau_subsetmean_vs_tierA=tau_nt,
                      order_legacy=">".join(x.sort_values("legacy_value", ascending=False).index),
                      order_tierA=">".join(x.sort_values("tierA_value", ascending=False).index),
                      order_subset_mean=">".join(x.sort_values("null_mean", ascending=False).index)))
    SB = pd.concat([sb, pd.DataFrame(T)], ignore_index=True); SB.to_csv(H / "selection_bias.csv", index=False)

    # ---------------- results.json
    att = pd.read_csv(H / "attrition.csv")
    def pick(frame, emb, task, est="full", s=None):
        x = frame[(frame.embedding == emb) & (frame.task == task) & (frame.estimator == est)]
        if s is not None: x = x[x.set == s]
        if not len(x): return None
        r = x.iloc[0]
        return {k: (None if (isinstance(r[k], float) and np.isnan(r[k])) else r[k]) for k in
                ("value", "ci_lo", "ci_hi", "ci_lo_deleted", "ci_hi_deleted", "p", "p_N", "best_band", "n_pos", "n_neg") if k in r}
    res = dict(experiment="E6_crossmodel", convention3_sha256=L.FROZEN,
               tierA_size=int(att[(att.stage == "Tier A") & (att.set.str.startswith("intersection"))].n_genes.iloc[0]),
               attrition=att.to_dict("records"),
               legacy_reproduction=LRdf.to_dict("records"),
               tierA={k: {t: {e: pick(tierA, k, t, e) for e in ("full", "SV2-SV7", "SV1-SV6", "bestband")}
                          for t in ("node_lambert", "node_trrust", "ppi", "edge")} for k in L.STATIC_KEYS},
               tierB={s: {e: {t: pick(tierB, e, t, "full", s) for t in ("node_lambert", "node_trrust", "ppi", "edge")}
                          for e in sorted(set(tierB[tierB.set == s].embedding)) if "@all" not in e}
                      for s in sorted(set(tierB.set))},
               tierB_crosslayer_maxt=tierB[tierB.estimator == "bestband_alllayers"].to_dict("records"),
               selection_bias=dict(meta=smeta, rows=SB.to_dict("records")),
               projection_control=Pdf.to_dict("records"),
               n_pairwise_rows=int(len(Ddf)),
               models={k: dict(label=(ESM_LABEL if k == "esm2_uce" else v[1]), file=v[0], provenance=v[2], kind=v[3])
                       for k, v in L.STATIC.items()},
               completed_embeddings={s_: sorted(set(df[df.set == s_].embedding)) for s_ in sorted(set(df.set))},
               corrections_by_resume_agent=[
                   "estimator (iii): the 12 one-dimensional bands were scored through stats3.l2norm_rows (|x|+1e-12), "
                   "which breaks exact ties by magnitude; observed and all 1,000 null draws of these bands were recomputed "
                   "with exact sign normalisation on the same null draws (e6_patch1d.py); max-T recomputed. Node 1-D values change "
                   "only where TF and non-TF scores tie exactly (e.g. all TFs share one sign on a direction: 8 of 232 legacy-set "
                   "node rows); PPI and edge 1-D values change generally. Effect per row: p1d_fix_effect.csv (e6_fixreport.py).",
                   "label: the UCE protein file (5,120-d) is ESM2-15B per-protein, reported as 'ESM2-15B (UCE per-protein)'; "
                   "the submitted manuscript called it ESM2-3B.",
                   "CONVENTION3_ERRATA item 3: every vertex-subsampling interval also reported with delete-d scaling "
                   "sqrt(m/(n-m)) (columns ci_lo_deleted / ci_hi_deleted); pre-registered sqrt(m/n) kept in ci_lo / ci_hi."],
               resolution_note="stats3.assert_resolution(N, 0.05) recorded per row (resolution_ok_alpha05); node/edge N=2,000, "
                               "max-T and PPI band nulls N=1,000, selection-bias null N=10,000.")
    # supplement (e6_supp.py): detection-stratified node CIs/p and paired deltas vs scGPT for every task
    if (H / "node_detstrat_results.csv").exists() and (H / "pairwise_deltas_vs_scgpt.csv").exists():
        nd = pd.read_csv(H / "node_detstrat_results.csv"); pdv = pd.read_csv(H / "pairwise_deltas_vs_scgpt.csv")
        res["node_detstrat"] = nd.where(nd.notna(), None).to_dict("records")
        key = pdv[((pdv.set == "tierA") & (pdv.model_b == "scgpt_static")) |
                  (pdv.set.str.startswith("tierB_") & (pdv.model_b.str.contains("@L"))) |
                  (pdv.set.str.startswith("legacy") & (pdv.model_b == "kidney_A@L11"))]
        res["pairwise_deltas_vs_scgpt"] = key.where(key.notna(), None).to_dict("records")
        res["supplement_note"] = ("e6_supp.py regenerated, with e6_run's seeds and rng streams, the PPI vertex-subsampling and edge "
                                  "two-way bootstrap draws (stored CIs reproduced to 1e-9) and added a gene-bootstrap CI and a "
                                  "within-decile permutation p (N=2,000) to the detection-stratified Lambert AUROC.")
    (H / "results.json").write_text(json.dumps(res, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    print("assembled", len(df), "rows")


if __name__ == "__main__":
    main()
