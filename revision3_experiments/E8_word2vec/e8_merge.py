"""E8 step 4 - merge per-context evaluations into the deliverables and compute representation-similarity
diagnostics (linear CKA) for the window calibration and seed stability.

Writes: w2v_results.csv, paired_vs_scgpt.csv, h1_assortativity_scgpt.csv, calibration.csv, cka.csv, results.json
"""
from __future__ import annotations
import hashlib, json, re
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
R3 = HERE.parent
CTXS = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
EV = HERE / "eval"; MD = HERE / "models"
LOCAL_EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"


def cka(A, B):
    """Linear CKA (Kornblith et al. 2019) of two representations of the same genes (rows), centred."""
    A = A - A.mean(0); B = B - B.mean(0)
    hsic = np.linalg.norm(A.T @ B) ** 2
    return float(hsic / (np.linalg.norm(A.T @ A) * np.linalg.norm(B.T @ B)))


def model(ctx, corpus, w, e, s):
    p = MD / f"{ctx}__{corpus}__w{w}__e{e}__s{s}.npy"
    return np.load(p).astype(np.float64) if p.exists() else None


def main():
    w2v = pd.concat([pd.read_csv(EV / f"{c}_w2v.csv") for c in CTXS if (EV / f"{c}_w2v.csv").exists()])
    par = pd.concat([pd.read_csv(EV / f"{c}_paired.csv") for c in CTXS if (EV / f"{c}_paired.csv").exists()])
    h1s = pd.concat([pd.read_csv(EV / f"{c}_h1_scgpt.csv") for c in CTXS if (EV / f"{c}_h1_scgpt.csv").exists()])
    w2v = w2v.reset_index(drop=True); par = par.reset_index(drop=True)
    # CONVENTION3_ERRATA item 3: every vertex-subsampling interval (H1b) is also reported with the delete-d
    # scaling sqrt(m/(n-m)) = 1 at m = n/2, i.e. the pre-registered sqrt(m/n) deviations x sqrt(2) on each side.
    # The pre-registered columns (ci_lo / ci_hi, delta_ci_lo / delta_ci_hi, scgpt_ci_lo / scgpt_ci_hi) are kept.
    r2 = np.sqrt(2.0)
    m = w2v.task == "H1b_pair_cos_auroc"
    w2v.loc[m, "ci_lo_deleted"] = w2v.value[m] - r2 * (w2v.value[m] - w2v.ci_lo[m])
    w2v.loc[m, "ci_hi_deleted"] = w2v.value[m] + r2 * (w2v.ci_hi[m] - w2v.value[m])
    m = par.task == "H1b_pair_cos_auroc"
    par.loc[m, "delta_ci_lo_deleted"] = par.delta[m] - r2 * (par.delta[m] - par.delta_ci_lo[m])
    par.loc[m, "delta_ci_hi_deleted"] = par.delta[m] + r2 * (par.delta_ci_hi[m] - par.delta[m])
    par.loc[m, "scgpt_ci_lo_deleted"] = par.scgpt[m] - r2 * (par.scgpt[m] - par.scgpt_ci_lo[m])
    par.loc[m, "scgpt_ci_hi_deleted"] = par.scgpt[m] + r2 * (par.scgpt_ci_hi[m] - par.scgpt[m])
    w2v.to_csv(HERE / "w2v_results.csv", index=False); par.to_csv(HERE / "paired_vs_scgpt.csv", index=False)
    h1s.to_csv(HERE / "h1_assortativity_scgpt.csv", index=False)

    # ---------------------------------------------------------------- CKA diagnostics
    rows = []
    for ctx in CTXS:
        obs_idx = np.load(R3 / "cache" / ctx / "obs_idx.npy")
        emb = np.load(LOCAL_EMB / f"{ctx}.npy", mmap_mode="r")
        full = [m for m in (model(ctx, "full", 5, 10, s) for s in range(5)) if m is not None]
        extr = [m for m in (model(ctx, "extr", 5, 10, s) for s in range(5)) if m is not None]
        for i in range(len(full)):
            for j in range(i + 1, len(full)):
                rows.append(dict(context=ctx, a=f"full_w5_e10_s{i}", b=f"full_w5_e10_s{j}", kind="seed_pair_full", cka=cka(full[i], full[j])))
        for i in range(len(extr)):
            for j in range(i + 1, len(extr)):
                rows.append(dict(context=ctx, a=f"extr_w5_e10_s{i}", b=f"extr_w5_e10_s{j}", kind="seed_pair_extr", cka=cka(extr[i], extr[j])))
        for i in range(min(len(full), len(extr))):
            rows.append(dict(context=ctx, a=f"full_w5_e10_s{i}", b=f"extr_w5_e10_s{i}", kind="full_vs_extr", cka=cka(full[i], extr[i])))
        for l in range(12):
            X = np.asarray(emb[l], np.float64)[obs_idx]
            for i, m in enumerate(full):
                rows.append(dict(context=ctx, a=f"full_w5_e10_s{i}", b=f"scgpt_L{l}", kind="w2v_full_vs_scgpt", layer=l, cka=cka(m, X)))
            for i, m in enumerate(extr):
                rows.append(dict(context=ctx, a=f"extr_w5_e10_s{i}", b=f"scgpt_L{l}", kind="w2v_extr_vs_scgpt", layer=l, cka=cka(m, X)))
        w512 = model(ctx, "extr", 512, 10, 0)
        if w512 is not None:
            for p in sorted(MD.glob(f"{ctx}__extr__w5__e*__s0.npy")):
                e = int(re.search(r"__e(\d+)__", p.name).group(1))
                rows.append(dict(context=ctx, a="extr_w512_e10_s0", b=f"extr_w5_e{e}_s0", kind="calibration", cka=cka(w512, np.load(p).astype(float))))
            for l in range(12):
                X = np.asarray(emb[l], np.float64)[obs_idx]
                rows.append(dict(context=ctx, a="extr_w512_e10_s0", b=f"scgpt_L{l}", kind="w512_vs_scgpt", layer=l, cka=cka(w512, X)))
        # training-budget saturation (e8_train_sat.py): full corpus, window 5, 40 vs 10 epochs, seed 0
        e40 = model(ctx, "full", 5, 40, 0)
        if e40 is not None and len(full):
            rows.append(dict(context=ctx, a="full_w5_e40_s0", b="full_w5_e10_s0", kind="saturation_e40_vs_e10", cka=cka(e40, full[0])))
            for l in range(12):
                X = np.asarray(emb[l], np.float64)[obs_idx]
                rows.append(dict(context=ctx, a="full_w5_e40_s0", b=f"scgpt_L{l}", kind="w2v_full_e40_vs_scgpt", layer=l, cka=cka(e40, X)))
    ck = pd.DataFrame(rows); ck.to_csv(HERE / "cka.csv", index=False)

    # ---------------------------------------------------------------- saturation table (full corpus, seed 0)
    sat = []
    for ctx in CTXS:
        sub = w2v[(w2v.context == ctx) & (w2v.corpus == "full") & (w2v.seed.astype(str) == "0")]
        if not (sub.epochs == 40).any(): continue
        for e, g in sub.groupby("epochs"):
            rec = dict(context=ctx, model=f"full_w5_e{int(e)}_s0",
                       pair_updates_total=json.loads((MD / f"{ctx}__full__w5__e{int(e)}__s0.json").read_text())["pair_updates_total"])
            for t in ("H1b_pair_cos_auroc", "H2_strat_tf_auroc", "H3_withinTF_auroc"):
                r = g[g.task == t].iloc[0]; rec[t] = r.value; rec[t + "_ci_lo"] = r.ci_lo; rec[t + "_ci_hi"] = r.ci_hi
                if t.startswith("H1b"): rec[t + "_ci_lo_deleted"] = r.ci_lo_deleted; rec[t + "_ci_hi_deleted"] = r.ci_hi_deleted
                if t.startswith("H2"): rec["H2_null_mean"] = r.null_mean; rec["H2_p"] = r.p_raw
            a = g[g.task == "H1a_assortativity"]
            rec["H1a_max_r"] = a.value.max(); rec["H1a_n_sv_maxT_lt_005"] = int((a.p_maxT < 0.05).sum())
            sat.append(rec)
    sat = pd.DataFrame(sat); sat.to_csv(HERE / "saturation.csv", index=False)

    # ---------------------------------------------------------------- calibration table (seed 0, extraction corpus)
    cal = []
    for ctx in CTXS:
        sub = w2v[(w2v.context == ctx) & (w2v.corpus == "extr") & (w2v.seed.astype(str) == "0")]
        for (w, e), g in sub.groupby(["window", "epochs"]):
            if not ((w == 512) or (ctx in ("lung", "immune"))):
                continue
            rec = dict(context=ctx, model=f"extr_w{int(w)}_e{int(e)}_s0")
            for t in ("H1b_pair_cos_auroc", "H2_strat_tf_auroc", "H3_withinTF_auroc"):
                r = g[g.task == t].iloc[0]; rec[t] = r.value; rec[t + "_ci_lo"] = r.ci_lo; rec[t + "_ci_hi"] = r.ci_hi
                if t.startswith("H2"): rec["H2_null_mean"] = r.null_mean; rec["H2_p"] = r.p_raw
            a = g[g.task == "H1a_assortativity"]
            rec["H1a_max_r"] = a.value.max(); rec["H1a_argmax_sv"] = int(a.sv[a.value.idxmax()] if a.value.idxmax() in a.index else 0)
            rec["H1a_min_p_maxT"] = a.p_maxT.min(); rec["H1a_n_sv_maxT_lt_005"] = int((a.p_maxT < 0.05).sum())
            cal.append(rec)
    cal = pd.DataFrame(cal)
    if len(cal):
        cal = cal[cal.context.isin(cal.loc[cal.model.str.contains("w512"), "context"])]
    cal.to_csv(HERE / "calibration.csv", index=False)

    # ---------------------------------------------------------------- results.json
    txt = (R3 / "FROZEN.txt").read_text().split(); frozen = txt[txt.index("sha256") + 1]
    actual = hashlib.sha256((R3 / "CONVENTION3.md").read_bytes()).hexdigest()
    head = {}
    for ctx in CTXS:
        h = {}
        for g in ("full_w5_e10", "extr_w5_e10"):
            m = w2v[(w2v.context == ctx) & (w2v.corpus == g) & (w2v.seed.astype(str) == "mean")]
            if not len(m): continue
            d = {}
            for t in ("H1b_pair_cos_auroc", "H2_strat_tf_auroc", "H3_withinTF_auroc"):
                r = m[m.task == t].iloc[0]
                d[t] = dict(mean=r.value, sd_seeds=r.sd_seeds, ci=[r.ci_lo, r.ci_hi], n_seeds=int(r.n_seeds))
                if t.startswith("H1b"): d[t].update(ci_deleted=[r.ci_lo_deleted, r.ci_hi_deleted],
                                                    ci_note="ci = pre-registered sqrt(m/n); ci_deleted = errata item 3")
                if t.startswith("H2"): d[t].update(p_seedmean=r.p_raw, null_mean=r.null_mean)
                if t.startswith("H3"): d[t].update(upper95_1sided=r.upper95_1sided)
            a = m[m.task == "H1a_assortativity"]
            d["H1a"] = dict(sv_with_all_seeds_maxT_lt_005=[int(x) for x in a.sv[a.frac_seeds_maxT_lt_005 == 1]],
                            max_mean_r=float(a.value.max()), argmax_sv=int(a.sv[a.value.idxmax()]))
            pp = par[(par.context == ctx) & (par.w2v == g)]
            for t in ("H1b_pair_cos_auroc", "H2_strat_tf_auroc", "H3_withinTF_auroc"):
                q = pp[pp.task == t]
                d[t]["paired_vs_scgpt"] = dict(
                    delta_range=[float(q.delta.min()), float(q.delta.max())],
                    n_layers_w2v_higher=int((q.delta_ci_lo > 0).sum()), n_layers_scgpt_higher=int((q.delta_ci_hi < 0).sum()),
                    n_layers_ci_covers_0=int(((q.delta_ci_lo <= 0) & (q.delta_ci_hi >= 0)).sum()),
                    best_scgpt_layer=int(q.layer[q.scgpt.idxmax()]), best_scgpt=float(q.scgpt.max()),
                    delta_vs_best_layer=float(q.delta[q.scgpt.idxmax()]),
                    delta_vs_best_layer_ci=[float(q.delta_ci_lo[q.scgpt.idxmax()]), float(q.delta_ci_hi[q.scgpt.idxmax()])])
                if t.startswith("H2"):
                    d[t]["paired_vs_scgpt"].update(
                        n_layers_w2v_higher_nullcentred=int((q.delta_nullcentred_ci_lo > 0).sum()),
                        n_layers_scgpt_higher_nullcentred=int((q.delta_nullcentred_ci_hi < 0).sum()))
                if t.startswith("H1b"):
                    d[t]["paired_vs_scgpt"].update(
                        n_layers_w2v_higher_deleted=int((q.delta_ci_lo_deleted > 0).sum()),
                        n_layers_scgpt_higher_deleted=int((q.delta_ci_hi_deleted < 0).sum()),
                        delta_vs_best_layer_ci_deleted=[float(q.delta_ci_lo_deleted[q.scgpt.idxmax()]),
                                                        float(q.delta_ci_hi_deleted[q.scgpt.idxmax()])])
            h[g] = d
        js = EV / f"{ctx}.json"
        if js.exists():
            j = json.loads(js.read_text()); h["n"] = {k: j[k] for k in ("n_genes", "n_string_edges", "n_rewired_nonpairs",
                                                                        "n_lambert_tf", "h3_pos", "h3_neg", "h3_tfs")}
            h["ms_null"] = j["ms_null"]
        head[ctx] = h
    corpus_meta = {c: json.loads((HERE / "corpus" / f"{c}_meta.json").read_text()) for c in CTXS}
    model_meta = {p.stem: json.loads(p.read_text()) for p in sorted(MD.glob("*.json")) if not p.name.startswith("._")}
    res = dict(experiment="E8_word2vec", label="sensitivity (CONVENTION3 §5 'E8 Word2Vec'); no pass/fail",
               convention3_sha256_frozen=frozen, convention3_sha256_actual=actual, hash_match=frozen == actual,
               spec=dict(gensim="4.4.0", sg=1, vector_size=256, min_count=1, negative=10, workers=2, epochs=10, seeds=5,
                         window_specified=512, window_production=5,
                         shuffle="token order and cell order re-drawn per pass, default_rng([20261026, seed, pass])",
                         corpora=dict(full="every cell of the extraction h5ad", extr="A2-reproduced extraction cells"),
                         tokenisation="A2_axis_verification/a2_lib.tokenised_genes (max_genes=512), validated exactly",
                         H1a=dict(svs=10, null="Maslov-Sneppen burn 10m thin m", N=1000, correction="max-T over SV1-10",
                                  tail="greater"),
                         H1b="pair-cosine AUROC vs one MS-rewired non-pair set; vertex subsampling m=n/2, B=2000",
                         H2="detection-stratified LOO TF AUROC (Lambert); within-decile permutation N=2000; gene bootstrap B=2000",
                         H3="within-TF stratified AUROC, h3_pairs; two-way cluster bootstrap B=2000"),
               deviations=["window=5 (gensim default) with per-pass shuffling instead of window=512 for all production "
                           "models: compute (see SUMMARY.txt); window=512 run on the extraction corpus as calibration",
                           "gensim workers=2 is not bit-reproducible",
                           "training budget: production full-corpus models carry ~1/58 (lung) to ~1/(cells/256) of the "
                           "pair updates of the specified window-512 full-corpus model; saturation check with 40 epochs "
                           "(seed 0, lung and immune; saturation.csv)",
                           "H1a uses N = 1,000 Maslov-Sneppen draws (as specified for E8), not the N = 10,000 / 2,000 of the "
                           "primary H1; H2 permutation N = 2,000 (specified for E8), not 10,000",
                           "H1b intervals: pre-registered sqrt(m/n) columns kept; delete-d (x sqrt 2) columns added "
                           "(CONVENTION3_ERRATA item 3)",
                           "gensim 4.4.0 / Cython 3 prints 'Exception ignored in our_dot_float' when a float dot product equals "
                           "exactly -1.0 (sdot's error sentinel); that single dot product is then returned as 0. Benign "
                           "(~10^2 occurrences per ~10^9 dot products per run)"],
               scripts=["e8_corpus.py", "e8_train.py", "e8_train_sat.py", "e8_eval.py", "e8_merge.py", "e8_driver.sh"],
               corpus=corpus_meta, models=model_meta, headline=head,
               calibration=cal.to_dict(orient="records") if len(cal) else [],
               saturation=sat.to_dict(orient="records") if len(sat) else [],
               cka_summary=ck.groupby(["context", "kind"]).cka.agg(["mean", "min", "max", "count"]).reset_index().to_dict(orient="records"))
    (HERE / "results.json").write_text(json.dumps(res, indent=1, default=float))
    print("merged", len(w2v), len(par), len(ck))


if __name__ == "__main__":
    main()
