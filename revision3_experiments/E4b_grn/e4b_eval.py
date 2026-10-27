"""E4b evaluation: H3 metrics for GENIE3 / GRNBoost2, paired deltas, seed variability, direction test.

python e4b_eval.py [ctx ...]        evaluate contexts (default: all) from partial/*.npz, then merge.
python e4b_eval.py --merge          only merge eval/<ctx>_*.csv into the top-level outputs.

Metrics (exactly the H3 statistic of CONVENTION3 §4, computed with stats3):
  * within-TF stratified AUROC (stats3.StratifiedPairs.auroc; ties count 1/2 = mid-rank);
  * two-way (TF x target) pigeonhole cluster bootstrap, B = 2,000 (stats3.StratifiedPairs.bootstrap):
    95% percentile CI and one-sided 95% upper bound (stats3.one_sided_upper); TF-only bootstrap as a check;
  * paired deltas: every score in a context is bootstrapped with the SAME weight draws (fresh
    np.random.default_rng(BOOT_SEED + context index) per score), so delta draws = draws_a - draws_b;
  * comparators on identical pairs: |co-expression| (|Pearson r| of log-normalised expression) on the same
    cells as the GRN tier, and the scGPT primary scorer (signed cosine of V-centred embeddings, each layer);
  * within-TF label-permutation p (stats3.permute_within_tf, N = 2,000, one-sided) for the GRN seed
    ensembles and |co-expression| -- a sensitivity analysis (no pass/fail);
  * run-to-run SD of the AUROC over seeds.
Speed: FastPairs subclasses stats3.StratifiedPairs and only caches each TF's pair indices; the arithmetic is
the parent's, in the same order. Equality with the parent's auroc() is asserted at run time on every score
(unweighted and on 3 bootstrap weight vectors).
Direction test: see e4b_common.py.
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import e4b_common as C  # noqa: E402
import stats3 as S  # noqa: E402

PART = C.HERE / "partial"
EVAL = C.HERE / "eval"
CACHE_DIR = EVAL / "resample_cache"   # bootstrap / permutation draws (not results)
B = 2000
N_PERM = 2000
BOOT_SEED = 41_2026
PERM_SEED = 52_2026
DIR_BOOT_SEED = 63_2026
N_LAYERS = 12


class FastPairs(S.StratifiedPairs):
    """stats3.StratifiedPairs with each TF's pair indices cached (identical arithmetic, same order)."""

    def __init__(self, tf, tg, y):
        super().__init__(tf, tg, y)
        self._idx = [np.flatnonzero(self.tf_code == c) for c in range(self.tfs.size)]

    def auroc(self, s, w=None) -> float:
        s = np.asarray(s, float); w = np.ones(s.size) if w is None else np.asarray(w, float)
        num = den = 0.0
        for idx in self._idx:
            yy = self.y[idx]; ss = s[idx]; ww = w[idx]
            sp, wp = ss[yy], ww[yy]; sn, wn = ss[~yy], ww[~yy]
            if wp.sum() == 0 or wn.sum() == 0: continue
            order = np.argsort(sn, kind="mergesort"); sn, wn = sn[order], wn[order]
            cw = np.concatenate([[0.0], np.cumsum(wn)])
            lo = np.searchsorted(sn, sp, side="left"); hi = np.searchsorted(sn, sp, side="right")
            num += float((wp * (cw[lo] + 0.5 * (cw[hi] - cw[lo]))).sum())
            den += float(wp.sum() * wn.sum())
        return num / den if den > 0 else float("nan")

    def check_equal(self, s, rng):
        a = self.auroc(s); b = S.StratifiedPairs.auroc(self, s)
        assert abs(a - b) < 1e-12, (a, b)
        for _ in range(3):
            w = np.bincount(rng.integers(0, self.tfs.size, self.tfs.size), minlength=self.tfs.size)[self.tf_code] \
                * np.bincount(rng.integers(0, self.tgs.size, self.tgs.size), minlength=self.tgs.size)[self.tg_code]
            a = self.auroc(s, w.astype(float)); b = S.StratifiedPairs.auroc(self, s, w.astype(float))
            assert (np.isnan(a) and np.isnan(b)) or abs(a - b) < 1e-12, (a, b)

    def perm_null(self, s, n, rng):
        y0 = self.y.copy(); out = np.empty(n)
        try:
            for r in range(n):
                self.y = self.permute_within_tf(rng); out[r] = self.auroc(s)
        finally:
            self.y = y0
        return out


def tie_stats(sp: FastPairs, s):
    """Fraction of pairs scored exactly 0, and fraction of pairs tied with >= 1 other pair of the same TF."""
    s = np.asarray(s, float); tied = np.zeros(s.size, bool)
    for idx in sp._idx:
        _, inv, cnt = np.unique(s[idx], return_inverse=True, return_counts=True); tied[idx] = cnt[inv] > 1
    tp = np.zeros(s.size, bool)          # tied with >= 1 pair of the opposite label in the same TF
    for idx in sp._idx:
        v = s[idx]; yy = sp.y[idx]
        tp[idx] = np.isin(v, v[yy]) & np.isin(v, v[~yy])
    return float(np.mean(s == 0)), float(tied.mean()), float(tp.mean())


def summ(draws):
    d = draws[np.isfinite(draws)]
    return dict(ci_lo=float(np.quantile(d, 0.025)), ci_hi=float(np.quantile(d, 0.975)),
                upper95=S.one_sided_upper(d, 0.95), lower95=float(np.quantile(d, 0.05)),
                boot_mean=float(d.mean()), boot_sd=float(d.std(ddof=1)), boot_n=int(d.size))


def abs_coexpr(expr, a, b):
    X = np.asarray(expr, np.float64); mu = X.mean(0); sd = X.std(0)
    const = sd == 0; sd[const] = 1.0
    Z = (X - mu) / sd; Z[:, const] = 0.0
    r = np.einsum("ij,ij->j", Z[:, a], Z[:, b]) / X.shape[0]
    return np.abs(r), float(np.mean(const[a] | const[b]))     # |r|, fraction of pairs with a constant gene (r := 0)


def scgpt_cos(ctx, a, b):
    emb = np.load(C.EMB / f"{ctx}.npy", mmap_mode="r")
    obs_idx = np.load(C.CACHE / ctx / "obs_idx.npy")
    assert emb.shape[0] == N_LAYERS
    out = np.empty((N_LAYERS, a.size))
    for l in range(N_LAYERS):
        Xl = np.asarray(emb[l], dtype=np.float64)[obs_idx]           # observed genes V only
        assert np.all(np.linalg.norm(Xl, axis=1) > 0)
        Xl = Xl - Xl.mean(0, keepdims=True)                           # V-centred (CONVENTION3 §2)
        Xl = Xl / np.linalg.norm(Xl, axis=1, keepdims=True)
        out[l] = np.einsum("ij,ij->i", Xl[a], Xl[b])
    return out


def load_grn(ctx, tier, method, seed):
    f = PART / f"{ctx}_{tier}_{method}_s{seed}.npz"
    if not f.exists(): return None
    z = np.load(f)
    return dict(imp=z["imp"].astype(np.float64), targets=z["targets"], regs=z["regs"],
                n_iter=z["n_iter"], fitted=z["fitted"], n_cells=int(z["n_cells"]))


def direction_rows(ctx, label_method, tier, seed_label, fwd, rev, src, rng_seed):
    d = fwd - rev; nz = d != 0; k = int((d > 0).sum()); n = int(nz.sum())
    row = dict(context=ctx, method=label_method, tier=tier, seed=seed_label, n_pairs=int(d.size),
               n_tied=int((~nz).sum()), frac_tied=float(np.mean(~nz)), n_untied=n, n_forward_gt=k,
               frac_forward_gt=(k / n if n else np.nan))
    if n:
        bt = binomtest(k, n, 0.5, alternative="two-sided"); ci = bt.proportion_ci(0.95, method="exact")
        row.update(sign_test_p_two_sided=float(bt.pvalue), cp_lo=float(ci.low), cp_hi=float(ci.high))
        # cluster bootstrap over source regulators (pairs sharing a regulator are not independent)
        rng = np.random.default_rng(rng_seed); us, code = np.unique(src, return_inverse=True)
        fr = np.empty(B)
        for r in range(B):
            w = np.bincount(rng.integers(0, us.size, us.size), minlength=us.size)[code].astype(float)
            den = w[nz].sum(); fr[r] = w[d > 0].sum() / den if den > 0 else np.nan
        fr = fr[np.isfinite(fr)]
        row.update(src_cluster_lo=float(np.quantile(fr, 0.025)), src_cluster_hi=float(np.quantile(fr, 0.975)),
                   n_source_regs=int(us.size))
        row.update(median_log_ratio_untied=float(np.median(np.log((fwd[nz] + 1e-12) / (rev[nz] + 1e-12)))))
    return row


def eval_context(ctx, log):
    ci = C.CONTEXTS.index(ctx)
    sets = C.load_sets(ctx)
    sp = FastPairs(sets["tf"], sets["tg"], sets["y"])
    assert sp.tf.size == sets["tf"].size, "StratifiedPairs dropped pairs"
    tf, tg = sp.tf, sp.tg
    chk = np.random.default_rng(7)
    scores = {}     # key -> (score vector, meta dict)
    # comparators
    cells = {}
    for tier in ("T1", "T2"):
        ex = C.load_expr(ctx, tier); cells[tier] = ex.shape[0]
        r, fconst = abs_coexpr(ex, tf, tg)
        scores[("abs_coexpr", tier, "-", -1)] = (r, dict(n_cells=ex.shape[0], frac_pairs_constant_gene=fconst))
        del ex
    cos = scgpt_cos(ctx, tf, tg)
    for l in range(N_LAYERS):
        scores[("scgpt_cos", "emb", "-", l)] = (cos[l], dict(n_cells=256 if ctx != "kidney_B" else 192))
    # GRN
    grn_imp = {}
    for method in C.METHODS:
        for tier in ("T1", "T2"):
            per = []
            for seed in C.SEEDS[tier]:
                g = load_grn(ctx, tier, method, seed)
                if g is None: continue
                s, ok = C.importance_lookup(g["imp"], g["targets"], g["regs"], tf, tg)
                assert ok.all(), "an h3 pair has no importance entry"
                scores[(method, tier, str(seed), -1)] = (s, dict(n_cells=g["n_cells"],
                                                                  median_n_iter=float(np.median(g["n_iter"])) if method == "grnboost2" else np.nan,
                                                                  n_targets_unfitted=int((g["fitted"] == 0).sum())))
                per.append((seed, g))
            if per:
                ens = np.mean([g["imp"] for _, g in per], axis=0)
                s, _ = C.importance_lookup(ens, per[0][1]["targets"], per[0][1]["regs"], tf, tg)
                scores[(method, tier, "ens", -1)] = (s, dict(n_cells=per[0][1]["n_cells"], n_seeds=len(per),
                                                             seeds=",".join(str(x) for x, _ in per)))
                grn_imp[(method, tier)] = (per, ens)
    # point estimates, equality check, bootstrap (headline scores only: ensembles + comparators)
    rows = []; draws = {}
    for key, (s, meta) in scores.items():
        method, tier, seed, layer = key
        sp.check_equal(s, chk)
        z, tied, tied_pn = tie_stats(sp, s)
        row = dict(context=ctx, method=method, tier=tier, seed=seed, layer=layer, n_pairs=int(s.size),
                   n_pos=int(sp.y.sum()), n_tfs=int(sp.tfs.size), n_targets=int(sp.tgs.size),
                   auroc=sp.auroc(s), macro_auroc=sp.macro_auroc(s), frac_zero=z, frac_tied_within_tf=tied,
                   frac_tied_pos_neg=tied_pn, **{k: v for k, v in meta.items()})
        headline = seed in ("ens", "-")
        if headline:
            t0 = time.time()
            # resample cache (added on resume): bootstrap / permutation draws depend only on the score vector,
            # the pairs and the fixed seeds, so a score whose bytes are unchanged reuses its stored draws.
            cf = CACHE_DIR / ctx / ("_".join(map(str, key)) + ".npz"); cf.parent.mkdir(parents=True, exist_ok=True)
            hsh = hashlib.sha1(np.ascontiguousarray(s, dtype=np.float64).tobytes() + sp.y.tobytes()).hexdigest()
            cz = np.load(cf) if cf.exists() else None
            if cz is not None and str(cz["hash"]) != hsh: cz = None
            dr = cz["boot"] if cz is not None else sp.bootstrap(s, B, np.random.default_rng(BOOT_SEED + ci))
            draws[key] = dr
            row.update(summ(dr))
            if method != "scgpt_cos":
                tfo = cz["tfo"] if cz is not None else sp.bootstrap(s, B, np.random.default_rng(BOOT_SEED + 100 + ci), two_way=False)
                row.update(tf_only_ci_lo=float(np.nanquantile(tfo, 0.025)), tf_only_ci_hi=float(np.nanquantile(tfo, 0.975)),
                           tf_only_upper95=S.one_sided_upper(tfo, 0.95))
                S.assert_resolution(N_PERM, 0.05, 6)          # 4 GRN ensembles + 2 |co-expression| per context
                nul = cz["nul"] if cz is not None else sp.perm_null(s, N_PERM, np.random.default_rng(PERM_SEED + ci))
                if cz is None: np.savez(cf, hash=np.array(hsh), boot=dr, tfo=tfo, nul=nul)
            elif cz is None:
                np.savez(cf, hash=np.array(hsh), boot=dr)
            if method != "scgpt_cos":
                p, b_ge, _ = S.emp_p(row["auroc"], nul, "greater")
                row.update(perm_p_one_sided=float(p), perm_b=int(b_ge), perm_N=N_PERM, null_mean=float(nul.mean()), null_sd=float(nul.std(ddof=1)),
                           null_q95=float(np.quantile(nul, 0.95)), null_q99=float(np.quantile(nul, 0.99)))
            log(f"  {ctx} {key}: AUROC {row['auroc']:.4f} [{row['ci_lo']:.3f}, {row['ci_hi']:.3f}]  ({time.time() - t0:.0f}s)")
        rows.append(row)
    summary = pd.DataFrame(rows)
    # seed variability
    var_rows = []
    for method in C.METHODS:
        for tier in ("T1", "T2"):
            sub = summary[(summary.method == method) & (summary.tier == tier) & ~summary.seed.isin(["ens", "-"])]
            if len(sub):
                var_rows.append(dict(context=ctx, method=method, tier=tier, seed="seed_mean_sd", layer=-1,
                                     n_seeds=len(sub), auroc=float(sub.auroc.mean()),
                                     seed_sd=float(sub.auroc.std(ddof=1)) if len(sub) > 1 else np.nan,
                                     seed_min=float(sub.auroc.min()), seed_max=float(sub.auroc.max())))
    summary = pd.concat([summary, pd.DataFrame(var_rows)], ignore_index=True)
    # paired deltas
    drows = []
    def add_delta(ka, kb, comparison):
        a, b = draws[ka], draws[kb]; d = a - b; d = d[np.isfinite(d)]
        drows.append(dict(context=ctx, method=ka[0], tier=ka[1], layer=ka[3], comparator=kb[0], comparator_tier=kb[1],
                          comparator_layer=kb[3], comparison=comparison,
                          delta=summary_lookup(ka) - summary_lookup(kb),
                          ci_lo=float(np.quantile(d, 0.025)), ci_hi=float(np.quantile(d, 0.975)),
                          upper95=S.one_sided_upper(d, 0.95), lower95=float(np.quantile(d, 0.05)),
                          frac_boot_le0=float(np.mean(d <= 0))))
    auc_of = {k: r["auroc"] for k, r in zip(scores.keys(), rows)}
    summary_lookup = lambda k: auc_of[k]  # noqa: E731
    for method in C.METHODS:
        for tier in ("T1", "T2"):
            ka = (method, tier, "ens", -1)
            if ka not in draws: continue
            add_delta(ka, ("abs_coexpr", tier, "-", -1), "grn_minus_abs_coexpr_same_cells")
            for l in range(N_LAYERS):
                add_delta(ka, ("scgpt_cos", "emb", "-", l), "grn_minus_scgpt_cos")
        if (method, "T1", "ens", -1) in draws and (method, "T2", "ens", -1) in draws:
            add_delta((method, "T2", "ens", -1), (method, "T1", "ens", -1), "T2_minus_T1_same_method")
    for l in range(N_LAYERS):
        for tier in ("T1", "T2"):
            add_delta(("scgpt_cos", "emb", "-", l), ("abs_coexpr", tier, "-", -1), "scgpt_minus_abs_coexpr")
    deltas = pd.DataFrame(drows)
    # direction test
    dp = sets["dir_pairs"]; src, dst = dp[:, 0], dp[:, 1]
    dir_rows = []
    for (method, tier), (per, ens) in grn_imp.items():
        g0 = per[0][1]
        for lab, M in [("ens", ens)] + [(str(sd), g["imp"]) for sd, g in per]:
            f, okf = C.importance_lookup(M, g0["targets"], g0["regs"], src, dst)
            r, okr = C.importance_lookup(M, g0["targets"], g0["regs"], dst, src)
            assert okf.all() and okr.all()
            dir_rows.append(direction_rows(ctx, method, tier, lab, f, r, src, DIR_BOOT_SEED + ci))
    for tier in ("T1", "T2"):
        ex = C.load_expr(ctx, tier)
        f, _ = abs_coexpr(ex, src, dst); r, _ = abs_coexpr(ex, dst, src)
        dir_rows.append(direction_rows(ctx, "abs_coexpr", tier, "-", f, r, src, DIR_BOOT_SEED + ci))
    cf = scgpt_cos(ctx, src, dst); cr = scgpt_cos(ctx, dst, src)
    for l in range(N_LAYERS):
        dir_rows.append(direction_rows(ctx, f"scgpt_cos_L{l}", "emb", "-", cf[l], cr[l], src, DIR_BOOT_SEED + ci))
    direction = pd.DataFrame(dir_rows)
    # per-seed pair-score files
    for tier in ("T1", "T2"):
        for seed in C.SEEDS[tier]:
            out = {}
            for method in C.METHODS:
                if (method, tier, str(seed), -1) not in scores: continue      # file appeared after scoring
                g = load_grn(ctx, tier, method, seed)
                out[f"{method}_pair_score"] = scores[(method, tier, str(seed), -1)][0]
                out[f"{method}_importance"] = g["imp"].astype(np.float32)
                out[f"{method}_n_iter"] = g["n_iter"]
                out["targets"] = g["targets"]; out["regulators"] = g["regs"]
            if out:
                out.update(tf=tf, tg=tg, y=sp.y, abs_coexpr=scores[("abs_coexpr", tier, "-", -1)][0],
                           note=np.array("importance rows = targets, columns = regulators (indices into "
                                         "cache/<ctx>/obs_genes.json); pair scores aligned with tf/tg/y = h3_pairs"))
                np.savez_compressed(C.HERE / f"grn_scores_{ctx}_{tier}_{seed}.npz", **out)
    EVAL.mkdir(exist_ok=True)
    summary.to_csv(EVAL / f"{ctx}_summary.csv", index=False)
    deltas.to_csv(EVAL / f"{ctx}_deltas.csv", index=False)
    direction.to_csv(EVAL / f"{ctx}_direction.csv", index=False)
    (EVAL / f"{ctx}_meta.json").write_text(json.dumps(dict(
        context=ctx, n_regulators=int(sets["regs"].size), n_targets_fitted=int(sets["targets"].size),
        n_h3_targets=int(sets["h3_targets"].size), n_direction_pairs=int(dp.shape[0]), cells=cells,
        grn_runs_used=sorted("_".join(map(str, k[:3])) for k in scores if k[0] in C.METHODS and k[2] != "ens")), indent=1))
    log(f"{ctx}: done")


def merge():
    s = pd.concat([pd.read_csv(f) for f in sorted(EVAL.glob("*_summary.csv")) if not f.name.startswith("._")], ignore_index=True)
    d = pd.concat([pd.read_csv(f) for f in sorted(EVAL.glob("*_deltas.csv")) if not f.name.startswith("._")], ignore_index=True)
    r = pd.concat([pd.read_csv(f) for f in sorted(EVAL.glob("*_direction.csv")) if not f.name.startswith("._")], ignore_index=True)
    order = {c: i for i, c in enumerate(C.CONTEXTS)}
    for df in (s, d, r): df.sort_values("context", key=lambda x: x.map(order), kind="stable", inplace=True)
    s.to_csv(C.HERE / "grn_summary.csv", index=False)
    d.to_csv(C.HERE / "grn_deltas.csv", index=False)
    r.to_csv(C.HERE / "direction_test.csv", index=False)
    print("merged", len(s), len(d), len(r))


def main():
    log = lambda m: print(time.strftime("%H:%M:%S ") + m, flush=True)  # noqa: E731
    args = sys.argv[1:]
    if args != ["--merge"]:
        for ctx in (args or C.CONTEXTS):
            eval_context(ctx, log)
    merge()


if __name__ == "__main__":
    main()
