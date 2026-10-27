"""Verifier check 5: (a) GRNBoost2 early stopping at the 26th stage per fit -> v5_early_stop.csv;
(b) deltas whose ONE-sided 95% bound excludes 0 (two-sided CI does not) -> v5_onesided.csv;
(c) grn_scores_<ctx>_<tier>_<seed>.npz pair scores equal an independent lookup into partial/*.npz -> printed + json;
(d) descriptive application of the H3 margins to the GRN ensembles -> v5_h3_margins.csv. Writes v5_extra.json."""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
V = Path(__file__).resolve().parent; E = V.parent; sys.path.insert(0, str(E)); import e4b_common as C
rows = []
for ctx in C.CONTEXTS:
    for tier in ("T1", "T2"):
        for sd in C.SEEDS[tier]:
            z = np.load(E / "partial" / f"{ctx}_{tier}_grnboost2_s{sd}.npz"); n = z["n_iter"]; f = n >= 0
            rows.append(dict(context=ctx, tier=tier, seed=sd, n_fitted=int(f.sum()), n_stop_at_26=int((n == 26).sum()),
                             frac_stop_at_26=float((n == 26).sum() / f.sum()), n_reached_1000=int((n == 1000).sum()),
                             median_n_iter=float(np.median(n[f]))))
es = pd.DataFrame(rows); es.to_csv(V / "v5_early_stop.csv", index=False)
d = pd.read_csv(E / "grn_deltas.csv")
one = d[((d.lower95 > 0) | (d.upper95 < 0))]; one.to_csv(V / "v5_onesided.csv", index=False)
chk = {}
for ctx in C.CONTEXTS:
    h = np.load(C.CACHE / ctx / "h3_pairs.npz")
    for tier in ("T1", "T2"):
        for sd in C.SEEDS[tier]:
            g = np.load(E / f"grn_scores_{ctx}_{tier}_{sd}.npz")
            ok = np.array_equal(g["tf"], h["tf"]) and np.array_equal(g["tg"], h["tg"])
            for m in C.METHODS:
                z = np.load(E / "partial" / f"{ctx}_{tier}_{m}_s{sd}.npz")
                rp = {int(r): k for k, r in enumerate(z["regs"])}; tp = {int(t): k for k, t in enumerate(z["targets"])}
                mine = np.array([z["imp"][tp[int(b)], rp[int(a)]] for a, b in zip(h["tf"], h["tg"])], float)
                ok = ok and float(np.abs(mine - g[f"{m}_pair_score"]).max()) == 0.0
            chk[f"{ctx}_{tier}_{sd}"] = bool(ok)
s = pd.read_csv(E / "grn_summary.csv", dtype={"seed": str})
mr = []
for x in s[s.seed == "ens"].itertuples():
    dc = d[(d.context == x.context) & (d.method == x.method) & (d.tier == x.tier) & (d.comparison == "grn_minus_abs_coexpr_same_cells")].iloc[0]
    mr.append(dict(context=x.context, method=x.method, tier=x.tier, auroc=x.auroc, perm_p=x.perm_p_one_sided, upper95=x.upper95,
                   perm_p_lt_05=x.perm_p_one_sided < 0.05, auroc_ge_060_excluded=x.upper95 < 0.60,
                   delta_vs_coexpr_upper95=dc.upper95, not_above_coexpr_by_005=dc.upper95 < 0.05))
mr = pd.DataFrame(mr); mr.to_csv(V / "v5_h3_margins.csv", index=False)
out = dict(grn_scores_files_match_partial=chk, all_match=all(chk.values()), n_onesided_only=int(len(one)),
           early_stop_T1_frac_range=[float(es[es.tier == "T1"].frac_stop_at_26.min()), float(es[es.tier == "T1"].frac_stop_at_26.max())],
           early_stop_T2_frac_range=[float(es[es.tier == "T2"].frac_stop_at_26.min()), float(es[es.tier == "T2"].frac_stop_at_26.max())])
(V / "v5_extra.json").write_text(json.dumps(out, indent=1))
pd.set_option("display.width", 250); print(json.dumps(out)); print(es.groupby(["context", "tier"]).frac_stop_at_26.agg(["min", "max"]).round(3))
print(mr.round(4).to_string())
