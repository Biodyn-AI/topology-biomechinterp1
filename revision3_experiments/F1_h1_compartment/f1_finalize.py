"""Merge the post-hoc controls summary and headline numbers (all read from the CSVs) into results.json."""
import json, hashlib
import pandas as pd
import f1_common as F
H = F.HERE
res = json.loads((H / "results.json").read_text())
sm = pd.read_csv(H / "resid_controls_summary.csv")
res["posthoc_resid_controls"] = dict(note="post hoc, descriptive, no test; see f1_controls.py", rows=sm.to_dict(orient="records"))
st = pd.read_csv(H / "strat_null.csv"); ci = pd.read_csv(H / "intervals.csv"); rs = pd.read_csv(H / "residual_assort.csv")
wb = pd.read_csv(H / "within_between.csv")
head = {}
for ctx in res["contexts"]:
    L, k = F.E2_TOP[ctx][0]
    s = st[(st.context == ctx) & (st.layer == L) & (st.sv == k)].iloc[0]
    def g(df, **kw):
        q = df[(df.context == ctx) & (df.layer == L) & (df.sv == k)]
        for a, b in kw.items(): q = q[q[a] == b]
        return q.iloc[0]
    rr = g(ci, analysis="resid30cc", statistic="ratio_resid_over_raw"); rc = g(ci, analysis="residclass", statistic="ratio_resid_over_raw")
    head[ctx] = dict(cell=f"L{L} SV{k}", r=float(s.r), ms_excess=float(s.f1ms_excess_r),
                     strat_perm_null_mean=float(s.strat_null_mean), strat_excess=float(s.strat_excess_r),
                     strat_excess_frac_of_ms_excess=float(s.strat_excess_over_ms_excess),
                     cpms_null_mean=float(s.cpms_null_mean), cpms_excess=float(s.cpms_excess_r),
                     cpms_excess_frac_of_ms_excess=float(s.cpms_excess_over_ms_excess),
                     resid30_ratio=float(rr.est), resid30_ratio_ci_prereg=[float(rr.ci_lo_prereg), float(rr.ci_hi_prereg)],
                     resid30_ratio_ci_deleted=[float(rr.ci_lo_deleted), float(rr.ci_hi_deleted)],
                     residclass_ratio=float(rc.est), residclass_ratio_ci_deleted=[float(rc.ci_lo_deleted), float(rc.ci_hi_deleted)],
                     within_class_centred_r=float(g(wb, edges="within", coords="class_centred").r),
                     between_class_centred_r=float(g(wb, edges="between", coords="class_centred").r),
                     median_over_120_cpms_excess_frac=float(st[st.context == ctx].cpms_excess_over_ms_excess.median()),
                     median_over_120_resid30_ratio=float(rs[(rs.context == ctx) & (rs.design == "resid30cc")].ratio_resid_over_raw.median()),
                     median_over_120_residclass_ratio=float(rs[(rs.context == ctx) & (rs.design == "residclass")].ratio_resid_over_raw.median()))
res["headline_e2_top_cell"] = head
res["files"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()[:16] for p in sorted(H.glob("*.csv"))}
F.jdump(H / "results.json", res)
print(json.dumps(head, indent=1))
