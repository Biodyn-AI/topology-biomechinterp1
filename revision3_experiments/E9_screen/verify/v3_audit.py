"""Independent verifier (A6 / E9): multiplicity numbers from audit_409_records.csv recomputed with statsmodels,
cross-checked against round-2 master_p_round2.csv; scan of post-63 records for p-values the regex misses."""
import json, re
from pathlib import Path
import numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
P41 = HERE.parents[2]
a = pd.read_csv(HERE.parent / "audit_409_records.csv")
r2 = pd.read_csv(P41 / "revision2_experiments/master_p_round2.csv")
out = {}
out["round2_master"] = dict(n=len(r2), n_ext=int(r2.p_round2.notna().sum()), n_bonf=int(r2.bonferroni_pass.sum()),
                            n_bh=int(r2.bh_pass.sum()), n_res_lim=int(r2.resolution_limited.sum()))
# floored p, own
p = a.primary_p_raw.copy(); n = a.n_perm_method.fillna(1000)
p = np.where(p == 0, 1 / (n + 1), p)
ok = ~np.isnan(p)
for name, sel in (("all", np.ones(len(a), bool)), ("excl_incomplete", ~a.aborted_duplicate_run.astype(bool).values),
                  ("iter_1_63", (a.iter_num <= 63).values)):
    m = sel & ok; pv = p[m]
    rej_bh = multipletests(pv, alpha=0.05, method="fdr_bh")[0]
    out[name] = dict(n_records=int(sel.sum()), n_ext=int(m.sum()), n_p0=int(((a.primary_p_raw == 0).values & sel).sum()),
                     bonf_over_next=int((pv < 0.05 / m.sum()).sum()), bonf_over_183=int((pv < 0.05 / 183).sum()),
                     bonf_over_nrecords=int((pv < 0.05 / sel.sum()).sum()), bh=int(rej_bh.sum()))
# merge with round-2 for the 183 (order)
s = a[a.iter_num <= 63].reset_index(drop=True)
out["183_floored_p_match_round2"] = bool(np.allclose(np.nan_to_num(p[(a.iter_num <= 63).values], nan=-1),
                                                     np.nan_to_num(r2.p_round2.values, nan=-1)))
# scan for missed p-values after iteration 63
recs = []
for f in sorted((P41 / "iterations").glob("iter_*/executor_hypothesis_screen.json")):
    it = f.parent.name
    if int(it[5:9]) <= 63: continue
    for h in json.load(open(f))["hypotheses"]:
        recs.append((it, h.get("id"), str(h.get("result_value", ""))))
pat = re.compile(r"(?i)\bp(?:[_ -]?(?:val(?:ue)?|emp|perm|adj|min|max|mean|median|fdr|bh|q))?\s*[=<≤:]\s*[0-9]")
hits = [(it, hid) for it, hid, rv in recs if pat.search(rv)]
ext = set(map(tuple, a[(a.iter_num > 63) & a.primary_p_raw.notna()][["iteration", "hypothesis_id"]].values))
missed = [h for h in hits if h not in ext]
out["post63_records"] = len(recs); out["post63_with_broad_p_pattern"] = len(hits)
out["post63_extracted_by_audit"] = len(ext); out["post63_broad_hits_not_extracted"] = len(missed)
out["post63_missed_examples"] = [dict(iteration=i, id=h, text=[rv for it, hh, rv in recs if it == i and hh == h][0][:250])
                                  for i, h in missed[:40]]
json.dump(out, open(HERE / "v3_audit.json", "w"), indent=1, ensure_ascii=False)
print(json.dumps(out, indent=1, ensure_ascii=False))
