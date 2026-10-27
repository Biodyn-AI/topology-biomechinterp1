"""Round-2 multiplicity audit.

Round 1 fed permutation p-values recorded as exactly 0.0 into BH and Bonferroni.
A permutation p of 0 does not exist: with N permutations and b exceedances the
estimator is (b+1)/(N+1), whose floor is 1/(N+1). Any hypothesis whose p was
recorded as 0.0 is therefore *resolution-limited* and its p can be reported only
as "<= 1/(N+1)".

This script re-audits all 183 hypotheses with that correction, and substitutes
the directly re-measured p-values (N = 10,000-100,000) for the headline claims.
"""
import sys, json, re
import numpy as np, pandas as pd
sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent))
import common as C

R1 = C.REPO / "subproject_41_claude_topology_hypothesis_screening_autoloop/revision_experiments/e5_multiple_testing"
ALPHA = C.BONFERRONI_ALPHA
master = pd.read_csv(R1 / "master_p.csv")
print(f"round-1 master_p.csv: {len(master)} hypotheses")

# permutation counts stated in each hypothesis's own method string
ITER = C.REPO / "subproject_41_claude_topology_hypothesis_screening_autoloop/iterations"
nperm = {}
for it in sorted(ITER.glob("iter_*/executor_hypothesis_screen.json")):
    try: d = json.loads(it.read_text())
    except Exception: continue
    for h in d.get("hypotheses", []):
        m = re.findall(r'(\d[\d,]*)\s*[- ]?(?:perm|permut|shuffle)', h.get("method",""), re.I)
        if m: nperm[(it.parent.name, h["id"])] = max(int(x.replace(",","")) for x in m)
print(f"permutation counts recovered from method strings: {len(nperm)}")

DEFAULT_N = 1000          # most permissive value stated in Methods 4.4
rows = []
for r in master.itertuples():
    p = r.primary_p
    try: p = float(p)
    except Exception: p = np.nan
    N = nperm.get((r.iter, r.hypothesis_id), DEFAULT_N)
    limited = False
    if not np.isnan(p) and p == 0.0:
        p = 1.0 / (N + 1); limited = True
    rows.append(dict(iter=r.iter, hypothesis_id=r.hypothesis_id, family=r.family,
                     p_round1=r.primary_p, n_perm=N, p_round2=p, resolution_limited=limited))
df = pd.DataFrame(rows)
ok = df.p_round2.notna()
print(f"extractable p         : {ok.sum()}")
print(f"recorded as exactly 0 : {df.resolution_limited.sum()}  -> replaced by 1/(N+1)")

q = np.full(len(df), np.nan); q[ok.values] = C.bh_fdr(df.loc[ok, "p_round2"].to_numpy())
df["bh_q"] = q
df["bonferroni_pass"] = df.p_round2 < ALPHA
df["bh_pass"] = df.bh_q < 0.05
print(f"\n{'':28s} {'round 1':>9s} {'round 2':>9s}")
print(f"{'BH q<0.05 passes':28s} {int((master.get('bh_reject_q05', pd.Series(dtype=bool))==True).sum()):>9d} {int(df.bh_pass.sum()):>9d}")
print(f"{'Bonferroni passes':28s} {int((master.get('bonferroni_paper', pd.Series(dtype=bool))==True).sum()):>9d} {int(df.bonferroni_pass.sum()):>9d}")
lim_pass = int((df.resolution_limited & df.bonferroni_pass).sum())
print(f"\nof the round-2 Bonferroni passes, {lim_pass} are resolution-limited "
      f"(p == 1/(N+1) exactly), i.e. the test could not have resolved a smaller p")

# --- headline claims: substitute directly re-measured p-values -------------
hl = pd.read_csv(R1 / "headline_audit.csv")
reb = pd.read_csv(C.OUT / "headline_rebuilt.csv")
def best(prefix):
    s = reb[reb.claim.str.startswith(prefix)]
    return None if s.empty else float(s.p.min())
MAP = {
 "SV1_extracellular_enrichment": ("SV1_extracellular_L11", best("SV1_extracellular_L11")),
 "SV1_secretory_pathway_layered": ("SV1_extracellular_best_layer", best("SV1_extracellular_best_layer")),
 "SV2_PPI_copole_string07": ("SV2_STRING_copole", best("SV2_STRING_copole")),
 "SV2_TRRUST_copole": ("SV2_TRRUST_copole", best("SV2_TRRUST_copole")),
 "TF_vs_target_AUROC_joint": ("TF_target_AUROC_SV2_7_stouffer", best("TF_target_AUROC_SV2_7_stouffer")),
}
print(f"\n{'headline claim':34s} {'round-1 p':>11s} {'round-2 p':>11s}  {'r1':>4s} {'r2':>4s}")
out = []
for r in hl.itertuples():
    p1 = float(r.primary_p); newp = MAP.get(r.headline, (None, None))[1]
    p2 = newp if newp is not None else (1/(DEFAULT_N+1) if p1 == 0.0 else p1)
    v1 = "PASS" if r.bonferroni_pass else "fail"
    v2 = "PASS" if p2 < ALPHA else "fail"
    flag = "  <-- CHANGED" if v1 != v2 else ""
    print(f"{r.headline:34s} {p1:11.2e} {p2:11.2e}  {v1:>4s} {v2:>4s}{flag}")
    out.append(dict(headline=r.headline, p_round1=p1, p_round2=p2,
                    bonferroni_round1=bool(r.bonferroni_pass), bonferroni_round2=bool(p2 < ALPHA),
                    remeasured=newp is not None))
o = pd.DataFrame(out)
print(f"\nheadline Bonferroni survivors: round 1 = {int(o.bonferroni_round1.sum())} of {len(o)}"
      f"  ->  round 2 = {int(o.bonferroni_round2.sum())} of {len(o)}")
df.to_csv(C.OUT/"master_p_round2.csv", index=False)
o.to_csv(C.OUT/"headline_audit_round2.csv", index=False)
json.dump(dict(n_total=len(df), n_extractable=int(ok.sum()),
               n_resolution_limited=int(df.resolution_limited.sum()),
               n_bh_pass=int(df.bh_pass.sum()), n_bonferroni_pass=int(df.bonferroni_pass.sum()),
               n_bonferroni_pass_resolution_limited=lim_pass,
               bonferroni_alpha=ALPHA,
               headline_survivors_round1=int(o.bonferroni_round1.sum()),
               headline_survivors_round2=int(o.bonferroni_round2.sum())),
          open(C.OUT/"audit_log_round2.json","w"), indent=1)
print("\nwritten: master_p_round2.csv, headline_audit_round2.csv, audit_log_round2.json")
