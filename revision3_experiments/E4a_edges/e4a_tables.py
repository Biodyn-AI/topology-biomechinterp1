"""Print the markdown tables embedded in SUMMARY.md (generated from the CSVs, so every number is traceable)."""
from pathlib import Path
import pandas as pd

H = Path(__file__).resolve().parent
s = pd.read_csv(H / "h3_summary.csv")
pw = pd.read_csv(H / "power.csv")


def f(x, d=3):
    return "" if pd.isna(x) else f"{x:.{d}f}"


print("### Table 1 - primary H3 per context (h3_summary.csv; per-layer values in h3_primary_<ctx>.csv)\n")
print("| context | n pos / neg / TF | AUROC range over 12 layers | best layer: AUROC [95% CI] | max one-sided UB | min max-T p | abs co-expr AUROC (extraction / 3k) | max UB of delta vs co-expr (extraction / 3k) | verdict (extraction comparator) | verdict (3k comparator) |")
print("|---|---|---|---|---|---|---|---|---|---|")
for _, r in s.iterrows():
    print(f"| {r.context} | {r.n_pos} / {r.n_neg} / {r.n_tf} | {f(r.auroc_min)}-{f(r.auroc_max)} | L{int(r.best_layer)}: {f(r.best_auroc)} [{f(r.best_ci_lo)}, {f(r.best_ci_hi)}] | "
          f"{f(r.max_ub95)} | {r.min_p_maxT:.1e} | {f(r.auroc_cx_ext)} / {f(r.auroc_cx_3k)} | {f(r.max_delta_cxext_ub95)} / {f(r.max_delta_cx3k_ub95)} | {r.verdict} | {r.verdict_3k_comparator} |")

print("\n### Table 2 - a priori minimum detectable effect (power.csv; power 0.8, two-sided alpha 0.05)\n")
print("| design | statistic | test | MDE AUROC (sim) | size at null (sim) | MDE delta vs co-expr (sim) | oracle-SD MDE AUROC | Hanley-McNeil MDE AUROC | observed AUROC range |")
print("|---|---|---|---|---|---|---|---|---|")
for _, r in pw.iterrows():
    print(f"| {r.design} | {r.statistic} | {r.test} | {f(r.mde_auroc_sim)} | {f(r.size_auroc_sim, 3)} | {f(r.mde_delta_sim)} | {f(r.mde_auroc_oracle_sd)} | "
          f"{f(r.mde_auroc_hanley_mcneil)} | {f(r.observed_auroc_min)}-{f(r.observed_auroc_max)} |")
