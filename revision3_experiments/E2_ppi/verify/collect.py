"""Collect the verifier's reproduced numbers next to the implementer's values -> verify/reproduced_numbers.csv."""
import json
from pathlib import Path
import pandas as pd
H = Path(__file__).resolve().parent
A = json.loads((H / "verify_a.json").read_text()); C = json.loads((H / "verify_c.json").read_text())
rows = []
def add(item, mine, theirs, src, method):
    rows.append(dict(item=item, verifier=mine, implementer=theirs, abs_diff=abs(mine - theirs), implementer_source=src, method=method))
hs = pd.read_csv(H.parent / "h1_summary.csv").set_index("context")
for ctx, a in A.items():
    if not isinstance(a, dict): continue
    add(f"{ctx} edges STRING-phys>=700 within V", a["edges_mine"], a["edges_cache"], "cache/edges_string.npy", "own parse of raw STRING files + h5ad var index")
    add(f"{ctx} top r ({a['top_cell']})", a["top_r_mine"], hs.loc[ctx, "top_r"], "h1_summary.csv top_r", "own eigh-based coordinates + corrcoef")
    add(f"{ctx} max |r diff| over 120 cells", a["max_abs_diff_r_vs_csv"], 0.0, f"h1_primary_{ctx}.csv r", "own")
    add(f"{ctx} family min p_adj (stored MS draws)", a["family_min_padj_mine"], hs.loc[ctx, "family_min_p_adj"], "h1_summary.csv", "own max-T on stored null")
    add(f"{ctx} n cells p_adj<0.05 (stored MS draws)", a["n_cells_padj_lt_05_mine"], hs.loc[ctx, "n_cells_p_adj_lt_05"], "h1_summary.csv", "own max-T on stored null")
for ctx in ("lung", "lung_external"):
    b = json.loads((H / f"verify_b_{ctx}.json").read_text()); ms = b["ms"]
    add(f"{ctx} MS null mean top cell (own chain N=1000)", ms["top_null_mean"], ms["csv_top_null_mean"], "h1_summary.csv top_null_mean", "own pure-Python double-edge swap chain")
    add(f"{ctx} MS null sd top cell (own chain N=1000)", ms["top_null_sd"], ms["csv_top_null_sd"], "h1_summary.csv top_null_sd", "own chain")
    add(f"{ctx} MS z top cell (own chain N=1000)", ms["top_z"], ms["csv_top_z"], "h1_summary.csv top_z", "own chain")
    add(f"{ctx} n cells max-T p_adj<0.05 (own chain N=1000)", ms["n_cells_padj_lt_05_mine"], 120, "h1_summary.csv", "own chain")
    add(f"{ctx} MS acceptance (own chain)", ms["acc_mean"], hs.loc[ctx, "ms_acc_mean"], "h1_summary.csv ms_acc_mean", "own chain")
    for nm, key in (("cosine top layer", [k for k in b["coexpr_extraction"] if k.startswith("cosine")][0]), ("-|dx| top cell", "top_cell_neg_absdiff"), ("|coexpr|", "abs_coexpr")):
        ce = pd.read_csv(H.parent / "sens_coexpr_matched.csv")
        sc = key
        t = ce[(ce.context == ctx) & (ce.coexpr_source == "extraction") & (ce.scorer == sc)].auroc_matched_stratified.iloc[0]
        add(f"{ctx} decile-stratified AUROC {nm} (extraction cells)", b["coexpr_extraction"][key]["strat"], t, "sens_coexpr_matched.csv", "own pool, own matching, sklearn per stratum")
for ctx in hs.index:
    c = C[ctx]
    add(f"{ctx} vertex CI lo (sqrt(m/n))", c["ci_sqrt_mn_mine"][0], c["ci_csv"][0], "h1_summary.csv top_r_ci_lo", "own subsampling, B=500")
    add(f"{ctx} vertex CI hi (sqrt(m/n))", c["ci_sqrt_mn_mine"][1], c["ci_csv"][1], "h1_summary.csv top_r_ci_hi", "own subsampling, B=500")
pd.DataFrame(rows).to_csv(H / "reproduced_numbers.csv", index=False)
print(pd.DataFrame(rows).to_string())
