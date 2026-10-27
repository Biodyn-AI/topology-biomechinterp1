"""Numbers behind the verifier's wording corrections, read from the implementer's CSVs (and gene2go for the MF
term counts). Writes verify/wording_numbers.json."""
import json, pickle, sys
from collections import Counter
from pathlib import Path
import pandas as pd
V = Path(__file__).resolve().parent; F1 = V.parent
sys.path.insert(0, str(V)); import verify_f1 as VF
st = pd.read_csv(F1 / "strat_null.csv"); rs = pd.read_csv(F1 / "residual_assort.csv"); rc = pd.read_csv(F1 / "resid_controls_summary.csv")
TOP3 = {"lung": [(3, 1), (6, 2), (8, 3)], "lung_external": [(8, 2), (3, 1), (7, 2)], "kidney_B": [(9, 3), (8, 3), (7, 2)]}
out = {}
for ctx, cells in TOP3.items():
    s = st[st.context == ctx]; q = pd.concat([s[(s.layer == L) & (s.sv == k)] for L, k in cells])
    o = dict(top3_strat_null_mean_range=[q.strat_null_mean.min(), q.strat_null_mean.max()],
             top3_cpms_null_mean_range=[q.cpms_null_mean.min(), q.cpms_null_mean.max()],
             top3_strat_excess_ratio_range=[q.strat_excess_over_ms_excess.min(), q.strat_excess_over_ms_excess.max()],
             top3_cpms_excess_ratio_range=[q.cpms_excess_over_ms_excess.min(), q.cpms_excess_over_ms_excess.max()])
    big = s[s.r > 0.3]
    o["n_cells_r_gt_0.3"] = int(len(big))
    o["lowest_cpms_excess_ratio_among_r_gt_0.3"] = big.nsmallest(4, "cpms_excess_over_ms_excess")[["layer", "sv", "r", "cpms_excess_over_ms_excess"]].to_dict("records")
    o["lowest_strat_excess_ratio_among_r_gt_0.3"] = big.nsmallest(3, "strat_excess_over_ms_excess")[["layer", "sv", "r", "strat_excess_over_ms_excess"]].to_dict("records")
    r30 = rs[(rs.context == ctx) & (rs.design == "resid30cc")]
    b30 = r30[r30.r_raw > 0.3]
    o["lowest_resid30_ratio_among_r_gt_0.3"] = b30.nsmallest(3, "ratio_resid_over_raw")[["layer", "sv", "r_raw", "ratio_resid_over_raw"]].to_dict("records")
    o["resid30_max_r_resid_cell"] = r30.loc[r30.r_resid.idxmax(), ["layer", "sv", "r_resid", "z"]].to_dict()
    o["resid30_max_absz_cell"] = r30.loc[r30.z.abs().idxmax(), ["layer", "sv", "r_resid", "z"]].to_dict()
    rr = rc[(rc.context == ctx) & rc.design.str.startswith("random30")]
    o["random30_median_ratio_range"] = [rr.median_ratio.min(), rr.median_ratio.max()]
    out[ctx] = o
s = st[(st.context == "lung") & (st.layer == 7) & (st.sv == 6)].iloc[0]
q = rs[(rs.context == "lung") & (rs.layer == 7) & (rs.sv == 6)].set_index("design")
out["lung_L7_SV6"] = dict(r=s.r, strat_null_mean=s.strat_null_mean, strat_excess_ratio=s.strat_excess_over_ms_excess,
                          cpms_excess_ratio=s.cpms_excess_over_ms_excess, resid30_ratio=q.loc["resid30cc", "ratio_resid_over_raw"],
                          resid30_R2=q.loc["resid30cc", "R2_coord_on_design"], residclass_ratio=q.loc["residclass", "ratio_resid_over_raw"])
g2go = pickle.load(open(VF.G2GO, "rb")); ns = VF.parse_obo()
for ctx in TOP3:
    genes = json.loads((VF.R3 / "cache" / ctx / "obs_genes.json").read_text())
    cnt = Counter(t for g in genes for t in g2go.get(g, ()) if ns.get(t) == "molecular_function")
    out[ctx]["mf_top3_terms_gene_counts"] = cnt.most_common(3); out[ctx]["n_V"] = len(genes)
(V / "wording_numbers.json").write_text(json.dumps(out, indent=1, default=float))
print(json.dumps(out, indent=1, default=float))
