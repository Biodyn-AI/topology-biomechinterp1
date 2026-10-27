"""E11: collect headline numbers into results.json (every value read from a CSV/JSON written by E11 scripts)."""
import json
import numpy as np, pandas as pd
from e11_common import HERE, CTXS, convention_hash


def rec(df):
    return json.loads(df.to_json(orient="records"))


def main():
    frozen, actual, ok = convention_hash()
    mt = pd.read_csv(HERE / "maxT.csv"); ct = pd.read_csv(HERE / "compartment_tests.csv")
    ax = pd.read_csv(HERE / "axis_summary.csv"); led = json.loads((HERE / "ledger_iter0006.json").read_text())
    fc = pd.read_csv(HERE / "further_candidates.csv"); p1 = json.loads((HERE / "part1_meta.json").read_text())
    p3 = json.loads((HERE / "part3_meta.json").read_text()); ffl = pd.read_csv(HERE / "cand_ii_ffl.csv")
    gobp = pd.read_csv(HERE / "cand_i_gobp_sv2.csv")
    out = dict(experiment="E11_compartment", convention_sha256=frozen, convention_hash_verified=ok,
               analysis_type="sensitivity (CONVENTION3 sec. 5; no pass/fail)")
    # ---- part 1
    v = ct[ct.universe == "V"]
    l11 = v[(v.layer == 11) & (v.direction == "SV1") & (v.term.isin(["GO:0005615", "GO:0005829"]))][
        ["context", "term_name", "auroc_signfree", "p_adj_auroc", "OR_best_pole", "p_best_pole", "best_pole"]]
    out["part1_compartment_maxT"] = dict(
        spec="12 layers x SV1..SV12 (observed-gene SVD) x {GO:0005615, GO:0005829}; poles K=floor(0.27n); two-sided Fisher; "
             "sign-invariant AUROC; max-T over N=2000 joint label permutations; p floor 1/2001",
        n_perm=int(mt.N.iloc[0]), p_floor=float(mt.p_floor.iloc[0]),
        best_cells=rec(mt[["context", "universe", "family", "statistic", "best_term_name", "best_layer", "best_direction",
                           "best_observed", "best_z", "best_p_adj", "n_cells_padj_lt_0p05", "family_size",
                           "maxnull_z_q95"]]),
        axis_summary=rec(ax), L11_SV1=rec(l11), context_meta=p1)
    # ---- part 2
    out["part2_ledger"] = dict(
        own_convention=dict(reproduced=led["own_convention"]["reproduced"],
                            OR=led["own_convention"]["orientations"]["numpy_sign"]["OR"],
                            p=led["own_convention"]["orientations"]["numpy_sign"]["p"],
                            zero_rows_in_extracellular_quartile=led["own_convention"]["orientations"]["numpy_sign"]["zero_rows_in_top"],
                            r_SV1_vs_observed_indicator=led["own_convention"]["r_SV1_vs_observed_indicator_L11"],
                            extracellular_rate_zero_row_named=led["own_convention"]["extracellular_rate_zero_row_named"],
                            extracellular_rate_observed_named=led["own_convention"]["extracellular_rate_observed_named"]),
        submitted_retraction=led["R2_retraction_convention"],
        paper_convention_named_genes=led["paper_convention_named_genes"],
        paper_convention_full_V_SV1=led["paper_convention_full_V_SV1"], source="ledger_iter0006.json; ledger_tests.csv")
    # ---- part 3
    gb = gobp.assign(sig=(gobp.q_top_ctx < 0.05) | (gobp.q_bottom_ctx < 0.05))
    out["part3_candidates"] = dict(
        i_GO_SV2=dict(own_convention=p3["i_own"], paper_convention=p3["i_paper"],
                      distinct_BP_terms_sig_ctx={c: int(gb[(gb.context == c) & gb.sig].term.nunique()) for c in gb.context.unique()}),
        iv_effrank=p3["iv"],
        ii_FFL={c: dict(n_ffl=int(d.n_ffl.iloc[0]),
                        best_label=rec(d.sort_values("p_adj_frac_label").head(1)[["layer", "subspace", "frac_between", "frac_null_mean_label", "p_adj_frac_label"]])[0],
                        best_role=rec(d.sort_values("p_adj_frac_role").head(1)[["layer", "subspace", "frac_between", "frac_null_mean_role", "p_adj_frac_role"]])[0],
                        min_p_adj_mean_t_2s_role=float(d.p_adj_mean_t_2s_role.min()))
                for c, d in ffl.groupby("context")},
        iii_PH=p3["iii"], table="further_candidates.csv")
    # ---- independent verification (verify/*.py, written by the adversarial verifier; own SVD / AUROC / max-T)
    vd = HERE / "verify"
    if vd.exists():
        out["verification"] = {p.stem: json.loads(p.read_text()) for p in sorted(vd.glob("verify_*.json"))}
        out["verification"]["fixes"] = [
            "cand_iv: single-context partial Spearman (ER vs AUROC | layer) set to NaN (not estimable): ER is strictly "
            "monotone in layer (rho = -1.00), so the rank residuals are ~1e-15 floating-point noise; the previously "
            "reported partial rho 0.048 / 0.012 / -0.068 / 0.095 were correlations of noise. Immune-triplet (n = 36) "
            "values unchanged."]
    (HERE / "results.json").write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps(out["part1_compartment_maxT"]["best_cells"][:2], indent=1))


if __name__ == "__main__":
    main()
