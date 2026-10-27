"""Collect verifier results into verify/verification.json (all numbers read from the v1-v5 outputs and the E4b CSVs)."""
import json, hashlib
from pathlib import Path
import pandas as pd
V = Path(__file__).resolve().parent; E = V.parent; R3 = E.parent
h = [l.split()[-1] for l in (R3 / "FROZEN.txt").read_text().splitlines() if "sha256" in l][0]
z = pd.read_csv(V / "v4_zero_rows.csv")
out = dict(
    convention3_sha256=h, convention3_sha256_verified=hashlib.sha256((R3 / "CONVENTION3.md").read_bytes()).hexdigest() == h,
    v1_inputs=json.loads((V / "v1_inputs.json").read_text()), v2_refit=json.loads((V / "v2_refit.json").read_text()),
    v3_headlines=json.loads((V / "v3_headlines.json").read_text()), v5_extra=json.loads((V / "v5_extra.json").read_text()),
    v4_fitted_zero_importance_rows=dict(total=int(z.n_fitted_zero_importance.sum()),
                                        genie3=int(z[z.method == "genie3"].n_fitted_zero_importance.sum()),
                                        grnboost2=int(z[z.method == "grnboost2"].n_fitted_zero_importance.sum()),
                                        max_per_fit=int(z.n_fitted_zero_importance.max()),
                                        max_frac_h3_pairs_in_zero_rows=float(z.frac_h3_pairs_in_zero_rows.max())),
    bug_changing_results=False,
    note="No re-run of E4b outputs was needed; corrections are to SUMMARY wording/rounding and added caveats.")
(V / "verification.json").write_text(json.dumps(out, indent=1, default=float))
print(out["convention3_sha256_verified"], out["v4_fitted_zero_importance_rows"])
