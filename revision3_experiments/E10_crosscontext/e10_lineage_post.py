"""Post-process lineage_tfs.csv: centred bootstrap intervals are clipped to [0, 1] and flagged unreliable when
fewer than 10 lineage TFs are observed (the leave-gene-out bootstrap of a 5-8 gene reference set is strongly
skewed; see SUMMARY). Adds columns lo_clip, hi_clip, ci_reliable. Idempotent."""
from pathlib import Path
import numpy as np, pandas as pd
H = Path(__file__).resolve().parent
d = pd.read_csv(H / "lineage_tfs.csv", low_memory=False)
ev = d.part == "eval"
d.loc[ev, "lo_clip"] = d.loc[ev, "lo"].clip(0, 1); d.loc[ev, "hi_clip"] = d.loc[ev, "hi"].clip(0, 1)
d.loc[ev, "ci_reliable"] = d.loc[ev, "n_observed"] >= 10
d.to_csv(H / "lineage_tfs.csv", index=False)
print(d[ev & (d.layer.astype(str) == "mean")][["context", "tissue_set", "n_observed", "auroc", "lo_clip", "hi_clip", "ci_reliable"]].to_string())
