"""Merge token-mask and vocab summaries into results.json (run order: a2_axis_verification.py,
a2_supplement.py, a2_token_masks.py, a2_vocab_check.py, a2_merge.py, a2_cross_file_overlap.py)."""
import json
from pathlib import Path
H = Path(__file__).resolve().parent
r = json.loads((H / "results.json").read_text())
r["token_mask_summary"] = json.loads((H / "token_mask_summary.json").read_text())
r["vocab_check"] = json.loads((H / "vocab_check.json").read_text())
r["run_order"] = ["a2_axis_verification.py", "a2_supplement.py", "a2_token_masks.py", "a2_vocab_check.py", "a2_merge.py", "a2_cross_file_overlap.py"]
(H / "results.json").write_text(json.dumps(r, indent=1, default=str))
print(sorted(r.keys()))
