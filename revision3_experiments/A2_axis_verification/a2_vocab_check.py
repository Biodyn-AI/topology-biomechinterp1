"""ScGPTDataset raises if any var gene is missing from the vocab; confirm every context's var genes are in the
vocab of the checkpoint it used (so no vocab filtering of the axis is possible)."""
import json, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent)); sys.path.insert(0, str(HERE))
import contexts as X, a2_lib as L
CK = X.REPO / "single_cell_mechinterp/external/scGPT_checkpoints"
use = {k: "whole-human" for k in X.CONTEXTS}; use["kidney_ckpt_kidney"] = "kidney"; use["kidney_ckpt_brain"] = "brain"
out = {}
for ctx, c in X.CONTEXTS.items():
    p = CK / use[ctx] / "vocab.json"
    if not p.exists():
        out[ctx] = dict(checkpoint=use[ctx], vocab_found=False); continue
    v = json.loads(p.read_text())
    keys = set(v.keys()) if isinstance(v, dict) else set(v)
    g = L.var_names(c["h5ad"])
    out[ctx] = dict(checkpoint=use[ctx], vocab_found=True, vocab_size=len(keys), G=len(g), n_missing=sum(x not in keys for x in g))
    print(ctx, out[ctx])
(HERE / "vocab_check.json").write_text(json.dumps(out, indent=1))
