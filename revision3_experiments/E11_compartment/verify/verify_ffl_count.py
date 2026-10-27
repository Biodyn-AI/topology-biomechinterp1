"""Independent FFL count from the raw TRRUST file (dedup to distinct TF->target pairs) restricted to observed genes."""
import json
from pathlib import Path
from collections import defaultdict
import pandas as pd
HERE = Path(__file__).resolve().parent; R3 = HERE.parents[1]; REPO = R3.parents[1]
t = pd.read_csv(REPO / "single_cell_mechinterp/external/networks/trrust_human.tsv", sep="\t", header=None)
pairs = {(a, b) for a, b in zip(t[0].astype(str), t[1].astype(str)) if a != b}
out = {}
for ctx in ("immune", "lung", "lung_external"):
    V = set(json.loads((R3 / "cache" / ctx / "obs_genes.json").read_text()))
    o = defaultdict(set)
    for a, b in pairs:
        if a in V and b in V: o[a].add(b)
    n = sum(1 for a in o for b in o[a] if b in o for c in (o[a] & o[b]) if c not in (a, b))
    out[ctx] = dict(n_trrust_pairs_in_V=sum(len(v) for v in o.values()), n_ffl=n)
print(out); (HERE / "verify_ffl_count.json").write_text(json.dumps(out, indent=1))
