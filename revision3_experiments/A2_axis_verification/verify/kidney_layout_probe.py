"""How do kidney_B file rows map to kidney_A file rows? (prefix identity, Smart-seq2 placement)"""
import sys, json
from pathlib import Path
import numpy as np, h5py
R3 = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(R3))
import contexts as X
def col(p, c):
    with h5py.File(p, "r") as f:
        g = f["obs"][c]
        if isinstance(g, h5py.Group):
            cats = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]]); return cats[g["codes"][:]]
        return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g[:]])
def names(p):
    with h5py.File(p, "r") as f:
        k = f["obs"].attrs.get("_index", "_index"); r = f["obs"][k][:]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in r])
nA, nB = names(X.H5_KIDNEY_A), names(X.H5_KIDNEY_B)
posA = {n: i for i, n in enumerate(nA)}; mapB = np.array([posA[n] for n in nB])
eq = mapB == np.arange(8000)
first_neq = int(np.argmin(eq)) if not eq.all() else 8000
aA = col(X.H5_KIDNEY_A, "assay")
out = dict(n_rows_B_equal_A_position=int(eq.sum()), first_B_row_not_equal=first_neq,
           mapB_head_after_break=[int(x) for x in mapB[first_neq:first_neq + 5]],
           A_smartseq_rows_min_max=[int(np.flatnonzero(aA == "Smart-seq2").min()), int(np.flatnonzero(aA == "Smart-seq2").max())],
           expected_positional_coincidences=float(eq.sum() * 256 / 11376 * 192 / 8000))
print(json.dumps(out, indent=1)); (Path(__file__).parent / "kidney_layout_probe.json").write_text(json.dumps(out, indent=1))
