"""Probe the kidney_A / kidney_B drawn-cell overlap excess (12 vs 4.39 expected)."""
import sys, json
from pathlib import Path
import numpy as np, h5py
from scipy.stats import hypergeom
R3 = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(R3))
import contexts as X
def names(p):
    with h5py.File(p, "r") as f:
        k = f["obs"].attrs.get("_index", "_index"); r = f["obs"][k][:]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in r])
nA, nB = names(X.H5_KIDNEY_A), names(X.H5_KIDNEY_B)
posA = {n: i for i, n in enumerate(nA)}
mapB = np.array([posA[n] for n in nB])
out = dict(B_rows_monotone_in_A=bool(np.all(np.diff(mapB) > 0)))
cA = np.load(R3/"cache/kidney_A/extraction_cells.npy"); cB = np.load(R3/"cache/kidney_B/extraction_cells.npy")
shared = sorted(set(cA) & set(mapB[cB]))
out["shared_A_rows"] = [int(x) for x in shared]
out["shared_B_rows"] = [int(np.flatnonzero(mapB == s)[0]) for s in shared]
# null by simulation: other seed pairs, same sizes
sims = []
for s in range(1000):
    a = np.sort(np.random.default_rng(10_000 + s).choice(11376, 256, replace=False))
    b = mapB[np.sort(np.random.default_rng(20_000 + s).choice(8000, 192, replace=False))]
    sims.append(len(set(a) & set(b)))
sims = np.array(sims)
out["sim_mean_overlap"] = float(sims.mean()); out["sim_frac_ge_12"] = float((sims >= 12).mean())
# does seed 42 x seed 67 specifically correlate? overlap for seed-42 A draw vs B draws with seeds 0..2000
ov = np.array([len(set(cA) & set(mapB[np.random.default_rng(s).choice(8000, 192, replace=False)])) for s in range(2001)])
out["A42_vs_Bseed_0_2000_mean"] = float(ov.mean()); out["rank_of_seed67_overlap"] = int((ov >= ov[67]).sum())
# lower-tail p for the zero overlaps A2 did not flag
out["p_lower_seed43_seed44_0_shared"] = float(hypergeom.cdf(0, 20000, 256, 256))
out["p_lower_immune_Bcell_0_shared"] = float(hypergeom.cdf(0, 3762, 50, 256))
print(json.dumps(out, indent=1)); (Path(__file__).parent/"kidney_overlap_probe.json").write_text(json.dumps(out, indent=1))
