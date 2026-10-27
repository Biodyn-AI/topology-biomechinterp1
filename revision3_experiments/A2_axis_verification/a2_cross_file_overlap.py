"""Do different-tissue h5ads share physical cells (TS immune subset vs TS lung / TS kidney)? obs_names and
observation_joinid overlaps, whole file and drawn cells. Merged into results.json['cross_file_overlap']."""
import json, sys, itertools
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent)); sys.path.insert(0, str(HERE))
import contexts as X, a2_lib as L
CACHE = HERE.parent / "cache"
F = {"immune": X.H5_IMMUNE, "kidney_A": X.H5_KIDNEY_A, "kidney_B": X.H5_KIDNEY_B, "lung": X.H5_LUNG, "lung_external": X.H5_EXT_LUNG}
names = {k: L.obs_names(h) for k, h in F.items()}
joinid = {k: L.obs_col(h, "observation_joinid") for k, h in F.items()}
drawn = {k: np.load(CACHE / k / "extraction_cells.npy") for k in F}
out = []
for a, b in itertools.combinations(F, 2):
    r = dict(a=a, b=b, file_shared_obs_names=len(set(names[a]) & set(names[b])),
             file_shared_joinid=len(set(joinid[a]) & set(joinid[b])),
             drawn_shared_obs_names=len(set(names[a][drawn[a]]) & set(names[b][drawn[b]])),
             drawn_shared_joinid=len(set(joinid[a][drawn[a]]) & set(joinid[b][drawn[b]])))
    out.append(r); print(r)
res = json.loads((HERE / "results.json").read_text()); res["cross_file_overlap"] = out
(HERE / "results.json").write_text(json.dumps(res, indent=1, default=str))
