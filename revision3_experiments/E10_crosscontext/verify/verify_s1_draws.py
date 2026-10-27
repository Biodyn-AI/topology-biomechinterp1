"""V11: rewiring-draw distribution of S1 (layer mean) with the stats3 Maslov-Sneppen chain (E10's generator) over 20
seeds, to separate draw noise from implementation differences against the networkx draws of V3. Writes ck_v11.json."""
import json, sys
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score
HERE = Path(__file__).resolve().parent; E10 = HERE.parent
sys.path.insert(0, str(E10))
import e10_common as K
out = {}
for ctx in ["immune", "kidney_A", "lung"]:
    c = K.load_cache(ctx); emb = K.load_emb(ctx); n = len(c["genes"]); ed = np.load(c["d"] / "edges_string.npy")
    Ens = [K.unit(K.centred(emb, L, c["obs_idx"])) for L in range(12)]
    vals = []
    for s in range(20 if ctx != "lung" else 8):
        epos, eneg, _ = K.rewired_nonpairs(ed, n, seed=777 + s)
        pa = np.r_[epos[:, 0], eneg[:, 0]]; pb = np.r_[epos[:, 1], eneg[:, 1]]; y = np.r_[np.ones(len(epos)), np.zeros(len(eneg))]
        vals.append(float(np.mean([roc_auc_score(y, (E[pa] * E[pb]).sum(1)) for E in Ens])))
    e10 = float(np.load(E10 / "partial" / f"stats_{ctx}.npz")["point"][0].mean())
    out[ctx] = dict(draws=vals, mean=float(np.mean(vals)), sd=float(np.std(vals, ddof=1)), e10_single=e10,
                    z_e10=(e10 - float(np.mean(vals))) / float(np.std(vals, ddof=1)))
    print(ctx, {k: v for k, v in out[ctx].items() if k != "draws"}, flush=True)
(HERE / "ck_v11.json").write_text(json.dumps(out, indent=1))
