"""V10: is the E10 Poisson node-bootstrap SE of S1 (layer mean) consistent with CONVENTION3 §6's network-statistic
interval (vertex subsampling m = n/2) under the delete-d scaling of the errata (sqrt(m/(n-m)) = 1)?
Uses E10's own non-pair set (same seed) so only the resampling scheme differs. Writes verify/ck_v10.json."""
import json, sys
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score
HERE = Path(__file__).resolve().parent; E10 = HERE.parent
sys.path.insert(0, str(E10))
import e10_common as K

out = {}
for ctx in ["immune", "kidney_A", "lung"]:
    c = K.load_cache(ctx); emb = K.load_emb(ctx); n = len(c["genes"])
    epos, eneg, _ = K.rewired_nonpairs(np.load(c["d"] / "edges_string.npy"), n, seed=K.SEED + K.CTXS.index(ctx))
    pa = np.r_[epos[:, 0], eneg[:, 0]]; pb = np.r_[epos[:, 1], eneg[:, 1]]; y = np.r_[np.ones(len(epos)), np.zeros(len(eneg))]
    S = np.array([(K.unit(K.centred(emb, L, c["obs_idx"]))[pa] * K.unit(K.centred(emb, L, c["obs_idx"]))[pb]).sum(1) for L in range(12)])
    rng = np.random.default_rng(5); m = n // 2; draws = []
    for r in range(300):
        keep = np.zeros(n, bool); keep[rng.choice(n, m, replace=False)] = True; e = keep[pa] & keep[pb]
        draws.append(np.mean([roc_auc_score(y[e], S[L, e]) for L in range(12)]))
    z = np.load(E10 / "partial" / f"stats_{ctx}.npz")
    se_boot = float(np.nanstd(z["boot"][0].mean(0), ddof=1))
    sd = float(np.std(draws, ddof=1))
    out[ctx] = dict(se_subsample_deleted=sd * np.sqrt(m / (n - m)), se_subsample_sqrt_m_over_n=sd * np.sqrt(m / n),
                    se_e10_poisson_node_bootstrap=se_boot, ratio_e10_over_deleted=se_boot / (sd * np.sqrt(m / (n - m))))
    print(ctx, out[ctx], flush=True)
(HERE / "ck_v10.json").write_text(json.dumps(out, indent=1))
