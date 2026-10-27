"""Independent re-computation of candidate (iv) in lung and immune: effective rank (exp entropy of S/sum S) of the
observed-gene centred matrix, H2 statistic (LOO mean cosine to Lambert TFs, AUROC pooled within detection deciles),
raw Spearman and partial Spearman on layer. Labels / detection from cache node_labels (shared, read-only).
Output: verify/verify_iv.json"""
import json
from pathlib import Path
import numpy as np
from scipy.stats import rankdata, spearmanr, pearsonr

HERE = Path(__file__).resolve().parent; R3 = HERE.parents[1]
EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"


def strat_auc(s, y, st):
    num = den = 0.0
    for v in np.unique(st):
        m = st == v; p, q = s[m & y], s[m & ~y]
        if p.size and q.size:
            r = rankdata(np.r_[p, q]); num += r[:p.size].sum() - p.size * (p.size + 1) / 2; den += p.size * q.size
    return num / den


out = {}
for ctx in ("lung", "immune"):
    E = np.load(EMB / f"{ctx}.npy", mmap_mode="r"); idx = np.load(R3 / "cache" / ctx / "obs_idx.npy")
    nl = np.load(R3 / "cache" / ctx / "node_labels.npz"); y = nl["is_tf_lambert"].astype(bool); det = nl["det_decile"]
    er, h2 = [], []
    for L in range(12):
        Xl = np.asarray(E[L], np.float64)[idx]; Xl -= Xl.mean(0)
        s = np.linalg.svd(Xl, compute_uv=False); p = s / s.sum(); er.append(float(np.exp(-(p * np.log(p)).sum())))
        En = Xl / np.linalg.norm(Xl, axis=1, keepdims=True); tot = En[y].sum(0); sc = En @ tot
        sc = np.where(y, (sc - 1.0) / (y.sum() - 1), sc / y.sum())
        h2.append(strat_auc(sc, y, det))
    lay = np.arange(12.)
    rx, ry, rz = rankdata(er), rankdata(h2), rankdata(lay); A = np.c_[np.ones(12), rz]
    ex = rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0]; ey = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]
    pr = pearsonr(ex, ey)
    out[ctx] = dict(eff_rank=er, h2=h2, rho_layer_er=float(spearmanr(lay, er).correlation),
                    rho_raw=float(spearmanr(er, h2).correlation), partial_rho=float(pr[0]), partial_p=float(pr[1]))
    print(ctx, out[ctx]["rho_raw"], out[ctx]["partial_rho"], out[ctx]["partial_p"], flush=True)
(HERE / "verify_iv.json").write_text(json.dumps(out, indent=1))
