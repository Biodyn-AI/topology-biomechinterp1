"""Independent verifier, part 3 (own code; no stats3/e4a imports): the borderline 'not above co-expression' label in
kidney_ckpt_kidney. Recomputes, from raw inputs, the L2 signed-cosine AUROC, the tie-corrected |co-expression|
(extraction cells) AUROC and the one-sided 95% upper bound of their paired delta with a two-way TF x target
pigeonhole bootstrap, B = 2,000: (a) replaying the implementer's weight stream (seed [seed_of(ctx), 2]) and
(b) two independent seeds -> Monte-Carlo spread of the bound. Writes verify/verify_part3.json."""
import json, os, time, zlib
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; E4 = HERE.parent; R3 = E4.parent; CACHE = R3 / "cache"
EMB = Path(os.path.expanduser("~/biomechinterp-local/subproject41_r3/emb"))
ctx = "kidney_ckpt_kidney"; LAYER = 2


def seed_of(c): return 20261026 + zlib.crc32(c.encode()) % 100000


def main():
    t0 = time.time(); d = CACHE / ctx; h = np.load(d / "h3_pairs.npz"); tf, tg, y = h["tf"], h["tg"], h["y"].astype(bool)
    ok = np.array([y[tf == t].any() and (~y[tf == t]).any() for t in np.unique(tf)]); k = np.isin(tf, np.unique(tf)[ok])
    tf, tg, y = tf[k], tg[k], y[k]
    emb = np.load(EMB / f"{ctx}.npy", mmap_mode="r"); oi = np.load(d / "obs_idx.npy")
    X = np.asarray(emb[LAYER], np.float64)[oi]; X = X - X.mean(0); X /= np.linalg.norm(X, axis=1, keepdims=True)
    s = np.einsum("ij,ij->i", X[tf], X[tg])
    M = np.load(d / "expr_extraction.npy").astype(np.float64); mu = M.mean(0); sd = M.std(0)
    c = np.round(np.abs(((M[:, tf] - mu[tf]) * (M[:, tg] - mu[tg])).mean(0) / (sd[tf] * sd[tg])), 12)
    tfs, tfc = np.unique(tf, return_inverse=True); tgs, tgc = np.unique(tg, return_inverse=True)
    blocks = [np.flatnonzero(tfc == j) for j in range(tfs.size)]
    pre = []
    for b in blocks:   # per TF: positives, negatives sorted per scorer
        yp = y[b]
        per = []
        for sc in (s, c):
            sn = sc[b][~yp]; o = np.argsort(sn, kind="mergesort"); sp_ = sc[b][yp]
            lo = np.searchsorted(sn[o], sp_, "left"); hi = np.searchsorted(sn[o], sp_, "right")
            per.append((o, lo, hi))
        pre.append((b[yp], b[~yp], per))

    def auc_pair(w):
        out = []
        for j in range(2):
            num = den = 0.0
            for bp, bn, per in pre:
                o, lo, hi = per[j]; wn = w[bn][o]; wp = w[bp]
                if wp.sum() == 0 or wn.sum() == 0: continue
                cw = np.r_[0.0, np.cumsum(wn)]
                num += (wp * (cw[lo] + 0.5 * (cw[hi] - cw[lo]))).sum(); den += wp.sum() * wn.sum()
            out.append(num / den)
        return out
    pt = auc_pair(np.ones(tf.size)); R = dict(context=ctx, layer=LAYER, auroc_L=pt[0], auroc_cx_ext_tiecorrected=pt[1], delta=pt[0] - pt[1])
    for tag, rng in (("replay_seed", np.random.default_rng([seed_of(ctx), 2])), ("indep_seed_1", np.random.default_rng(111)),
                     ("indep_seed_2", np.random.default_rng(222))):
        bs = np.empty((2000, 2))
        for r in range(2000):
            wt = np.bincount(rng.integers(0, tfs.size, tfs.size), minlength=tfs.size)
            wg = np.bincount(rng.integers(0, tgs.size, tgs.size), minlength=tgs.size)
            bs[r] = auc_pair((wt[tfc] * wg[tgc]).astype(float))
        dl = bs[:, 0] - bs[:, 1]
        R[tag] = dict(ub95_auroc=float(np.quantile(bs[:, 0], .95)), delta_ub95=float(np.quantile(dl, .95)), delta_sd=float(dl.std(ddof=1)),
                      cx_ci=[float(np.quantile(bs[:, 1], .025)), float(np.quantile(bs[:, 1], .975))])
        print(tag, R[tag], f"{time.time()-t0:.0f}s", flush=True)
    (HERE / "verify_part3.json").write_text(json.dumps(R, indent=1, default=float))


if __name__ == "__main__":
    main()
