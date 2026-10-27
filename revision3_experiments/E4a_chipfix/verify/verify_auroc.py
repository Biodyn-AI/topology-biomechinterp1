"""Independent verifier, part 2: lung ChIP AUROCs recomputed from scratch.

No import of stats3 / e4a_lib / e4a_main / e4a_chipfix. Scores, the within-TF stratified AUROC (Mann-Whitney
pooled over TFs, ties 1/2, via scipy rankdata), the two-way TF x target pigeonhole bootstrap (same integer
stream as stats3.StratifiedPairs.bootstrap: per replicate rng.integers(0,nTF,nTF) then rng.integers(0,nTG,nTG),
rng = default_rng([seed_of(ctx), 40]), seed_of = 20261026 + crc32(ctx) % 100000) and percentile CI/UB are written
here. Targets: lung L0 and L11 signed cosine (V-centred full 512-d), |co-expression| on extraction cells and on the
3,000-cell tier (|r| rounded to 12 decimals), on the FIXED pairs and, as code-path check, on the OLD pairs against
E4a_edges/sens_chip.csv. Also lung_external L0 / |co-expr| extraction point estimates.
"""
from __future__ import annotations
import json, os, sys, time, zlib
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent; FIX = HERE.parent; R3 = FIX.parent
EMB = Path(os.path.expanduser("~/biomechinterp-local/subproject41_r3/emb"))
B = 2000


def seed_of(ctx): return 20261026 + zlib.crc32(ctx.encode()) % 100000


def filt(tf, tg, y):
    """keep TFs with >=1 positive and >=1 candidate (the stats3 container rule)."""
    y = y.astype(bool); df = pd.DataFrame(dict(tf=tf, y=y))
    g = df.groupby("tf").y.agg(["any", "all"]); ok = g.index[g["any"] & ~g["all"]]
    m = np.isin(tf, ok); return m


def auroc_strat(s, tf, y):
    y = y.astype(bool); num = den = 0.0
    order = np.argsort(tf, kind="stable"); tfs, start = np.unique(tf[order], return_index=True)
    bounds = list(start) + [tf.size]
    for i in range(tfs.size):
        ix = order[bounds[i]:bounds[i + 1]]; ss = s[ix]; yy = y[ix]
        n1 = yy.sum(); n0 = yy.size - n1
        if n1 == 0 or n0 == 0: continue
        r = rankdata(ss)                       # mid-ranks: ties count 1/2
        num += r[yy].sum() - n1 * (n1 + 1) / 2; den += n1 * n0
    return num / den


class Boot:
    def __init__(self, tf, tg, y, s_list):
        self.y = y.astype(bool)
        _, self.tfc = np.unique(tf, return_inverse=True); _, self.tgc = np.unique(tg, return_inverse=True)
        self.ntf = self.tfc.max() + 1; self.ntg = self.tgc.max() + 1
        self.prep = [self._prep(s) for s in s_list]

    def _prep(self, s):
        o = np.lexsort((s, self.tfc)); st = self.tfc[o]; ss = s[o]; n = s.size
        newg = np.r_[True, (st[1:] != st[:-1]) | (ss[1:] != ss[:-1])]; gid = np.cumsum(newg) - 1
        gs = np.flatnonzero(newg); ge = np.r_[gs[1:], n]
        news = np.r_[True, st[1:] != st[:-1]]; segs = np.flatnonzero(news); sid = np.cumsum(news) - 1
        return dict(o=o, gs=gs[gid], ge=ge[gid], s0=segs[sid], segs=segs, ys=self.y[o])

    def auc(self, k, w):
        p = self.prep[k]; ws = w[p["o"]]; ys = p["ys"]; nw = np.where(ys, 0.0, ws)
        C = np.r_[0.0, np.cumsum(nw)]
        below = C[p["gs"]] - C[p["s0"]]; tied = C[p["ge"]] - C[p["gs"]]
        num = np.sum(np.where(ys, ws * (below + 0.5 * tied), 0.0))
        P = np.add.reduceat(np.where(ys, ws, 0.0), p["segs"]); N = np.add.reduceat(nw, p["segs"])
        den = np.sum(P * N)
        return num / den if den > 0 else np.nan

    def run(self, B, rng):
        out = np.empty((B, len(self.prep)))
        for r in range(B):
            wtf = np.bincount(rng.integers(0, self.ntf, self.ntf), minlength=self.ntf)
            wtg = np.bincount(rng.integers(0, self.ntg, self.ntg), minlength=self.ntg)
            w = wtf[self.tfc].astype(float) * wtg[self.tgc]
            for k in range(len(self.prep)): out[r, k] = self.auc(k, w)
        return out


def brute_weighted(s, tf, y, w):
    """slow reference for a weighted stratified AUROC (pairwise per TF, ties 1/2)."""
    y = y.astype(bool); num = den = 0.0
    for t in np.unique(tf):
        m = tf == t; sp, wp = s[m & y], w[m & y]; sn, wn = s[m & ~y], w[m & ~y]
        if wp.sum() == 0 or wn.sum() == 0: continue
        o = np.argsort(sn); sn, wn = sn[o], wn[o]; cw = np.r_[0.0, np.cumsum(wn)]
        lo = np.searchsorted(sn, sp, "left"); hi = np.searchsorted(sn, sp, "right")
        num += np.sum(wp * (cw[lo] + 0.5 * (cw[hi] - cw[lo]))); den += wp.sum() * wn.sum()
    return num / den


def scores(ctx, tf, tg, layers):
    d = R3 / "cache" / ctx; obs = np.load(d / "obs_idx.npy"); E = np.load(EMB / f"{ctx}.npy", mmap_mode="r"); out = {}
    for l in layers:
        Xl = np.asarray(E[l], np.float64)[obs]; Xl = Xl - Xl.mean(0); nrm = np.sqrt((Xl ** 2).sum(1))
        out[f"L{l:02d}"] = (Xl[tf] * Xl[tg]).sum(1) / (nrm[tf] * nrm[tg])
    for key, f in (("cx_ext", "expr_extraction.npy"), ("cx_3k", "expr.npy")):
        Mx = np.load(d / f).astype(np.float64); Z = Mx - Mx.mean(0); sd = np.sqrt((Z ** 2).mean(0))
        Z = Z / np.where(sd > 0, sd, np.inf)
        r = np.empty(tf.size)
        for i in range(0, tf.size, 50000):
            a, b = tf[i:i + 50000], tg[i:i + 50000]; r[i:i + 50000] = (Z[:, a] * Z[:, b]).mean(0)
        out[key + "_raw"] = np.abs(r); out[key] = np.round(np.abs(r), 12)
    return out


def ci(b): return dict(ci_lo=float(np.nanquantile(b, .025)), ci_hi=float(np.nanquantile(b, .975)), ub95=float(np.nanquantile(b, .95)), sd=float(np.nanstd(b, ddof=1)))


def evaluate(ctx, pairs_file, ref_csv, layers, do_boot, tag):
    z = dict(np.load(pairs_file)); m = filt(z["tf"], z["tg"], z["y"])
    tf, tg, y = z["tf"][m], z["tg"][m], z["y"][m].astype(bool)
    sc = scores(ctx, tf, tg, layers)
    keys = ["cx_ext", "cx_3k"] + [f"L{l:02d}" for l in layers]
    est = {k: auroc_strat(sc[k], tf, y) for k in keys}
    est_raw = {k: auroc_strat(sc[k + "_raw"], tf, y) for k in ("cx_ext", "cx_3k")}
    ref = pd.read_csv(ref_csv); ref = ref[ref.context == ctx].reset_index(drop=True)
    def refrow(k):
        if k == "cx_ext": return ref.iloc[0]
        if k == "cx_3k": return ref.iloc[1]
        return ref[(ref.scorer == "signed cosine full-512") & (ref.layer == int(k[1:]))].iloc[0]
    rows = []
    for k in keys:
        rr = refrow(k)
        rows.append(dict(context=ctx, pairs=tag, scorer=k, n_pos=int(y.sum()), n_neg=int((~y).sum()), n_tf=int(np.unique(tf).size),
                         ref_n_pos=int(rr.n_pos), ref_n_neg=int(rr.n_neg), mine_auroc=est[k], ref_auroc=float(rr.auroc),
                         absdiff_auroc=abs(est[k] - float(rr.auroc)),
                         mine_auroc_unrounded_cx=est_raw.get(k, np.nan)))
    if do_boot:
        bt = Boot(tf, tg, y, [sc[k] for k in keys]); t0 = time.time()
        # brute check of the weighted evaluator on 2 replicates (fresh stream; same draw rule)
        rr_ = np.random.default_rng(7)
        for _ in range(2):
            w = np.bincount(rr_.integers(0, bt.ntf, bt.ntf), minlength=bt.ntf)[bt.tfc].astype(float) * \
                np.bincount(rr_.integers(0, bt.ntg, bt.ntg), minlength=bt.ntg)[bt.tgc]
            for j in (0, 2):
                assert abs(bt.auc(j, w) - brute_weighted(sc[keys[j]], tf, y, w)) < 1e-12
        bb = bt.run(B, np.random.default_rng([seed_of(ctx), 40])); print(ctx, tag, f"boot {time.time()-t0:.0f}s", flush=True)
        for i, k in enumerate(keys):
            rr = refrow(k); c = ci(bb[:, i]); row = rows[i]
            row.update({f"mine_{q}": v for q, v in c.items()})
            row.update(ref_ci_lo=float(rr.ci_lo), ref_ci_hi=float(rr.ci_hi), ref_ub95=float(rr.ub95), ref_sd=float(rr.boot_sd))
            row["maxdiff_ci"] = max(abs(c["ci_lo"] - rr.ci_lo), abs(c["ci_hi"] - rr.ci_hi), abs(c["ub95"] - rr.ub95), abs(c["sd"] - rr.boot_sd))
            if k.startswith("L"):
                for ck, rk in (("cx_ext", "delta_cxext"), ("cx_3k", "delta_cx3k")):
                    dl = bb[:, i] - bb[:, keys.index(ck)]; cd = ci(dl)
                    row[f"mine_{rk}"] = est[k] - est[ck]; row[f"ref_{rk}"] = float(rr[rk])
                    row[f"mine_{rk}_ci_lo"] = cd["ci_lo"]; row[f"mine_{rk}_ub95"] = cd["ub95"]
                    row[f"ref_{rk}_ci_lo"] = float(rr[rk + "_ci_lo"]); row[f"ref_{rk}_ub95"] = float(rr[rk + "_ub95"])
                    row[f"maxdiff_{rk}"] = max(abs(row[f"mine_{rk}"] - row[f"ref_{rk}"]), abs(cd["ci_lo"] - rr[rk + "_ci_lo"]), abs(cd["ub95"] - rr[rk + "_ub95"]))
    return rows


if __name__ == "__main__":
    E4 = R3 / "E4a_edges"; rows = []
    rows += evaluate("lung", FIX / "cache/lung/h3_chip_pairs_fixed.npz", FIX / "sens_chip_fixed.csv", [0, 11], True, "fixed")
    rows += evaluate("lung", R3 / "cache/lung/h3_chip_pairs.npz", E4 / "sens_chip.csv", [0, 11], True, "old")
    rows += evaluate("lung_external", FIX / "cache/lung_external/h3_chip_pairs_fixed.npz", FIX / "sens_chip_fixed.csv", [0, 11], False, "fixed")
    df = pd.DataFrame(rows); df.to_csv(HERE / "verify_auroc.csv", index=False)
    pd.set_option("display.width", 250); print(df.T.to_string())
