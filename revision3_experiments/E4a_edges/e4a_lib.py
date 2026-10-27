"""E4a helpers: vectorised, output-identical evaluators of the stats3 within-TF stratified AUROC.

Why this file exists. stats3.StratifiedPairs.auroc loops over TFs in Python (~32 ms per call on lung's
52,426 pairs), stats3.StratifiedPairs.bootstrap calls it once per replicate, and permute_within_tf scans
the whole pair list once per TF. The pre-registered H3 analysis needs 10,000 joint permutations x 12
layers and 2,000 two-way bootstrap replicates for ~200 scorers per context; with the stats3 loops that is
days of CPU on a shared 10-core machine. Nothing here defines a new statistic:

  * FastStrat.auc(s, W, y)    == StratifiedPairs.auroc(s, w) for every row w of W (same definition:
                                 Mann-Whitney pooled over TFs, pair weights, ties 1/2). Exact in float64
                                 because bootstrap weights are integers.
  * FastStrat.boot_weights()  draws the SAME integers in the SAME order as StratifiedPairs.bootstrap
                                 (rng.integers(0, nTF, nTF) then rng.integers(0, nTG, nTG) per replicate),
                                 so with the same seed it reproduces stats3's bootstrap replicate for
                                 replicate. Reusing one weight matrix across scorers = "same resamples".
  * FastStrat.perm_null()     within-TF label permutation (counts per TF preserved), the null of
                                 StratifiedPairs.permute_within_tf, drawn by a vectorised generator
                                 (random keys sorted within TF blocks). The permutation distribution is
                                 identical; the random stream differs from permute_within_tf's.
  * verify_against_stats3()   checks all three claims numerically on the real pairs of each context;
                                 the result is written to every results.json.

Observed point estimates, macro AUROC, pooled AUROC and max-T adjustment always call stats3 directly.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import stats3 as S  # noqa: E402


class FastStrat:
    """Wraps a stats3.StratifiedPairs (so filtering and codes are stats3's own).

    strata: 'tf' (within-TF stratified, the primary statistic) or 'pooled' (one stratum: ordinary
    pooled AUROC, still with TF x target cluster weights in the bootstrap)."""

    def __init__(self, sp: S.StratifiedPairs, strata: str = "tf"):
        self.sp = sp
        self.y = sp.y.astype(bool)
        self.n = self.y.size
        self.tf_code = sp.tf_code
        self.tg_code = sp.tg_code
        self.ntf = sp.tfs.size
        self.ntg = sp.tgs.size
        self.strata = sp.tf_code.copy() if strata == "tf" else np.zeros(self.n, np.int64)
        self._prep_cache = {}

    # ------------------------------------------------------------------ ordering / ties
    def _prep(self, s):
        s = np.asarray(s, np.float64)
        order = np.lexsort((s, self.strata))
        ss = s[order]; st = self.strata[order]
        newg = np.r_[True, (st[1:] != st[:-1]) | (ss[1:] != ss[:-1])]
        gid = np.cumsum(newg) - 1
        gstart = np.flatnonzero(newg); gend = np.r_[gstart[1:], self.n]
        newseg = np.r_[True, st[1:] != st[:-1]]
        segs = np.flatnonzero(newseg); seg_of = np.cumsum(newseg) - 1
        return dict(order=order, gs=gstart[gid], ge=gend[gid], s0=segs[seg_of], segs=segs)

    # ------------------------------------------------------------------ weighted AUROC, many weight rows
    def auc(self, s, wts=None, y=None, batch=None):
        """Stratified AUROC for each bootstrap replicate. wts = (WTF (B x nTF), WTG (B x nTG) or None) as
        returned by boot_weights; pair weight = WTF[tf] * WTG[target]. wts=None -> unit weights (point
        estimate). Weights are built per batch, never as a full B x n matrix."""
        p = self._prep(s)
        y = self.y if y is None else np.asarray(y, bool)
        ys = y[p["order"]]
        posi = np.flatnonzero(ys)
        gs, ge, s0 = p["gs"][posi], p["ge"][posi], p["s0"][posi]
        tfo = self.tf_code[p["order"]]; tgo = self.tg_code[p["order"]]
        if wts is None:
            WTF = np.ones((1, self.ntf)); WTG = None
        else:
            WTF, WTG = wts
        Bn = WTF.shape[0]
        out = np.empty(Bn)
        if batch is None:
            batch = max(1, int(4e6 // self.n))
        for a in range(0, Bn, batch):
            Ws = WTF[a:a + batch][:, tfo].astype(np.float64)
            if WTG is not None:
                Ws *= WTG[a:a + batch][:, tgo]
            Nw = Ws * (~ys)
            C = np.zeros((Ws.shape[0], self.n + 1)); np.cumsum(Nw, axis=1, out=C[:, 1:])
            val = (C[:, gs] - C[:, s0]) + 0.5 * (C[:, ge] - C[:, gs])
            num = (Ws[:, posi] * val).sum(1)
            P = np.add.reduceat(Ws * ys, p["segs"], axis=1)
            N = np.add.reduceat(Nw, p["segs"], axis=1)
            den = (P * N).sum(1)
            with np.errstate(invalid="ignore", divide="ignore"):
                out[a:a + batch] = np.where(den > 0, num / den, np.nan)
        return out

    # ------------------------------------------------------------------ bootstrap weights (stats3 stream)
    def boot_weights(self, B, rng, two_way=True):
        """(WTF, WTG): TF and target multiplicities drawn exactly as StratifiedPairs.bootstrap draws them
        (per replicate: rng.integers(0, nTF, nTF), then rng.integers(0, nTG, nTG) if two_way)."""
        WTF = np.empty((B, self.ntf), np.int32); WTG = np.empty((B, self.ntg), np.int32) if two_way else None
        for r in range(B):
            WTF[r] = np.bincount(rng.integers(0, self.ntf, self.ntf), minlength=self.ntf)
            if two_way:
                WTG[r] = np.bincount(rng.integers(0, self.ntg, self.ntg), minlength=self.ntg)
        return WTF, WTG

    # ------------------------------------------------------------------ within-TF permutation null
    def midranks(self, s):
        """Within-stratum mid-ranks (1-based) of s, in the original pair order."""
        p = self._prep(s)
        r = np.empty(self.n)
        r[p["order"]] = (p["gs"] + 1 + p["ge"]) / 2.0 - p["s0"]
        return r

    def perm_null(self, score_cols, n_perm, rng, batch=200, return_labels_first=0):
        """Joint within-TF label permutation applied to every column of score_cols (n x m).
        Returns (null (n_perm x m), first `return_labels_first` permuted label vectors for verification)."""
        assert np.array_equal(self.strata, self.tf_code), "perm_null is defined for the within-TF statistic"
        R = np.column_stack([self.midranks(score_cols[:, j]) for j in range(score_cols.shape[1])])
        o0 = np.argsort(self.tf_code, kind="stable"); tfc0 = self.tf_code[o0].astype(np.float64)
        y0 = self.y[o0]; R0 = R[o0]
        n1 = np.bincount(self.tf_code, weights=self.y.astype(float), minlength=self.ntf)
        n0 = np.bincount(self.tf_code, minlength=self.ntf) - n1
        const = float((n1 * (n1 + 1) / 2).sum()); den = float((n1 * n0).sum())
        out = np.empty((n_perm, R.shape[1])); keep = []; done = 0
        while done < n_perm:
            b = min(batch, n_perm - done)
            idx = np.argsort(tfc0[None, :] + rng.random((b, self.n)), axis=1)
            Y0 = y0[idx]                                   # permuted labels at TF-sorted positions
            out[done:done + b] = (Y0.astype(np.float64) @ R0 - const) / den
            if len(keep) < return_labels_first:
                for r in range(min(b, return_labels_first - len(keep))):
                    yo = np.empty(self.n, bool); yo[o0] = Y0[r]; keep.append(yo)
            done += b
        return out, keep


def verify_against_stats3(sp: S.StratifiedPairs, s: np.ndarray, seed: int, n_boot=10, n_perm=5):
    """Numerical equality checks of the fast evaluators against stats3 on real pairs/scores."""
    fs = FastStrat(sp)
    a_s3 = sp.auroc(s); a_f = fs.auc(s)[0]
    b_s3 = sp.bootstrap(s, n_boot, np.random.default_rng(seed), two_way=True)
    b_f = fs.auc(s, fs.boot_weights(n_boot, np.random.default_rng(seed), True))
    b1_s3 = sp.bootstrap(s, n_boot, np.random.default_rng(seed), two_way=False)
    b1_f = fs.auc(s, fs.boot_weights(n_boot, np.random.default_rng(seed), False))
    nul, labs = fs.perm_null(s[:, None], n_perm, np.random.default_rng(seed), return_labels_first=n_perm)
    p_s3 = np.array([S.StratifiedPairs(sp.tf, sp.tg, yy).auroc(s) for yy in labs])
    cnt_ok = all(np.array_equal(np.bincount(sp.tf_code, weights=yy.astype(float), minlength=sp.tfs.size),
                                np.bincount(sp.tf_code, weights=sp.y.astype(float), minlength=sp.tfs.size))
                 for yy in labs)
    return dict(point_absdiff=float(abs(a_s3 - a_f)),
                boot2way_max_absdiff=float(np.nanmax(np.abs(b_s3 - b_f))),
                boot_tf_only_max_absdiff=float(np.nanmax(np.abs(b1_s3 - b1_f))),
                perm_max_absdiff=float(np.max(np.abs(p_s3 - nul[:, 0]))),
                perm_counts_per_tf_preserved=bool(cnt_ok), n_boot_checked=n_boot, n_perm_checked=n_perm)


def abs_corr_pairs(Xc, a, b, chunk=40000):
    """|Pearson correlation| of columns a[i], b[i] of a cells x genes matrix (zero-variance genes -> 0)."""
    X = np.asarray(Xc, np.float64)
    mu = X.mean(0); sd = X.std(0)
    Z = (X - mu) / np.where(sd > 0, sd, np.inf)
    out = np.empty(a.size)
    for i in range(0, a.size, chunk):
        aa, bb = a[i:i + chunk], b[i:i + chunk]
        out[i:i + chunk] = np.abs(np.einsum("ij,ij->j", Z[:, aa], Z[:, bb]) / X.shape[0])
    return out


def rowcos(En, a, b, chunk=100000):
    out = np.empty(a.size)
    for i in range(0, a.size, chunk):
        out[i:i + chunk] = np.einsum("ij,ij->i", En[a[i:i + chunk]], En[b[i:i + chunk]])
    return out


def summarise_null(x):
    x = np.asarray(x, float)
    q = np.quantile(x, [0.01, 0.025, 0.05, 0.5, 0.95, 0.975, 0.99])
    return dict(N=int(x.size), mean=float(x.mean()), sd=float(x.std(ddof=1)),
                q01=float(q[0]), q025=float(q[1]), q05=float(q[2]), q50=float(q[3]),
                q95=float(q[4]), q975=float(q[5]), q99=float(q[6]), max=float(x.max()))


def ci_row(est, boot, prefix=""):
    boot = np.asarray(boot, float)
    return {f"{prefix}est": float(est),
            f"{prefix}ci_lo": float(np.nanquantile(boot, 0.025)),
            f"{prefix}ci_hi": float(np.nanquantile(boot, 0.975)),
            f"{prefix}ub95": S.one_sided_upper(boot, 0.95),
            f"{prefix}boot_sd": float(np.nanstd(boot, ddof=1)),
            f"{prefix}boot_mean": float(np.nanmean(boot))}
