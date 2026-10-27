"""E6 supplement (third, post-reboot agent): intervals and paired deltas that e6_run.py did not store.

e6_run.py stored bootstrap draws only for the node tasks, and the detection-stratified Lambert AUROC only as a
point estimate. This script regenerates, with the SAME seeds and rng streams as e6_run.py (so every draw is
paired across embeddings and reproduces the stored intervals), for every (gene set, embedding), estimator (i):
  * node_lambert detection-stratified AUROC (stats3.stratified_auroc over gs['det_dec']):
      - gene-bootstrap CI: the exact gene draws of e6lib.node_boot (seed e6_run.sd(set,'node_lambert')), each
        replicate scored leave-gene-out and pooled within detection deciles with bootstrap-count weights;
        percentile re-centred on the estimate (same convention as the node rows);
      - within-decile label-permutation p (stats3.permute_within_strata, N=2,000).
  * ppi: the B=200 vertex-subsampling draws of stats3.vertex_subsample_ci (rng seed+2), node subsets identical
    across embeddings.
  * edge: the B=2,000 two-way TF x target bootstrap draws (rng seed+5, e6lib.two_way_weights, chunks of 250).
Checks (abort on mismatch): stored value/CI of every ppi and edge row, the stored node bootstrap draws (float32),
and the stored detection-stratified point estimate are reproduced.

usage: python e6_supp.py worker            (claims jobs via work/supp/<set>/<emb>.lock; resumable)
       python e6_supp.py aggregate         (writes node_detstrat_results.csv, pairwise_deltas_vs_scgpt.csv)
"""
import json, os, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import e6lib as L
import e6_run as R
S = L.S
SUP = L.WORK / "supp"; SUP.mkdir(exist_ok=True)
SETS = ["tierA", "tierB_lung", "tierB_lung_external", "tierB_kidney_B", "tierB_immune", "tierB_kidney_A",
        "legacy186", "legacy178"]
N_DPERM = 2000


def embs(setname):
    return [e for e in R.embeddings_for(setname) if "@rp" not in e]


def jobs():
    return [(s, e) for s in SETS for e in embs(s)]


def stored_rows(setname, emb):
    return {(r["task"], r["estimator"]): r for r in json.loads((L.WORK / setname / f"{emb}.json").read_text())["rows"]}


def wauroc(sp, wp, sn, wn):
    """weighted tie-aware Mann-Whitney numerator and denominator."""
    o = np.argsort(sn, kind="mergesort"); sn, wn = sn[o], wn[o]
    cw = np.concatenate([[0.0], np.cumsum(wn)])
    lo = np.searchsorted(sn, sp, side="left"); hi = np.searchsorted(sn, sp, side="right")
    return float((wp * (cw[lo] + 0.5 * (cw[hi] - cw[lo]))).sum()), float(wp.sum() * wn.sum())


def strat_w(score_tf, w_tf, st_tf, score_ng, w_ng, st_ng):
    num = den = 0.0
    for st in np.unique(np.r_[st_tf, st_ng]):
        a = (st_tf == st) & (w_tf > 0); b = (st_ng == st) & (w_ng > 0)
        if a.any() and b.any():
            n_, d_ = wauroc(score_tf[a], w_tf[a], score_ng[b], w_ng[b]); num += n_; den += d_
    return num / den if den else np.nan


def node_det(E, gs, setname, stored):
    lamb = gs["lambert"]; tp = np.flatnonzero(lamb); ng = np.flatnonzero(~lamb); dd = np.asarray(gs["det_dec"])
    n = lamb.size; isp = lamb.copy()
    seed = R.sd(setname, "node_lambert")
    s = S.tf_scores_loo(E, tp); val = S.stratified_auroc(s, isp, dd)
    ref = stored[("node_lambert", "full")].get("value_detection_stratified")
    assert ref is not None and abs(ref - val) < 1e-9, ("detstrat point mismatch", ref, val)
    # bootstrap: exact draws of e6lib.node_boot
    lab = np.concatenate([tp, ng]); ntf = tp.size; nlab = lab.size
    En = S.l2norm_rows(np.asarray(E, dtype=np.float64))
    Gm = En[lab] @ En[tp].T; selfsim = (En[lab] * En[lab]).sum(1)
    st_tf, st_ng = dd[lab[:ntf]], dd[lab[ntf:]]
    rng = np.random.default_rng(seed); B = R.B_NODE
    boot = np.empty(B); boot_u = np.empty(B)
    for r in range(B):
        c = np.bincount(rng.integers(0, nlab, nlab), minlength=nlab)
        ctf = c[:ntf]; K = ctf.sum(); tot = Gm @ ctf
        sc_tf = (tot[:ntf] - ctf * selfsim[:ntf]) / np.maximum(K - ctf, 1); sc_ng = tot[ntf:] / K
        boot[r] = strat_w(sc_tf, ctf.astype(float), st_tf, sc_ng, c[ntf:].astype(float), st_ng)
        if r < 20:
            boot_u[r] = S.auroc(np.repeat(sc_tf, ctf), np.repeat(sc_ng, c[ntf:]))
        if r == 0:   # explicit check of the weighted stratified statistic on the repeated arrays
            sr = np.r_[np.repeat(sc_tf, ctf), np.repeat(sc_ng, c[ntf:])]
            yr = np.r_[np.ones(ctf.sum(), bool), np.zeros(c[ntf:].sum(), bool)]
            stv = np.r_[np.repeat(st_tf, ctf), np.repeat(st_ng, c[ntf:])]
            assert abs(S.stratified_auroc(sr, yr, stv) - boot[r]) < 1e-9, "weighted strat mismatch"
    zb = np.load(L.WORK / setname / f"{emb_g[0]}.npz")["boot_node_lambert_full"].astype(float)
    assert np.max(np.abs(zb[:20] - boot_u[:20])) < 1e-5, ("node bootstrap draws not reproduced", zb[:3], boot_u[:3])
    # within-decile label permutation null (N_DPERM), batched LOO scoring
    prng = np.random.default_rng(seed + 77); null = np.empty(N_DPERM); k = tp.size
    strata = [np.flatnonzero(dd == st) for st in np.unique(dd)]
    from scipy.stats import rankdata
    for c0 in range(0, N_DPERM, 250):
        b = min(250, N_DPERM - c0)
        Y = np.stack([S.permute_within_strata(isp, dd, prng) for _ in range(b)], 1)
        M = Y.astype(float); tot = En @ (En.T @ M)
        ss = (En * En).sum(1)[:, None]
        sc = np.where(Y, (tot - ss) / max(k - 1, 1), tot / k)
        num = np.zeros(b); den = np.zeros(b)
        for idx in strata:
            y = Y[idx]; n1 = y.sum(0); n0 = idx.size - n1
            rk = rankdata(sc[idx], axis=0)
            u = (rk * y).sum(0) - n1 * (n1 + 1) / 2
            ok = (n1 > 0) & (n0 > 0); num[ok] += u[ok]; den[ok] += (n1 * n0)[ok]
        null[c0:c0 + b] = num / den
    p, bb, N = S.emp_p(val, null)
    return val, boot, null, (p, bb, N)


def ppi_draws(E, gs, setname, stored):
    pi, pj, ni, nj = gs["pi"], gs["pj"], gs["ni"], gs["nj"]; n = len(gs["genes"])
    if pi.size < 3: return None
    seed = R.sd(setname, "ppi")
    En = S.l2norm_rows(np.asarray(E, float)); cp = L.pair_cos(En, pi, pj); cn = L.pair_cos(En, ni, nj)
    full = S.auroc(cp, cn); rng = np.random.default_rng(seed + 2); m = int(0.5 * n); B = R.B_PPI
    draws = np.empty(B)
    for r in range(B):
        nodes = np.sort(rng.choice(n, m, replace=False)); keep = np.zeros(n, bool); keep[nodes] = True
        e = keep[pi] & keep[pj]
        if e.sum() > 2:
            x = cp[keep[pi] & keep[pj]]; y = cn[keep[ni] & keep[nj]]
            draws[r] = S.auroc(x, y) if x.size > 2 and y.size > 2 else np.nan
        else:
            draws[r] = np.nan
    dev = (draws - np.nanmean(draws)) * np.sqrt(m / n)
    lo, hi = full - np.nanquantile(dev, 0.975), full - np.nanquantile(dev, 0.025)
    st = stored[("ppi", "full")]
    assert abs(st["value"] - full) < 1e-9 and abs(st["ci_lo"] - lo) < 1e-9 and abs(st["ci_hi"] - hi) < 1e-9, \
        ("ppi CI not reproduced", st["value"], full, st["ci_lo"], lo, st["ci_hi"], hi)
    return full, draws, m, n


def edge_draws(E, gs, setname, stored):
    if gs["e_y"].size == 0 or gs["e_y"].sum() < 3: return None
    seed = R.sd(setname, "edge")
    sp = S.StratifiedPairs(gs["e_tf"], gs["e_tg"], gs["e_y"])
    En = S.l2norm_rows(np.asarray(E, float)); s = L.pair_cos(En, sp.tf, sp.tg)
    auc = sp.auroc(s); fs = L.FastStrat(sp, s); rngw = np.random.default_rng(seed + 5); boot = []
    for c in range(0, R.B_EDGE_FULL, 250):
        W = L.two_way_weights(sp, min(250, R.B_EDGE_FULL - c), rngw); boot.append(fs.auroc_w(W)); del W
    boot = np.concatenate(boot); lo, hi = L.ci(boot)
    st = stored[("edge", "full")]
    assert abs(st["value"] - auc) < 1e-9 and abs(st["ci_lo"] - lo) < 1e-9 and abs(st["ci_hi"] - hi) < 1e-9, \
        ("edge CI not reproduced", st["value"], auc, st["ci_lo"], lo)
    return auc, boot


emb_g = [None]


def run_one(setname, emb):
    out = SUP / setname / f"{emb}.npz"
    if out.exists(): return
    t0 = time.time(); emb_g[0] = emb
    gs = L.load_geneset(setname); P = R.provide(emb, gs["genes"]); E = P["E"]; stored = stored_rows(setname, emb)
    res = {}
    v, b, nul, (p, bb, N) = node_det(E, gs, setname, stored)
    res.update(det_value=v, det_boot=b, det_null=nul, det_p=p, det_pb=bb, det_pN=N)
    pp = ppi_draws(E, gs, setname, stored)
    if pp is not None: res.update(ppi_value=pp[0], ppi_draws=pp[1], ppi_m=pp[2], ppi_n=pp[3])
    ee = edge_draws(E, gs, setname, stored)
    if ee is not None: res.update(edge_value=ee[0], edge_boot=ee[1])
    tmp = out.with_name(out.stem + f".{os.getpid()}.tmp.npz"); np.savez(tmp, **res); tmp.rename(out)
    print(f"{setname} {emb} done {time.time()-t0:.0f}s det={v:.4f} p={p:.4g}", flush=True)


def worker():
    for s, e in jobs():
        d = SUP / s; d.mkdir(exist_ok=True)
        if (d / f"{e}.npz").exists(): continue
        try:
            fd = os.open(d / f"{e}.lock", os.O_CREAT | os.O_EXCL | os.O_WRONLY); os.write(fd, str(os.getpid()).encode()); os.close(fd)
        except FileExistsError:
            continue
        run_one(s, e)


# ------------------------------------------------------------------------------------------- aggregate
def aggregate():
    A = pd.concat([pd.read_csv(L.HERE / "tierA_results.csv"), pd.read_csv(L.HERE / "tierB_results.csv", low_memory=False),
                   pd.read_csv(L.HERE / "legacy_reproduction.csv").head(0)], ignore_index=True)
    allrows = []
    for d in sorted(L.WORK.iterdir()):
        if d.is_dir() and not d.name.startswith(("_", ".")) and d.name in SETS:
            for f in sorted(x for x in d.glob("*.json") if not x.name.startswith("._")):
                allrows += json.loads(f.read_text())["rows"]
    df = pd.DataFrame(allrows); df = df[df.estimator == "full"]
    lab = {(r.set, r.embedding): r.label for r in df.itertuples()}
    lab = {k: v.replace("ESM2 (UCE protein file; 'ESM2-3B' in submitted version)", "ESM2-15B (UCE per-protein)") for k, v in lab.items()}
    Z = {}; missing = []
    for s, e in jobs():
        f = SUP / s / f"{e}.npz"
        if f.exists(): Z[(s, e)] = dict(np.load(f))
        else: missing.append((s, e))
    if missing: print("MISSING supp jobs:", len(missing), missing[:5])
    # ---- detection-stratified node results with CI and p
    D = []
    for (s, e), z in Z.items():
        v = float(z["det_value"]); b = z["det_boot"]; plo, phi = L.ci(b); sh = v - float(np.mean(b))
        x = df[(df.set == s) & (df.embedding == e) & (df.task == "node_lambert")].iloc[0]
        ns = L.null_summary(z["det_null"])
        D.append(dict(set=s, embedding=e, label=lab[(s, e)], n_genes=int(x.n_genes), n_pos=int(x.n_pos), n_neg=int(x.n_neg),
                      task="node_lambert", estimator="full", metric="detection-stratified AUROC (LOO mean cosine to Lambert TF set, "
                      "pooled within detection deciles)", value=v, ci_lo=max(plo + sh, 0), ci_hi=min(phi + sh, 1),
                      ci_method="gene bootstrap B=1000 (draws of e6lib.node_boot), leave-gene-out, count-weighted within-decile pooling, "
                      "percentile re-centred", p=float(z["det_p"]), p_b=int(z["det_pb"]), p_N=int(z["det_pN"]),
                      p_method="within-decile label permutation (stats3.permute_within_strata)", null_mean=ns["mean"], null_q975=ns["q975"],
                      unstratified_value=float(x.value), unstratified_ci_lo=float(x.ci_lo), unstratified_ci_hi=float(x.ci_hi)))
    Ddf = pd.DataFrame(D); Ddf.to_csv(L.HERE / "node_detstrat_results.csv", index=False)
    # ---- paired deltas vs scGPT, all tasks
    def refs_for(s):
        if s == "tierA": return ["scgpt_static"]
        if s.startswith("tierB_"): return ["scgpt_static"] + [f"{s[6:]}@L{l}" for l in range(12)]
        return ["kidney_A@L11", "scgpt_static"]
    def nodeboot(s, e, t): return np.load(L.WORK / s / f"{e}.npz")[f"boot_{t}_full"].astype(float)
    def val(s, e, t):
        x = df[(df.set == s) & (df.embedding == e) & (df.task == t)]; return float(x.value.iloc[0])
    O = []
    for s in SETS:
        for ref in refs_for(s):
            if (s, ref) not in Z: continue
            for e in embs(s):
                if e == ref or (s, e) not in Z: continue
                if "@L" in e and "@L" in ref: continue           # contextual vs contextual layers not needed here
                za, zb = Z[(s, e)], Z[(s, ref)]
                base = dict(set=s, model_a=e, label_a=lab[(s, e)], model_b=ref, label_b=lab[(s, ref)])
                for t in ("node_lambert", "node_trrust"):
                    va, vb = val(s, e, t), val(s, ref, t); d = va - vb; db = nodeboot(s, e, t) - nodeboot(s, ref, t)
                    lo, hi = L.ci(db); sh = d - float(np.mean(db))
                    O.append(dict(**base, task=t, value_a=va, value_b=vb, delta=d, ci_lo=lo + sh, ci_hi=hi + sh, B=db.size,
                                  frac_le0=float(np.mean(db + sh <= 0)),
                                  ci_method="paired gene bootstrap (stored draws), percentile re-centred"))
                va, vb = float(za["det_value"]), float(zb["det_value"]); d = va - vb; db = za["det_boot"] - zb["det_boot"]
                lo, hi = L.ci(db); sh = d - float(np.mean(db))
                O.append(dict(**base, task="node_lambert_detstrat", value_a=va, value_b=vb, delta=d, ci_lo=lo + sh, ci_hi=hi + sh,
                              B=db.size, frac_le0=float(np.mean(db + sh <= 0)),
                              ci_method="paired gene bootstrap, detection-stratified, percentile re-centred"))
                if "ppi_draws" in za and "ppi_draws" in zb:
                    va, vb = float(za["ppi_value"]), float(zb["ppi_value"]); d = va - vb
                    db = za["ppi_draws"] - zb["ppi_draws"]; m, n = int(za["ppi_m"]), int(za["ppi_n"])
                    dev = db - np.nanmean(db)
                    r1 = dict(ci_lo=d - np.sqrt(m / n) * np.nanquantile(dev, .975), ci_hi=d - np.sqrt(m / n) * np.nanquantile(dev, .025))
                    r2 = dict(ci_lo_deleted=d - np.sqrt(m / (n - m)) * np.nanquantile(dev, .975),
                              ci_hi_deleted=d - np.sqrt(m / (n - m)) * np.nanquantile(dev, .025))
                    O.append(dict(**base, task="ppi", value_a=va, value_b=vb, delta=d, **r1, **r2, B=int(np.isfinite(db).sum()),
                                  frac_le0=float(np.nanmean(d + np.sqrt(m / n) * dev <= 0)),
                                  ci_method="paired vertex subsampling m=n/2 (identical node subsets), ci = sqrt(m/n) pre-registered; "
                                            "ci_*_deleted = delete-d sqrt(m/(n-m))"))
                if "edge_boot" in za and "edge_boot" in zb:
                    va, vb = float(za["edge_value"]), float(zb["edge_value"]); d = va - vb; db = za["edge_boot"] - zb["edge_boot"]
                    lo, hi = L.ci(db)
                    O.append(dict(**base, task="edge", value_a=va, value_b=vb, delta=d, ci_lo=lo, ci_hi=hi, B=db.size,
                                  frac_le0=float(np.mean(db <= 0)),
                                  ci_method="paired two-way TF x target bootstrap (identical weights), percentile (as edge rows)"))
    Odf = pd.DataFrame(O)
    Odf["ci_excludes_0"] = (Odf.ci_lo > 0) | (Odf.ci_hi < 0)
    # verifier addition: for PPI rows, ci_excludes_0 uses the pre-registered sqrt(m/n) interval; the manuscript uses
    # delete-d (CONVENTION3_ERRATA item 3), so record the delete-d decision too (NaN for non-PPI rows).
    Odf["ci_excludes_0_deleted"] = np.where(Odf.task == "ppi", (Odf.ci_lo_deleted > 0) | (Odf.ci_hi_deleted < 0), np.nan)
    Odf.to_csv(L.HERE / "pairwise_deltas_vs_scgpt.csv", index=False)
    print("node_detstrat rows", len(Ddf), "delta rows", len(Odf), "missing", len(missing))


if __name__ == "__main__":
    worker() if sys.argv[1] == "worker" else aggregate()
