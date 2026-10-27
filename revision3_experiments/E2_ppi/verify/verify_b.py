"""Independent verifier for E2_ppi, part B (stochastic re-derivations with my own code and seeds).

No import of e2_common / stats3. Per context (lung, lung_external):
  1. my own pure-Python Maslov-Sneppen double-edge-swap chain (different swap orientation rule, different RNG):
     burn-in 10m accepted swaps, thinning m accepted swaps, N_MS joint draws; all 120 assortativities per draw;
     degree sequence and simplicity asserted on every draw. -> null mean / sd / z / max-T family.
  2. vertex-subsampling interval for the top cell's r (my implementation, m = n/2, B = 500), with the
     sqrt(m/n) rescaling of CONVENTION3 s6 AND the delete-d (d = n/2) jackknife scaling (factor 1) for comparison.
  3. co-expression-matched AUROC (extraction cells): own non-edge pool (no physical STRING link at any score),
     own decile matching (5 per positive), stratified AUROC via sklearn per stratum, for full-dim cosine at the
     top layer, -|dx| of the top coordinate, and |co-expression|.
Writes verify/verify_b_<ctx>.json (resumable per context).
"""
import gzip, json, random, sys, time
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_a import var_names, coords_eig, H5, EMB, CACHE, SDIR, E2   # my own helpers only

HERE = Path(__file__).resolve().parent
N_MS = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
TOP = {"lung": (3, 1), "lung_external": (8, 2)}


def all_coords(E, obs):
    return np.concatenate([coords_eig(E, L, obs) for L in range(12)], axis=1)


def assort_all(Xc, a, b):
    A = np.concatenate([Xc[a], Xc[b]]); B = np.concatenate([Xc[b], Xc[a]])
    A = A - A.mean(0); B = B - B.mean(0)
    return (A * B).sum(0) / np.sqrt((A * A).sum(0) * (B * B).sum(0))


def ms_chain(a, b, n, N, seed, Xc, log):
    rnd = random.Random(seed)
    a = list(map(int, a)); b = list(map(int, b)); m = len(a)
    keys = {min(x, y) * n + max(x, y) for x, y in zip(a, b)}
    assert len(keys) == m
    deg0 = np.bincount(np.r_[a, b], minlength=n)

    def swaps(k):
        acc = tried = 0
        while acc < k:
            tried += 1
            i = rnd.randrange(m); j = rnd.randrange(m)
            if i == j: continue
            u, v, x, y = a[i], b[i], a[j], b[j]
            if rnd.random() < 0.5:     # (u,x),(v,y)
                p1, p2, q1, q2 = u, x, v, y
            else:                      # (u,y),(v,x)
                p1, p2, q1, q2 = u, y, v, x
            if p1 == p2 or q1 == q2: continue
            k1 = min(p1, p2) * n + max(p1, p2); k2 = min(q1, q2) * n + max(q1, q2)
            if k1 == k2 or k1 in keys or k2 in keys: continue
            keys.discard(min(u, v) * n + max(u, v)); keys.discard(min(x, y) * n + max(x, y))
            keys.add(k1); keys.add(k2)
            a[i], b[i], a[j], b[j] = p1, p2, q1, q2
            acc += 1
        return acc / tried
    swaps(10 * m)
    R = np.empty((N, Xc.shape[1])); accs = np.empty(N)
    for d in range(N):
        accs[d] = swaps(m)
        aa = np.array(a); bb = np.array(b)
        if d % 100 == 0:
            assert np.array_equal(np.bincount(np.r_[aa, bb], minlength=n), deg0)
            assert len(keys) == m and np.all(aa != bb)
        R[d] = assort_all(Xc, aa, bb)
        if (d + 1) % 100 == 0: log(f"  MS {d+1}/{N}")
    return R, accs


def strat_auc(s, y, st):
    num = den = 0.0
    for q in np.unique(st):
        k = st == q; yp = y[k]
        n1 = yp.sum(); n0 = (~yp).sum()
        if n1 and n0:
            num += roc_auc_score(yp, s[k]) * n1 * n0; den += n1 * n0
    return num / den


def main(ctx):
    out = HERE / f"verify_b_{ctx}.json"
    if out.exists(): print("exists", ctx); return
    logf = open(HERE / f"verify_b_{ctx}.log", "a")
    def log(s): print(s, flush=True); logf.write(s + "\n"); logf.flush()
    t0 = time.time()
    genes_all = var_names(H5[ctx]); E = np.load(EMB / f"{ctx}.npy", mmap_mode="r")
    obs = np.flatnonzero(np.linalg.norm(np.asarray(E[0], np.float64), axis=1) > 0); n = obs.size
    V = [genes_all[i] for i in obs]; pos = {g: k for k, g in enumerate(V)}
    e2s = {}
    with gzip.open(SDIR / "string_protein_info_v12.txt.gz", "rt") as f:
        next(f)
        for line in f:
            p = line.split("\t", 2); e2s[p[0]] = p[1]
    e700, eany = set(), set()
    with gzip.open(SDIR / "string_links_v12.txt.gz", "rt") as f:
        next(f)
        for line in f:
            x, y, s = line.split(); gx, gy = e2s.get(x), e2s.get(y)
            if gx in pos and gy in pos and gx != gy:
                k = (min(pos[gx], pos[gy]), max(pos[gx], pos[gy])); eany.add(k)
                if int(s) >= 700: e700.add(k)
    ed = np.array(sorted(e700)); a, b = ed[:, 0], ed[:, 1]; m = a.size
    Xc = all_coords(E, obs)
    obsr = assort_all(Xc, a, b)
    log(f"[{ctx}] n={n} m={m} coords {time.time()-t0:.0f}s")
    res = dict(context=ctx, n=int(n), m=int(m), N_ms=N_MS)
    # 1. MS
    R, accs = ms_chain(a, b, n, N_MS, seed=777 + len(ctx), Xc=Xc, log=log)
    mu = R.mean(0); sd = R.std(0, ddof=1); z = (obsr - mu) / sd
    zn = np.abs((R - mu) / sd); mx = zn.max(1)
    padj = np.array([(np.count_nonzero(mx >= abs(q)) + 1) / (N_MS + 1) for q in z])
    L, k = TOP[ctx]; c = L * 10 + k - 1
    import pandas as pd
    df = pd.read_csv(E2 / f"h1_primary_{ctx}.csv")
    res["ms"] = dict(acc_mean=float(accs.mean()), top_cell=f"L{L} SV{k}", top_r=float(obsr[c]), top_null_mean=float(mu[c]),
                     top_null_sd=float(sd[c]), top_z=float(z[c]), csv_top_null_mean=float(df.null_mean[c]),
                     csv_top_null_sd=float(df.null_sd[c]), csv_top_z=float(df.z[c]),
                     max_abs_diff_null_mean_120=float(np.abs(mu - df.null_mean.values).max()),
                     median_ratio_sd_mine_over_csv=float(np.median(sd / df.null_sd.values)),
                     min_ratio_sd=float((sd / df.null_sd.values).min()), max_ratio_sd=float((sd / df.null_sd.values).max()),
                     corr_z_mine_csv=float(np.corrcoef(z, df.z.values)[0, 1]), min_abs_z_mine=float(np.abs(z).min()),
                     family_min_padj_mine=float(padj.min()), n_cells_padj_lt_05_mine=int((padj < 0.05).sum()),
                     maxT_null_q95=float(np.quantile(mx, 0.95)),
                     lag1_top=float(np.corrcoef(R[:-1, c], R[1:, c])[0, 1]))
    log(f"[{ctx}] MS done {time.time()-t0:.0f}s {res['ms']}")
    # 2. vertex subsampling
    rng = np.random.default_rng(99); x = Xc[:, c]; B = 500; mm = n // 2; dr = np.empty(B)
    for r_ in range(B):
        keep = np.zeros(n, bool); keep[rng.choice(n, mm, replace=False)] = True
        e = keep[a] & keep[b]
        dr[r_] = np.corrcoef(np.r_[x[a[e]], x[b[e]]], np.r_[x[b[e]], x[a[e]]])[0, 1]
    dev = dr - dr.mean()
    res["vertex"] = dict(B=B, est=float(obsr[c]),
                         ci_sqrt_mn=[float(obsr[c] - np.quantile(dev * np.sqrt(mm / n), 0.975)),
                                     float(obsr[c] - np.quantile(dev * np.sqrt(mm / n), 0.025))],
                         ci_delete_d=[float(obsr[c] - np.quantile(dev, 0.975)), float(obsr[c] - np.quantile(dev, 0.025))],
                         sd_halfsample=float(dr.std(ddof=1)))
    # delete-1 jackknife SE for comparison
    jk = []
    deg = np.bincount(np.r_[a, b], minlength=n)
    for g in np.flatnonzero(deg > 0):
        e = (a != g) & (b != g)
        jk.append(np.corrcoef(np.r_[x[a[e]], x[b[e]]], np.r_[x[b[e]], x[a[e]]])[0, 1])
    jk = np.array(jk); res["vertex"]["jackknife_se_delete1"] = float(np.sqrt((jk.size - 1) / jk.size * ((jk - jk.mean()) ** 2).sum()))
    log(f"[{ctx}] vertex {res['vertex']}")
    # 3. co-expression matched (extraction cells)
    M = np.load(CACHE / ctx / "expr_extraction.npy").astype(np.float64)
    M = M - M.mean(0); sdv = M.std(0); Z = np.where(sdv > 0, M / np.where(sdv > 0, sdv, 1), 0).T.copy()
    def corr(p, q): return np.einsum("ij,ij->i", Z[p], Z[q]) / Z.shape[1]
    rng = np.random.default_rng(1234)
    pa, pb = [], []; seen = set()
    while len(pa) < 30 * m:
        u, v = rng.integers(0, n, 2)
        if u == v: continue
        k_ = (min(u, v), max(u, v))
        if k_ in eany or k_ in seen: continue
        seen.add(k_); pa.append(k_[0]); pb.append(k_[1])
    pa = np.array(pa); pb = np.array(pb)
    cp, cq = np.abs(corr(a, b)), np.abs(corr(pa, pb))
    edges_q = np.quantile(np.r_[cp, cq], np.linspace(0, 1, 11))
    dpos = np.clip(np.searchsorted(edges_q, cp, side="right") - 1, 0, 9); dq = np.clip(np.searchsorted(edges_q, cq, side="right") - 1, 0, 9)
    sel = []
    for d in range(10):
        need = 5 * int((dpos == d).sum()); cand = np.flatnonzero(dq == d)
        sel.append(rng.choice(cand, min(need, cand.size), replace=False))
    sel = np.concatenate(sel)
    A = np.r_[a, pa[sel]]; Bq = np.r_[b, pb[sel]]; y = np.r_[np.ones(m, bool), np.zeros(sel.size, bool)]; st = np.r_[dpos, dq[sel]]
    El = np.asarray(E[L], np.float64)[obs]; El = El - El.mean(0); El /= np.linalg.norm(El, axis=1, keepdims=True)
    sc = {f"cosine_L{L}": np.einsum("ij,ij->i", El[A], El[Bq]), "top_cell_neg_absdiff": -np.abs(x[A] - x[Bq]),
          "abs_coexpr": np.abs(corr(A, Bq))}
    res["coexpr_extraction"] = {k_: dict(strat=float(strat_auc(s, y, st)), pooled=float(roc_auc_score(y, s))) for k_, s in sc.items()}
    res["coexpr_extraction"]["n_matched_neg"] = int(sel.size)
    log(f"[{ctx}] coexpr {res['coexpr_extraction']}")
    out.write_text(json.dumps(res, indent=1))
    log(f"[{ctx}] done {time.time()-t0:.0f}s")


if __name__ == "__main__":
    for c_ in sys.argv[1].split(","):
        main(c_)
