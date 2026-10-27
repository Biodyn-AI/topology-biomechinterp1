"""Independent verifier, part 2: an independent exact conditional-logit implementation (log-space elementary-symmetric
DP in torch, autograd gradient, scipy BFGS + Newton polish with autograd Hessian), not using e5_lib.

Checks (written to verify_clogit_<tag>.json):
  point:  adjusted and unadjusted score coefficients for selected (context, layer) cells vs ckpt beta_adj / beta_unadj;
          reduced model (covariates only) and the LRT of the score term.
  boot:   bootstrap refits with TF multiplicity implemented by PHYSICALLY duplicating strata (E5 used stratum weights)
          and target multiplicity by duplicating rows; the draw b is regenerated from the documented RNG recipe
          default_rng([20261026, 5, ctx_id, b]); compared draw-by-draw with ckpt boot_adj; if all 500 are run, the
          percentile CI is recomputed.
  perm:   within-TF permutation refits (stats3.StratifiedPairs.permute_within_tf, default_rng([20261026, 55, ctx_id, r])).
usage: python verify_clogit.py <tag> <ctx> <layer> <n_boot> <n_perm>
"""
import json, sys, time
from pathlib import Path
import numpy as np, torch
from scipy.optimize import minimize
from scipy.stats import chi2

torch.set_num_threads(2); torch.set_default_dtype(torch.float64)
V = Path(__file__).resolve().parent; E5 = V.parent; R3 = E5.parent
sys.path.insert(0, str(R3))
import stats3 as S  # noqa: shared estimator (permutation recipe only)
CTXS = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"


class CL:
    def __init__(self, X, y, g):
        order = np.argsort(g, kind="stable"); X, y, g = X[order], y[order], g[order]
        ug, st, cnt = np.unique(g, return_index=True, return_counts=True)
        keep = [i for i in range(ug.size) if 0 < y[st[i]:st[i] + cnt[i]].sum() < cnt[i]]
        G = len(keep); nmax = int(cnt[keep].max()); p = X.shape[1]
        ks = np.array([int(y[st[i]:st[i] + cnt[i]].sum()) for i in keep]); kmax = int(ks.max())
        Xp = np.zeros((G, nmax, p)); M = np.zeros((G, nmax), bool); sx = np.zeros((G, p))
        for j, i in enumerate(keep):
            a, n = st[i], cnt[i]; Xp[j, :n] = X[a:a + n]; M[j, :n] = True; sx[j] = X[a:a + n][y[a:a + n] == 1].sum(0)
        self.X = torch.tensor(Xp); self.M = torch.tensor(M); self.k = torch.tensor(ks); self.sx = torch.tensor(sx)
        self.kmax = kmax; self.nmax = nmax; self.G = G; self.p = p

    def ll(self, beta):
        eta = self.X @ beta                                   # G x nmax
        c = torch.where(self.M, eta, torch.full_like(eta, -1e300)).max(1).values.detach()
        eta = eta - c[:, None]
        lf = torch.full((self.G, self.kmax + 1), -1e100); lf[:, 0] = 0.0
        for i in range(self.nmax):
            new_tail = torch.logaddexp(lf[:, 1:], lf[:, :-1] + eta[:, i:i + 1])
            new = torch.cat([lf[:, :1], new_tail], 1)
            lf = torch.where(self.M[:, i:i + 1], new, lf)
        loge = lf.gather(1, self.k[:, None]).squeeze(1) + self.k * c
        return (self.sx @ beta).sum() - loge.sum()

    def fit(self, beta0=None):
        b0 = np.zeros(self.p) if beta0 is None else np.asarray(beta0, float)

        def f(b):
            bt = torch.tensor(b, requires_grad=True); v = -self.ll(bt); v.backward()
            return float(v.detach()), bt.grad.numpy().copy()
        r = minimize(f, b0, jac=True, method="BFGS", options=dict(gtol=1e-8, maxiter=500))
        b = r.x

        def hess(b, h=1e-5):                                 # central differences of the exact autograd gradient
            H = np.zeros((self.p, self.p))
            for j in range(self.p):
                e = np.zeros(self.p); e[j] = h
                H[:, j] = (f(b + e)[1] - f(b - e)[1]) / (2 * h)
            return (H + H.T) / 2
        for _ in range(3):                                   # Newton polish
            _, g = f(b); step = np.linalg.solve(hess(b), g); b = b - step
            if np.max(np.abs(step)) < 1e-10: break
        H = hess(b)
        val, g = f(b)
        return dict(beta=b, ll=-val, se=np.sqrt(np.diag(np.linalg.inv(H))), max_grad=float(np.max(np.abs(g))))


def z(x):
    x = np.asarray(x, float); return (x - x.mean(0)) / x.std(0)


def main(tag, ctx, layer, n_boot, n_perm):
    t0 = time.time()
    d = R3 / "cache" / ctx
    obs_idx = np.load(d / "obs_idx.npy")
    P = np.load(d / "h3_pairs.npz"); tf, tg, y = P["tf"].astype(int), P["tg"].astype(int), P["y"].astype(int)
    cov = np.load(V / f"cov_{ctx}.npy")                    # verifier's own covariates (verify_part1.py)
    emb = np.load(EMB / f"{ctx}.npy", mmap_mode="r")
    M = np.asarray(emb[layer], np.float64)[obs_idx]; M = M - M.mean(0); M = M / np.linalg.norm(M, axis=1, keepdims=True)
    s = z((M[tf] * M[tg]).sum(1)); cz = z(cov)
    Xa = np.column_stack([s, cz]); tfu, tfc = np.unique(tf, return_inverse=True); tgu, tgc = np.unique(tg, return_inverse=True)
    ck = np.load(E5 / "ckpt" / f"{ctx}_L{layer:02d}.npz")
    res = dict(tag=tag, context=ctx, layer=layer)
    fa = CL(Xa, y, tfc).fit(); fu = CL(s[:, None], y, tfc).fit(); fr = CL(cz, y, tfc).fit()
    lrt = 2 * (fa["ll"] - fr["ll"])
    res.update(beta_adj_verify=fa["beta"].tolist(), beta_adj_e5=ck["beta_adj"].tolist(),
               max_abs_diff_beta_adj=float(np.max(np.abs(fa["beta"] - ck["beta_adj"]))),
               se_score_verify=float(fa["se"][0]), se_score_e5=float(ck["se_adj"][0]),
               beta_unadj_verify=float(fu["beta"][0]), beta_unadj_e5=float(ck["beta_unadj"][0]),
               ll_adj_verify=fa["ll"], ll_adj_e5=float(ck["ll_adj"]), lrt_verify=lrt, lrt_e5=float(ck["lrt_adj"]),
               p_lrt_verify=float(chi2.sf(lrt, 1)), p_lrt_e5=float(ck["p_lrt_adj"]), max_grad=fa["max_grad"],
               or_adj_verify=float(np.exp(fa["beta"][0])))
    print(json.dumps(res), flush=True)
    bd = []
    for b in range(n_boot):
        rng = np.random.default_rng([20261026, 5, CTXS.index(ctx), b])
        wtf = np.bincount(rng.integers(0, tfu.size, tfu.size), minlength=tfu.size)
        wtg = np.bincount(rng.integers(0, tgu.size, tgu.size), minlength=tgu.size)
        Xs, ys, gs = [], [], []; gid = 0
        for t in np.flatnonzero(wtf):
            rows = np.flatnonzero(tfc == t); rows = np.repeat(rows, wtg[tgc[rows]])
            for _ in range(wtf[t]):
                Xs.append(Xa[rows]); ys.append(y[rows]); gs.append(np.full(rows.size, gid)); gid += 1
        fb = CL(np.vstack(Xs), np.concatenate(ys), np.concatenate(gs)).fit(beta0=ck["beta_adj"])
        bd.append(float(fb["beta"][0]))
        print(tag, "boot", b, bd[-1], float(ck["boot_adj"][b, 0]), f"{time.time()-t0:.0f}s", flush=True)
    if n_boot:
        bd = np.array(bd); e = ck["boot_adj"][:n_boot, 0]
        res.update(n_boot=n_boot, max_abs_diff_boot=float(np.max(np.abs(bd - e))))
        if n_boot == 500:
            res.update(or_boot_ci_verify=np.exp(np.quantile(bd, [0.025, 0.975])).tolist(),
                       or_boot_ci_e5=np.exp(np.quantile(e, [0.025, 0.975])).tolist())
    sp = S.StratifiedPairs(tf, tg, y); pdiff = []
    for r in range(n_perm):
        yp = sp.permute_within_tf(np.random.default_rng([20261026, 55, CTXS.index(ctx), r])).astype(int)
        fp = CL(Xa, yp, tfc).fit()
        pdiff.append(abs(float(fp["beta"][0]) - float(ck["perm_adj"][r, 0])))
        print(tag, "perm", r, fp["beta"][0], ck["perm_adj"][r, 0], flush=True)
    if n_perm: res.update(n_perm=n_perm, max_abs_diff_perm=float(max(pdiff)))
    res["seconds"] = time.time() - t0
    (V / f"verify_clogit_{tag}.json").write_text(json.dumps(res, indent=1))
    print("DONE", json.dumps(res), flush=True)


if __name__ == "__main__":
    a = sys.argv[1:]; main(a[0], a[1], int(a[2]), int(a[3]), int(a[4]))
