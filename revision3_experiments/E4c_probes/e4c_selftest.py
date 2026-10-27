"""Self-tests for e4c_lib: solver vs sklearn, fast metric vs stats3, fold disjointness, rewiring invariants."""
import numpy as np
from sklearn.linear_model import LogisticRegression
import e4c_lib as L

rng = np.random.default_rng(0)
# 1. batched solver == sklearn LogisticRegression (C = 1/(n alpha)) for several alphas
n, d = 3000, 40
X = rng.normal(size=(n, d)).astype(np.float32)
wtrue = rng.normal(size=d) * 0.3
Y = np.column_stack([(rng.random(n) < 1 / (1 + np.exp(-(X @ wtrue - 2)))) for _ in range(3)]).astype(np.float32)
for a in (1e-4, 1e-2, 1e-1):
    W, b, it = L.fit_logistic_batched(X, Y, a, max_iter=100, tol=1e-12)
    for k in range(3):
        sk = LogisticRegression(C=1 / (n * a), tol=1e-10, max_iter=5000).fit(X, Y[:, k])
        err = np.abs(sk.coef_[0] - W[:, k]).max()
        assert err < 2e-3 * max(1, np.abs(sk.coef_).max()), (a, k, err)
    print(f"solver alpha={a}: iters={it}, max|w-w_sklearn|={err:.2e}")
# 2. fast within-TF AUROC == stats3.StratifiedPairs (incl. ties, fallback)
tf = rng.integers(0, 15, 600); tg = rng.integers(0, 300, 600)
Ys = rng.random((600, 4)) < 0.2; Sc = np.round(rng.normal(size=(600, 4)), 1)
fast = L.fast_within_tf_auroc(tf, Ys, Sc)
for k in range(4):
    ref = L.fold_metric(tf, tg, Ys[:, k], Sc[:, k])[0]
    assert abs(ref - fast[k]) < 1e-9, (ref, fast[k])
tf1 = np.zeros(10, int); Ys1 = np.zeros((10, 1), bool); Ys1[:2] = True; Sc1 = rng.normal(size=(10, 1))
assert L.fold_metric(tf1, np.arange(10), Ys1[:, 0], Sc1[:, 0])[1] == "pooled"
assert abs(L.fast_within_tf_auroc(tf1, Ys1, Sc1)[0] - L.fold_metric(tf1, np.arange(10), Ys1[:, 0], Sc1[:, 0])[0]) < 1e-12
print("fast within-TF AUROC == stats3.StratifiedPairs: OK")
# 3. folds on a real context: C2 TF-disjoint, C3 gene-disjoint, C1 partition
c = L.Ctx("immune")
for reg in ("C1", "C2", "C3"):
    fs = L.make_folds(c, reg, 0)
    allte = np.concatenate([te for _, te in fs])
    if reg != "C3": assert np.array_equal(np.sort(allte), np.arange(c.n))
    for tr, te in fs:
        assert np.intersect1d(tr, te).size == 0
        if reg == "C2": assert np.intersect1d(c.tf[tr], c.tf[te]).size == 0
        if reg == "C3":
            gtr = np.union1d(c.tf[tr], c.tg[tr]); gte = np.union1d(c.tf[te], c.tg[te])
            assert np.intersect1d(gtr, gte).size == 0
        itr, iva = L.inner_split(c, reg, tr, 1)
        assert np.intersect1d(itr, iva).size == 0
    assert np.array_equal(np.concatenate([te for _, te in L.make_folds(c, reg, 0)]), allte)   # deterministic
print("fold regimes: OK")
# 4. rewiring preserves every TF's and target's positive count; stays inside the universe
Yr, st = L.rewire_draws(c, 5, seed=7)
for yr in Yr:
    assert yr.sum() == c.y.sum()
    assert np.array_equal(np.bincount(c.tf[yr], minlength=c.G), np.bincount(c.tf[c.y], minlength=c.G))
    assert np.array_equal(np.bincount(c.tg[yr], minlength=c.G), np.bincount(c.tg[c.y], minlength=c.G))
print("rewiring invariants: OK; retained fraction", [round(float((yr & c.y).sum() / c.y.sum()), 3) for yr in Yr], st[:2])
print("E4c self-tests: ALL PASS")
# 5. Gather design == Dense design (concat) for solver and gram
E1 = rng.normal(size=(50, 6)).astype(np.float32); E2 = rng.normal(size=(50, 6)).astype(np.float32)
ti = rng.integers(0, 50, 2000); gi = rng.integers(0, 50, 2000)
Xc = np.hstack([E1[ti], E2[gi]]); Yc = (rng.random((2000, 3)) < 0.2).astype(np.float32)
Wd, bd, _ = L.fit_logistic_batched(L.Dense(Xc), Yc, 1e-3, max_iter=60, tol=1e-12)
Wg, bg, _ = L.fit_logistic_batched(L.Gather([E1, E2], [ti, gi], 50), Yc, 1e-3, max_iter=60, tol=1e-12)
assert np.abs(Wd - Wg).max() < 1e-4 and np.abs(bd - bg).max() < 1e-4, np.abs(Wd - Wg).max()
dv = rng.random(2000)
Hd, cd = L.Dense(Xc).gram(dv); Hg, cg = L.Gather([E1, E2], [ti, gi], 50).gram(dv)
assert np.abs(Hd - Hg).max() < 1e-2 and np.abs(cd - cg).max() < 1e-2
print("Gather design == Dense design: OK")
print("E4c self-tests (incl. Gather): ALL PASS")
