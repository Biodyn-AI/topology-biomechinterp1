"""Independent pure-Python class-pair-preserving double-edge-swap chain (no numba, no f1_common), run on the
verifier's own coordinates/classes (verify/_v_<ctx>.npz from verify_f1.py). A swap of (u,v),(x,y) to (u,y),(x,v)
is accepted iff the multiset of endpoint-class pairs is unchanged and the graph stays simple. Burn 10m accepted,
thin m accepted, N_DRAWS draws. Reports the null mean / sd of r at the E2 top cell vs strat_null.csv cpms columns,
checks invariants on every draw. Also a plain MS chain (no class constraint) as a reference.
Usage: python verify_cpms.py --ctx lung --draws 150
"""
import argparse, json, random, time
from pathlib import Path
import numpy as np, pandas as pd

V = Path(__file__).resolve().parent
TOP = {"lung": (3, 1), "lung_external": (8, 2), "kidney_B": (9, 3)}


def chain(a, b, cls, n_draws, seed, constrained):
    rnd = random.Random(seed); A = list(map(int, a)); B = list(map(int, b)); m = len(A)
    es = {(min(u, v), max(u, v)) for u, v in zip(A, B)}
    def pair(p, q): return (min(cls[p], cls[q]), max(cls[p], cls[q]))
    def step(k):
        acc = 0; tr = 0
        while acc < k:
            tr += 1
            i = rnd.randrange(m); j = rnd.randrange(m)
            if i == j: continue
            u, v = A[i], B[i]; x, y = A[j], B[j]
            if rnd.random() < 0.5: u, v = v, u
            if rnd.random() < 0.5: x, y = y, x
            if u == y or x == v: continue
            e1 = (min(u, y), max(u, y)); e2 = (min(x, v), max(x, v))
            if e1 == e2 or e1 in es or e2 in es: continue
            if constrained and sorted([pair(u, v), pair(x, y)]) != sorted([pair(u, y), pair(x, v)]): continue
            es.discard((min(u, v), max(u, v))); es.discard((min(x, y), max(x, y))); es.add(e1); es.add(e2)
            A[i], B[i], A[j], B[j] = u, y, x, v; acc += 1
        return acc / tr
    step(10 * m)
    for _ in range(n_draws):
        ar = step(m)
        yield np.array(A), np.array(B), ar


def r_assort(x, a, b):
    return float(np.corrcoef(np.r_[x[a], x[b]], np.r_[x[b], x[a]])[0, 1])


def main(ctx, n_draws):
    z = np.load(V / f"_v_{ctx}.npz"); X, cls, a, b = z["X"], z["cls"], z["a"], z["b"]
    L, k = TOP[ctx]; x = X[:, L * 10 + k - 1]; n = len(x); cl = [int(c) for c in cls]
    deg0 = np.bincount(np.r_[a, b], minlength=n)
    cp0 = sorted((min(cl[p], cl[q]), max(cl[p], cl[q])) for p, q in zip(a, b))
    res = {}
    for nm, con in (("cpms", True), ("ms", False)):
        t0 = time.time(); vals = []; accs = []
        for A, B, ar in chain(a, b, cl, n_draws, seed=4242 + con, constrained=con):
            assert np.array_equal(np.bincount(np.r_[A, B], minlength=n), deg0)
            if con: assert sorted((min(cl[p], cl[q]), max(cl[p], cl[q])) for p, q in zip(A, B)) == cp0
            vals.append(r_assort(x, A, B)); accs.append(ar)
        vals = np.array(vals)
        res[nm] = dict(null_mean=float(vals.mean()), null_sd=float(vals.std(ddof=1)), N=n_draws, acc=float(np.mean(accs)),
                       secs=time.time() - t0)
        print(ctx, nm, res[nm], flush=True)
    st = pd.read_csv(V.parent / "strat_null.csv"); q = st[(st.context == ctx) & (st.layer == L) & (st.sv == k)].iloc[0]
    res["implementer"] = dict(cpms_null_mean=float(q.cpms_null_mean), cpms_null_sd=float(q.cpms_null_sd),
                              ms_null_mean=float(q.f1ms_null_mean), ms_null_sd=float(q.f1ms_null_sd))
    (V / f"verify_cpms_{ctx}.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--ctx", default="lung"); ap.add_argument("--draws", type=int, default=150)
    p = ap.parse_args()
    for c in p.ctx.split(","): main(c, p.draws)
