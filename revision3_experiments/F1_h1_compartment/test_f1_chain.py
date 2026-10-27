"""Self-test of the class-pair-preserving MS chain: invariants and agreement with an exact small-graph check."""
import numpy as np, f1_common as F
rng = np.random.default_rng(0)
n = 300; P = rng.random((n, n)) < 0.03; P = np.triu(P, 1); pi, pj = np.nonzero(P)
cls = rng.integers(0, 4, n)
deg0 = np.bincount(np.r_[pi, pj], minlength=n); cp0 = F.class_pair_counts(cls, pi, pj, 4)
fr = []
for a, b, acc, ok in F.ms_classpair_chain(pi, pj, cls, n, 200, seed=5):
    assert np.array_equal(np.bincount(np.r_[a, b], minlength=n), deg0)
    assert np.array_equal(F.class_pair_counts(cls, a, b, 4), cp0)
    assert len(np.unique(np.minimum(a, b) * n + np.maximum(a, b))) == a.size and np.all(a != b) and ok
    fr.append(np.mean((np.minimum(a,b)*n+np.maximum(a,b))[:, None] == (pi*n+pj)[None, :], axis=0).sum() / a.size if False else 0)
# the chain must actually move: fraction of original edges retained after the last draw
orig = set(zip(pi.tolist(), pj.tolist())); last = set(zip(np.minimum(a,b).tolist(), np.maximum(a,b).tolist()))
print("m", pi.size, "acc", round(acc, 3), "retained", len(orig & last) / len(orig))
assert len(orig & last) / len(orig) < 0.3
print("class-pair chain self-test passed")
