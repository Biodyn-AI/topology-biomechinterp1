"""E2 shared helpers: context loading, SVk coordinates, vectorised family statistics, null runners.

All estimators come from ../stats3.py. The only local code is a vectorised evaluation of
stats3.assortativity / stats3.copole_rate over many attribute columns at once (one gather of the edge
endpoints per null draw); `check_vectorised` asserts equality with the stats3 scalar functions.
"""
from __future__ import annotations

import json, os, sys, hashlib
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
R3 = HERE.parent
sys.path.insert(0, str(R3))
import stats3 as S          # noqa: E402
import contexts as X        # noqa: E402

CACHE = R3 / "cache"
CKPT = HERE / "_ckpt"
EMB_LOCAL = Path.home() / "biomechinterp-local/subproject41_r3/emb"
N_LAYERS, K_SV = 12, 10
POLE_FRACS = (0.10, 0.20, 0.27, 0.33, 0.50)
FROZEN_HASH = "3e3249d1d0a8ce1bc3cf64639fe233b016535456f82432a8bbfcacdd64a05c01"
PRIMARY_CTX = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
N_PRIMARY = {"lung": 10000, "lung_external": 10000, "kidney_B": 2000, "immune": 2000, "kidney_A": 2000}
ROBUST_CTX = ["immune_seed43", "immune_seed44", "immune_TSP14", "immune_TSP21", "immune_TSP25",
              "immune_Bcell", "immune_CD4T", "immune_CD8T",
              "kidney_ckpt_whole_human", "kidney_ckpt_kidney", "kidney_ckpt_brain"]
SEED = 20261026


def seed_for(*parts) -> int:
    """Deterministic 31-bit seed from a tuple of labels."""
    h = hashlib.sha256("|".join(map(str, (SEED,) + parts)).encode()).hexdigest()
    return int(h[:7], 16)


def frozen_hash_check():
    txt = (R3 / "FROZEN.txt").read_text()
    assert FROZEN_HASH in txt
    actual = hashlib.sha256((R3 / "CONVENTION3.md").read_bytes()).hexdigest()
    assert actual == FROZEN_HASH, "CONVENTION3.md changed since freezing"
    return FROZEN_HASH


def load_emb(ctx):
    return np.load(EMB_LOCAL / f"{ctx}.npy", mmap_mode="r")


def load_ctx(ctx):
    d = CACHE / ctx
    genes = json.loads((d / "obs_genes.json").read_text())
    obs_idx = np.load(d / "obs_idx.npy")
    edges = np.load(d / "edges_string.npy").astype(np.int64)
    lab = dict(np.load(d / "node_labels.npz"))
    return dict(genes=genes, obs_idx=obs_idx, edges=edges, labels=lab, n=len(genes),
                ncells=np.load(d / "ncells.npy"))


def coords(ctx, obs_idx, want_vt=False):
    """Observed-gene SVD (stats3.svd_observed) at every layer. Returns Xc (n x 120) with column
    c = layer*10 + (k-1) holding the SVk coordinate (U_k S_k) at that layer; optionally Vt[:, :10]."""
    emb = load_emb(ctx)
    assert emb.shape[0] == N_LAYERS
    # V check: every obs_idx row is non-zero, every other row is zero
    assert np.array_equal(np.flatnonzero(S.observed_mask(emb)), obs_idx), f"{ctx}: V mismatch with cache"
    cols, vts, svals = [], [], []
    for L in range(N_LAYERS):
        US, Sv, Vt = S.svd_observed(emb, L, obs_idx)
        cols.append(US[:, :K_SV]); vts.append(Vt[:K_SV]); svals.append(Sv[:K_SV + 1])
    Xc = np.ascontiguousarray(np.concatenate(cols, axis=1))
    return (Xc, np.stack(vts), np.stack(svals)) if want_vt else Xc


def cell_name(c):
    return int(c // K_SV), int(c % K_SV) + 1        # (layer, k)


def pole_matrix(Xc, fracs=POLE_FRACS):
    """int8 (n x (ncol*len(fracs))) pole labels; column = f_index*ncol + c."""
    n, p = Xc.shape
    L = np.empty((n, p * len(fracs)), dtype=np.int8)
    for fi, f in enumerate(fracs):
        for c in range(p):
            L[:, fi * p + c] = S.pole_labels(Xc[:, c], f)
    return L


# ------------------------------------------------------------------ vectorised family statistics
def assort_cols(Xc, a, b):
    """stats3.assortativity for every column of Xc over undirected edges (a, b): Pearson correlation over
    both orientations. r = (2 sum x_a x_b - M mu^2) / (sum x_a^2 + sum x_b^2 - M mu^2), M = 2m."""
    A = Xc[a]; B = Xc[b]; M = 2.0 * a.size
    sab = np.einsum("ij,ij->j", A, B)
    s1 = A.sum(0) + B.sum(0); s2 = np.einsum("ij,ij->j", A, A) + np.einsum("ij,ij->j", B, B)
    mu = s1 / M
    return (2.0 * sab - M * mu * mu) / (s2 - M * mu * mu)


def copole_cols(Lm, a, b):
    """stats3.copole_rate for every column of the pole-label matrix: same non-zero pole <=> product == 1."""
    P = Lm[a] * Lm[b]
    return (P == 1).sum(0) / a.size


def check_vectorised(Xc, Lm, a, b, cols=(0, 7, 55, 119)):
    for c in cols:
        assert abs(assort_cols(Xc[:, [c]], a, b)[0] - S.assortativity(Xc[:, c], a, b)) < 1e-10
    for c in (0, Lm.shape[1] - 1, Lm.shape[1] // 2):
        assert abs(copole_cols(Lm[:, [c]], a, b)[0] - S.copole_rate(Lm[:, c], a, b)) < 1e-12


# ------------------------------------------------------------------ null runners
def run_ms(Xc, Lm, pi, pj, n, N, seed, log=None, every=1000):
    """Maslov-Sneppen chain (stats3.maslov_sneppen_chain, burn 10m, thin m): ALL statistics on each draw."""
    R = np.empty((N, Xc.shape[1])); C = None if Lm is None else np.empty((N, Lm.shape[1]), np.float32)
    acc = np.empty(N)
    for d, (a, b, ar) in enumerate(S.maslov_sneppen_chain(pi, pj, n, N, seed=seed, burn_factor=10, thin_factor=1)):
        R[d] = assort_cols(Xc, a, b)
        if C is not None: C[d] = copole_cols(Lm, a, b)
        acc[d] = ar
        if log and (d + 1) % every == 0: log(f"    MS draw {d+1}/{N} acc={ar:.3f}")
    return R, C, acc


def perm_within(groups, n, rng):
    sig = np.arange(n)
    for g in groups:
        if g.size > 1: sig[g] = rng.permutation(sig[g])
    return sig


def run_perm(Xc, Lm, pi, pj, n, N, rng, strata=None, log=None, every=2000):
    """Gene-label permutation (optionally within strata), ONE permutation applied to every statistic per
    draw. Draw sigma; permuted attribute x'_i = x_sigma(i), so x'[pi] = x[sigma[pi]]."""
    groups = [np.arange(n)] if strata is None else [np.flatnonzero(strata == s) for s in np.unique(strata)]
    R = np.empty((N, Xc.shape[1])); C = None if Lm is None else np.empty((N, Lm.shape[1]), np.float32)
    for d in range(N):
        sig = perm_within(groups, n, rng)
        a, b = sig[pi], sig[pj]
        R[d] = assort_cols(Xc, a, b)
        if C is not None: C[d] = copole_cols(Lm, a, b)
        if log and (d + 1) % every == 0: log(f"    perm draw {d+1}/{N}")
    return R, C


def null_summary(obs, null, tail="two-sided"):
    """Per-statistic summary + max-T adjustment (stats3.max_t_adjust)."""
    z_abs, p_adj, p_raw = S.max_t_adjust(obs, null, tail=tail)
    mu = null.mean(0); sd = null.std(0, ddof=1)
    z = (obs - mu) / sd
    q = np.quantile(null, [0.025, 0.5, 0.975], axis=0)
    return dict(null_mean=mu, null_sd=sd, z=z, p_raw=p_raw, p_adj=p_adj, q025=q[0], q50=q[1], q975=q[2])


def lag1(x):
    x = np.asarray(x, float) - np.mean(x)
    return float((x[:-1] * x[1:]).sum() / (x * x).sum())


def logger(path):
    f = open(path, "a")
    def log(msg):
        print(msg, flush=True); f.write(msg + "\n"); f.flush()
    return log
