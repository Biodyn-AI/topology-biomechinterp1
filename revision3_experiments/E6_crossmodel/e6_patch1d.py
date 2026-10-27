"""E6 fix (resume agent): exact sign normalisation for the 12 one-dimensional bands of estimator (iii).

Problem found on resume: e6_run.node_maxt / band_cos score 1-D bands (SVk alone) through
stats3.l2norm_rows / (|x| + 1e-12) normalisation, so a 1-D 'cosine' is sign(x) * (1 - O(1e-12/|x|)) and exact
ties are broken by magnitude (known pitfall in the task rules). Observed and null were computed the same way, so
the old max-T tests were valid tests of a slightly different statistic, but not of the cosine statistic.

This script recomputes, for every finished (set, embedding) that carries band statistics, the observed value and
all N_MAXT joint null draws of the 12 1-D bands with exact normalisation x -> sign(x), using the SAME null draws
as the original run (node: label permutations regenerated from the same rng seed; PPI: the cached Maslov-Sneppen
draws; edge: the cached within-TF permutations). The 66 multi-dimensional bands are copied unchanged.
Output: work/<set>/<emb>.p1d.npz (tobs_<task>, tnull_<task>, old_tobs1d_<task>). e6_assemble uses these.
"""
import json, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import e6lib as L
import e6_run as R
S = L.S
ONE_D = [j for j, (a, b) in enumerate(L.BANDS) if a == b]          # 12 band indices
N_MAXT = R.N_MAXT


def node_perms(n_lab, k, seed):
    rng = np.random.default_rng(seed); perms = []; done = 0
    while done < N_MAXT:
        b = min(500, N_MAXT - done); perms.append(np.argsort(rng.random((n_lab, b)), axis=0)[:k]); done += b
    return np.concatenate(perms, 1)


def node_scores_sign(sg, tf_mask):
    """Exact LOO mean cosine for a sign vector (cos(i,j) = sg_i*sg_j): integer numerators, exact ties."""
    k = int(tf_mask.sum()); T = float(sg[tf_mask].sum())
    return np.where(tf_mask, (sg * T - sg * sg) / max(k - 1, 1), sg * T / k)


def patch_node(C, tp, ng, seed, z, task):
    lab = np.concatenate([tp, ng]); k = tp.size; n = lab.size
    perms = node_perms(n, k, seed)
    to = z[f"tobs_{task}"].astype(float).copy(); tn = z[f"tnull_{task}"].astype(float).copy(); old = to[ONE_D].copy()
    for j in ONE_D:
        a = L.BANDS[j][0]; sg = np.sign(C[lab, a - 1])
        m0 = np.zeros(n, bool); m0[:k] = True
        s = node_scores_sign(sg, m0); to[j] = S.auroc(s[m0], s[~m0])
        # check against stats3.tf_scores_loo on the sign vector (rounded to remove fp noise)
        s_ref = np.round(S.tf_scores_loo(sg[:, None].astype(float), np.arange(k)), 9)
        assert abs(S.auroc(s_ref[m0], s_ref[~m0]) - to[j]) < 1e-12, "node sign check"
        for r in range(N_MAXT):
            m = np.zeros(n, bool); m[perms[:, r]] = True
            s = node_scores_sign(sg, m); tn[r, j] = S.auroc(s[m], s[~m])
    return to, tn, old


def auc_counts(pos, neg):
    """AUROC for values in {-1, 0, 1} (exact, tie-aware)."""
    pc = np.array([(pos == v).sum() for v in (-1, 0, 1)], float); nc = np.array([(neg == v).sum() for v in (-1, 0, 1)], float)
    num = pc[2] * (nc[0] + nc[1]) + pc[1] * nc[0] + 0.5 * (pc * nc).sum()
    return float(num / (pc.sum() * nc.sum()))


def patch_ppi(C, gs, setname, z):
    pi, pj, ni, nj = gs["pi"], gs["pj"], gs["ni"], gs["nj"]
    seed = R.sd(setname, "ppi") + 3
    f = L.WORK / f"msdraws_{setname}_{seed}_{N_MAXT}.npz"
    D = np.load(f); a, b, off = D["a"], D["b"], D["off"]
    to = z["tobs_ppi"].astype(float).copy(); tn = z["tnull_ppi"].astype(float).copy(); old = to[ONE_D].copy()
    SG = np.sign(C).astype(np.int8)
    for j in ONE_D:
        c = L.BANDS[j][0] - 1; sg = SG[:, c]
        cp = (sg[pi] * sg[pj]).astype(float); cn = (sg[ni] * sg[nj]).astype(float)
        to[j] = S.auroc(cp, cn); assert abs(auc_counts(cp, cn) - to[j]) < 1e-12
        for r in range(off.size - 1):
            x = (sg[a[off[r]:off[r + 1]]] * sg[b[off[r]:off[r + 1]]]).astype(float)
            tn[r, j] = auc_counts(x, cn)
            if r < 2: assert abs(S.auroc(x, cn) - tn[r, j]) < 1e-12
    return to, tn, old


def patch_edge(C, gs, setname, z):
    sp = S.StratifiedPairs(gs["e_tf"], gs["e_tg"], gs["e_y"])
    Y = L.edge_perms_cached(setname, sp, R.N_EDGE, R.sd(setname, "edge") + 4)[:, :N_MAXT]
    to = z["tobs_edge"].astype(float).copy(); tn = z["tnull_edge"].astype(float).copy(); old = to[ONE_D].copy()
    for j in ONE_D:
        c = L.BANDS[j][0] - 1; s = np.sign(C[sp.tf, c]) * np.sign(C[sp.tg, c])
        to[j] = sp.auroc(s); fs = L.FastStrat(sp, s)
        tn[:, j] = np.concatenate([fs.perm_null(Y[:, q:q + 250]) for q in range(0, N_MAXT, 250)])
        for r in range(2):
            assert abs(S.StratifiedPairs(sp.tf, sp.tg, Y[:, r]).auroc(s) - tn[r, j]) < 1e-9, "edge sign check"
    return to, tn, old


def patch_one(setname, emb):
    d = L.WORK / setname; fz = d / f"{emb}.npz"; fo = d / f"{emb}.p1d.npz"
    if fo.exists() and fo.stat().st_mtime > fz.stat().st_mtime: return False
    z = np.load(fz)
    tasks = [t for t in ("node_lambert", "node_trrust", "ppi", "edge") if f"tobs_{t}" in z.files]
    if not tasks: return False
    t0 = time.time()
    gs = L.load_geneset(setname); P = R.provide(emb, gs["genes"]); C = np.asarray(P["C"], float)
    out = {}
    for t in tasks:
        if t.startswith("node"):
            msk = "lambert" if t == "node_lambert" else None
            tp = np.flatnonzero(gs["lambert"]) if msk else np.flatnonzero(gs["trrust_tf"])
            ng = np.flatnonzero(~gs["lambert"]) if msk else np.flatnonzero(gs["trrust_tgt"])
            to, tn, old = patch_node(C, tp, ng, R.sd(setname, t) + 10, z, t)
        elif t == "ppi":
            to, tn, old = patch_ppi(C, gs, setname, z)
        else:
            to, tn, old = patch_edge(C, gs, setname, z)
        # untouched multi-dimensional bands must be identical
        md = np.setdiff1d(np.arange(len(L.BANDS)), ONE_D)
        assert np.array_equal(to[md], z[f"tobs_{t}"][md])
        out[f"tobs_{t}"] = to; out[f"tnull_{t}"] = tn.astype(np.float32); out[f"old_tobs1d_{t}"] = old
    tmp = d / f"{emb}.p1d.tmp.npz"; np.savez(tmp, **out); tmp.rename(fo)
    print(f"patched {setname} {emb} {tasks} {time.time()-t0:.0f}s", flush=True)
    return True


def main():
    for d in sorted(L.WORK.iterdir()):
        if not d.is_dir() or d.name.startswith("._"): continue
        for fj in sorted(d.glob("*.json")):
            if fj.name.startswith("._"): continue
            emb = fj.stem
            if (d / f"{emb}.npz").exists():
                patch_one(d.name, emb)


if __name__ == "__main__":
    if len(sys.argv) > 2: patch_one(sys.argv[1], sys.argv[2])
    else: main()
