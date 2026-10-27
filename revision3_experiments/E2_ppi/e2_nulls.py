"""E2 pass 1: the 120-statistic assortativity family (SV1-SV10 x layers 0-11, observed-gene SVD) on the
STRING-physical >= 700 edges within V, under joint null draws.

kinds: ms      Maslov-Sneppen degree-preserving rewiring (PRIMARY null; burn 10m, thin m)
       label   gene-label shuffle (one permutation applied to all 120 statistics)
       deg_v   shuffle within within-V degree bins (stats3.degree_bins, 10 bins)
       deg_g   shuffle within genome-wide STRING-physical(>=700) degree bins
       det     shuffle within detection deciles (cache node_labels det_decile)
Writes _ckpt/<ctx>_<kind>_<N>.npz (obs, null, acc); skips files that exist (resumable per context x null).
Usage: python e2_nulls.py --ctx lung --kinds ms,label,deg_v,deg_g,det --N 10000
"""
import argparse, json, time
import numpy as np
import e2_common as E

S = E.S


def strata_for(kind, c, pi, pj):
    n = c["n"]
    if kind == "deg_v":
        return S.degree_bins(np.bincount(np.r_[pi, pj], minlength=n), 10)
    if kind == "deg_g":
        gd = json.loads((E.CKPT / "genome_degree700.json").read_text())
        return S.degree_bins(np.array([gd.get(g, 0) for g in c["genes"]]), 10)
    if kind == "det":
        return c["labels"]["det_decile"]
    return None


def run(ctx, kinds, N, edges_override=None, tag=""):
    log = E.logger(E.HERE / f"e2_nulls_{ctx}{tag}.log")
    c = E.load_ctx(ctx)
    edges = c["edges"] if edges_override is None else edges_override
    pi, pj = edges[:, 0].copy(), edges[:, 1].copy()
    todo = [k for k in kinds if not (E.CKPT / f"{ctx}{tag}_{k}_{N}.npz").exists()]
    if not todo:
        log(f"[{ctx}{tag}] all done"); return
    t0 = time.time()
    Xc = E.coords(ctx, c["obs_idx"])
    E.check_vectorised(Xc, E.pole_matrix(Xc[:, :4], (0.27,)), pi, pj, cols=(0, 1, 2, 3))
    obs = E.assort_cols(Xc, pi, pj)
    log(f"[{ctx}{tag}] n={c['n']} m={pi.size} N={N} kinds={todo} svd {time.time()-t0:.0f}s; max|r|={np.abs(obs).max():.3f}")
    for k in todo:
        t1 = time.time()
        if k == "ms":
            R, _, acc = E.run_ms(Xc, None, pi, pj, c["n"], N, seed=E.seed_for(ctx + tag, "ms"), log=log)
        else:
            rng = np.random.default_rng(E.seed_for(ctx + tag, k))
            R, _ = E.run_perm(Xc, None, pi, pj, c["n"], N, rng, strata=strata_for(k, c, pi, pj), log=log)
            acc = np.array([])
        np.savez(E.CKPT / f"{ctx}{tag}_{k}_{N}.npz", obs=obs, null=R, acc=acc, n=c["n"], m=pi.size)
        log(f"[{ctx}{tag}] {k} done in {time.time()-t1:.0f}s")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ctx", required=True); ap.add_argument("--kinds", required=True)
    ap.add_argument("--N", type=int, required=True)
    a = ap.parse_args()
    for ctx in a.ctx.split(","):
        run(ctx, a.kinds.split(","), a.N)
