"""F1 null runs (resumable per context x step; each step writes _ckpt/<ctx>_<step>.npz and is skipped if present).

steps:
  setup  observed-gene SVD coordinates (SV1-SV10 x 12 layers), compartment classes, 30 GO-CC indicators,
         residualised coordinates (30 CC + log n_obs; and, as an extra, 7 class dummies + log n_obs)
  perm   compartment-stratified null: one within-class permutation of genes applied to all 120 statistics, N=2000
  ms     Maslov-Sneppen (stats3.maslov_sneppen_chain, burn 10m, thin m), N=2000; on each draw: raw, resid30 and
         residclass assortativity (360 statistics) and the same-class edge fraction
  cpms   class-pair-preserving MS (f1_common.ms_classpair_chain, burn 10m, thin m accepted swaps), N=2000, raw;
         degree sequence and class-pair counts checked on every draw
  wb     within-class edges and between-class edges, each with its own MS chain (N=1000); raw and class-centred
  nobs   V restricted to n_obs >= 3, 5, 10: SVD recomputed on restricted V, STRING edges within it, MS N=2000
Usage: python f1_run.py --ctx lung --steps setup,perm,ms,cpms,wb,nobs
"""
import argparse, time
import numpy as np
import f1_common as F

E, S = F.E, F.S
N_MAIN, N_WB = 2000, 1000
NOBS_T = (3, 5, 10)


def ck(ctx, step):
    return F.CKPT / f"{ctx}_{step}.npz"


def setup(ctx, log):
    p = ck(ctx, "setup")
    if p.exists():
        return dict(np.load(p, allow_pickle=True))
    t0 = time.time()
    c = E.load_ctx(ctx)
    Xc = E.coords(ctx, c["obs_idx"])
    g2go = F.load_gene2go(); info = F.go_terms_info()
    cls = F.assign_classes(c["genes"], g2go)
    ids, names, counts, M, n_cc = F.cc_design(c["genes"], g2go, info)
    lnobs = c["labels"]["log_nobs"].astype(np.float64)
    assert np.allclose(lnobs, np.log(c["ncells"]))
    R30, r2_30 = F.residualise(Xc, np.column_stack([M, lnobs]))
    dummies = np.column_stack([(cls == k).astype(float) for k in np.unique(cls)[1:]])
    Rcl, r2_cl = F.residualise(Xc, np.column_stack([dummies, lnobs]))
    out = dict(Xc=Xc, R30=R30, Rcl=Rcl, r2_30=r2_30, r2_cl=r2_cl, cls=cls, cc_ids=np.array(ids),
               cc_names=np.array(names), cc_counts=np.array(counts), n_cc_annot=n_cc,
               n_any_annot=int(sum(len(g2go.get(g, ())) > 0 for g in c["genes"])),
               edges=c["edges"], ncells=c["ncells"], n=c["n"])
    np.savez(p, **out)
    log(f"[{ctx}] setup {time.time()-t0:.0f}s n={c['n']} m={len(c['edges'])} classes={np.bincount(cls).tolist()}")
    return dict(np.load(p, allow_pickle=True))


def run(ctx, steps):
    log = E.logger(F.HERE / f"f1_run_{ctx}.log")
    E.frozen_hash_check()
    S.assert_resolution(N_MAIN, 0.05, 1); S.assert_resolution(N_WB, 0.05, 1)
    d = setup(ctx, log)
    Xc, R30, Rcl, cls = d["Xc"], d["R30"], d["Rcl"], d["cls"]
    n = int(d["n"]); edges = d["edges"].astype(np.int64); pi, pj = edges[:, 0].copy(), edges[:, 1].copy()
    obs = E.assort_cols(Xc, pi, pj)
    for cc in (0, 21, 119):
        assert abs(obs[cc] - S.assortativity(Xc[:, cc], pi, pj)) < 1e-10
        assert abs(E.assort_cols(R30[:, [cc]], pi, pj)[0] - S.assortativity(R30[:, cc], pi, pj)) < 1e-10

    if "perm" in steps and not ck(ctx, "perm").exists():
        t0 = time.time(); rng = np.random.default_rng(E.seed_for("F1", ctx, "perm"))
        Rn, _ = E.run_perm(Xc, None, pi, pj, n, N_MAIN, rng, strata=cls, log=log, every=500)
        np.savez(ck(ctx, "perm"), obs=obs, null=Rn, N=N_MAIN)
        log(f"[{ctx}] perm done {time.time()-t0:.0f}s")

    if "ms" in steps and not ck(ctx, "ms").exists():
        t0 = time.time()
        ch = S.maslov_sneppen_chain(pi, pj, n, N_MAIN, seed=E.seed_for("F1", ctx, "ms"), burn_factor=10, thin_factor=1)
        (Rr, R3_, Rc), ex, acc, _ = F.run_chain(ch, [Xc, R30, Rcl], extra=lambda a, b: F.same_class_frac(cls, a, b),
                                               log=log, tag="ms")
        np.savez(ck(ctx, "ms"), obs_raw=obs, obs_r30=E.assort_cols(R30, pi, pj), obs_rcl=E.assort_cols(Rcl, pi, pj),
                 null_raw=Rr, null_r30=R3_, null_rcl=Rc, same_class=ex, acc=acc, N=N_MAIN)
        log(f"[{ctx}] ms done {time.time()-t0:.0f}s acc={acc.mean():.3f}")

    if "cpms" in steps and not ck(ctx, "cpms").exists():
        t0 = time.time()
        deg0 = np.bincount(np.r_[pi, pj], minlength=n); ncl = int(cls.max()) + 1
        cp0 = F.class_pair_counts(cls, pi, pj, ncl)

        def chk(a, b):
            return np.array([F.same_class_frac(cls, a, b), np.array_equal(np.bincount(np.r_[a, b], minlength=n), deg0),
                             np.array_equal(F.class_pair_counts(cls, a, b, ncl), cp0),
                             len(np.unique(np.minimum(a, b) * n + np.maximum(a, b))) == a.size and bool(np.all(a != b))], float)
        ch = F.ms_classpair_chain(pi, pj, cls, n, N_MAIN, seed=E.seed_for("F1", ctx, "cpms"))
        (Rr,), ex, acc, ok = F.run_chain(ch, [Xc], extra=chk, log=log, tag="cpms")
        assert ex[:, 1:].all(), "class-pair chain broke an invariant"
        np.savez(ck(ctx, "cpms"), obs=obs, null=Rr, same_class=ex[:, 0], acc=acc, completed=ok, N=N_MAIN)
        log(f"[{ctx}] cpms done {time.time()-t0:.0f}s acc={acc.mean():.3f} completed={ok.mean():.3f}")

    if "wb" in steps and not ck(ctx, "wb").exists():
        t0 = time.time(); Xcc = F.class_centre(Xc, cls); out = {}
        same = cls[pi] == cls[pj]
        for nm, msk in (("within", same), ("between", ~same)):
            a0, b0 = pi[msk].copy(), pj[msk].copy()
            ch = S.maslov_sneppen_chain(a0, b0, n, N_WB, seed=E.seed_for("F1", ctx, "wb", nm), burn_factor=10, thin_factor=1)
            (Rr, Rcc), _, acc, _ = F.run_chain(ch, [Xc, Xcc], log=log, tag=f"wb-{nm}")
            out.update({f"{nm}_obs_raw": E.assort_cols(Xc, a0, b0), f"{nm}_obs_cc": E.assort_cols(Xcc, a0, b0),
                        f"{nm}_null_raw": Rr, f"{nm}_null_cc": Rcc, f"{nm}_acc": acc, f"{nm}_m": a0.size})
            log(f"[{ctx}] wb {nm} m={a0.size} done {time.time()-t0:.0f}s acc={acc.mean():.3f}")
        np.savez(ck(ctx, "wb"), N=N_WB, **out)

    if "nobs" in steps:
        emb = E.load_emb(ctx); obs_idx = E.load_ctx(ctx)["obs_idx"]; nc = d["ncells"]
        for t in NOBS_T:
            p = ck(ctx, f"nobs{t}")
            if p.exists(): continue
            t0 = time.time()
            keep = nc >= t; sub = obs_idx[keep]; nn = int(keep.sum())
            Xs = np.ascontiguousarray(np.concatenate([S.svd_observed(emb, L, sub)[0][:, :E.K_SV]
                                                      for L in range(E.N_LAYERS)], axis=1))
            remap = np.full(n, -1); remap[np.flatnonzero(keep)] = np.arange(nn)
            ek = keep[pi] & keep[pj]; a0, b0 = remap[pi[ek]], remap[pj[ek]]
            o = E.assort_cols(Xs, a0, b0)
            ch = S.maslov_sneppen_chain(a0, b0, nn, N_MAIN, seed=E.seed_for("F1", ctx, "nobs", t), burn_factor=10, thin_factor=1)
            (Rr,), _, acc, _ = F.run_chain(ch, [Xs], log=log, tag=f"nobs{t}")
            np.savez(p, obs=o, null=Rr, acc=acc, n=nn, m=a0.size, Xs=Xs, a=a0, b=b0, N=N_MAIN, keep=keep)
            log(f"[{ctx}] nobs>={t} n={nn} m={a0.size} done {time.time()-t0:.0f}s max r={o.max():.3f}")
    log(f"[{ctx}] steps {steps} complete")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ctx", required=True); ap.add_argument("--steps", default="setup,perm,ms,cpms,wb,nobs")
    a = ap.parse_args()
    for cx in a.ctx.split(","):
        run(cx, a.steps.split(","))
