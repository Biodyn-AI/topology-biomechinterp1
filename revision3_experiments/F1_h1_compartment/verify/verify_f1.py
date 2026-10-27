"""Independent verifier recomputation for F1 (does NOT import stats3, e2_common or f1_common).

From raw inputs (local embedding copy, cache/<ctx>/obs_genes.json, obs_idx.npy, ncells.npy, edges_string.npy,
gene2go_all.pkl, go-basic.obo) recomputes, with independently written code:
  V1 raw assortativity of every SVk x layer cell (own centred SVD; sign-invariant statistic) vs strat_null.csv r
  V2 compartment classes, class sizes, same-class edge fraction and its stub-matching expectation (classes.csv)
  V3 class-stratified permutation null for the E2 top cell (own permutation code, N = 1000, own seed)
  V4 residualisation on the 30 most frequent direct GO-CC terms + log n_obs (own obo parser, own lstsq) and on the
     8-class factor + log n_obs: residual r and ratio for all cells, medians (residual_assort.csv)
  V5 within-/between-class edges: raw and class-centred r at the E2 top cell (within_between.csv)
  V6 n_obs >= 3/5/10: own SVD on restricted V, top |r| cell r (nobs_threshold_cells.csv)
  V7 vertex-subsampling interval of the lung top-cell r with both scalings (own subsampler, B = 300)
  V8 30 GO-BP + log n_obs median ratio (resid_controls_summary.csv)
Writes verify/verify_f1.json and verify/verify_f1_compare.csv.
Usage: python verify_f1.py --ctx lung,lung_external,kidney_B
"""
import argparse, json, pickle
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

V = Path(__file__).resolve().parent
F1 = V.parent
R3 = F1.parent
REPO = R3.parents[1]
EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"
G2GO = REPO / "single_cell_mechinterp/data/perturb/gene2go_all.pkl"
OBO = REPO / "subproject_61_two_axes_revision/data/go-basic.obo"
TOP = {"lung": (3, 1), "lung_external": (8, 2), "kidney_B": (9, 3)}
PRIO = [["GO:0005615"], ["GO:0005886"], ["GO:0005783", "GO:0005788", "GO:0005794"], ["GO:0005739"],
        ["GO:0005634"], ["GO:0005829"]]


def r_assort(x, a, b):
    u = np.r_[x[a], x[b]]; v = np.r_[x[b], x[a]]
    return float(np.corrcoef(u, v)[0, 1])


def r_cols(X, a, b):
    U = np.vstack([X[a], X[b]]); W = np.vstack([X[b], X[a]])
    U = U - U.mean(0); W = W - W.mean(0)
    return (U * W).sum(0) / np.sqrt((U * U).sum(0) * (W * W).sum(0))


def coords(emb, rows, layers=range(12), k=10):
    out = []
    for L in layers:
        M = np.array(emb[L][rows], dtype=np.float64); M -= M.mean(0)
        U, s, _ = np.linalg.svd(M, full_matrices=False)
        out.append(U[:, :k] * s[:k])
    return np.hstack(out)


def parse_obo():
    ns, cur, alts = {}, None, []
    with open(OBO) as f:
        for line in f:
            line = line.strip()
            if line.startswith("["):
                cur = None; alts = []
                if line == "[Term]": cur = {}
            elif cur is None: continue
            elif line.startswith("id: "): cur["id"] = line[4:]
            elif line.startswith("namespace: "):
                cur["ns"] = line[11:]; ns[cur["id"]] = cur["ns"]
                for a_ in alts: ns.setdefault(a_, cur["ns"])
            elif line.startswith("alt_id: "):
                alts.append(line[8:])
                if "ns" in cur: ns.setdefault(line[8:], cur["ns"])
    return ns


def resid(X, D):
    Z = np.column_stack([np.ones(len(X)), D])
    beta, *_ = np.linalg.lstsq(Z, X, rcond=None)
    return X - Z @ beta


def top_terms(genes, g2go, ns, space, k=30):
    per = [{t for t in g2go.get(g, ()) if ns.get(t) == space} for g in genes]
    cnt = Counter(t for s in per for t in s)
    ids = [t for t, _ in sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0]))[:k]]
    return np.array([[t in s for t in ids] for s in per], float), ids


def main(ctxs):
    g2go = pickle.load(open(G2GO, "rb")); ns = parse_obo()
    out, cmp = {}, []
    st = pd.read_csv(F1 / "strat_null.csv"); cl = pd.read_csv(F1 / "classes.csv")
    rs = pd.read_csv(F1 / "residual_assort.csv"); wb = pd.read_csv(F1 / "within_between.csv")
    nc_ = pd.read_csv(F1 / "nobs_threshold_cells.csv"); rcs = pd.read_csv(F1 / "resid_controls_summary.csv")
    ci = pd.read_csv(F1 / "intervals.csv")

    def add(ctx, what, mine, theirs, tol):
        cmp.append(dict(context=ctx, quantity=what, verifier=mine, implementer=theirs, abs_diff=abs(mine - theirs),
                        tol=tol, ok=bool(abs(mine - theirs) <= tol)))

    for ctx in ctxs:
        d = R3 / "cache" / ctx
        genes = json.loads((d / "obs_genes.json").read_text()); oi = np.load(d / "obs_idx.npy")
        nc = np.load(d / "ncells.npy").astype(float); E = np.load(d / "edges_string.npy").astype(np.int64)
        emb = np.load(EMB / f"{ctx}.npy", mmap_mode="r")
        nz = np.flatnonzero(np.abs(np.asarray(emb[0])).sum(1) > 0)
        assert np.array_equal(nz, oi), "obs_idx != non-zero rows"
        assert len(genes) == oi.size == nc.size and nc.min() >= 1 and E.max() < oi.size and np.all(E[:, 0] != E[:, 1])
        und = np.unique(np.sort(E, 1), axis=0); assert len(und) == len(E), "duplicate edges"
        a, b = E[:, 0], E[:, 1]; n = len(genes)
        X = coords(emb, oi)
        r = r_cols(X, a, b)
        s_ = st[st.context == ctx].sort_values(["layer", "sv"])
        add(ctx, "V1 max |r diff| over 120 cells", float(np.abs(r - s_.r.values).max()), 0.0, 1e-9)
        L0, k0 = TOP[ctx]; c0 = L0 * 10 + k0 - 1
        add(ctx, f"V1 r L{L0} SV{k0}", float(r[c0]), float(s_.r.values[c0]), 1e-9)

        # V2 classes
        cls = np.full(n, -1)
        for i, g in enumerate(genes):
            ann = set(g2go.get(g, ()))
            for ci_, terms in enumerate(PRIO):
                if ann.intersection(terms): cls[i] = ci_; break
            if cls[i] < 0: cls[i] = 6 if ann else 7
        sizes = np.bincount(cls, minlength=8)
        theirs = cl[(cl.context == ctx) & (cl.cls != "ALL_same_class")].n_genes.values
        add(ctx, "V2 class sizes max |diff|", float(np.abs(sizes - theirs).max()), 0.0, 0)
        same = cls[a] == cls[b]; stub = np.bincount(cls[np.r_[a, b]], minlength=8) / (2 * len(a))
        allrow = cl[(cl.context == ctx) & (cl.cls == "ALL_same_class")].iloc[0]
        add(ctx, "V2 same-class edge fraction", float(same.mean()), float(allrow.frac_edges_within_this_class), 1e-12)
        add(ctx, "V2 stub-matching expectation", float((stub ** 2).sum()), float(allrow.config_expected_frac_edges_within_this_class), 1e-12)

        # V3 class-stratified permutation, E2 top cell, own code
        rng = np.random.default_rng(777001 + len(ctx)); x = X[:, c0]; groups = [np.flatnonzero(cls == c) for c in range(8)]
        nul = np.empty(1000)
        for t in range(1000):
            xp = x.copy()
            for g in groups: xp[g] = x[rng.permutation(g)]
            nul[t] = r_assort(xp, a, b)
        row = s_.iloc[c0]
        add(ctx, "V3 strat-perm null mean (top cell; N=1000 vs 2000)", float(nul.mean()), float(row.strat_null_mean), 0.003)
        add(ctx, "V3 strat-perm null sd (top cell)", float(nul.std(ddof=1)), float(row.strat_null_sd), 0.0015)
        add(ctx, "V3 strat-perm excess / MS excess (top cell)", float((r[c0] - nul.mean()) / (r[c0] - row.f1ms_null_mean)),
            float(row.strat_excess_over_ms_excess), 0.01)

        # V4 residualisation
        CC, ids = top_terms(genes, g2go, ns, "cellular_component")
        ln = np.log(nc)
        R30 = resid(X, np.column_stack([CC, ln])); r30 = r_cols(R30, a, b)
        dm = np.column_stack([(cls == c).astype(float) for c in range(1, 8)])
        Rcl = resid(X, np.column_stack([dm, ln])); rcl = r_cols(Rcl, a, b)
        for des, rr in (("resid30cc", r30), ("residclass", rcl)):
            q = rs[(rs.context == ctx) & (rs.design == des)].sort_values(["layer", "sv"])
            add(ctx, f"V4 {des} max |r_resid diff| over 120", float(np.abs(rr - q.r_resid.values).max()), 0.0, 1e-8)
            add(ctx, f"V4 {des} ratio top cell", float(rr[c0] / r[c0]), float(q.ratio_resid_over_raw.values[c0]), 1e-8)
            add(ctx, f"V4 {des} median ratio over 120", float(np.median(rr / r)), float(q.ratio_resid_over_raw.median()), 1e-8)
        BP, _ = top_terms(genes, g2go, ns, "biological_process")
        rbp = r_cols(resid(X, np.column_stack([BP, ln])), a, b)
        q = rcs[(rcs.context == ctx) & (rcs.design == "bp30+lnobs")].iloc[0]
        add(ctx, "V8 bp30+lnobs median ratio", float(np.median(rbp / r)), float(q.median_ratio), 1e-8)
        add(ctx, "V8 frac edges sharing >=1 cc30 term", float(((CC[a] * CC[b]).sum(1) > 0).mean()),
            float(rcs[(rcs.context == ctx) & (rcs.design == "cc30+lnobs")].iloc[0].frac_edges_sharing_a_term), 1e-12)

        # V5 within / between
        Xcc = X.copy()
        for c in range(8): Xcc[cls == c] -= Xcc[cls == c].mean(0)
        for part, m_ in (("within", same), ("between", ~same)):
            for nm, M in (("raw", X), ("class_centred", Xcc)):
                q = wb[(wb.context == ctx) & (wb.edges == part) & (wb.coords == nm)].sort_values(["layer", "sv"])
                add(ctx, f"V5 {part} {nm} r top cell", r_assort(M[:, c0], a[m_], b[m_]), float(q.r.values[c0]), 1e-9)

        # V6 n_obs thresholds
        for t_ in (3, 5, 10):
            keep = nc >= t_; Xs = coords(emb, oi[keep]); mp = -np.ones(n, int); mp[keep] = np.arange(keep.sum())
            ek = keep[a] & keep[b]; rr = r_cols(Xs, mp[a[ek]], mp[b[ek]])
            q = nc_[(nc_.context == ctx) & (nc_.nobs_min == t_)].sort_values(["layer", "sv"])
            add(ctx, f"V6 nobs>={t_} n_V", float(keep.sum()), float(q.n_V.iloc[0]), 0)
            add(ctx, f"V6 nobs>={t_} m", float(ek.sum()), float(q.m.iloc[0]), 0)
            add(ctx, f"V6 nobs>={t_} max |r diff| over 120", float(np.abs(rr - q.r.values).max()), 0.0, 1e-9)

        # V7 interval (lung top cell only; own subsampler)
        if ctx == "lung":
            rng = np.random.default_rng(99); mm = n // 2; dr = []
            for _ in range(300):
                kp = np.zeros(n, bool); kp[rng.choice(n, mm, replace=False)] = True; e = kp[a] & kp[b]
                dr.append(r_assort(x, a[e], b[e]))
            dr = np.array(dr); dev = dr - dr.mean()
            for nm, sc in (("prereg", np.sqrt(mm / n)), ("deleted", np.sqrt(mm / (n - mm)))):
                lo, hi = r[c0] - np.quantile(dev * sc, 0.975), r[c0] - np.quantile(dev * sc, 0.025)
                q = ci[(ci.context == ctx) & (ci.analysis == "raw_all_edges") & (ci.layer == L0) & (ci.sv == k0)].iloc[0]
                add(ctx, f"V7 CI lo {nm} (B=300, own draws)", float(lo), float(q[f"ci_lo_{nm}"]), 0.01)
                add(ctx, f"V7 CI hi {nm} (B=300, own draws)", float(hi), float(q[f"ci_hi_{nm}"]), 0.01)
        out[ctx] = dict(cc30_ids=ids, n=n, m=int(len(a)))
        np.savez(V / f"_v_{ctx}.npz", X=X, cls=cls, a=a, b=b)
        print(ctx, "done", flush=True)
    df = pd.DataFrame(cmp); df.to_csv(V / "verify_f1_compare.csv", index=False)
    out["n_checks"] = len(df); out["n_ok"] = int(df.ok.sum())
    (V / "verify_f1.json").write_text(json.dumps(out, indent=1))
    pd.set_option("display.width", 250); print(df.to_string())


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--ctx", default="lung,lung_external,kidney_B")
    main(ap.parse_args().ctx.split(","))
