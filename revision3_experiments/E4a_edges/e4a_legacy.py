"""E4a sensitivity (1): the submitted (round-2) held-out regulatory-edge benchmark, reproduced and re-split.

Answers R1.1 literally. Step 1 reproduces revision2_experiments/run_trees.py's split exactly (same RNG
seed 20260819, same order of RNG consumption: 3,000-cell draw, permutation of the sorted TRRUST pairs,
negative draw), on kidney_A (the matrix of the submitted version; local byte-identical embedding copy),
and checks the reproduced AUROCs and bootstrap CIs against revision2_experiments/GRN_benchmark_final.json
(GENIE3/GRNBoost2 are not refitted: they do not consume the RNG before the split, and their bootstrap draws
are replayed because the draws depend only on the pair counts).
Step 2: 500 random 80/20 splits with the same construction. Scorers: full-512 L0 cosine, |co-expression|,
spectral SV2-SV7 and SV5-SV7 at L0 under the SUBMITTED convention (SVD of the full 4,803-row matrix incl.
never-observed zero rows; uncentred full-512 cosine) and under the PAPER (round-3) convention
(observed-gene SVD on V, V-centred embedding), and random.
NOTE: nothing is trained on the 80%; the spread over splits is the variability of a 20% subsample of
pairs, NOT a confidence interval for the benchmark's AUROC.
Also writes work/legacy_design.npz (the reproduced 49/147 design) for the power simulation.
"""
from __future__ import annotations
import json, os, sys, time
from pathlib import Path
import numpy as np, pandas as pd, h5py

HERE = Path(__file__).resolve().parent; R3 = HERE.parent
R2 = R3.parent / "revision2_experiments"
sys.path.insert(0, str(R3)); sys.path.insert(0, str(HERE)); sys.path.insert(0, str(R2))
import stats3 as S        # noqa: E402
import common as C        # noqa: E402  (round-2 module, read-only)
import e4a_lib as L       # noqa: E402

EMB = Path(os.path.expanduser("~/biomechinterp-local/subproject41_r3/emb/kidney_A.npy"))
N_SPLITS = 500


def main():
    t0 = time.time()
    rng = np.random.default_rng(20260819)
    genes = C.h5ad_genes(); emb = np.load(EMB, mmap_mode="r"); nz0 = C.nonzero_mask(emb, 0)
    sub = np.array([i for i in C.curated_genes(genes) if nz0[i]])
    names = [genes[i] for i in sub]; loc = {g: k for k, g in enumerate(names)}
    t = pd.read_csv(C.TRRUST_PATH, sep="\t", header=None, names=["TF", "target", "mode", "pmid"])
    tr = sorted({(a, b) for a, b in zip(t.TF.astype(str), t.target.astype(str)) if a in loc and b in loc and a != b})
    # ---- expression exactly as run_trees.py (same rng call)
    with h5py.File(C.H5AD, "r") as f:
        gs = [g.decode() if isinstance(g, bytes) else str(g) for g in f["var"]["_index"][:]]
        gi = {g: i for i, g in enumerate(gs)}
        keep = [g for g in names if g in gi]; cols = [gi[g] for g in keep]
        Xg = f["X"]
        if isinstance(Xg, h5py.Group):
            data, indices, indptr = Xg["data"][:], Xg["indices"][:], Xg["indptr"][:]
            ncell = indptr.shape[0] - 1
            sel = np.sort(rng.choice(ncell, min(3000, ncell), replace=False))
            colpos = np.full(len(gs), -1); colpos[cols] = np.arange(len(cols))
            E = np.zeros((sel.size, len(cols)), dtype=np.float32)
            for r, ci in enumerate(sel):
                a, b = indptr[ci], indptr[ci + 1]; cp = colpos[indices[a:b]]; k = cp >= 0
                E[r, cp[k]] = data[a:b][k]
        else:
            ncell = Xg.shape[0]
            sel = np.sort(rng.choice(ncell, min(3000, ncell), replace=False))
            E = np.asarray(Xg[sel, :][:, cols], dtype=np.float32)
    assert keep == names
    kloc = {g: k for k, g in enumerate(keep)}
    Zs = (E - E.mean(0)) / (E.std(0) + 1e-9)
    CX = np.abs(Zs.T @ Zs / Zs.shape[0])                            # |co-expression| matrix, run_trees definition
    # ---- geometry: submitted and paper conventions
    U0_sub = C.layer_us(emb, 0)[sub]                                  # zero-row SVD (submitted)
    full0_sub = np.asarray(emb[0], np.float64)[sub]                   # uncentred (submitted)
    obs_idx = np.load(R3 / "cache/kidney_A/obs_idx.npy")
    V = {int(i): r for r, i in enumerate(obs_idx)}; rows_V = np.array([V[int(i)] for i in sub])
    US_V, _, _ = S.svd_observed(emb, 0, obs_idx)
    XV = np.asarray(emb[0], np.float64)[obs_idx]; XV = XV - XV.mean(0)
    U0_pap = US_V[rows_V]; full0_pap = XV[rows_V]
    def cosmat(M):
        Mn = S.l2norm_rows(M); return Mn @ Mn.T
    SC = {"full-512 L0 cosine [submitted: uncentred]": cosmat(full0_sub),
          "full-512 L0 cosine [paper: V-centred]": cosmat(full0_pap),
          "SV2-SV7 L0 cosine [submitted: zero-row SVD]": cosmat(U0_sub[:, 1:7]),
          "SV5-SV7 L0 cosine [submitted: zero-row SVD]": cosmat(U0_sub[:, 4:7]),
          "SV2-SV7 L0 cosine [paper: observed-gene SVD]": cosmat(U0_pap[:, 1:7]),
          "SV5-SV7 L0 cosine [paper: observed-gene SVD]": cosmat(U0_pap[:, 4:7]),
          "|co-expression| (3,000 cells, run_trees definition)": CX}
    # ---- the exact submitted split
    pos_list = sorted(tr); perm = rng.permutation(len(pos_list))
    n_test = max(1, int(0.20 * len(pos_list)))
    held_pos = [pos_list[i] for i in perm[:n_test]]
    train_pos = {pos_list[i] for i in perm[n_test:]}; hp = set(held_pos)
    held_tfs = sorted({a for a, _ in held_pos})
    cand = [(tf, g) for tf in held_tfs for g in names if g != tf and (tf, g) not in train_pos and (tf, g) not in hp]
    held_neg = [cand[i] for i in rng.choice(len(cand), min(len(held_pos) * 3, len(cand)), replace=False)]
    ref = json.loads((R2 / "GRN_benchmark_final.json").read_text())
    refm = {m["method"]: m for m in ref["methods"]}
    def sc_pairs(M, P): return np.array([M[loc[a], loc[b]] for a, b in P])
    # replay the bootstrap stream: GRNBoost2 and GENIE3 bootstraps consume 2 x B integer draws each
    Bb = 2000
    repro = []
    for nm_ref, nm in (("GRNBoost2", None), ("GENIE3", None), ("scGPT full-512 (L0)", "full-512 L0 cosine [submitted: uncentred]"),
                       ("scGPT spectral SV2-SV7 (L0)", "SV2-SV7 L0 cosine [submitted: zero-row SVD]"),
                       ("scGPT spectral SV5-SV7 (L0)", "SV5-SV7 L0 cosine [submitted: zero-row SVD]"),
                       ("|co-expression|", "|co-expression| (3,000 cells, run_trees definition)")):
        if nm is None:
            for _ in range(Bb): rng.integers(0, len(held_pos), len(held_pos)); rng.integers(0, len(held_neg), len(held_neg))
            continue
        p = sc_pairs(SC[nm], held_pos); q = sc_pairs(SC[nm], held_neg); a = S.auroc(p, q)
        bs = np.array([S.auroc(p[rng.integers(0, p.size, p.size)], q[rng.integers(0, q.size, q.size)]) for _ in range(Bb)])
        lo, hi = np.percentile(bs, [2.5, 97.5])
        repro.append(dict(method=nm_ref, scorer=nm, auroc=a, ci_lo=lo, ci_hi=hi, ref_auroc=refm[nm_ref]["auroc"],
                          ref_ci_lo=refm[nm_ref]["ci_lo"], ref_ci_hi=refm[nm_ref]["ci_hi"],
                          exact=bool(abs(a - refm[nm_ref]["auroc"]) < 1e-12 and abs(lo - refm[nm_ref]["ci_lo"]) < 1e-12
                                     and abs(hi - refm[nm_ref]["ci_hi"]) < 1e-12)))
        print(repro[-1], flush=True)
    p = np.array([float(rng.random()) for _ in held_pos]); q = np.array([float(rng.random()) for _ in held_neg])
    repro.append(dict(method="random", scorer="random", auroc=S.auroc(p, q), ref_auroc=refm["random"]["auroc"],
                      exact=bool(abs(S.auroc(p, q) - refm["random"]["auroc"]) < 1e-12)))
    print(repro[-1], flush=True)
    pd.DataFrame(repro).to_csv(HERE / "legacy_split_reproduction.csv", index=False)
    # ---- legacy design for the power simulation (+ all-candidate scores for ICC estimation)
    P_all = held_pos + held_neg
    allc = [(tf, g) for tf in held_tfs for g in names if g != tf]
    np.savez(HERE / "work" / "legacy_design.npz",
             tf=np.array([a for a, _ in P_all]), tg=np.array([b for _, b in P_all]),
             y=np.r_[np.ones(len(held_pos)), np.zeros(len(held_neg))].astype(np.int8),
             s_full=sc_pairs(SC["full-512 L0 cosine [submitted: uncentred]"], P_all),
             s_cx=sc_pairs(SC["|co-expression| (3,000 cells, run_trees definition)"], P_all),
             all_tf=np.array([a for a, _ in allc]), all_tg=np.array([b for _, b in allc]),
             all_s_full=sc_pairs(SC["full-512 L0 cosine [submitted: uncentred]"], allc),
             all_s_cx=sc_pairs(SC["|co-expression| (3,000 cells, run_trees definition)"], allc))
    # ---- 500 random splits, same construction
    cx_nm = "|co-expression| (3,000 cells, run_trees definition)"
    draws = []
    for i in range(N_SPLITS):
        r2 = np.random.default_rng([20261026, 500, i])
        perm = r2.permutation(len(pos_list)); hpos = [pos_list[j] for j in perm[:n_test]]
        trn = {pos_list[j] for j in perm[n_test:]}; hps = set(hpos); htf = sorted({a for a, _ in hpos})
        cnd = [(tf, g) for tf in htf for g in names if g != tf and (tf, g) not in trn and (tf, g) not in hps]
        hneg = [cnd[j] for j in r2.choice(len(cnd), min(len(hpos) * 3, len(cnd)), replace=False)]
        rec = dict(split=i, n_pos=len(hpos), n_neg=len(hneg), n_tf=len(htf))
        au = {}
        for nm, M in SC.items():
            p = sc_pairs(M, hpos); q = sc_pairs(M, hneg); au[nm] = S.auroc(p, q)
        p = r2.random(len(hpos)); q = r2.random(len(hneg)); au["random"] = S.auroc(p, q)
        # within-TF stratified AUROC of the submitted full-512 scorer on the same split (extra column)
        tfs = [a for a, _ in hpos] + [a for a, _ in hneg]; tgs = [b for _, b in hpos] + [b for _, b in hneg]
        spx = S.StratifiedPairs(np.array(tfs), np.array(tgs), np.r_[np.ones(len(hpos)), np.zeros(len(hneg))])
        sfull = np.r_[sc_pairs(SC["full-512 L0 cosine [submitted: uncentred]"], hpos), sc_pairs(SC["full-512 L0 cosine [submitted: uncentred]"], hneg)]
        kk = np.isin(np.array(tfs), spx.tfs)
        rec["stratified_full512_submitted"] = spx.auroc(sfull[kk]); rec["stratified_n_pos"] = int(spx.y.sum())
        for nm, a in au.items():
            rec[f"auroc::{nm}"] = a
            if nm != cx_nm: rec[f"delta::{nm}"] = a - au[cx_nm]
        draws.append(rec)
    dfd = pd.DataFrame(draws); dfd.to_csv(HERE / "legacy_splits_draws.csv", index=False)
    out = []
    for c in dfd.columns:
        if not (c.startswith("auroc::") or c.startswith("delta::") or c == "stratified_full512_submitted"): continue
        v = dfd[c].to_numpy(float)
        kind = "AUROC" if c.startswith("auroc::") or c.startswith("strat") else "delta vs |co-expression|"
        nm = c.split("::", 1)[-1]
        r = dict(scorer=nm, quantity=kind, n_splits=v.size, mean=v.mean(), sd=v.std(ddof=1), p2_5=np.percentile(v, 2.5),
                 p97_5=np.percentile(v, 97.5), min=v.min(), max=v.max(), frac_above_0_5=np.mean(v > 0.5) if kind == "AUROC" else np.mean(v > 0))
        if kind == "AUROC":
            r["percentile_of_0.591"] = 100 * (np.mean(v < 0.591) + 0.5 * np.mean(v == 0.591))
        out.append(r)
    pd.DataFrame(out).to_csv(HERE / "legacy_splits.csv", index=False)
    json.dump(dict(n_trrust_pairs=len(pos_list), n_test=n_test, n_neg=len(held_neg), n_held_tfs=len(held_tfs), n_genes=len(names),
                   n_splits=N_SPLITS, seconds=time.time() - t0, split_reproduced_exactly=all(r["exact"] for r in repro)),
              open(HERE / "work" / "legacy_meta.json", "w"), indent=1, default=float)
    print("done", time.time() - t0, flush=True)


if __name__ == "__main__":
    main()
