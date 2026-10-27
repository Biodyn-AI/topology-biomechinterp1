"""E2 sensitivity (i): kidney_A bridge from the submitted version's PPI statistic to the round-3 primary H1.

Cells: SV2, SV3, SV4 x layers 0-11 (36 statistics) at every step; factors changed one at a time, cumulatively:
  S0 submitted  : zero-row full-matrix SVD (all 4,803 rows, centred over all rows, round-1 `layer_us`), restricted
                  to the 187 legacy genes; functional STRING (iter_0015 cache) score >= 0.7 among them; co-pole
                  rate at K = floor(0.27 n); gene-label shuffle over the 187 genes (one-sided, greater)
  S1 +obsSVD    : observed-gene SVD (stats3.svd_observed on V), coordinates restricted to the 187 genes
  S2 +V         : universe = V (1,972 genes): poles over V, labels shuffled over V; same functional pairs
  S3 +physical  : STRING v12 physical >= 700 edges within V (cache edges_string.npy)
  S4 +assort    : scalar assortativity r instead of co-pole (label shuffle, two-sided)
  S5 +MS        : Maslov-Sneppen null (= the primary H1 cells SV2-SV4 of kidney_A; read from the primary run)
Off-path check: S1p = S1 with STRING physical >= 700 among the 187 legacy genes (network type at the legacy
universe). Label shuffles: N = 10,000, one permutation applied to all 36 statistics per draw (joint); max-T over
the 36 per step is descriptive. Writes bridge_kidneyA.csv, bridge_kidneyA_cells.csv.
"""
import json
import numpy as np, pandas as pd
import e2_common as E

S, X = E.S, E.X
N = 10000
CELLS = [(L, k) for L in range(12) for k in (2, 3, 4)]
BONF183 = 0.05 / 183


def us_full_zero_rows(emb, L):
    """Round-1 / round-2 `common.layer_us`: SVD of the FULL gene x 512 matrix (zero rows included)."""
    Xf = np.asarray(emb[L], dtype=np.float64); Xf = Xf - Xf.mean(axis=0, keepdims=True)
    U, Sv, _ = np.linalg.svd(Xf, full_matrices=False)
    return U * Sv


def family(name, Xs, pi, pj, stat, tail, seed, n):
    """Observed + joint label-shuffle null for 36 columns. stat: 'copole' or 'assort'."""
    Lm = E.pole_matrix(Xs, (0.27,))
    rng = np.random.default_rng(seed)
    R, C = E.run_perm(Xs, Lm if stat == "copole" else None, pi, pj, n, N, rng)
    if stat == "copole":
        obs = E.copole_cols(Lm, pi, pj); nul = C.astype(np.float64)
        for j in range(3):
            assert abs(obs[j] - S.copole_rate(Lm[:, j], pi, pj)) < 1e-12
    else:
        obs = E.assort_cols(Xs, pi, pj); nul = R
    return obs, nul


def summarise(step, desc, obs, nul, tail, n, m, N_):
    sm = E.null_summary(obs, nul, tail=tail)
    if tail == "greater":
        pr = np.array([S.emp_p(obs[j], nul[:, j], "greater")[0] for j in range(obs.size)])
        assert np.allclose(pr, sm["p_raw"])
    rows = [dict(step=step, layer=L, sv=k, stat=obs[j], null_mean=sm["null_mean"][j], null_sd=sm["null_sd"][j],
                 z=sm["z"][j], p_raw=sm["p_raw"][j], p_adj36=sm["p_adj"][j], n_genes=n, n_pairs=m, N=N_, tail=tail)
            for j, (L, k) in enumerate(CELLS)]
    j = int(np.argmax(sm["z"] if tail == "greater" else np.abs(sm["z"])))
    summ = dict(step=step, description=desc, n_genes=n, n_pairs=m, N=N_, tail=tail,
                n_cells_p_lt_05=int((sm["p_raw"] < 0.05).sum()), n_cells_p_lt_bonf183=(int((sm["p_raw"] < BONF183).sum()) if 1 / (N_ + 1) < BONF183 else np.nan),
                bonf183_resolvable=bool(1 / (N_ + 1) < BONF183), p_floor=1 / (N_ + 1),
                best_layer=CELLS[j][0], best_sv=CELLS[j][1], best_stat=obs[j], best_null_mean=sm["null_mean"][j],
                best_ratio=obs[j] / sm["null_mean"][j] if tail == "greater" else np.nan,
                best_z=sm["z"][j], best_p_raw=sm["p_raw"][j], min_p_adj36=sm["p_adj"].min(),
                SV3_L5_stat=obs[CELLS.index((5, 3))], SV3_L5_z=sm["z"][CELLS.index((5, 3))],
                mean_stat_36=float(np.mean(obs)), median_z_36=float(np.median(sm["z"])))
    return rows, summ


def main():
    log = E.logger(E.HERE / "e2_bridge.log")
    ctx = "kidney_A"; emb = E.load_emb(ctx)
    genes_all = X.h5ad_genes(X.CONTEXTS[ctx]["h5ad"])
    c = E.load_ctx(ctx); V = c["genes"]; posV = {g: i for i, g in enumerate(V)}; nV = c["n"]
    legacy = json.loads((X.REPO / "subproject_41_claude_topology_hypothesis_screening_autoloop/revision3_experiments/"
                         "E6_crossmodel/legacy_genes.json").read_text())["curated_nonzero_kidneyA"]
    assert len(legacy) == 187 and all(g in posV for g in legacy)
    gidx_all = np.array([genes_all.index(g) for g in legacy])
    loc = {g: k for k, g in enumerate(legacy)}
    obj = json.loads((X.S41 / "iter_0015/string_ppi_score04_cache.json").read_text())
    st = sorted({(p["g1"], p["g2"]) for p in obj["pairs"] if p["score"] >= 0.7 and p["g1"] in loc and p["g2"] in loc
                 and p["g1"] != p["g2"]})                        # exactly as revision2 rebuild_headline.py
    n_unordered = len({tuple(sorted(e)) for e in st})
    fi = np.array([loc[a] for a, b in st]); fj = np.array([loc[b] for a, b in st])
    log(f"legacy genes 187; functional pairs (submitted construction) {len(st)}, unordered {n_unordered}")
    # coordinates
    USf = {L: us_full_zero_rows(emb, L)[gidx_all] for L in range(12)}
    Xfull = np.stack([USf[L][:, k - 1] for (L, k) in CELLS], 1)
    Xobs = E.coords(ctx, c["obs_idx"])                                     # n_V x 120
    colsV = [L * 10 + (k - 1) for (L, k) in CELLS]
    XobsV = Xobs[:, colsV]
    leg_in_V = np.array([posV[g] for g in legacy])
    Xobs187 = XobsV[leg_in_V]
    rows, summ = [], []
    # S0
    obs, nul = family("S0", Xfull, fi, fj, "copole", "greater", E.seed_for("bridge", "S0"), 187)
    ref = pd.read_csv(X.S41.parent / "revision2_experiments/headline_rebuilt.csv")
    chk = []
    for j, (L, k) in enumerate(CELLS):
        v = ref.loc[ref.claim == f"SV{k}_STRING_copole_L{L}", "stat"].values[0]; chk.append(abs(v - obs[j]))
    log(f"S0 reproduces submitted/round-2 observed co-pole: max |diff| = {max(chk):.2e}")
    assert max(chk) < 1e-12
    r_, s_ = summarise("S0_submitted", "zero-row full SVD; 187 legacy genes; functional STRING>=0.7; co-pole K=0.27n; label shuffle",
                       obs, nul, "greater", 187, len(st), N); rows += r_; summ.append(s_)
    # S1
    obs, nul = family("S1", Xobs187, fi, fj, "copole", "greater", E.seed_for("bridge", "S1"), 187)
    r_, s_ = summarise("S1_observed_SVD", "S0 + observed-gene SVD (V), restricted to 187", obs, nul, "greater", 187, len(st), N)
    rows += r_; summ.append(s_)
    # S1p off-path: physical STRING among legacy genes
    eV = c["edges"]; legset = set(leg_in_V.tolist()); invleg = {v: k for k, v in enumerate(leg_in_V)}
    pl = np.array([(invleg[a], invleg[b]) for a, b in eV if a in legset and b in legset]).reshape(-1, 2)
    obs, nul = family("S1p", Xobs187, pl[:, 0], pl[:, 1], "copole", "greater", E.seed_for("bridge", "S1p"), 187)
    r_, s_ = summarise("S1p_offpath_physical_on_187", "S1 + physical STRING>=700 among the 187 (off-path check)",
                       obs, nul, "greater", 187, len(pl), N); rows += r_; summ.append(s_)
    # S2
    vi = leg_in_V[fi]; vj = leg_in_V[fj]
    obs, nul = family("S2", XobsV, vi, vj, "copole", "greater", E.seed_for("bridge", "S2"), nV)
    r_, s_ = summarise("S2_V_universe", "S1 + universe V (poles and shuffle over V); same functional pairs",
                       obs, nul, "greater", nV, len(st), N); rows += r_; summ.append(s_)
    # S3
    pi, pj = eV[:, 0].copy(), eV[:, 1].copy()
    obs, nul = family("S3", XobsV, pi, pj, "copole", "greater", E.seed_for("bridge", "S3"), nV)
    r_, s_ = summarise("S3_physical_STRING", "S2 + STRING physical>=700 within V", obs, nul, "greater", nV, pi.size, N)
    rows += r_; summ.append(s_)
    # S4
    obs, nul = family("S4", XobsV, pi, pj, "assort", "two-sided", E.seed_for("bridge", "S4"), nV)
    r_, s_ = summarise("S4_assortativity", "S3 + scalar assortativity r (label shuffle, two-sided)",
                       obs, nul, "two-sided", nV, pi.size, N); rows += r_; summ.append(s_)
    # S5 (primary run)
    z = np.load(E.CKPT / f"{ctx}_ms_{E.N_PRIMARY[ctx]}.npz")
    assert np.allclose(z["obs"][colsV], obs, rtol=0, atol=1e-10)
    r_, s_ = summarise("S5_Maslov_Sneppen", "S4 + Maslov-Sneppen null (primary H1 kidney_A run, N=2,000; max-T here over SV2-SV4 only)",
                       z["obs"][colsV], z["null"][:, colsV], "two-sided", nV, pi.size, int(z["null"].shape[0]))
    rows += r_; summ.append(s_)
    pd.DataFrame(rows).to_csv(E.HERE / "bridge_kidneyA_cells.csv", index=False)
    pd.DataFrame(summ).to_csv(E.HERE / "bridge_kidneyA.csv", index=False)
    meta = dict(functional_pairs_submitted=len(st), functional_pairs_unordered=n_unordered, physical_pairs_187=int(len(pl)),
                s0_reproduction_max_abs_diff=float(max(chk)))
    (E.CKPT / "bridge_meta.json").write_text(json.dumps(meta, indent=1))
    log(json.dumps(meta)); log(pd.DataFrame(summ).to_string())


if __name__ == "__main__":
    main()
