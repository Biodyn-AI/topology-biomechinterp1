"""Independent verifier for E2_ppi, part C.
 1. Vertex-subsampling intervals for each primary context's top cell (my implementation, B = 500, m = n/2):
    CONVENTION3 s6 sqrt(m/n) rescaling vs delete-d jackknife scaling (factor 1 when d = n/2), and the
    delete-1 jackknife SE reported in sens_jackknife_summary.csv.
 2. kidney_A: does the zero-row full-matrix SVD shift every index by one? |cos| of left vectors restricted to V
    between full-SVD SV(k+1) and observed-gene SVk, k = 1..9, every layer.
Writes verify/verify_c.json.
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_a import var_names, coords_eig, H5, EMB, CACHE, E2

HERE = Path(__file__).resolve().parent


def main():
    res = {}
    hs = pd.read_csv(E2 / "h1_summary.csv").set_index("context")
    jk = pd.read_csv(E2 / "sens_jackknife_summary.csv").set_index("context")
    for ctx in hs.index:
        E = np.load(EMB / f"{ctx}.npy", mmap_mode="r")
        obs = np.load(CACHE / ctx / "obs_idx.npy"); n = obs.size
        ed = np.load(CACHE / ctx / "edges_string.npy"); a, b = ed[:, 0], ed[:, 1]
        L, k = int(hs.loc[ctx, "top_layer"]), int(hs.loc[ctx, "top_sv"])
        x = coords_eig(E, L, obs)[:, k - 1]
        r0 = np.corrcoef(np.r_[x[a], x[b]], np.r_[x[b], x[a]])[0, 1]
        rng = np.random.default_rng(2024); m = n // 2; dr = []
        for _ in range(500):
            keep = np.zeros(n, bool); keep[rng.choice(n, m, replace=False)] = True; e = keep[a] & keep[b]
            dr.append(np.corrcoef(np.r_[x[a[e]], x[b[e]]], np.r_[x[b[e]], x[a[e]]])[0, 1])
        dr = np.array(dr); dev = dr - dr.mean()
        res[ctx] = dict(cell=f"L{L} SV{k}", r=float(r0),
                        ci_sqrt_mn_mine=[float(r0 - np.quantile(dev * np.sqrt(m / n), q)) for q in (0.975, 0.025)],
                        ci_csv=[float(hs.loc[ctx, "top_r_ci_lo"]), float(hs.loc[ctx, "top_r_ci_hi"])],
                        ci_delete_d=[float(r0 - np.quantile(dev, q)) for q in (0.975, 0.025)],
                        se_implied_sqrt_mn=float(dev.std(ddof=1) * np.sqrt(m / n)), se_delete_d=float(dev.std(ddof=1)),
                        se_jackknife_delete1_csv=float(jk.loc[ctx, "jackknife_se"]),
                        null_mean_ms=float(hs.loc[ctx, "top_null_mean"]))
        print(ctx, res[ctx], flush=True)
    # 2. kidney_A index shift
    ctx = "kidney_A"; E = np.load(EMB / f"{ctx}.npy", mmap_mode="r"); obs = np.load(CACHE / ctx / "obs_idx.npy")
    sh = []
    for L in range(12):
        Xf = np.asarray(E[L], np.float64); Xf = Xf - Xf.mean(0)
        U, S, _ = np.linalg.svd(Xf, full_matrices=False)
        Uf = U[obs, :11]
        Xo = np.asarray(E[L], np.float64)[obs]; Xo = Xo - Xo.mean(0)
        Uo, So, _ = np.linalg.svd(Xo, full_matrices=False)
        cs = [float(abs(np.corrcoef(Uf[:, kk + 1], Uo[:, kk])[0, 1])) for kk in range(9)]
        sh.append(dict(layer=L, corr_full_k1_vs_obs_k=cs, min_corr=min(cs)))
    res["kidney_A_index_shift"] = sh
    res["kidney_A_index_shift_min_over_layers_SV1_4"] = float(min(min(s["corr_full_k1_vs_obs_k"][:4]) for s in sh))
    res["kidney_A_index_shift_min_over_layers_SV1_9"] = float(min(s["min_corr"] for s in sh))
    print({k: v for k, v in res.items() if k.startswith("kidney_A_index")}, flush=True)
    (HERE / "verify_c.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
