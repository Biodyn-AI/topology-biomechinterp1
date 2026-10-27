"""Re-bootstrap (B = 1000, independent seed) the 9 core-context near-degenerate flags among SV1..SV4
(gaps.csv, k <= 3) to test whether the B = 200 flags are Monte Carlo stable. Own Gram-spectrum code.
Writes verify/gap_B1000.csv."""
from pathlib import Path
import numpy as np, pandas as pd
from joblib import Parallel, delayed
H = Path(__file__).resolve().parent
R3 = H.parent.parent
EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"
B = 1000
c = pd.read_csv(H / "gap_flag_stability_core_k1to3.csv")

def spec(X):
    X = X - X.mean(0)
    return np.sqrt(np.clip(np.linalg.eigvalsh(X.T @ X)[::-1], 0, None))

def one(i, ctx, L, k):
    emb = np.load(EMB / f"{ctx}.npy", mmap_mode="r"); oi = np.load(R3 / "cache" / ctx / "obs_idx.npy")
    X = np.asarray(emb[L], np.float64)[oi]; s = spec(X); r = s[k - 1] / s[k]
    rng = np.random.default_rng([4242, i])
    rb = np.empty(B)
    for b in range(B):
        sb = spec(X[rng.integers(0, oi.size, oi.size)]); rb[b] = sb[k - 1] / sb[k]
    q = np.quantile(rb, [0.025, 0.975])
    return dict(context=ctx, layer=L, k=k, ratio=r, B=B, boot_mean=rb.mean(), boot_sd=rb.std(ddof=1),
                basic_lo_B1000=2 * r - q[1], basic_hi_B1000=2 * r - q[0], flag_B1000=bool(2 * r - q[1] <= 1))

rows = Parallel(n_jobs=3)(delayed(one)(i, r.context, int(r.layer), int(r.k)) for i, r in enumerate(c.itertuples()))
out = pd.DataFrame(rows).merge(c[["context", "layer", "k", "basic_lo", "near_degenerate"]], on=["context", "layer", "k"])
out.to_csv(H / "gap_B1000.csv", index=False); print(out.round(4).to_string())
