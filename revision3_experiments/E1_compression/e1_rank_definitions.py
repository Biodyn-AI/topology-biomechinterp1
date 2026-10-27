"""Three spectral-concentration measures on the observed-gene SVD, L0..L11, all 16 contexts.
er_sv  : exp(entropy of S/sum S)            (Roy & Vetterli 2007; the paper's effective rank)
er_var : exp(entropy of S^2/sum S^2)        (entropy of variance shares; the screen's 'effective rank')
pr     : (sum S^2)^2 / sum S^4              (participation ratio)
Writes rank_definitions.csv."""
import os, sys, json
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import contexts as X, stats3 as S
EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"
C = Path(__file__).resolve().parents[1] / "cache"
rows = []
for ctx in X.CONTEXTS:
    E = np.load(EMB / f"{ctx}.npy", mmap_mode="r") if (EMB / f"{ctx}.npy").exists() else X.load_context(ctx)[0]
    idx = np.load(C / ctx / "obs_idx.npy")
    for L in range(12):
        M = np.asarray(E[L][idx], dtype=np.float64); M = M - M.mean(0)
        s = np.linalg.svd(M, compute_uv=False)
        v = s**2; q = v / v.sum(); q = q[q > 0]
        rows.append(dict(context=ctx, layer=L, n_V=len(idx), er_sv=S.effective_rank(s),
                         er_var=float(np.exp(-(q*np.log(q)).sum())), pr=S.participation_ratio(s)))
    print(ctx, rows[-12]["er_sv"], rows[-1]["er_sv"], rows[-12]["er_var"], rows[-1]["er_var"], flush=True)
df = pd.DataFrame(rows); df.to_csv(Path(__file__).parent / "rank_definitions.csv", index=False)
f = df[df.layer.isin([0, 11])].pivot(index="context", columns="layer", values=["er_sv", "er_var", "pr"])
for m in ["er_sv", "er_var", "pr"]: f[(m, "fold")] = f[(m, 0)] / f[(m, 11)]
print(f.round(2).to_string())
