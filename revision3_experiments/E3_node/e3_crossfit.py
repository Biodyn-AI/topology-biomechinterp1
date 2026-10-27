"""E3 sensitivity (2): cross-fitted (selection-unbiased) band estimate.

200 random gene halves. On half A the best (band, layer) of the 936-cell 'cos' scan (width 1: raw coordinate,
sign-invariant; width >= 2: LOO cosine within the band; TF centroid from A's TFs only) is chosen by detection-
stratified AUROC; on half B the chosen cell is re-scored with B's own TF centroid (width 1: the sign chosen on A)
and its stratified AUROC is the held-out estimate. The same is done for the full 512-d score (layer chosen on A).
The SVD basis is fitted on all of V (unsupervised; no labels), as in the scan.
usage: python e3_crossfit.py <context> <task: lambert|trrust>
"""
import sys, time
import numpy as np, pandas as pd
import e3_lib as L
S = L.S

ctx, task = sys.argv[1], sys.argv[2]
CI = L.CONTEXTS.index(ctx); NSPLIT = 200
SEED = 20261026 + 1000 * CI + (0 if task == "lambert" else 300) + 60
t0 = time.time()
D = L.load(ctx)
if task == "lambert":
    U = np.arange(len(D["genes"])); y = D["is_tf"]
else:
    U = np.flatnonzero(D["is_trrust_tf"] | D["is_trrust_tgt"]); y = D["is_trrust_tf"][U]
det = D["det"][U]; n = U.size
bands = L.bands12()
USs = [S.svd_observed(D["emb"], l, D["obs_idx"])[0][U, :12].copy() for l in range(12)]
Xs = [L.centred(D["emb"], l, D["obs_idx"])[U] for l in range(12)]
names = [(l, i, j) for l in range(12) for (i, j) in bands]

def band_scores(rows, l, i, j, pos_local):
    Ub = USs[l][rows, i:j + 1]
    return Ub[:, 0] if i == j else S.tf_scores_loo(Ub, pos_local)

out = []
for s in range(NSPLIT):
    rng = np.random.default_rng([SEED, s]); p = rng.permutation(n); A = np.sort(p[: n // 2]); B = np.sort(p[n // 2:])
    yA, yB = y[A], y[B]; dA, dB = det[A], det[B]; pA, pB = np.flatnonzero(yA), np.flatnonzero(yB)
    fsA = L.FastStrata(dA, yA)
    SC = np.column_stack([band_scores(A, l, i, j, pA) for (l, i, j) in names])
    a_raw = fsA.cols(SC, yA[:, None].astype(float))
    width1 = np.array([i == j for (_, i, j) in names])
    a_sel = np.where(width1, 0.5 + np.abs(a_raw - 0.5), a_raw)
    c = int(np.argmax(a_sel)); l, i, j = names[c]
    sB = band_scores(B, l, i, j, pB)
    if i == j and a_raw[c] < 0.5: sB = -sB
    heldout = S.stratified_auroc(sB, yB, dB)
    fA = [S.stratified_auroc(S.tf_scores_loo(Xs[ll][A], pA), yA, dA) for ll in range(12)]
    lf = int(np.argmax(fA)); fB = S.stratified_auroc(S.tf_scores_loo(Xs[lf][B], pB), yB, dB)
    out.append(dict(context=ctx, task=task, split=s, nA=A.size, nB=B.size, tfA=int(yA.sum()), tfB=int(yB.sum()),
                    chosen_layer=l, chosen_band=f"SV{i+1}-{j+1}" if j > i else f"SV{i+1}", insample_A=float(a_sel[c]),
                    heldout_B=float(heldout), full_chosen_layer=lf, full_insample_A=float(fA[lf]), full_heldout_B=float(fB)))
    if s % 20 == 0: print(f"[{time.time()-t0:6.1f}s] split {s}", out[-1], flush=True)
pd.DataFrame(out).to_csv(L.HERE / "partial" / f"crossfit_{ctx}_{task}.csv", index=False)
print("DONE", time.time() - t0, flush=True)
