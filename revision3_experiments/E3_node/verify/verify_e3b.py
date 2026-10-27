"""Independent verification, part 2 (no stats3 / e3_lib imports):
 (a) ledger 0.744: label set from the raw edge dataset; screen procedure (SV2-7, scaler + logistic C=1, 5-fold
     StratifiedKFold seed 42, fold-mean AUROC, mean over layers); LOO cosine SV2-7 stratified / unstratified
     layer means; full 512 stratified layer mean;
 (b) cross-fit split 0 and 1 in lung (same split RNG as e3_crossfit.py): chosen cell, in-sample A, held-out B;
 (c) own N=300 within-stratum null for the lung SV3-12 @ L2 cell -> z (reported 'cos' z from bandscan_lung.csv).
Writes verify/verify_e3b.csv.
"""
import csv, json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import rankdata
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent; E3 = HERE.parent; R3 = E3.parent
EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"
sys.path.insert(0, str(HERE)); from verify_common import mw, strat, deciles, loo, ctx_data  # own implementations
sys.path.insert(0, str(R3)); import contexts as X

rows = []
def add(q, mine, rep): rows.append(dict(quantity=q, mine=float(mine), reported=float(rep))); print(q, mine, rep, flush=True)

def svd_fixed(Xc):
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False); US = U * S
    sg = np.sign(US[np.abs(US).argmax(0), np.arange(US.shape[1])]); sg[sg == 0] = 1
    return US * sg

# ---------------- (a) ledger
emb, obs, genes, cnt, ytf, det_all = ctx_data("immune")
er = list(csv.DictReader(open(X.S38 / "cycle4_immune_main/cycle1_edge_dataset.tsv"), delimiter="\t"))
named = {}
for r in er: named.setdefault(r["source"], int(r["source_idx"])); named.setdefault(r["target"], int(r["target_idx"]))
src = {r["source"] for r in er}; rp = {int(o): i for i, o in enumerate(obs)}
og = sorted(g for g in named if named[g] in rp)
tf = [g for g in og if g in src]; tg = [g for g in og if g not in src]
lab = np.array([rp[named[g]] for g in tf + tg]); y = np.r_[np.ones(len(tf)), np.zeros(len(tg))].astype(bool); det = det_all[lab]
L = json.loads((E3 / "ledger_0744.json").read_text())
add("ledger n_tf", len(tf), L["n_tf"]); add("ledger n_target_only", len(tg), L["n_target_only"])
cvm, loo_u, loo_s, f_s = [], [], [], []
for l in range(12):
    Xc = np.asarray(emb[l], float)[obs]; Xc = Xc - Xc.mean(0)
    F = svd_fixed(Xc)[lab, 1:7]
    skf = StratifiedKFold(5, shuffle=True, random_state=42); fold = []
    for a, b in skf.split(F, y.astype(int)):
        m = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=1000)).fit(F[a], y[a].astype(int))
        fold.append(roc_auc_score(y[b].astype(int), m.predict_proba(F[b])[:, 1]))
    cvm.append(np.mean(fold))
    s = loo(F, y); loo_u.append(mw(s, y)); loo_s.append(strat(s, y, det))
    f_s.append(strat(loo(Xc[lab], y), y, det))
sp = L["screen_procedure"]
add("ledger screen procedure mean over layers", np.mean(cvm), sp["mean_over_layers"])
add("ledger LOO SV2-7 unstrat mean", np.mean(loo_u), L["loo_cos_SV2_7"]["unstratified"]["mean_over_layers"])
add("ledger LOO SV2-7 strat mean", np.mean(loo_s), L["loo_cos_SV2_7"]["stratified"]["mean_over_layers"])
add("ledger LOO full512 strat mean", np.mean(f_s), L["loo_cos_full512"]["stratified"]["mean_over_layers"])

# ---------------- (b) cross-fit lung splits 0, 1
emb, obs, genes, cnt, y, det = ctx_data("lung")
Xs = []; USs = []
for l in range(12):
    Xc = np.asarray(emb[l], float)[obs]; Xc = Xc - Xc.mean(0); Xs.append(Xc); USs.append(svd_fixed(Xc)[:, :12])
bands = [(i, j) for i in range(12) for j in range(i, 12)]; names = [(l, i, j) for l in range(12) for (i, j) in bands]
cf = pd.read_csv(E3 / "partial/crossfit_lung_lambert.csv")
SEED = 20261026 + 1000 * 0 + 60; n = y.size
for s_ in (0, 1):
    rng = np.random.default_rng([SEED, s_]); p = rng.permutation(n); A = np.sort(p[: n // 2]); B = np.sort(p[n // 2:])
    best = (-1, None, None)
    for (l, i, j) in names:
        Ub = USs[l][A, i:j + 1]; sc = Ub[:, 0] if i == j else loo(Ub, y[A]); a = strat(sc, y[A], det[A])
        sel = 0.5 + abs(a - 0.5) if i == j else a
        if sel > best[0]: best = (sel, (l, i, j), a)
    sel, (l, i, j), araw = best
    Ub = USs[l][B, i:j + 1]; sB = Ub[:, 0] if i == j else loo(Ub, y[B])
    if i == j and araw < 0.5: sB = -sB
    r = cf[cf.split == s_].iloc[0]
    band = f"SV{i+1}-{j+1}" if j > i else f"SV{i+1}"
    print("split", s_, "mine", band, l, "reported", r.chosen_band, r.chosen_layer)
    add(f"crossfit lung split{s_} chosen cell match", float(band == r.chosen_band and l == r.chosen_layer), 1.0)
    add(f"crossfit lung split{s_} insample A", sel, r.insample_A)
    add(f"crossfit lung split{s_} heldout B", strat(sB, y[B], det[B]), r.heldout_B)

# ---------------- (c) own null for lung SV3-12 @ L2
Ub = USs[2][:, 2:12]; En = Ub / np.linalg.norm(Ub, axis=1)[:, None]; obs_a = strat(loo(Ub, y), y, det)
rng = np.random.default_rng(5); nul = []
for b in range(300):
    yp = y.copy()
    for st in np.unique(det):
        idx = np.flatnonzero(det == st); yp[idx] = yp[rng.permutation(idx)]
    v = En[yp].sum(0); s = En @ v; k = yp.sum(); s = np.where(yp, (s - 1) / (k - 1), s / k); nul.append(strat(s, yp, det))
nul = np.array(nul)
bs = pd.read_csv(E3 / "bandscan_lung.csv"); bs = bs[(bs.row_type == "band_cell") & (bs.layer == 2) & (bs.band == "SV3-12")].iloc[0]
add("band lung SV3-12 L2 null mean", nul.mean(), bs.null_mean_cos); add("band lung SV3-12 L2 null sd", nul.std(ddof=1), bs.null_sd_cos)
add("band lung SV3-12 L2 z", (obs_a - nul.mean()) / nul.std(ddof=1), bs.z_cos)
df = pd.DataFrame(rows); df["absdiff"] = (df.mine - df.reported).abs(); df.to_csv(HERE / "verify_e3b.csv", index=False)
print(df.to_string()); print("DONE")
