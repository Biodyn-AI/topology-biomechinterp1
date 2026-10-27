"""E4c independent verifier, part 1 (no model fitting).

1. Gene alignment: obs_genes == h5ad var index at obs_idx; observed rows non-zero, never-observed rows zero;
   h3 pair labels agree with TRRUST (y=1 <=> pair in TRRUST; y=0 pairs absent from any_regulatory both ways).
2. CV folds: C3 gene-disjoint (train genes and test genes disjoint) for all 5 repeats; C2 TF-disjoint.
3. Rewired labels: every TF's and every target's positive count preserved, all 100 draws.
4. Recompute the unsupervised reference (signed cosine, full 512-d, V-centred) within-TF AUROC on every
   C3 fold with an independent within-TF Mann-Whitney implementation (scipy rankdata), compare with probe_folds.csv.
5. Recompute Nadeau-Bengio intervals (lung scGPT bilinear vs target-only, C3) from probe_folds.csv.
6. Recompute empirical p and max-T from the raw per-draw checkpoint AUROCs (lung bilinear, 12 layers).
7. Added test: difference-in-differences under rewiring. Statistic D = mean_fold AUROC(bilinear) - AUROC(concat)
   on true labels; null = the same difference on each of the 100 rewired draws (same draws, same folds). concat
   is additive (f(x_i)+g(x_j)) and therefore carries no within-TF relational information, so D isolates the
   interaction term. Also reports concat true vs rewired (is the rewired null 'relation-clean'?).
Outputs: verify/v1_*.csv, verify/v1_summary.json
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd, h5py
from scipy.stats import rankdata, t as tdist

HERE = Path(__file__).resolve().parent; E4C = HERE.parent; R3 = E4C.parent
sys.path.insert(0, str(E4C)); sys.path.insert(0, str(R3))
import e4c_lib as L     # used ONLY for make_folds (fold identity must match) and Ctx paths
import contexts as X
CTX = ["lung", "lung_external", "kidney_B"]
out = {}


def within_tf_auc(tf, y, s, min_pos=3, min_neg=3):
    num = den = 0.0
    for t in np.unique(tf):
        m = tf == t; yy = y[m]; n1 = yy.sum(); n0 = m.sum() - n1
        if n1 < min_pos or n0 < min_neg: continue
        r = rankdata(s[m]); num += r[yy].sum() - n1 * (n1 + 1) / 2; den += n1 * n0
    return num / den if den else np.nan


# ---------------------------------------------------------------- 1. alignment
align = []
t, trrust = X.load_trrust(); trs = set(trrust)
for ctx in CTX:
    d = R3 / "cache" / ctx
    genes = json.loads((d / "obs_genes.json").read_text()); oi = np.load(d / "obs_idx.npy")
    hg = X.h5ad_genes(X.CONTEXTS[ctx]["h5ad"])
    emb = np.load(Path.home() / f"biomechinterp-local/subproject41_r3/emb/{ctx}.npy", mmap_mode="r")
    nrm = np.linalg.norm(np.asarray(emb[0], np.float64), axis=1)
    z = np.load(d / "h3_pairs.npz"); tf, tg, y = z["tf"], z["tg"], z["y"].astype(bool)
    anyr = {(int(a), int(b)) for a, b in np.load(d / "any_regulatory.npy")}
    lab_ok = all(((genes[a], genes[b]) in trs) == bool(yy) for a, b, yy in zip(tf, tg, y))
    neg_ok = all((int(a), int(b)) not in anyr and (int(b), int(a)) not in anyr for a, b, yy in zip(tf, tg, y) if not yy)
    align.append(dict(context=ctx, genes_match_h5ad=[hg[i] for i in oi] == genes, n_V=len(genes),
                      observed_rows_nonzero=bool((nrm[oi] > 0).all()),
                      unobserved_rows_zero=bool((np.delete(nrm, oi) == 0).all()),
                      emb_rows=int(emb.shape[1]), h5ad_genes=len(hg), labels_match_trrust=lab_ok,
                      negatives_absent_from_any_regulatory_both_directions=neg_ok, n_pairs=int(y.size), n_pos=int(y.sum())))
    print(align[-1], flush=True)
pd.DataFrame(align).to_csv(HERE / "v1_alignment.csv", index=False)

# ---------------------------------------------------------------- 2-3. folds and rewiring
fchk, rchk = [], []
for ctx in CTX:
    c = L.Ctx(ctx)
    for reg in ("C2", "C3"):
        for rep in range(5):
            for k, (tr, te) in enumerate(L.make_folds(c, reg, rep)):
                gtr = set(c.tf[tr]) | set(c.tg[tr]); gte = set(c.tf[te]) | set(c.tg[te])
                fchk.append(dict(context=ctx, regime=reg, repeat=rep, fold=k, n_train=tr.size, n_test=te.size,
                                 gene_overlap=len(gtr & gte) if reg == "C3" else np.nan,
                                 tf_overlap=len(set(c.tf[tr]) & set(c.tf[te])),
                                 pair_overlap=len(set(tr) & set(te))))
    zr = np.load(E4C / "rewired" / f"{ctx}.npz"); Yr = np.unpackbits(zr["Y"], axis=1, count=int(zr["n"])).astype(bool)
    dtf = np.bincount(c.tf[c.y], minlength=c.G); dtg = np.bincount(c.tg[c.y], minlength=c.G)
    okt = all(np.array_equal(np.bincount(c.tf[v], minlength=c.G), dtf) for v in Yr)
    okg = all(np.array_equal(np.bincount(c.tg[v], minlength=c.G), dtg) for v in Yr)
    rchk.append(dict(context=ctx, n_draws=len(Yr), tf_counts_preserved=okt, target_counts_preserved=okg,
                     retained_mean=float((Yr & c.y).sum(1).mean() / c.y.sum())))
F = pd.DataFrame(fchk); F.to_csv(HERE / "v1_folds.csv", index=False)
out["c3_max_gene_overlap"] = int(F[F.regime == "C3"].gene_overlap.max())
out["c2_max_tf_overlap"] = int(F[F.regime == "C2"].tf_overlap.max())
out["max_pair_overlap"] = int(F.pair_overlap.max())
out["rewiring"] = rchk
print(out, flush=True)

# ---------------------------------------------------------------- 4. reference cosine recompute (C3 all repeats)
PF = pd.read_csv(E4C / "probe_folds.csv")
rc = []
for ctx in CTX:
    c = L.Ctx(ctx)
    emb = np.load(Path.home() / f"biomechinterp-local/subproject41_r3/emb/{ctx}.npy", mmap_mode="r")
    oi = np.load(R3 / "cache" / ctx / "obs_idx.npy")
    for l in (0, 2, 6, 11):
        Xl = np.asarray(emb[l], np.float64)[oi]; Xl = Xl - Xl.mean(0)
        Xn = Xl / np.linalg.norm(Xl, axis=1, keepdims=True)
        s = (Xn[c.tf] * Xn[c.tg]).sum(1)
        for rep in range(5):
            for k, (tr, te) in enumerate(L.make_folds(c, "C3", rep)):
                mine = within_tf_auc(c.tf[te], c.y[te], s[te])
                theirs = PF.query("context==@ctx and model=='ref' and feature=='cosine' and layer==@l and regime=='C3' and repeat==@rep and fold==@k").auc.values[0]
                rc.append(dict(context=ctx, layer=l, repeat=rep, fold=k, mine=mine, theirs=theirs))
RC = pd.DataFrame(rc); RC["absdiff"] = (RC.mine - RC.theirs).abs(); RC.to_csv(HERE / "v1_ref_cosine.csv", index=False)
out["ref_cosine_max_absdiff"] = float(RC.absdiff.max())
out["ref_cosine_layer_means"] = RC.groupby(["context", "layer"]).mine.mean().round(4).reset_index().to_dict("records")
print("ref cosine max abs diff", RC.absdiff.max(), flush=True)


# ---------------------------------------------------------------- 5. NB recompute
def nb(d, ratio):
    d = np.asarray(d, float); J = d.size; m = d.mean(); se = np.sqrt((1 / J + ratio) * d.var(ddof=1))
    q = tdist.ppf(0.975, J - 1); return m, m - q * se, m + q * se


nbr = []
for ctx in CTX:
    for l in range(12):
        a = PF.query("context==@ctx and model=='M1' and source=='scgpt' and feature=='bilinear' and layer==@l and regime=='C3'").sort_values(["repeat", "fold"])
        b = PF.query("context==@ctx and model=='M1' and source=='scgpt' and feature=='tg_only' and layer==@l and regime=='C3'").sort_values(["repeat", "fold"])
        assert (a[["repeat", "fold"]].values == b[["repeat", "fold"]].values).all()
        ratio = a.n_test.mean() / a.n_train.mean()
        m, lo, hi = nb(a.auc.values - np.fmax(b.auc.values, 0.5), ratio)
        pm, plo, phi = nb(a.auc.values - 0.5, ratio)
        nbr.append(dict(context=ctx, layer=l, bil_auc=a.auc.mean(), bil_lo=0.5 + plo, bil_hi=0.5 + phi,
                        tg_auc=b.auc.mean(), delta=m, lo=lo, hi=hi))
NBR = pd.DataFrame(nbr); NBR.to_csv(HERE / "v1_nb_bilinear_vs_tg.csv", index=False)
print(NBR.round(4).to_string(), flush=True)

# ---------------------------------------------------------------- 6-7. empirical p, max-T, DiD
def null_matrix(ctx, src, feat, l):
    d = pd.read_csv(E4C / "partial" / "m1" / ctx / f"{src}_{feat}_L{l}_C3_r0_null.csv")
    piv = d.pivot_table(index="fold", columns="label_col", values="auc")
    return piv.mean(0).values      # (101,) col 0 = true


did = []
for ctx in CTX:
    Tb, Tc, Tt = [], [], []
    for l in range(12):
        b = null_matrix(ctx, "scgpt", "bilinear", l); cc = null_matrix(ctx, "scgpt", "concat", l)
        Tb.append(b); Tc.append(cc)
    Tb = np.array(Tb).T; Tc = np.array(Tc).T                          # (101, 12)
    D = Tb - Tc
    for name, M in (("bilinear_vs_rewired", Tb), ("concat_vs_rewired", Tc), ("DiD_bilinear_minus_concat", D)):
        obs, nul = M[0], M[1:]
        p = (np.sum(nul >= obs[None, :], 0) + 1) / (nul.shape[0] + 1)
        mu = nul.mean(0); sd = nul.std(0, ddof=1); z = (obs - mu) / sd; zn = (nul - mu) / sd; mx = zn.max(1)
        padj = np.array([(np.sum(mx >= zz) + 1) / (nul.shape[0] + 1) for zz in z])
        for l in range(12):
            did.append(dict(context=ctx, test=name, layer=l, obs=obs[l], null_mean=mu[l], null_sd=sd[l],
                            null_q975=np.quantile(nul[:, l], .975), z=z[l], p_emp=p[l], p_maxT=padj[l], N=nul.shape[0]))
DID = pd.DataFrame(did); DID.to_csv(HERE / "v1_rewired_did.csv", index=False)
NR = pd.read_csv(E4C / "null_results.csv")
chk = DID[DID.test == "bilinear_vs_rewired"].merge(
    NR[(NR.model == "M1") & (NR.source == "scgpt") & (NR.feature == "bilinear")][["context", "layer", "obs_auc_mean", "p_emp", "p_maxT_layers"]],
    on=["context", "layer"])
out["null_recompute_max_absdiff_obs"] = float((chk.obs - chk.obs_auc_mean).abs().max())
out["null_recompute_max_absdiff_p"] = float((chk.p_emp_x - chk.p_emp_y).abs().max())
out["null_recompute_max_absdiff_pmaxT"] = float((chk.p_maxT - chk.p_maxT_layers).abs().max())
summ = DID.groupby(["context", "test"]).agg(obs_min=("obs", "min"), obs_max=("obs", "max"),
                                             null_mean_min=("null_mean", "min"), null_mean_max=("null_mean", "max"),
                                             n_p_lt05=("p_emp", lambda x: int((x < 0.05).sum())),
                                             n_maxT_lt05=("p_maxT", lambda x: int((x < 0.05).sum()))).round(4)
print(summ.to_string(), flush=True)
out["did_summary"] = summ.reset_index().to_dict("records")
json.dump(out, open(HERE / "v1_summary.json", "w"), indent=1, default=float)
print("done")
