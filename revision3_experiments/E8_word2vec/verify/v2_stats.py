"""Verifier step 2: independent recomputation of headline statistics (no stats3, no e8_eval code).

Own implementations: centring, cosine, sklearn roc_auc_score, per-stratum Mann-Whitney, own LOO, own SVD and
assortativity (np.corrcoef over both edge orientations), own vertex subsampling (delete-d scaling sqrt(m/(n-m))),
own within-TF pooled AUROC, own two-way cluster bootstrap, own within-decile permutation, networkx double-edge
swap MS null. Inputs: cache/<ctx>/ (genes, labels, edges, pairs), models/*.npy, local scGPT embeddings,
eval/<ctx>_rewired_nonpairs.npy (the stored non-pair set; checked for overlap with edges).
Writes verify/v2_stats.csv (mine vs stored) and verify/v2_extra.json.
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd, networkx as nx
from sklearn.metrics import roc_auc_score
HERE = Path(__file__).resolve().parent; E8 = HERE.parent; R3 = E8.parent
EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"
W = pd.read_csv(E8 / "w2v_results.csv"); P = pd.read_csv(E8 / "paired_vs_scgpt.csv"); H1S = pd.read_csv(E8 / "h1_assortativity_scgpt.csv")
rows = []; extra = {}

def rec(name, mine, stored): rows.append(dict(check=name, mine=float(mine), stored=float(stored), absdiff=abs(float(mine) - float(stored)))); print(name, round(float(mine), 5), round(float(stored), 5), flush=True)

def load(ctx):
    d = R3 / "cache" / ctx
    obs = np.load(d / "obs_idx.npy"); nl = np.load(d / "node_labels.npz"); e = np.load(d / "edges_string.npy"); h3 = np.load(d / "h3_pairs.npz")
    return obs, nl, e, h3

def cen(M): M = np.asarray(M, np.float64); return M - M.mean(0)
def unit(M): return M / np.linalg.norm(M, axis=1, keepdims=True)
def w2v(ctx, s, corpus="full", w=5, e=10): return cen(np.load(E8 / "models" / f"{ctx}__{corpus}__w{w}__e{e}__s{s}.npy"))
def scg(ctx, l, obs): return cen(np.load(EMB / f"{ctx}.npy", mmap_mode="r")[l][obs])
def cos(U, a, b): return (U[a] * U[b]).sum(1)

def h2_stat(E, is_tf, dec):
    U = unit(E); tot = U @ U[is_tf].sum(0); k = is_tf.sum()
    s = np.where(is_tf, (tot - 1.0) / (k - 1), tot / k)
    num = den = 0.0
    for st in np.unique(dec):
        m = dec == st; y = is_tf[m]
        if 0 < y.sum() < m.sum():
            n1, n0 = y.sum(), (~y).sum(); num += roc_auc_score(y, s[m]) * n1 * n0; den += n1 * n0
    return num / den, s

def h3_stat(sc, tf, y, w=None):
    w = np.ones(sc.size) if w is None else w; num = den = 0.0
    for t in np.unique(tf):
        m = tf == t; p, q = sc[m & y], sc[m & ~y]; wp, wq = w[m & y], w[m & ~y]
        if wp.sum() == 0 or wq.sum() == 0: continue
        G = (p[:, None] > q[None, :]) + 0.5 * (p[:, None] == q[None, :])
        num += (wp[:, None] * wq[None, :] * G).sum(); den += wp.sum() * wq.sum()
    return num / den

def stored_w(ctx, corpus, task, seed, sv=None, w=5, e=10):
    q = W[(W.context == ctx) & (W.task == task) & (W.seed.astype(str) == str(seed))]
    q = q[q.corpus == corpus] if seed == "mean" else q[(q.corpus == corpus) & (q.window == w) & (q.epochs == e)]
    if sv is not None: q = q[q.sv == sv]
    assert len(q) == 1, (ctx, corpus, task, seed, len(q)); return q.iloc[0]

def stored_p(ctx, task, layer, g="full_w5_e10"):
    q = P[(P.context == ctx) & (P.task == task) & (P.layer == layer) & (P.w2v == g)]; assert len(q) == 1; return q.iloc[0]

# ================================================================ H1b (lung, immune): AUROC + delete-d CI of delta
for ctx, L in (("lung", 1), ("immune", 2), ("immune", 11)):
    obs, nl, e, h3 = load(ctx); n = obs.size
    nonp = np.load(E8 / "eval" / f"{ctx}_rewired_nonpairs.npy")
    ek = {tuple(sorted(x)) for x in e.tolist()}
    extra[f"{ctx}_nonpairs_overlap_edges"] = int(sum(tuple(sorted(x)) in ek for x in nonp.tolist()))
    extra[f"{ctx}_nonpairs_selfloops"] = int((nonp[:, 0] == nonp[:, 1]).sum())
    Us = [unit(w2v(ctx, s)) for s in range(5)]; Uc = unit(scg(ctx, L, obs))
    y = np.r_[np.ones(len(e)), np.zeros(len(nonp))]
    def aucs(mask_nodes=None):
        if mask_nodes is None: ep = np.ones(len(e), bool); en = np.ones(len(nonp), bool)
        else: ep = mask_nodes[e[:, 0]] & mask_nodes[e[:, 1]]; en = mask_nodes[nonp[:, 0]] & mask_nodes[nonp[:, 1]]
        yy = np.r_[np.ones(ep.sum()), np.zeros(en.sum())]
        f = lambda U: roc_auc_score(yy, np.r_[cos(U, e[ep, 0], e[ep, 1]), cos(U, nonp[en, 0], nonp[en, 1])])
        return np.mean([f(U) for U in Us]), f(Uc)
    a_w, a_c = aucs()
    if L in (1, 2):
        rec(f"{ctx} H1b W2V full seed-mean AUROC", a_w, stored_w(ctx, "full_w5_e10", "H1b_pair_cos_auroc", "mean").value)
    rec(f"{ctx} H1b scGPT L{L} AUROC", a_c, stored_p(ctx, "H1b_pair_cos_auroc", L).scgpt)
    rec(f"{ctx} H1b delta W2V-scGPT L{L}", a_w - a_c, stored_p(ctx, "H1b_pair_cos_auroc", L).delta)
    rng = np.random.default_rng(4242); m = n // 2; dv = []
    for r in range(400):
        mk = np.zeros(n, bool); mk[rng.choice(n, m, replace=False)] = True; x, z = aucs(mk); dv.append(x - z)
    dv = np.array(dv); dev = (dv - dv.mean()) * np.sqrt(m / (n - m)); d0 = a_w - a_c
    lo, hi = d0 - np.quantile(dev, 0.975), d0 - np.quantile(dev, 0.025)
    sp = stored_p(ctx, "H1b_pair_cos_auroc", L)
    rec(f"{ctx} H1b delta L{L} delete-d CI lo (own draws B=400)", lo, sp.delta_ci_lo_deleted)
    rec(f"{ctx} H1b delta L{L} delete-d CI hi (own draws B=400)", hi, sp.delta_ci_hi_deleted)

# ================================================================ H2 (lung W2V, immune scGPT L2, kidney_A W2V) + permutation p
for ctx in ("lung", "immune", "kidney_A"):
    obs, nl, e, h3 = load(ctx); is_tf = nl["is_tf_lambert"].astype(bool); dec = nl["det_decile"]
    vals = [h2_stat(w2v(ctx, s), is_tf, dec)[0] for s in range(5)]
    rec(f"{ctx} H2 W2V full seed-mean", np.mean(vals), stored_w(ctx, "full_w5_e10", "H2_strat_tf_auroc", "mean").value)
    best = int(P[(P.context == ctx) & (P.task == "H2_strat_tf_auroc") & (P.w2v == "full_w5_e10")].sort_values("scgpt").layer.iloc[-1])
    sc_c = h2_stat(scg(ctx, best, obs), is_tf, dec)[0]
    rec(f"{ctx} H2 scGPT best layer L{best}", sc_c, stored_p(ctx, "H2_strat_tf_auroc", best).scgpt)
    # own within-decile permutation for W2V seed 0 and the scGPT best layer (N=2000, different seed)
    rng = np.random.default_rng(99); E0 = w2v(ctx, 0); Ec = scg(ctx, best, obs); nw, ncs = [], []
    for b in range(2000 if ctx != "lung" else 600):
        yp = is_tf.copy()
        for st in np.unique(dec):
            ix = np.flatnonzero(dec == st); yp[ix] = rng.permutation(yp[ix])
        nw.append(h2_stat(E0, yp, dec)[0]); ncs.append(h2_stat(Ec, yp, dec)[0])
    nw, ncs = np.array(nw), np.array(ncs); v0 = h2_stat(E0, is_tf, dec)[0]
    st0 = stored_w(ctx, "full", "H2_strat_tf_auroc", 0)
    rec(f"{ctx} H2 W2V s0 perm p (own N={nw.size})", (np.sum(nw >= v0) + 1) / (nw.size + 1), st0.p_raw)
    rec(f"{ctx} H2 W2V s0 null mean (own)", nw.mean(), st0.null_mean)
    rec(f"{ctx} H2 scGPT L{best} null mean (own)", ncs.mean(), stored_p(ctx, "H2_strat_tf_auroc", best).scgpt_null_mean)

# ================================================================ H3 (kidney_A, immune): stat + two-way bootstrap CI of delta
for ctx, L in (("kidney_A", 1), ("immune", 0)):
    obs, nl, e, h3 = load(ctx); tf, tg, y = h3["tf"], h3["tg"], h3["y"].astype(bool)
    cw = [cos(unit(w2v(ctx, s)), tf, tg) for s in range(5)]; cc = cos(unit(scg(ctx, L, obs)), tf, tg)
    vw = np.mean([h3_stat(c, tf, y) for c in cw]); vc = h3_stat(cc, tf, y)
    rec(f"{ctx} H3 W2V full seed-mean", vw, stored_w(ctx, "full_w5_e10", "H3_withinTF_auroc", "mean").value)
    rec(f"{ctx} H3 scGPT L{L}", vc, stored_p(ctx, "H3_withinTF_auroc", L).scgpt)
    utf, itf = np.unique(tf, return_inverse=True); utg, itg = np.unique(tg, return_inverse=True)
    rng = np.random.default_rng(5150); dd = []
    for r in range(300):
        w = np.bincount(rng.integers(0, utf.size, utf.size), minlength=utf.size)[itf] * \
            np.bincount(rng.integers(0, utg.size, utg.size), minlength=utg.size)[itg]
        w = w.astype(float); dd.append(np.mean([h3_stat(c, tf, y, w) for c in cw]) - h3_stat(cc, tf, y, w))
    sp = stored_p(ctx, "H3_withinTF_auroc", L)
    rec(f"{ctx} H3 delta L{L} CI lo (own two-way B=300)", np.nanquantile(dd, 0.025), sp.delta_ci_lo)
    rec(f"{ctx} H3 delta L{L} CI hi (own two-way B=300)", np.nanquantile(dd, 0.975), sp.delta_ci_hi)

# ================================================================ H1a: own SVD + assortativity; own MS null (immune W2V s0)
def assort(x, a, b): return np.corrcoef(np.r_[x[a], x[b]], np.r_[x[b], x[a]])[0, 1]
for ctx in ("lung", "immune"):
    obs, nl, e, h3 = load(ctx); E0 = w2v(ctx, 0)
    U, S, Vt = np.linalg.svd(E0, full_matrices=False); X10 = U[:, :10] * S[:10]
    for j in range(10):
        rec(f"{ctx} H1a W2V s0 SV{j+1} r", assort(X10[:, j], e[:, 0], e[:, 1]), stored_w(ctx, "full", "H1a_assortativity", 0, sv=j + 1).value)
    Uc, Sc, _ = np.linalg.svd(scg(ctx, 2, obs), full_matrices=False)
    rec(f"{ctx} H1a scGPT L2 SV2 r", assort(Uc[:, 1] * Sc[1], e[:, 0], e[:, 1]), H1S[(H1S.context == ctx) & (H1S.layer == 2) & (H1S.sv == 2)].r.iloc[0])
obs, nl, e, h3 = load("immune"); E0 = w2v("immune", 0)
U, S, Vt = np.linalg.svd(E0, full_matrices=False); X10 = U[:, :10] * S[:10]
G = nx.Graph(); G.add_edges_from(e.tolist()); mE = G.number_of_edges()
nx.double_edge_swap(G, nswap=10 * mE, max_tries=10 ** 8, seed=11)
nullr = []
for t in range(300):
    nx.double_edge_swap(G, nswap=mE, max_tries=10 ** 8, seed=100 + t)
    ed = np.array(G.edges()); nullr.append([assort(X10[:, j], ed[:, 0], ed[:, 1]) for j in range(10)])
nullr = np.array(nullr); js = json.loads((E8 / "eval" / "immune.json").read_text())["null_summaries"]["h1a"]["full_w5_e10_s0"]
for j in (0, 1, 4, 9):
    rec(f"immune H1a W2V s0 SV{j+1} MS null mean (own nx, 300)", nullr[:, j].mean(), js[f"SV{j+1}"]["mean"])
    rec(f"immune H1a W2V s0 SV{j+1} MS null sd (own nx, 300)", nullr[:, j].std(ddof=1), js[f"SV{j+1}"]["sd"])
z = (np.array([assort(X10[:, j], e[:, 0], e[:, 1]) for j in range(10)]) - nullr.mean(0)) / nullr.std(0, ddof=1)
extra["immune_W2V_s0_H1a_own_z_SV1_10"] = z.round(2).tolist()
extra["immune_W2V_s0_H1a_own_max_null_z_quantile95"] = float(np.quantile(((nullr - nullr.mean(0)) / nullr.std(0, ddof=1)).max(1), 0.95))
# CKA spot check
def cka(A, B):
    A = A - A.mean(0); B = B - B.mean(0); return np.linalg.norm(A.T @ B) ** 2 / (np.linalg.norm(A.T @ A) * np.linalg.norm(B.T @ B))
ck = pd.read_csv(E8 / "cka.csv"); obs, *_ = load("lung")
rec("lung CKA w512 vs w5_e10 (calibration)", cka(w2v("lung", 0, "extr", 512, 10), w2v("lung", 0, "extr", 5, 10)),
    ck[(ck.context == "lung") & (ck.kind == "calibration") & (ck.b == "extr_w5_e10_s0")].cka.iloc[0])
pd.DataFrame(rows).to_csv(HERE / "v2_stats.csv", index=False)
(HERE / "v2_extra.json").write_text(json.dumps(extra, indent=1))
print(json.dumps(extra)); print("DONE")
