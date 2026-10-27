"""E6 selection-bias test: is the legacy 186-gene set an unusual sample of Tier A for the TRRUST TF/target task?

Null: 10,000 random gene subsets of Tier A, each matched to the legacy set cell-by-cell on
(global TRRUST label: TF / target-only) x (decile of lung mean expression; 'NA' for genes outside lung V).
Lung mean expression = column mean of cache/lung/expr.npy (3,000 lung cells, log-normalised). Decile edges are
those of the Tier A TRRUST-labelled genes inside lung V; legacy genes are placed with the same edges.
Statistic per subset and static model = estimator (i): rows centred over the subset, full dimension,
LOO mean cosine to the subset's TFs (stats3.tf_scores_loo semantics), AUROC TF vs target-only (stats3.auroc).
Computed through each model's uncentred Gram matrix on the pool (double-centring inside the subset); checked
against a direct computation on the legacy set.
Variants: 'legacy186' (size 186; tGPT lacks 8 HLA genes, so its legacy value uses the 178 it has),
'legacy178' (legacy & Tier A, 178 genes, every model on identical genes), and the added sensitivity
'legacy186_degmatched': cells = label x quartile of genome-wide TRRUST degree (out-degree for TFs, in-degree for
target-only genes; quartiles within label over the pool) x quintile of lung mean expression (+ 'NA').
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import e6lib as L
S = L.S
N_SUB = 10000


def main():
    rf = L.refs()
    A = L.load_geneset("tierA")["genes"]; Aset = set(A)
    leg = json.loads((L.HERE / "legacy_genes.json").read_text())["legacy186"]
    lab = lambda g: "TF" if g in rf["trrust_tf"] else ("TG" if g in rf["trrust_target_only"] else None)
    pool = [g for g in A if lab(g) is not None]
    lg, _ = L.ctx_obs("lung"); ex = np.load(L.CACHE / "lung/expr.npy", mmap_mode="r")
    mexp = dict(zip(lg, np.asarray(ex, dtype=np.float64).mean(0)))
    edges_src = np.array([mexp[g] for g in pool if g in mexp])
    q = np.unique(np.quantile(edges_src, np.linspace(0, 1, 11)))
    def dec(g):
        if g not in mexp: return "NA"
        return int(np.clip(np.searchsorted(q, mexp[g], side="right") - 1, 0, len(q) - 2))
    cell = {g: (lab(g), dec(g)) for g in set(pool) | set(leg)}
    outd, ind = {}, {}
    for a_, b_ in rf["trrust"]: outd[a_] = outd.get(a_, 0) + 1; ind[b_] = ind.get(b_, 0) + 1
    degv = lambda g: outd.get(g, 0) if lab(g) == "TF" else ind.get(g, 0)
    qdeg = {L_: np.unique(np.quantile([degv(g) for g in pool if lab(g) == L_], [0, .25, .5, .75, 1])) for L_ in ("TF", "TG")}
    q5 = np.unique(np.quantile(edges_src, np.linspace(0, 1, 6)))
    def cell2(g):
        dq = int(np.clip(np.searchsorted(qdeg[lab(g)], degv(g), side="right") - 1, 0, len(qdeg[lab(g)]) - 2))
        eq = "NA" if g not in mexp else int(np.clip(np.searchsorted(q5, mexp[g], side="right") - 1, 0, len(q5) - 2))
        return (lab(g), dq, eq)
    cellD = {g: cell2(g) for g in set(pool) | set(leg)}
    universe = sorted(set(pool) | set(leg)); upos = {g: i for i, g in enumerate(universe)}
    rows_meta = []
    variants = {"legacy186": sorted(leg), "legacy178": sorted(set(leg) & Aset), "legacy186_degmatched": sorted(leg)}
    # matched draws (shared by every model)
    rng = np.random.default_rng(L.SEED + 606)
    draws = {}
    comp = {}
    for v, gl in variants.items():
        cf = cellD if v.endswith("degmatched") else cell
        pool_by_cell = {}
        for g in pool: pool_by_cell.setdefault(cf[g], []).append(g)
        need = {}
        for g in gl: need[cf[g]] = need.get(cf[g], 0) + 1
        comp[v] = {f"{k[0]}|{k[1]}": int(n) for k, n in sorted(need.items(), key=lambda x: str(x))}
        for c, n in need.items():
            assert len(pool_by_cell.get(c, [])) >= n, (v, c, n, len(pool_by_cell.get(c, [])))
        D = np.empty((N_SUB, len(gl)), dtype=np.int64)
        for r in range(N_SUB):
            sel = []
            for c in sorted(need, key=str):
                pc = pool_by_cell[c]; sel += [pc[i] for i in rng.choice(len(pc), need[c], replace=False)]
            D[r] = [upos[g] for g in sel]
        draws[v] = D
    istf_u = np.array([cell[g][0] == "TF" for g in universe])
    rows = []
    null_store = {}
    for key in L.STATIC_KEYS:
        g, E = L.load_static(key); gpos = {x: i for i, x in enumerate(g)}
        have = np.array([u in gpos for u in universe])
        Xu = np.zeros((len(universe), E.shape[1])); Xu[have] = E[[gpos[u] for u in np.array(universe)[have]]]
        K = Xu @ Xu.T

        def stat(idx):
            Ks = K[np.ix_(idx, idx)]; r = Ks.mean(1); c = r.mean()
            Kc = Ks - r[:, None] - r[None, :] + c
            nr = np.sqrt(np.maximum(np.diag(Kc), 0)) + 1e-12
            Cs = Kc / nr[:, None] / nr[None, :]
            t = istf_u[idx]; k = t.sum()
            tot = Cs[:, t].sum(1); sc = np.where(t, (tot - np.diag(Cs)) / max(k - 1, 1), tot / k)
            return S.auroc(sc[t], sc[~t])
        for v, gl in variants.items():
            gl_m = [x for x in gl if x in gpos]
            idx_leg = np.array([upos[x] for x in gl_m])
            obs = stat(idx_leg)
            # direct check (estimator (i) as implemented in e6_run)
            Xl = E[[gpos[x] for x in gl_m]].astype(np.float64); Xl -= Xl.mean(0)
            tfp = np.flatnonzero(istf_u[idx_leg]); tgp = np.flatnonzero(~istf_u[idx_leg])
            sdir = S.tf_scores_loo(Xl, tfp); assert abs(S.auroc(sdir[tfp], sdir[tgp]) - obs) < 1e-6, (key, v)
            nul = np.array([stat(d) for d in draws[v]])
            ns = L.null_summary(nul)
            pct = float((np.sum(nul < obs) + 0.5 * np.sum(nul == obs)) / nul.size)
            p2 = min(1.0, 2 * min((np.sum(nul >= obs) + 1) / (nul.size + 1), (np.sum(nul <= obs) + 1) / (nul.size + 1)))
            rows.append(dict(variant=v, model=key, label=L.STATIC[key][1], n_legacy_genes_used=len(gl_m),
                             legacy_value=obs, null_N=ns["N"], null_mean=ns["mean"], null_sd=ns["sd"],
                             null_q025=ns["q025"], null_q50=ns["q50"], null_q975=ns["q975"],
                             legacy_percentile=pct, p_two_sided=p2, legacy_minus_null_mean=obs - ns["mean"]))
            null_store[f"{v}|{key}"] = nul.astype(np.float32)
            print(v, key, round(obs, 4), round(ns["mean"], 4), round(pct, 4), flush=True)
    df = pd.DataFrame(rows)
    np.savez(L.WORK / "selection_null_draws.npz", **null_store)
    df.to_csv(L.WORK / "selection_bias_core.csv", index=False)
    (L.WORK / "selection_meta.json").write_text(json.dumps(dict(
        pool_size=len(pool), n_TF_pool=int(sum(cell[g][0] == "TF" for g in pool)),
        decile_edges=q.tolist(), composition=comp, degree_quartile_edges={k: v.tolist() for k, v in qdeg.items()},
        legacy_mean_trrust_degree={L_: float(np.mean([degv(g) for g in leg if lab(g) == L_])) for L_ in ("TF", "TG")},
        pool_mean_trrust_degree={L_: float(np.mean([degv(g) for g in pool if lab(g) == L_])) for L_ in ("TF", "TG")}, n_subsets=N_SUB,
        legacy_not_in_lungV=[g for g in leg if g not in mexp]), indent=1))


if __name__ == "__main__":
    main()
