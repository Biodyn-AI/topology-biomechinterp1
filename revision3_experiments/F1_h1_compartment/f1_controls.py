"""F1 post-hoc descriptive controls for part (c) (NOT pre-specified; no test attached): which part of the
30-GO-CC + log n_obs design removes the assortativity, and is GO cellular component special?
Residualise SVk coordinates (OLS, f1_common.residualise) on alternative designs and report the observed residual
assortativity over the STRING edges. Under the Maslov-Sneppen null the residual-r null means in part (c) were
-0.002 to 0.000 (residual_assort.csv), so observed residual r is ~ the excess over MS.
designs: lnobs | cc30 | cc30+lnobs (= part c) | class8_dummies+lnobs (= residclass) | compartment8_multilabel+lnobs
(the 8 GO terms defining the classes, multi-label) | bp30+lnobs | mf30+lnobs (30 most frequent direct GO BP / MF terms
among V) | random30+lnobs (30 random gene sets with the sizes of the cc30 terms; 5 seeds).
Writes resid_controls.csv (per context x design x cell) and resid_controls_summary.csv.
"""
import numpy as np, pandas as pd
from collections import Counter
import f1_common as F
E, S = F.E, F.S

def top_terms(genes, g2go, info, ns, k=30):
    cnt = Counter(); per = []
    for g in genes:
        s = {t for t in g2go.get(g, ()) if info.get(t, ("", ""))[0] == ns}; per.append(s); cnt.update(s)
    ids = [t for t, _ in sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0]))[:k]]
    return np.array([[t in s for t in ids] for s in per], float), ids

E.frozen_hash_check()
g2go = F.load_gene2go(); info = F.go_terms_info()
rows, summ = [], []
COMP8 = ["GO:0005615", "GO:0005886", "GO:0005783", "GO:0005788", "GO:0005794", "GO:0005739", "GO:0005634", "GO:0005829"]
for ctx in ["lung", "lung_external", "kidney_B"]:
    d = dict(np.load(F.CKPT / f"{ctx}_setup.npz", allow_pickle=True)); Xc = d["Xc"]; cls = d["cls"]
    pi, pj = d["edges"][:, 0], d["edges"][:, 1]; genes = E.load_ctx(ctx)["genes"]; n = len(genes)
    ln = np.log(d["ncells"].astype(float))[:, None]
    CC, _ = top_terms(genes, g2go, info, "cellular_component"); BP, _ = top_terms(genes, g2go, info, "biological_process")
    MF, _ = top_terms(genes, g2go, info, "molecular_function")
    C8 = np.array([[t in set(g2go.get(g, ())) for t in COMP8] for g in genes], float)
    dummies = np.column_stack([(cls == k).astype(float) for k in np.unique(cls)[1:]])
    designs = {"lnobs": ln, "cc30": CC, "cc30+lnobs": np.column_stack([CC, ln]),
               "class8_dummies+lnobs": np.column_stack([dummies, ln]), "compartment8_multilabel+lnobs": np.column_stack([C8, ln]),
               "bp30+lnobs": np.column_stack([BP, ln]), "mf30+lnobs": np.column_stack([MF, ln])}
    sizes = CC.sum(0).astype(int)
    for s in range(5):
        rng = np.random.default_rng(E.seed_for("F1", ctx, "random30", s))
        Rn = np.zeros((n, 30))
        for j, sz in enumerate(sizes): Rn[rng.choice(n, sz, replace=False), j] = 1
        designs[f"random30+lnobs_seed{s}"] = np.column_stack([Rn, ln])
    r = E.assort_cols(Xc, pi, pj)
    shared = {"cc30": CC, "bp30": BP, "mf30": MF}
    # fraction of edges whose endpoints share >= 1 of the design's terms, vs the same for MS-free expectation (random pairs)
    rp = np.random.default_rng(E.seed_for("F1", ctx, "randpairs")); ra, rb = rp.integers(0, n, 200000), rp.integers(0, n, 200000)
    for nm, D in designs.items():
        R, r2 = F.residualise(Xc, D); o = E.assort_cols(R, pi, pj)
        for c in range(120):
            L, k = E.cell_name(c)
            rows.append(dict(context=ctx, design=nm, layer=L, sv=k, r_raw=r[c], r_resid=o[c], ratio=o[c] / r[c], R2=r2[c]))
        top = [F.col_of(L, k) for L, k in F.E2_TOP[ctx]]
        key = nm.split("+")[0]
        sh = shared.get(key)
        summ.append(dict(context=ctx, design=nm, n_regressors=D.shape[1], median_ratio=float(np.median(o / r)),
                         median_R2=float(np.median(r2)), top1_cell=f"L{F.E2_TOP[ctx][0][0]} SV{F.E2_TOP[ctx][0][1]}",
                         top1_r_raw=float(r[top[0]]), top1_r_resid=float(o[top[0]]), top1_ratio=float(o[top[0]] / r[top[0]]),
                         top1_R2=float(r2[top[0]]), top3_mean_ratio=float(np.mean(o[top] / r[top])),
                         frac_edges_sharing_a_term=(float(((D[pi] * D[pj])[:, :sh.shape[1]].sum(1) > 0).mean()) if sh is not None else np.nan),
                         frac_random_pairs_sharing_a_term=(float(((sh[ra] * sh[rb]).sum(1) > 0).mean()) if sh is not None else np.nan)))
    print(ctx, "done", flush=True)
pd.DataFrame(rows).to_csv(F.HERE / "resid_controls.csv", index=False)
sm = pd.DataFrame(summ); sm.to_csv(F.HERE / "resid_controls_summary.csv", index=False)
pd.set_option("display.width", 250); print(sm.round(3).to_string())
