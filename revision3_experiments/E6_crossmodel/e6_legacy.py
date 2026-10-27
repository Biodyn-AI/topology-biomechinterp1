"""E6 legacy bridge: submitted-style values for every model on the 186-gene set, and the functional-STRING
PPI bridge under estimator (i). Writes work/legacy_bridge.csv (assembled into legacy_reproduction.csv)."""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import e6lib as L
S = L.S
sys.path.insert(0, str(L.R3.parent / "revision2_experiments"))
import common as C  # noqa: E402  (submitted round-2 estimators, read-only)


def main():
    gs = L.load_geneset("legacy186"); genes = gs["genes"]; n = len(genes); loc = {g: k for k, g in enumerate(genes)}
    obj = json.loads(C.STRING_CACHE.read_text())
    st = sorted({(p["g1"], p["g2"]) if p["g1"] < p["g2"] else (p["g2"], p["g1"]) for p in obj["pairs"]
                 if p["score"] >= 0.7 and p["g1"] in loc and p["g2"] in loc and p["g1"] != p["g2"]})
    pi = np.array([loc[a] for a, b in st]); pj = np.array([loc[b] for a, b in st])
    a, b, acc = next(S.maslov_sneppen_chain(pi, pj, n, 1, seed=L.SEED + 77, burn_factor=10))
    true = set(zip(pi.tolist(), pj.tolist()))
    k = np.array([(int(x), int(y)) not in true for x, y in zip(np.minimum(a, b), np.maximum(a, b))]); ni, nj = a[k], b[k]
    tfs = set(L.refs()["trrust_tf"]); tf_pos = np.array([i for i, g in enumerate(genes) if g in tfs])
    tg_pos = np.array([i for i, g in enumerate(genes) if g in L.refs()["trrust_target_only"]])
    rows = []
    embs = [k_ for k_ in L.STATIC_KEYS if k_ != "tgpt"] + [f"kidney_A@L{l}" for l in range(12)]
    import e6_run as R
    for e in embs:
        P = R.provide(e, genes)
        # submitted estimator: SVD on the 186-gene subset, SV2..SV7, LOO cosine (round-2 code)
        Xr = P["E"]           # already centred over the 186 genes; the round-2 us_of centres the same way
        U, Sv, _ = np.linalg.svd(Xr, full_matrices=False); US = U * Sv
        au_sub = C.tf_target_auroc(US, slice(1, 7), tf_pos, tg_pos)
        lab = C.pole_labels(US[:, 1]); cop = C.copole_rate(lab, pi, pj)
        nul = C.copole_null(lab, pi, pj, 20000, np.random.default_rng(L.SEED + 78))
        z = (cop - nul.mean()) / (nul.std() + 1e-12)
        En = S.l2norm_rows(P["E"]); cp = L.pair_cos(En, pi, pj); cn = L.pair_cos(En, ni, nj)

        def stat(nodes, a_, b_, cp=cp, cn=cn):
            kp = np.zeros(n, bool); kp[nodes] = True
            x = cp[kp[pi] & kp[pj]]; y = cn[kp[ni] & kp[nj]]
            return S.auroc(x, y) if x.size > 2 and y.size > 2 else np.nan
        est, lo, hi, _ = S.vertex_subsample_ci(stat, n, pi, pj, 200, np.random.default_rng(L.SEED + 79))
        rows.append(dict(embedding=e, label=P["label"], d=P["d"],
                         submitted_style_tf_target_auroc_subsetSVD_SV2_SV7=au_sub,
                         submitted_style_copole_SV2_functional=cop, submitted_style_copole_z_label_shuffle=z,
                         functional_string_pairs=int(pi.size), functional_rewired_nonpairs=int(ni.size),
                         ppi_functional_auroc_full=est, ppi_functional_ci_lo=lo, ppi_functional_ci_hi=hi))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(L.WORK / "legacy_bridge.csv", index=False)


if __name__ == "__main__":
    main()
