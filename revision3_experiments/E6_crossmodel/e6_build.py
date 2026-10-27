"""E6 step 1: vocabularies, spectral caches, gene sets (Tier A, Tier B x 5 contexts, legacy) and attrition.csv."""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import e6lib as L

S = L.S


def main():
    rf = L.refs(); stg = L.string700()
    print("STRING>=700 pairs genome-wide:", len(stg), flush=True)
    voc = {k: set(L.vocab(k)) for k in L.STATIC_KEYS}
    raw_n = {}
    for k in L.STATIC_KEYS:
        z = np.load(L.REPO / L.STATIC[k][0]); raw_n[k] = int(z["gene_names"].shape[0])
    for k in L.STATIC_KEYS:
        L.static_spectral(k); print("spectral cached", k, flush=True)
    A = sorted(set.intersection(*voc.values()))
    print("Tier A", len(A), flush=True)
    att = []
    for k in L.STATIC_KEYS:
        att.append(dict(stage="vocabulary", set=k, n_genes=len(voc[k]),
                        note=f"{raw_n[k]} rows in npz; {raw_n[k]-len(voc[k])} all-zero rows dropped"))
    core8 = [k for k in L.STATIC_KEYS if k not in ("esm2_650m", "seq_features", "hpa_coexp")]
    att.append(dict(stage="intersection", set="8 GPL-rollout models (scGPT..ESM2-UCE)",
                    n_genes=len(set.intersection(*[voc[k] for k in core8])), note=""))
    att.append(dict(stage="intersection", set="Geneformer V1-10M & V2-104M & ESM2-UCE (submitted models)",
                    n_genes=len(voc["gf_v1_10m"] & voc["gf_v2_104m"] & voc["esm2_uce"]), note=""))
    att.append(dict(stage="Tier A", set="intersection of all 11 static vocabularies", n_genes=len(A), note=""))
    Aset = set(A)
    TR = rf["trrust"]

    def trrust_part(G):
        return sorted({g for p in TR if p[0] in G and p[1] in G for g in p})
    att.append(dict(stage="Tier A", set="TRRUST-participating within Tier A", n_genes=len(trrust_part(Aset)), note=""))
    att.append(dict(stage="Tier A", set="Lambert TFs in Tier A", n_genes=int(sum(g in rf["lambert"] for g in A)), note=""))

    # Tier A detection deciles: lung deciles where the gene is in lung V, else the median lung decile
    lg, _ = L.ctx_obs("lung"); nl = np.load(L.CACHE / "lung/node_labels.npz")
    dl = dict(zip(lg, nl["det_decile"].astype(int)))
    inl = [dl[g] for g in A if g in dl]; med = int(np.median(inl))
    detA = np.array([dl.get(g, med) for g in A])
    gsA = L.build_geneset("tierA", A, detA, rf, stg); L.save_geneset(gsA)
    print("tierA built", len(A), "STRING", gsA["pi"].size, "TRRUST edges", int(gsA["e_y"].sum()),
          "pairs", gsA["e_y"].size, flush=True)
    meta = dict(tierA=dict(n=len(A), n_in_lung_V=len(inl), median_decile_imputed=med,
                           n_imputed=len(A) - len(inl)))

    for ctx in L.CORE_CTX:
        og, oi = L.ctx_obs(ctx); nlab = np.load(L.CACHE / ctx / "node_labels.npz")
        dd = dict(zip(og, nlab["det_decile"].astype(int)))
        Gc = len(X_genes(ctx))
        B_ = sorted(set(og) & Aset)
        gs = L.build_geneset(f"tierB_{ctx}", B_, np.array([dd[g] for g in B_]), rf, stg); L.save_geneset(gs)
        att += [dict(stage=f"context {ctx}", set="gene axis (h5ad var)", n_genes=Gc, note=""),
                dict(stage=f"context {ctx}", set="V (observed, non-zero embedding)", n_genes=len(og), note=""),
                dict(stage=f"context {ctx}", set="TRRUST-participating within V", n_genes=len(trrust_part(set(og))), note=""),
                dict(stage=f"context {ctx}", set="Tier B = V & Tier A", n_genes=len(B_), note=""),
                dict(stage=f"context {ctx}", set="TRRUST-participating within Tier B", n_genes=gs["n_trrust_participating"], note=""),
                dict(stage=f"context {ctx}", set="Lambert TFs in Tier B", n_genes=int(gs["lambert"].sum()), note="")]
        meta[f"tierB_{ctx}"] = dict(n=len(B_), string_pairs=int(gs["pi"].size), trrust_within=gs["n_trrust_within"],
                                   h3_pos=int(gs["e_y"].sum()), h3_neg=int((1 - gs["e_y"]).sum()))
        print(ctx, len(B_), "STRING", gs["pi"].size, "TRRUST pos", int(gs["e_y"].sum()), flush=True)

    # legacy chain
    leg = json.loads((L.HERE / "legacy_genes.json").read_text())
    tsv = L.REPO / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main/cycle1_edge_dataset.tsv"
    e = pd.read_csv(tsv, sep="\t"); named = set(e.source.astype(str)) | set(e.target.astype(str))
    kg = X_genes("kidney_A"); in_axis = sorted(g for g in kg if g in named)
    k_og, _ = L.ctx_obs("kidney_A")
    lost_nz = sorted(set(in_axis) - set(leg["curated_nonzero_kidneyA"]))
    lost_vocab = sorted(set(leg["curated_nonzero_kidneyA"]) - set(leg["legacy186"]))
    l186 = leg["legacy186"]; l178 = sorted(set(l186) & Aset)
    att += [dict(stage="legacy", set="genes named in cycle1_edge_dataset.tsv", n_genes=len(named),
                 note=f"{len(e)} rows (label 1/0) of the round-1 TRRUST-filtered kidney_A edge list"),
            dict(stage="legacy", set="... in the kidney_A gene axis", n_genes=len(in_axis), note=""),
            dict(stage="legacy", set="... non-zero in kidney_A at L0 (curated set)", n_genes=len(leg["curated_nonzero_kidneyA"]),
                 note=f"{len(lost_nz)} never-observed genes dropped"),
            dict(stage="legacy", set="... & Geneformer V1-10M & V2-104M & ESM2-UCE vocabularies (submitted Table 1)",
                 n_genes=len(l186), note="lost: " + ",".join(lost_vocab)),
            dict(stage="legacy", set="... & Tier A (all 11 vocabularies)", n_genes=len(l178),
                 note="lost (absent from tGPT): " + ",".join(sorted(set(l186) - Aset))),
            dict(stage="legacy", set="kidney_A V & Tier A, TRRUST-participating (what the vocabularies would have allowed)",
                 n_genes=len(trrust_part(set(k_og) & Aset)), note="")]
    for g in l186: assert g in set(k_og)
    nk = np.load(L.CACHE / "kidney_A/node_labels.npz"); dk = dict(zip(k_og, nk["det_decile"].astype(int)))
    for nm, gl in (("legacy186", l186), ("legacy178", l178)):
        gs = L.build_geneset(nm, sorted(gl), np.array([dk[g] for g in sorted(gl)]), rf, stg); L.save_geneset(gs)
        meta[nm] = dict(n=len(gl), string_pairs=int(gs["pi"].size), trrust_within=gs["n_trrust_within"],
                        h3_pos=int(gs["e_y"].sum()), h3_tfs=int(np.unique(gs["e_tf"]).size),
                        n_trrust_tf=int(gs["trrust_tf"].sum()), n_trrust_tgt=int(gs["trrust_tgt"].sum()),
                        n_lambert=int(gs["lambert"].sum()))
        print(nm, meta[nm], flush=True)
    pd.DataFrame(att).to_csv(L.HERE / "attrition.csv", index=False)
    (L.WORK / "build_meta.json").write_text(json.dumps(meta, indent=1))


def X_genes(ctx):
    return L.X.h5ad_genes(L.X.CONTEXTS[ctx]["h5ad"])


if __name__ == "__main__":
    main()
