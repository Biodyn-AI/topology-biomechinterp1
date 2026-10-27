"""E9 / A6 step 6 -- crosswalk between the screen's ("legacy") SVD directions and the observed-gene SVD.

Geometry only: no outcome (TRRUST / STRING / label) data are read here, so the band mapping written to
legacy_band_mapping.json is fixed before any outcome analysis.

Conventions compared (right-singular vectors, 512-d feature space; sign-free |cos|):
  observed            stats3.svd_observed on the observed (non-zero) rows, mean-centred over them (CONVENTION3)
  legacy_full_zero    SVD of the FULL matrix incl. all-zero rows, mean-centred over all rows
                      (kidney_A: 4,803 rows; the submitted manuscript's convention -- defect D2; also the
                       convention of screen iters 6-7 SV1 analyses on kidney_A and, on a 2,000-row random
                       subsample, of iter_0046's effective-rank/SV1 claims on immune)
  screen_209_true     SVD of the 209 edge-dataset genes at their TRUE rows (22 of them all-zero), mean-centred
                      over the 209 (kidney_A; screen iters 8-26: SV1/SV2/SV3/SV4 axes, SV2-SV4 PPI axes)
  screen_195_misaligned  SVD of rows 0..208 minus the 14 zero rows among them (the '195 in-vocab genes' of
                      iters 29-40 under the misaligned gene_list.txt mapping) (kidney_A)
  screen_nonzero      the screen's own phrase "SVD of centered nonzero embeddings" (immune iters 47-63:
                      [2039,512]; kidney_A in phase-2 iterations) -- computed independently here and compared
  screen_nonzero_sum  kidney_A observed rows of layer_gene_embedding_sum.npy (iters 0059-0060 used sums, not means)
"""
from __future__ import annotations
import json, time
import numpy as np, pandas as pd
import e9_lib as L

OUT = L.HERE
NSV = 12
BANDS = {"SV2-4": (2, 4), "SV5-7": (5, 7), "SV2-7": (2, 7)}


def top_vt(M, k=40):
    _, s, vt = L.svd_rows(M)
    return vt[:k], s


def compare(ctx, conv, layer, Vleg, Vobs, rows_v, rows_b):
    C = np.abs(Vleg[:NSV] @ Vobs[:40].T)                 # 12 x 40
    for k in range(NSV):
        j = int(C[k].argmax())
        rows_v.append(dict(kind="vector", context=ctx, legacy_convention=conv, layer=layer, legacy_sv=f"SV{k+1}",
                           best_observed_sv=f"SV{j+1}", abs_cos_best=float(C[k, j]),
                           abs_cos_same_index=float(C[k, k]),
                           second_best_abs_cos=float(np.sort(C[k])[-2]),
                           norm_in_observed_SV1_12=float(np.linalg.norm(Vobs[:NSV] @ Vleg[k]))))
    for lname, (a, b) in [("SV2-7", (2, 7)), ("SV2-4", (2, 4)), ("SV5-7", (5, 7))]:
        w = b - a + 1
        targets = {"SV1-6": (1, 6), "SV2-7": (2, 7)} if lname == "SV2-7" else {}
        for tname, (c, d) in targets.items():
            cs = L.principal_cosines(Vleg[a - 1:b], Vobs[c - 1:d])
            rows_b.append(dict(kind="band_principal_angles", context=ctx, legacy_convention=conv, layer=layer,
                               legacy_band=lname, observed_band=tname,
                               principal_angles_deg=";".join(f"{x:.2f}" for x in np.degrees(np.arccos(cs))),
                               mean_cos2=float((cs ** 2).mean()), max_angle_deg=float(np.degrees(np.arccos(cs.min())))))
        # alignment with every contiguous observed band of the same width ending at <= SV12
        for start in range(1, NSV - w + 2):
            cs = L.principal_cosines(Vleg[a - 1:b], Vobs[start - 1:start - 1 + w])
            rows_b.append(dict(kind="band_scan", context=ctx, legacy_convention=conv, layer=layer, legacy_band=lname,
                               observed_band=f"SV{start}-{start + w - 1}",
                               principal_angles_deg=";".join(f"{x:.2f}" for x in np.degrees(np.arccos(cs))),
                               mean_cos2=float((cs ** 2).mean()), max_angle_deg=float(np.degrees(np.arccos(cs.min())))))


def main():
    t0 = time.time()
    rows_v, rows_b, meta = [], [], {}
    # ---------------- kidney_A ----------------
    emb, genes = X_load("kidney_A")
    nz = L.nonzero_rows(emb); obs_idx = np.flatnonzero(nz)
    true209 = L.named_genes_true(); idx209 = np.array(sorted(true209.values()))
    mis = L.named_genes_misaligned(); rows_mis = np.array(sorted(mis.values()))
    rows195 = rows_mis[nz[rows_mis]]
    emb_sum = np.load(L.CYC1 / "layer_gene_embedding_sum.npy", mmap_mode="r")
    nz_sum = np.linalg.norm(np.asarray(emb_sum[0], dtype=np.float64), axis=1) > 0
    meta["kidney_A"] = dict(G=int(emb.shape[1]), n_observed=int(nz.sum()), n_209=int(idx209.size),
                            n_209_zero_rows=int((~nz[idx209]).sum()), misaligned_rows=f"{rows_mis.min()}..{rows_mis.max()}",
                            n_misaligned_nonzero=int(rows195.size),
                            sum_matrix_observed_rows_equal=bool(np.array_equal(nz_sum, nz)))
    for l in range(emb.shape[0]):
        _, _, Vobs = L.S.svd_observed(emb, l, obs_idx)
        E = np.asarray(emb[l], dtype=np.float64)
        Vleg, _ = top_vt(E)
        compare("kidney_A", "legacy_full_zero", l, Vleg, Vobs, rows_v, rows_b)
        V209, _ = top_vt(E[idx209])
        compare("kidney_A", "screen_209_true", l, V209, Vobs, rows_v, rows_b)
        V195, _ = top_vt(E[rows195])
        compare("kidney_A", "screen_195_misaligned", l, V195, Vobs, rows_v, rows_b)
        Vnz, _ = top_vt(E[nz])
        compare("kidney_A", "screen_nonzero", l, Vnz, Vobs, rows_v, rows_b)
        Vsum, _ = top_vt(np.asarray(emb_sum[l], dtype=np.float64)[nz_sum])
        compare("kidney_A", "screen_nonzero_sum", l, Vsum, Vobs, rows_v, rows_b)
    # ---------------- immune ----------------
    emb, genes = X_load("immune")
    nz = L.nonzero_rows(emb); obs_idx = np.flatnonzero(nz)
    meta["immune"] = dict(G=int(emb.shape[1]), n_observed=int(nz.sum()))
    for l in range(emb.shape[0]):
        _, _, Vobs = L.S.svd_observed(emb, l, obs_idx)
        E = np.asarray(emb[l], dtype=np.float64)
        Vnz, _ = top_vt(E[nz])
        compare("immune", "screen_nonzero", l, Vnz, Vobs, rows_v, rows_b)
        Vleg, _ = top_vt(E)
        compare("immune", "legacy_full_zero", l, Vleg, Vobs, rows_v, rows_b)
    dv = pd.DataFrame(rows_v); db = pd.DataFrame(rows_b)
    pd.concat([dv, db[db.kind == "band_principal_angles"]], ignore_index=True).to_csv(OUT / "crosswalk.csv", index=False)
    db[db.kind == "band_scan"].to_csv(OUT / "crosswalk_band_scan.csv", index=False)

    # ---------------- fixed legacy-band mapping ----------------
    # which legacy band was defined on which context/convention in the screen (see band_provenance.md)
    DEFINED_ON = [
        ("SV2-4", "kidney_A", "screen_209_true", "iters 13-19: SV2, SV3, SV4 PPI co-pole axes on the 209-gene subset"),
        ("SV2-4", "immune", "screen_nonzero", "iter_0049 onwards: SV2-SV4 TRRUST proximity / TF-target tests"),
        ("SV5-7", "immune", "screen_nonzero", "iter_0050 onwards"),
        ("SV2-7", "immune", "screen_nonzero", "iter_0056 onwards (joint 6-D classifier)"),
        ("SV5-7", "kidney_A", "screen_nonzero_sum", "iters 0059-0060 (layer_gene_embedding_sum.npy, observed rows)"),
        ("SV5-7", "kidney_A", "screen_nonzero", "phase-2 iterations that load cycle1_main as 'kidney' (mean embeddings, nonzero rows)"),
        ("SV2-4", "kidney_A", "legacy_full_zero", "submitted manuscript (D2 convention)"),
        ("SV5-7", "kidney_A", "legacy_full_zero", "submitted manuscript (D2 convention)"),
        ("SV2-7", "kidney_A", "legacy_full_zero", "submitted manuscript (D2 convention)"),
    ]
    mapping = {"_note": ("Fixed from geometry alone (no outcome data read). For each legacy band, the contiguous "
                         "observed-basis band of the same width (ending at or before SV12) with maximal mean "
                         "cos^2 of principal angles, averaged over the 12 layers; per-layer winners also given. "
                         "Observed basis = stats3.svd_observed (CONVENTION3)."),
               "created_unix": time.time(), "meta": meta, "bands": []}
    for band, ctx, conv, where in DEFINED_ON:
        s = db[(db.kind == "band_scan") & (db.context == ctx) & (db.legacy_convention == conv) & (db.legacy_band == band)]
        avg = s.groupby("observed_band").mean_cos2.mean().sort_values(ascending=False)
        per_layer = {int(l): g.sort_values("mean_cos2", ascending=False).iloc[0][["observed_band", "mean_cos2"]].to_dict()
                     for l, g in s.groupby("layer")}
        winners = pd.Series([v["observed_band"] for v in per_layer.values()]).value_counts().to_dict()
        mapping["bands"].append(dict(legacy_band=band, context=ctx, legacy_convention=conv, where_used=where,
                                     mapped_observed_band=avg.index[0], mean_cos2_layer_avg=float(avg.iloc[0]),
                                     runner_up=avg.index[1], runner_up_mean_cos2=float(avg.iloc[1]),
                                     identity_band_mean_cos2=float(avg.get(band, np.nan)),
                                     per_layer_winner=per_layer, n_layers_won=winners))
    (OUT / "legacy_band_mapping.json").write_text(json.dumps(mapping, indent=1, default=float))
    print(json.dumps(mapping["bands"], indent=1, default=float)[:6000])
    print("elapsed", round(time.time() - t0, 1), "s")


def X_load(ctx):
    return L.X.load_context(ctx)


if __name__ == "__main__":
    main()
