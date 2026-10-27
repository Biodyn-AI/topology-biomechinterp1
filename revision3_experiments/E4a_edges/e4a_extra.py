"""E4a extra sensitivity rows (added after the kidney_A pipeline test showed that the pre-specified shortcut
row 'target TRRUST in-degree' is not ~0.5; the PRIMARY specification is unchanged).

Why: TRRUST in-degree of a positive's target counts the pair's own edge, so 'total in-degree' gives every
positive +1 relative to an otherwise identical candidate; and the in-degree deciles used for matching
collapse when most targets have in-degree 1-2 (np.unique of tied quantiles). Rows written here:
  (a) shortcut 'other-TF in-degree' = total in-degree minus the pair's own documented edge (y);
  (b) primary scorer (signed cosine, 12 layers) and |co-expression| restricted to the fully matched
      candidates (match_level 0) plus all positives;
  (c) primary scorer and |co-expression| stratified within (TF x exact other-TF in-degree) cells, with the
      same two-way TF x target cluster bootstrap (clusters still TF and target) -- removes any residual
      target-popularity contribution;
  (d) the legacy-mapped bands (E9_screen/legacy_band_mapping.json) if that file exists.
All are sensitivity analyses: effect sizes with 95% CIs, no pass/fail.
Usage: python e4a_extra.py <ctx> ...
"""
from __future__ import annotations
import json, sys, time, zlib
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; R3 = HERE.parent
sys.path.insert(0, str(R3)); sys.path.insert(0, str(HERE))
import stats3 as S   # noqa: E402
import e4a_lib as L  # noqa: E402
from e4a_main import rowdict, scores_layer, seed_of, jdump, log, NL, B, EMB_LOCAL, CACHE  # noqa: E402

MAPPING = R3 / "E9_screen" / "legacy_band_mapping.json"


def run(ctx):
    W = HERE / "work" / ctx; f = W / "sens_extra.json"
    if f.exists(): log(ctx, "extra exists"); return
    h3 = dict(np.load(CACHE / ctx / "h3_pairs.npz"))
    sp = S.StratifiedPairs(h3["tf"], h3["tg"], h3["y"]); keep = np.isin(h3["tf"], sp.tfs)
    lev = h3["match_level"][keep]
    cx = dict(np.load(W / "coexp.npz")); seed = seed_of(ctx)
    indeg_tot = cx["indeg"][sp.tg]; indeg_oth = indeg_tot - sp.y
    out = dict(context=ctx, rows=[])
    fs = L.FastStrat(sp); wts = fs.boot_weights(B, np.random.default_rng([seed, 2]))   # same resamples as primary
    for nm, s in (("shortcut: target TRRUST in-degree (total, incl. own edge)", indeg_tot),
                  ("shortcut: target other-TF in-degree (excl. own edge)", indeg_oth)):
        r = rowdict(nm, -1, sp.auroc(s), fs.auc(s, wts)); r["block"] = "shortcut"; out["rows"].append(r)
    # (b) fully matched candidates only
    m = sp.y | (lev == 0)
    sub = S.StratifiedPairs(sp.tf[m], sp.tg[m], sp.y[m]); ii = np.flatnonzero(m)[np.isin(sp.tf[m], sub.tfs)]
    fsb = L.FastStrat(sub); wb = fsb.boot_weights(B, np.random.default_rng([seed, 50]))
    ce = sub.auroc(cx["cx_ext"][ii]); cb = fsb.auc(cx["cx_ext"][ii], wb)
    r = rowdict("|co-expression| extraction cells", -1, ce, cb); r.update(block="match_level 0 only", n_pos=int(sub.y.sum()), n_neg=int((~sub.y).sum())); out["rows"].append(r)
    r = rowdict("shortcut: target other-TF in-degree (excl. own edge)", -1, sub.auroc(indeg_oth[ii]), fsb.auc(indeg_oth[ii], wb)); r["block"] = "match_level 0 only"; out["rows"].append(r)
    for l in range(NL):
        s = scores_layer(W, l)["prim"][ii]
        r = rowdict("signed cosine full-512", l, sub.auroc(s), fsb.auc(s, wb), ce, cb); r.update(block="match_level 0 only", n_pos=int(sub.y.sum()), n_neg=int((~sub.y).sum())); out["rows"].append(r)
    # (c) strata = TF x exact other-TF in-degree (clusters unchanged)
    key = sp.tf_code.astype(np.int64) * 100000 + indeg_oth.astype(np.int64)
    fsc = L.FastStrat(sp); fsc.strata = np.unique(key, return_inverse=True)[1]
    n_eff = 0
    for c in np.unique(fsc.strata):
        mm = fsc.strata == c; n_eff += int(sp.y[mm].sum() > 0 and (~sp.y[mm]).sum() > 0) * int(sp.y[mm].sum())
    ce = fsc.auc(cx["cx_ext"])[0]; cb = fsc.auc(cx["cx_ext"], wts)
    r = rowdict("|co-expression| extraction cells", -1, ce, cb); r.update(block="strata TF x other-TF in-degree", n_pos_in_informative_strata=n_eff); out["rows"].append(r)
    for l in range(NL):
        s = scores_layer(W, l)["prim"]
        r = rowdict("signed cosine full-512", l, fsc.auc(s)[0], fsc.auc(s, wts), ce, cb)
        r.update(block="strata TF x other-TF in-degree", n_pos_in_informative_strata=n_eff); out["rows"].append(r)
    # (d) legacy-mapped bands
    out["mapping_file_present"] = MAPPING.exists()
    if MAPPING.exists():
        mp = json.loads(MAPPING.read_text()); out["mapping_raw_keys"] = list(mp)[:50]
        out["mapped_bands"] = mapped_bands(ctx, sp, fs, wts, mp, W, cx)
    jdump(out, f); log(ctx, "extra done")


def _cols(band):
    a, b = band.replace("SV", "").split("-"); return np.arange(int(a) - 1, int(b))


def mapped_bands(ctx, sp, fs, wts, mp, W, cx):
    """Cosine on the observed-gene SVD band that E9_screen/legacy_band_mapping.json maps each legacy band to
    (entries whose 'context' is this context). Two variants: the layer-averaged mapped band at every layer
    ('mapped_observed_band') and the per-layer winner ('per_layer_winner'). Hypothesis-generating only."""
    ents = [e for e in mp.get("bands", []) if e.get("context") == ctx]
    if not ents:
        return dict(status="no mapping entry for this context")
    emb = np.load(EMB_LOCAL / f"{ctx}.npy", mmap_mode="r"); oi = np.load(CACHE / ctx / "obs_idx.npy")
    ce = sp.auroc(cx["cx_ext"]); cb = fs.auc(cx["cx_ext"], wts); rows = []; done = {}
    for l in range(NL):
        US, _, _ = S.svd_observed(emb, l, oi)
        for e in ents:
            for variant, band in (("layer-averaged mapping", e["mapped_observed_band"]),
                                  ("per-layer winner", e["per_layer_winner"][str(l)]["observed_band"])):
                key = (l, band)
                if key not in done:
                    s = np.round(L.rowcos(S.l2norm_rows(US[:, _cols(band)]), sp.tf, sp.tg), 10)
                    done[key] = (sp.auroc(s), fs.auc(s, wts))
                est, bt = done[key]
                r = rowdict(f"mapped band: legacy {e['legacy_band']} [{e['legacy_convention']}] -> observed {band} ({variant})", l, est, bt, ce, cb)
                r.update(legacy_band=e["legacy_band"], legacy_convention=e["legacy_convention"], observed_band=band, variant=variant,
                         mapping_mean_cos2=e["mean_cos2_layer_avg"]); rows.append(r)
    return dict(status="ok", rows=rows)


if __name__ == "__main__":
    for c in sys.argv[1:]:
        run(c)
