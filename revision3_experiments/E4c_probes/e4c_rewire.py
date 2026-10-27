"""Generate the 100 degree-preserving rewired label vectors per context (shared by every layer, feature,
regime and model; M2/M3 use draws 0-9). Writes rewired/<ctx>.npz and rewire_stats.csv."""
import sys
import numpy as np, pandas as pd
import e4c_lib as L

OUT = L.HERE / "rewired"; OUT.mkdir(exist_ok=True)
rows = []
for ci, ctx in enumerate(["lung", "lung_external", "kidney_B"]):
    c = L.Ctx(ctx)
    Yr, st = L.rewire_draws(c, 100, seed=L.SEED + 1000 * ci, burn=10, thin=5)
    deg_tf = np.bincount(c.tf[c.y], minlength=c.G); deg_tg = np.bincount(c.tg[c.y], minlength=c.G)
    for yr in Yr:
        assert yr.sum() == c.y.sum()
        assert np.array_equal(np.bincount(c.tf[yr], minlength=c.G), deg_tf)
        assert np.array_equal(np.bincount(c.tg[yr], minlength=c.G), deg_tg)
    ret = (Yr & c.y).sum(1) / c.y.sum()
    succ = np.array([(Yr[i] & Yr[i + 1]).sum() / c.y.sum() for i in range(len(Yr) - 1)])
    # mixing diagnostic: one long independent chain (burn-in 50 x n_pos)
    Ylong, stl = L.rewire_draws(c, 1, seed=L.SEED + 1000 * ci + 777, burn=50)
    ret_long = float((Ylong[0] & c.y).sum() / c.y.sum())
    # expected retention if each TF's positives were redrawn uniformly among its own candidate pairs
    exp_unif = float(sum(c.y[c.tf == t].sum() ** 2 / (c.tf == t).sum() for t in np.unique(c.tf)) / c.y.sum())
    acc_rate = np.mean([a / max(t, 1) for _, a, t in st[1:]])
    np.savez_compressed(OUT / f"{ctx}.npz", Y=np.packbits(Yr, axis=1), n=c.n)
    rows.append(dict(context=ctx, n_pairs=c.n, n_pos=int(c.y.sum()), n_draws=len(Yr),
                     burn_factor=10, thin_factor=5, burn_accepted=st[0][1], burn_tried=st[0][2], thin_accept_rate=acc_rate,
                     retained_frac_mean=ret.mean(), retained_frac_min=ret.min(), retained_frac_max=ret.max(),
                     successive_overlap_mean=succ.mean(), retained_frac_long_chain_burn50=ret_long,
                     retained_expected_uniform_within_tf=exp_unif, degrees_preserved=True))
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(L.HERE / "rewire_stats.csv", index=False)
