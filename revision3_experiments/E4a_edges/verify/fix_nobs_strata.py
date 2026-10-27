"""Verifier addition: TF-detection strata with the CONVENTION3 §2 definition of detection.

CONVENTION3 §2 defines detection n_obs as the number of sampled (extraction) cells whose scGPT INPUT contained the
gene (cache/<ctx>/ncells.npy). e4a_main.py's strata 'TF detected in >= / < 10% of extraction cells' use instead the
fraction of extraction cells with expression > 0 (expr_extraction.npy). The two coincide in the immune cell-type
contexts (immune_Bcell/CD4T/CD8T: identical counts) but not in lung / lung_external / kidney contexts, where the
input contains only part of the expressed genes (n_obs is a median 0.24-0.45 of the expression count). This script
appends the n_obs-defined strata to work/<ctx>/sens_strata.json (same estimator, bootstrap and seed scheme as
e4a_main.stage_sens, tie-corrected comparator). Runs only after verify/fixed/<ctx>.done (fix_ties.py).
Usage: python fix_nobs_strata.py <ctx> ...
"""
from __future__ import annotations
import json, sys, time, zlib
from pathlib import Path
import numpy as np

V = Path(__file__).resolve().parent; E4 = V.parent; R3 = E4.parent
sys.path.insert(0, str(R3)); sys.path.insert(0, str(E4))
import stats3 as S          # noqa: E402
import e4a_lib as L         # noqa: E402
import e4a_main as M        # noqa: E402

NAMES = ("TF n_obs (scGPT input) >= 10% of extraction cells [CONVENTION3 detection]",
         "TF n_obs (scGPT input) < 10% of extraction cells [CONVENTION3 detection]")


def run(ctx):
    W = E4 / "work" / ctx; done = V / "fixed" / f"{ctx}_nobs.done"
    if done.exists(): return
    while not (V / "fixed" / f"{ctx}.done").exists(): time.sleep(30)
    t0 = time.time(); seed = M.seed_of(ctx); d = M.CACHE / ctx
    h3 = dict(np.load(d / "h3_pairs.npz")); sp = S.StratifiedPairs(h3["tf"], h3["tg"], h3["y"])
    cx = dict(np.load(W / "coexp.npz"))                     # tie-corrected by fix_ties.py
    ncell = np.load(d / "expr_extraction.npy", mmap_mode="r").shape[0]
    det = np.load(d / "ncells.npy")[sp.tf] / ncell
    st = json.loads((W / "sens_strata.json").read_text()); st = [r for r in st if r.get("stratum") not in NAMES]
    for sname, m in ((NAMES[0], det >= 0.10), (NAMES[1], det < 0.10)):
        if m.sum() == 0: st.append(dict(stratum=sname, n_pos=0)); continue
        sub = S.StratifiedPairs(sp.tf[m], sp.tg[m], sp.y[m])
        if sub.tf.size == 0 or sub.y.sum() == 0: st.append(dict(stratum=sname, n_pos=0)); continue
        idx_m = np.flatnonzero(m); ii = idx_m[np.isin(sp.tf[m], sub.tfs)]; assert np.array_equal(sp.tg[ii], sub.tg)
        fs = L.FastStrat(sub); Bs = M.B if sub.tf.size < 20000 else M.B_SECONDARY
        wts = fs.boot_weights(Bs, np.random.default_rng([seed, 30, zlib.crc32(sname.encode())]))
        ce = sub.auroc(cx["cx_ext"][ii]); cb = fs.auc(cx["cx_ext"][ii], wts); c3 = sub.auroc(cx["cx_3k"][ii]); c3b = fs.auc(cx["cx_3k"][ii], wts)
        rws = [M.rowdict("|co-expression| extraction cells", -1, ce, cb)]
        for l in range(12):
            s = M.scores_layer(W, l)["prim"][ii]
            rws.append(M.rowdict("signed cosine full-512", l, sub.auroc(s), fs.auc(s, wts), ce, cb, c3, c3b))
        for r in rws: r.update(stratum=sname, n_pos=int(sub.y.sum()), n_neg=int((~sub.y).sum()), n_tf=int(sub.tfs.size),
                               delta_cx_basis="tie-corrected", added_by="verifier (CONVENTION3 detection definition)")
        st += rws; print(f"[{time.strftime('%H:%M:%S')}] {ctx} {sname} n_tf={sub.tfs.size} {time.time()-t0:.0f}s", flush=True)
    M.jdump(st, W / "sens_strata.json"); done.write_text("ok\n")


if __name__ == "__main__":
    for c in sys.argv[1:]: run(c)
