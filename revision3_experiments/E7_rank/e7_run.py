"""E7 runner (H4 primary, CONVENTION3 §4).

  python e7_run.py curves --ctx immune_seed43     # family T / R / N curves for one context
  python e7_run.py select                         # k* per task from immune_seed43 ONLY -> kstar.json
  python e7_run.py eval --ctx lung                # curves + Delta(k*) paired CIs + descriptive any-k test

Order enforced: `eval` refuses to run unless kstar.json exists; kstar.json is written from the selection
context's curves alone (the evaluation contexts' curves are not read by `select`).
Selection rule (fixed before any evaluation context was computed): k* = the top-k value in
{1,2,3,5,7,10,16,32,64,128,256,512} maximising the layer-mean metric in immune_seed43; ties -> smallest k.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

import e7lib as L

S = L.S
HERE = L.HERE


def log(msg, f=None):
    s = f"[{time.strftime('%H:%M:%S')}] {msg}"; print(s, flush=True)


# ------------------------------------------------------------------------------------------------- curves
def curves(ctx, D=None, svds=None):
    out = HERE / f"curves_{ctx}.csv"
    if out.exists():
        log(f"{ctx}: curves exist"); return pd.read_csv(out)
    t0 = time.time()
    D = D or L.load_ctx(ctx); emb = L.load_emb(ctx)
    svds = svds or L.layer_svds(emb, D["obs_idx"])
    chk = max(L.check_full_equals_centred(emb, D["obs_idx"], svds[l][0], l) for l in (0, 11))
    assert chk < 1e-8, chk
    reps = L.representations(ctx); rows = []
    for l in range(L.N_LAYERS):
        US = svds[l][0]
        for rid, fam, par, cols in reps:
            En = L.unit_rep(US, cols); dim = cols.stop - cols.start
            for t in L.TASKS:
                rows.append(dict(ctx=ctx, task=t, family=fam, rep=rid, param=par, dim=dim, layer=l,
                                 value=L.METRIC[t](En, D)))
        log(f"{ctx}: curves layer {l} done ({time.time()-t0:.0f}s)")
    df = pd.DataFrame(rows)
    full = df[df.rep == "T512"].set_index(["task", "layer"]).value
    df["full_value"] = [full[(t, l)] for t, l in zip(df.task, df.layer)]
    df["diff_vs_full"] = df.value - df.full_value
    # family N partners: same-width top-k
    part = {3: "T3", 6: "T6"}
    pv = df.set_index(["task", "layer", "rep"]).value
    df["partner_rep"] = [part[d] if f.startswith("N_") else "" for f, d in zip(df.family, df.dim)]
    df["partner_value"] = [pv[(t, l, p)] if p else np.nan for t, l, p in zip(df.task, df.layer, df.partner_rep)]
    df["diff_vs_partner"] = df.value - df.partner_value
    # layer-mean rows
    m = df.groupby(["ctx", "task", "family", "rep", "param", "dim"], sort=False)[
        ["value", "full_value", "diff_vs_full", "partner_value", "diff_vs_partner"]].mean().reset_index()
    m["layer"] = "mean"
    pr = df.drop_duplicates("rep").set_index("rep").partner_rep
    m["partner_rep"] = m.rep.map(pr)
    df = pd.concat([df, m[df.columns]], ignore_index=True)
    df.to_csv(out, index=False)
    (L.WORK / f"meta_{ctx}.json").write_text(json.dumps(dict(D["meta"], full_vs_centred_cos_maxabsdiff=chk,
                                                              singular_values_L0_top12=svds[0][1][:12].tolist(),
                                                              seconds=time.time() - t0), indent=1))
    log(f"{ctx}: curves written ({time.time()-t0:.0f}s)")
    return df


# ------------------------------------------------------------------------------------------------- select
def select():
    out = HERE / "kstar.json"
    assert not out.exists(), "kstar.json already written; selection is done once"
    df = pd.read_csv(HERE / f"curves_{L.SELECTION_CTX}.csv")
    df = df[(df.family == "T") & (df.layer == "mean")]
    res = dict(rule="k* = argmax over top-k values {1,2,3,5,7,10,16,32,64,128,256,512} of the layer-mean metric "
                    "in immune_seed43 (selection context); ties -> smallest k",
               selection_context=L.SELECTION_CTX, convention3=L.convention_hash(),
               written_unix=time.time(), written=time.strftime("%Y-%m-%d %H:%M:%S %Z"),
               evaluation_contexts_computed_before_selection=[c for c in L.EVAL_CTX
                                                              if (HERE / f"curves_{c}.csv").exists()],
               tasks={})
    for t in L.TASKS:
        d = df[df.task == t].copy(); d["k"] = d.param.astype(int); d = d.sort_values("k")
        best = d.value.max(); ks = d[d.value == best].k.tolist(); kstar = int(min(ks))
        res["tasks"][t] = dict(k_star=kstar, layer_mean_at_kstar=float(best),
                               layer_mean_full=float(d[d.k == 512].value.iloc[0]),
                               selection_gap=float(best - d[d.k == 512].value.iloc[0]),
                               layer_mean_curve={int(k): float(v) for k, v in zip(d.k, d.value)})
    out.write_text(json.dumps(res, indent=1))
    log(f"kstar: { {t: res['tasks'][t]['k_star'] for t in L.TASKS} }")


# ------------------------------------------------------------------------------------------------- eval
def scores_for(t, En, D):
    return L.METRIC[t](En, D, return_scores=True)


def delta(ctx, D, svds, kstar, FS, fs_sha):
    for t in L.TASKS:
        f = L.WORK / f"{ctx}_delta_{t}.npz"
        if f.exists():
            log(f"{ctx} {t}: delta checkpoint exists"); continue
        t0 = time.time(); ks = kstar[t]
        pts = np.empty((2, L.N_LAYERS)); draws = np.empty((2, L.N_LAYERS, L.B_BOOT)); ver = {}
        if t == "NODE":
            W = L.gene_boot_weights(D["n"], L.B_BOOT, L.SEED + 101)
            for l in range(L.N_LAYERS):
                for r, k in enumerate((ks, 512)):
                    En = L.unit_rep(svds[l][0], slice(0, k)); pts[r, l] = L.metric_node(En, D)
                    if l == 0 and r == 1: ver = L.verify_node_boot(En, D, W)
                    draws[r, l] = draws[0, l] if (r == 1 and ks == 512) else L.node_boot(En, D, W)
                log(f"{ctx} NODE boot layer {l} ({time.time()-t0:.0f}s)")
            ci_method = "gene bootstrap (all V genes resampled with replacement; leave-gene-out LOO scores), percentile"
        elif t == "EDGE":
            fs = FS.FastStrat(D["sp"]); wts = fs.boot_weights(L.B_BOOT, np.random.default_rng(L.SEED + 102), True)
            for l in range(L.N_LAYERS):
                for r, k in enumerate((ks, 512)):
                    En = L.unit_rep(svds[l][0], slice(0, k)); v, s = L.metric_edge(En, D, True); pts[r, l] = v
                    if l == 0 and r == 1: ver = FS.verify_against_stats3(D["sp"], s, L.SEED)
                    draws[r, l] = draws[0, l] if (r == 1 and ks == 512) else fs.auc(s, wts)
                log(f"{ctx} EDGE boot layer {l} ({time.time()-t0:.0f}s)")
            ci_method = ("two-way TF x target pigeonhole cluster bootstrap (StratifiedPairs.bootstrap draws via "
                         f"E4a FastStrat, e4a_lib.py sha256 {fs_sha[:12]}), percentile")
        else:
            sc = []
            for l in range(L.N_LAYERS):
                for k in (ks, 512):
                    v, s = L.metric_ppi(L.unit_rep(svds[l][0], slice(0, k)), D, True); sc.append(s)
            for l in range(L.N_LAYERS):
                pts[0, l] = S.auroc(sc[2 * l][D["py"]], sc[2 * l][~D["py"]])
                pts[1, l] = S.auroc(sc[2 * l + 1][D["py"]], sc[2 * l + 1][~D["py"]])
            ver = L.verify_ppi_draws(D, sc[0], sc[1], L.SEED + 103)
            dr, m = L.ppi_vertex_draws(D, sc, L.B_BOOT, L.SEED + 103)
            draws = dr.reshape(L.N_LAYERS, 2, L.B_BOOT).transpose(1, 0, 2)
            ci_method = "vertex subsampling without replacement (m = n/2, sqrt(m/n) rescaling) of the paired difference"
            log(f"{ctx} PPI vertex draws ({time.time()-t0:.0f}s)")
        np.savez(f, pts=pts, draws=draws, k_star=ks, ci_method=ci_method, verify=json.dumps(ver))
        log(f"{ctx} {t}: delta done ({time.time()-t0:.0f}s) verify={ver}")


def anyk(ctx, D, svds, FS):
    for t in L.TASKS:
        f = L.WORK / f"{ctx}_anyk_{t}.npz"
        if f.exists():
            log(f"{ctx} {t}: anyk checkpoint exists"); continue
        t0 = time.time(); nK = len(L.K_T); obs = np.empty((nK, L.N_LAYERS)); N = L.N_PERM_ANYK
        null = np.empty((N, nK, L.N_LAYERS)); ver = None
        if t == "NODE":
            Y = L.node_perm_labels(D, N, L.SEED + 201)
            for l in range(L.N_LAYERS):
                for i, k in enumerate(L.K_T):
                    En = L.unit_rep(svds[l][0], slice(0, k)); obs[i, l] = L.metric_node(En, D)
                    if l == 0 and k == 512: ver = L.verify_node_perm(En, D, Y)
                    null[:, i, l] = L.node_perm_null(En, D, Y)
                log(f"{ctx} NODE perm layer {l} ({time.time()-t0:.0f}s)")
        else:
            cols = []
            for l in range(L.N_LAYERS):
                for i, k in enumerate(L.K_T):
                    v, s = L.METRIC[t](L.unit_rep(svds[l][0], slice(0, k)), D, True); obs[i, l] = v; cols.append(s)
            Cm = np.column_stack(cols); del cols
            if t == "EDGE":
                fs = FS.FastStrat(D["sp"])
                nl, labs = fs.perm_null(Cm, N, np.random.default_rng(L.SEED + 202), return_labels_first=2)
                ref = [S.StratifiedPairs(D["sp"].tf, D["sp"].tg, yy).auroc(Cm[:, -1]) for yy in labs]
                ver = float(np.max(np.abs(np.array(ref) - nl[:2, -1]))); assert ver < 1e-9, ver
            else:
                Y = L.ppi_perm_labels(D, N, L.SEED + 203)
                nl = L.ppi_perm_null(Cm, Y)
                ref = [S.auroc(Cm[Y[:, r], -1], Cm[~Y[:, r], -1]) for r in range(2)]
                ver = float(np.max(np.abs(np.array(ref) - nl[:2, -1]))); assert ver < 1e-9, ver
            # column order: layer-major, k-minor
            null = nl.reshape(N, L.N_LAYERS, nK).transpose(0, 2, 1)
        np.savez(f, obs=obs, null=null, verify=ver)
        log(f"{ctx} {t}: anyk done ({time.time()-t0:.0f}s) verify={ver}")


def run_eval(ctx):
    kpath = HERE / "kstar.json"
    assert kpath.exists(), "run `select` first"
    kj = json.loads(kpath.read_text()); kstar = {t: kj["tasks"][t]["k_star"] for t in L.TASKS}
    FS, fs_sha = L.load_faststrat()
    D = L.load_ctx(ctx); emb = L.load_emb(ctx); svds = L.layer_svds(emb, D["obs_idx"])
    curves(ctx, D, svds)
    if ctx in L.EVAL_CTX:
        delta(ctx, D, svds, kstar, FS, fs_sha)
        anyk(ctx, D, svds, FS)
    log(f"{ctx}: eval complete")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("stage"); ap.add_argument("--ctx")
    a = ap.parse_args()
    assert L.convention_hash()["match"], "CONVENTION3 hash mismatch"
    L.check_mapping()
    if a.stage == "curves": curves(a.ctx)
    elif a.stage == "select": select()
    elif a.stage == "eval": run_eval(a.ctx)
    else: raise SystemExit("unknown stage")
