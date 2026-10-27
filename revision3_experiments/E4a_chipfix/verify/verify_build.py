"""Independent verifier, part 1: the corrected ChIP-only benchmark (all 16 contexts).

Independent of e4a_chipfix.py and of prep_pairs.build / prep_pairs.load_refs: the reference sets are parsed here
from the raw TSVs, and the candidate builder is re-written here. Only stats3.quantile_strata (the decile helper
whose output defines the match strata) and the context list are imported.

Checks per context:
  R1  my own builder, replaying the RNG stream of prep_pairs.main() (h3 build -> legacy draws -> ChIP draw),
      reproduces ../../cache/<ctx>/h3_pairs.npz, h3_legacy_pairs.npz and the OLD h3_chip_pairs.npz (old rule
      anyreg - set(ch_v)), and with the full anyreg reproduces cache/<ctx>/h3_chip_pairs_fixed.npz bit for bit.
  R2  fixed candidates: none in anyreg (either direction), none equal to any ChIP positive of ANY TF (either
      direction), none equal to a TRRUST / DoRothEA pair, no self pairs, no duplicate (tf, tg) rows, no candidate
      repeated across positive and negative label.
  R3  positives: fixed y==1 rows == old y==1 rows (same order); == the ChIP-independent pairs within V for the
      TFs kept; TFs dropped (none expected).
  R4  match levels: shares recomputed; consistency: every level-0 candidate matches (indeg decile, det decile) of
      at least one positive of its TF, every level-1 candidate matches an indeg decile of one positive; a level-2
      candidate is assigned only when its positive had no level-0/1 option (checked by replay = R1).
  R5  the contamination audit numbers in build_audit.csv (old_neg_is_own_chip_positive etc.) recomputed.
"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent; FIX = HERE.parent; R3 = FIX.parent
sys.path.insert(0, str(R3))
import contexts as X                       # noqa: E402  (context list, paths)
from stats3 import quantile_strata         # noqa: E402

SEED = 20261026; NPP = 10
NET = X.NET


def refs():
    t = pd.read_csv(NET / "trrust_human.tsv", sep="\t", header=None, names=["TF", "target", "mode", "pmid"])
    t = t[t.TF.astype(str) != t.target.astype(str)]
    trrust = set(zip(t.TF.astype(str), t.target.astype(str)))
    mn = pd.read_csv(NET / "dorothea_human.tsv", sep="\t"); ch = pd.read_csv(NET / "dorothea_chipseq_human.tsv", sep="\t")
    main_all = set(zip(mn.source.astype(str), mn.target.astype(str)))
    chip_all = set(zip(ch.source.astype(str), ch.target.astype(str)))
    anyreg = trrust | main_all | chip_all
    main_abc = {(a, b) for a, b, c in zip(mn.source.astype(str), mn.target.astype(str), mn.confidence.astype(str))
                if set(c.split(";")) & {"A", "B", "C"}}
    chip_d = {(a, b) for a, b, c in zip(ch.source.astype(str), ch.target.astype(str), ch.confidence.astype(str))
              if c.split(";") == ["D"]}
    indep = {p for p in chip_d if p not in trrust and p not in main_abc and p[0] != p[1]}
    return trrust, main_all, chip_all, anyreg, indep


def my_build(pos_pairs, genes, idx, excl, indeg, det_dec, rng, others):
    ind_dec = quantile_strata(np.log1p(np.array([indeg.get(g, 0) for g in genes], float)), 10)
    by = {}
    for a, b in pos_pairs: by.setdefault(a, []).append(b)
    out = []
    for tf in sorted(by):
        tg = by[tf]
        pool = sorted(g for g in others.get(tf, ()) if g != tf and (tf, g) not in excl and (g, tf) not in excl)
        if len(tg) < 3 or len(pool) < 3: continue
        pool = np.array(pool); pi = np.array([idx[g] for g in pool]); used = np.zeros(pool.size, bool); neg = []
        for t in tg:
            ti = idx[t]
            m0 = (ind_dec[pi] == ind_dec[ti]) & (det_dec[pi] == det_dec[ti]) & ~used
            m1 = (ind_dec[pi] == ind_dec[ti]) & ~used
            m2 = ~used
            for lev, m in ((0, m0), (1, m1), (2, m2)):
                c = np.flatnonzero(m)
                if c.size: break
            if c.size:
                k = rng.choice(c, min(NPP, c.size), replace=False); used[k] = True
                neg += [(idx[pool[j]], lev) for j in k]
        if len(neg) < 3: continue
        out += [(idx[tf], idx[t], 1, 0) for t in tg] + [(idx[tf], g, 0, lev) for g, lev in neg]
    if not out: return None
    a = np.array(out, np.int64)
    return dict(tf=a[:, 0], tg=a[:, 1], y=a[:, 2].astype(np.int8), match_level=a[:, 3].astype(np.int8)), ind_dec


def eq(a, b, keys): return all(np.array_equal(np.asarray(a[k]), np.asarray(b[k])) for k in keys)


def run(ctx, R):
    trrust, main_all, chip_all, anyreg, indep = R
    d = R3 / "cache" / ctx
    genes = json.loads((d / "obs_genes.json").read_text()); idx = {g: i for i, g in enumerate(genes)}; V = set(genes)
    det_dec = quantile_strata(np.log(np.maximum(np.load(d / "ncells.npy").astype(float), 1)), 10)
    rng = np.random.default_rng(SEED)
    tr = sorted(p for p in trrust if p[0] in V and p[1] in V)
    indeg = {}
    for _, b in tr: indeg[b] = indeg.get(b, 0) + 1
    oth = {t: {b for a, b in tr if a != t} for t in {a for a, _ in tr}}
    h3, _ = my_build(tr, genes, idx, anyreg, indeg, det_dec, rng, oth)
    r1_h3 = eq(h3, dict(np.load(d / "h3_pairs.npz")), ["tf", "tg", "y", "match_level"])
    leg = []
    for t in sorted({a for a, _ in tr}):
        tg = [b for a, b in tr if a == t]
        if len(tg) < 3: continue
        cand = np.array([idx[g] for g in genes if g != t and (t, g) not in anyreg and (g, t) not in anyreg])
        k = rng.choice(cand, min(NPP * len(tg), cand.size), replace=False)
        leg += [(idx[t], idx[b]) for b in tg] + [(idx[t], int(j)) for j in k]
    leg = np.array(leg); ol = np.load(d / "h3_legacy_pairs.npz")
    r1_leg = np.array_equal(leg[:, 0], ol["tf"]) and np.array_equal(leg[:, 1], ol["tg"])
    cv = sorted(p for p in indep if p[0] in V and p[1] in V)
    indeg_c = {}
    for _, b in cv: indeg_c[b] = indeg_c.get(b, 0) + 1
    oth_c = {t: {b for a, b in cv if a != t} | {b for a, b in tr if a != t} for t in {a for a, _ in cv}}
    st = rng.bit_generator.state
    rng_new = np.random.default_rng(); rng_new.bit_generator.state = st
    old_mine, _ = my_build(cv, genes, idx, anyreg - set(cv), indeg_c, det_dec, rng, oth_c)
    new_mine, ind_dec_c = my_build(cv, genes, idx, anyreg, indeg_c, det_dec, rng_new, oth_c)
    old = dict(np.load(d / "h3_chip_pairs.npz")); fx = dict(np.load(FIX / "cache" / ctx / "h3_chip_pairs_fixed.npz"))
    K = ["tf", "tg", "y", "match_level"]
    r1_old = eq(old_mine, old, K); r1_new = eq(new_mine, fx, K)
    # ---- R2 contamination of fixed candidates
    G = np.array(genes, dtype=object)
    neg = fx["y"] == 0; posm = fx["y"] == 1
    cand = list(zip(G[fx["tf"][neg]], G[fx["tg"][neg]]))
    cset = set(cand)
    posset_all = set(cv)                                        # every ChIP-independent pair within V (any TF)
    bench_pos = set(zip(G[fx["tf"][posm]], G[fx["tg"][posm]]))
    n_anyreg = sum((p in anyreg) or ((p[1], p[0]) in anyreg) for p in cand)
    n_pos_any = sum((p in posset_all) or ((p[1], p[0]) in posset_all) for p in cand)
    n_bench_pos = sum((p in bench_pos) or ((p[1], p[0]) in bench_pos) for p in cand)
    n_tr = sum((p in trrust) or ((p[1], p[0]) in trrust) for p in cand)
    n_main = sum((p in main_all) or ((p[1], p[0]) in main_all) for p in cand)
    n_chipall = sum((p in chip_all) or ((p[1], p[0]) in chip_all) for p in cand)
    n_self = int(np.sum(fx["tf"] == fx["tg"]))
    n_dup = len(cand) - len(cset)
    pair_ids = fx["tf"].astype(np.int64) * 100000 + fx["tg"]
    n_dup_all = pair_ids.size - np.unique(pair_ids).size
    # ---- R3 positives
    po = old["y"] == 1
    pos_same = np.array_equal(old["tf"][po], fx["tf"][posm]) and np.array_equal(old["tg"][po], fx["tg"][posm])
    kept_tfs = set(G[np.unique(fx["tf"])])
    exp_pos = sorted(p for p in cv if p[0] in kept_tfs)
    pos_eq_ref = sorted(bench_pos) == exp_pos and len(bench_pos) == int(posm.sum())
    old_tfs = set(G[np.unique(old["tf"])])
    # ---- R4 match levels and consistency
    lev = fx["match_level"][neg]; shares = [float(np.mean(lev == k)) for k in (0, 1, 2)]
    oshares = [float(np.mean(old["match_level"][old["y"] == 0] == k)) for k in (0, 1, 2)]
    bad0 = bad1 = 0
    pos_keys = {}
    for a, b in zip(fx["tf"][posm], fx["tg"][posm]):
        pos_keys.setdefault(a, [set(), set()]); pos_keys[a][0].add((ind_dec_c[b], det_dec[b])); pos_keys[a][1].add(ind_dec_c[b])
    for a, b, l in zip(fx["tf"][neg], fx["tg"][neg], lev):
        if l == 0 and (ind_dec_c[b], det_dec[b]) not in pos_keys[a][0]: bad0 += 1
        if l == 1 and ind_dec_c[b] not in pos_keys[a][1]: bad1 += 1
    # ---- R5 old contamination audit
    oneg = old["y"] == 0
    ocand = list(zip(G[old["tf"][oneg]], G[old["tg"][oneg]]))
    cvs = set(cv); o_own = sum(p in cvs for p in ocand)
    o_any = sum((p in anyreg) or ((p[1], p[0]) in anyreg) for p in ocand)
    so = set(ocand)
    row = dict(context=ctx, R1_h3_identical=r1_h3, R1_legacy_identical=r1_leg, R1_old_chip_identical=r1_old,
               R1_fixed_identical_to_my_rebuild=r1_new,
               n_pos=int(posm.sum()), n_neg=int(neg.sum()), n_tf=len(kept_tfs), old_n_tf=len(old_tfs),
               R2_cand_in_anyreg_either_dir=n_anyreg, R2_cand_is_any_chip_indep_pair_either_dir=n_pos_any,
               R2_cand_is_benchmark_positive_either_dir=n_bench_pos, R2_cand_in_trrust=n_tr, R2_cand_in_dorothea_main=n_main,
               R2_cand_in_dorothea_chip=n_chipall, R2_self_pairs=n_self, R2_duplicate_candidates=n_dup, R2_duplicate_rows_any_label=int(n_dup_all),
               R3_positives_identical_to_old=pos_same, R3_positives_equal_reference_for_kept_tfs=pos_eq_ref,
               R3_tfs_dropped=len(old_tfs - kept_tfs),
               R4_share0=shares[0], R4_share1=shares[1], R4_share2=shares[2],
               R4_old_share0=oshares[0], R4_old_share1=oshares[1], R4_old_share2=oshares[2],
               R4_level0_inconsistent=bad0, R4_level1_inconsistent=bad1,
               R5_old_neg_is_own_chip_positive=o_own, R5_old_neg_in_anyreg_either_dir=o_any,
               R5_neg_common=len(so & cset), R5_frac_old_retained=len(so & cset) / len(so))
    return row


if __name__ == "__main__":
    ctxs = sys.argv[1:] or list(X.CONTEXTS)
    R = refs(); rows = []
    for c in ctxs:
        t0 = time.time(); r = run(c, R); rows.append(r)
        print(f"[{time.strftime('%H:%M:%S')}] {c} {time.time()-t0:.0f}s", {k: v for k, v in r.items() if k.startswith(("R1", "R2_cand_in_anyreg", "R3_pos"))}, flush=True)
    df = pd.DataFrame(rows); df.to_csv(HERE / "verify_build.csv", index=False)
    aud = pd.read_csv(FIX / "build_audit.csv").set_index("context")
    cmp = []
    for r in rows:
        a = aud.loc[r["context"]]
        cmp.append(dict(context=r["context"],
                        own_pos_match=int(a.old_neg_is_own_chip_positive) == r["R5_old_neg_is_own_chip_positive"],
                        anyreg_match=int(a.old_neg_in_anyreg_either_direction) == r["R5_old_neg_in_anyreg_either_dir"],
                        new_anyreg_zero=int(a.new_neg_in_anyreg_either_direction) == 0 == r["R2_cand_in_anyreg_either_dir"],
                        common_match=int(a.neg_common) == r["R5_neg_common"],
                        share_match=max(abs(a[f"new_share_level{k}"] - r[f"R4_share{k}"]) for k in (0, 1, 2)) < 1e-12,
                        nneg_match=int(a.new_n_neg) == r["n_neg"]))
    cdf = pd.DataFrame(cmp); cdf.to_csv(HERE / "verify_build_vs_audit.csv", index=False)
    print(cdf.to_string()); print(df.T.to_string())
