"""Build the shared benchmark sets (CONVENTION3 §3-4) once, deterministically, for every context.

Writes cache/<context>/:
  node_labels.npz   is_tf_lambert (bool per observed gene), is_trrust_tf, is_trrust_target_only,
                    log_nobs (log of detection), det_decile
  h3_pairs.npz      primary H3 benchmark: tf, tg (indices into obs_genes), y (1 TRRUST edge / 0 candidate
                    non-edge), match_level (0 = matched on in-degree AND detection decile, 1 = in-degree only,
                    2 = unmatched fallback), mode (TRRUST mode of positives: Activation/Repression/Unknown),
                    n_pmid (positives), for TFs with >= 3 positives and >= 3 candidate non-edges
  h3_legacy_pairs.npz  round-2 style: same TF, any observed gene that is not a documented target of it
  h3_chip_pairs.npz    independent positives: DoRothEA ChIP-file level-D pairs not in TRRUST and not in the
                       main DoRothEA A-C; candidate non-edges built with the same rules as h3_pairs
  pairs_meta.json   counts
Candidate non-edges are "TRRUST-negative candidate pairs": undocumented, not confirmed absent.
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
import contexts as X
import stats3 as S

C = Path(__file__).parent / "cache"
SEED = 20261026
NEG_PER_POS = 10
LAMBERT = X.REPO / "subproject_61_two_axes_revision/data/humantfs_v1.01.txt"


def load_refs():
    t, _ = X.load_trrust()
    t = t[t.TF != t.target]
    agg = t.groupby(["TF", "target"]).agg(n_pmid=("pmid", "nunique"),
                                           mode=("mode", lambda m: "/".join(sorted(set(map(str, m))))))
    trrust = {k: dict(n_pmid=int(v.n_pmid), mode=v["mode"]) for k, v in agg.iterrows()}
    main = pd.read_csv(X.NET / "dorothea_human.tsv", sep="\t")
    chip = pd.read_csv(X.NET / "dorothea_chipseq_human.tsv", sep="\t")
    anyreg = set(trrust) | set(zip(main.source.astype(str), main.target.astype(str))) \
             | set(zip(chip.source.astype(str), chip.target.astype(str)))
    main_abc = {(a, b) for a, b, c in zip(main.source.astype(str), main.target.astype(str), main.confidence.astype(str))
                if any(l in c.split(";") for l in "ABC")}
    chip_d = {(a, b) for a, b, c in zip(chip.source.astype(str), chip.target.astype(str), chip.confidence.astype(str))
              if c.split(";") == ["D"]}
    chip_indep = {e for e in chip_d if e not in trrust and e not in main_abc and e[0] != e[1]}
    lam = pd.read_csv(LAMBERT, sep="\t")
    lambert = set(lam.loc[lam["Is TF?"].astype(str).str.strip() == "Yes", "HGNC symbol"].astype(str))
    return trrust, anyreg, chip_indep, lambert


def build(pos_pairs, genes, pos, anyreg, indeg, det_dec, rng, other_targets_by_tf):
    """Target-matched candidate non-edges for a list of positive (tf, tg) symbol pairs."""
    indeg_dec = S.quantile_strata(np.log1p(np.array([indeg.get(g, 0) for g in genes], float)), 10)
    rows = []
    by_tf = {}
    for a, b in pos_pairs: by_tf.setdefault(a, []).append(b)
    for tf, tgts in sorted(by_tf.items()):
        pool = [g for g in other_targets_by_tf.get(tf, ()) if g != tf and (tf, g) not in anyreg and (g, tf) not in anyreg]
        if len(tgts) < 3 or len(pool) < 3: continue
        pool = np.array(sorted(pool)); pi = np.array([pos[g] for g in pool])
        used = np.zeros(pool.size, bool); negs = []
        for t in tgts:
            ti = pos[t]
            for level, mask in ((0, (indeg_dec[pi] == indeg_dec[ti]) & (det_dec[pi] == det_dec[ti])),
                                (1, indeg_dec[pi] == indeg_dec[ti]), (2, np.ones(pool.size, bool))):
                cand = np.flatnonzero(mask & ~used)
                if cand.size >= 1: break
            take = rng.choice(cand, min(NEG_PER_POS, cand.size), replace=False) if cand.size else []
            used[take] = True
            negs += [(pool[k], level) for k in take]
        if len(negs) < 3: continue
        for t in tgts: rows.append((pos[tf], pos[t], 1, 0))
        for g, level in negs: rows.append((pos[tf], pos[g], 0, level))
    if not rows: return None
    a = np.array(rows, dtype=np.int64)
    return a[:, 0], a[:, 1], a[:, 2].astype(np.int8), a[:, 3].astype(np.int8)


def main():
    trrust, anyreg, chip_indep, lambert = load_refs()
    meta = {}
    for ctx in X.CONTEXTS:
        d = C / ctx; rng = np.random.default_rng(SEED)
        genes = json.loads((d / "obs_genes.json").read_text()); pos = {g: k for k, g in enumerate(genes)}
        V = set(genes); nobs = np.load(d / "ncells.npy").astype(float)
        log_nobs = np.log(np.maximum(nobs, 1)); det_dec = S.quantile_strata(log_nobs, 10)
        tr_v = sorted(p for p in trrust if p[0] in V and p[1] in V)
        tfs_tr = {a for a, _ in tr_v}; targets_tr = {b for _, b in tr_v}
        indeg = {}
        for a, b in tr_v: indeg[b] = indeg.get(b, 0) + 1
        other = {tf: {b for a, b in tr_v if a != tf} for tf in tfs_tr}
        np.savez(d / "node_labels.npz",
                 is_tf_lambert=np.array([g in lambert for g in genes]),
                 is_trrust_tf=np.array([g in tfs_tr for g in genes]),
                 is_trrust_target_only=np.array([(g in targets_tr) and (g not in tfs_tr) for g in genes]),
                 log_nobs=log_nobs, det_decile=det_dec)
        res = build(tr_v, genes, pos, anyreg, indeg, det_dec, rng, other)
        tf, tg, y, lev = res
        mode = np.array([trrust[(genes[a], genes[b])]["mode"] if yy else "" for a, b, yy in zip(tf, tg, y)])
        npm = np.array([trrust[(genes[a], genes[b])]["n_pmid"] if yy else 0 for a, b, yy in zip(tf, tg, y)])
        np.savez(d / "h3_pairs.npz", tf=tf, tg=tg, y=y, match_level=lev, mode=mode, n_pmid=npm)
        # legacy: any observed gene not documented for this TF (round-2 construction, all-V)
        leg = []
        for t in sorted(tfs_tr):
            tg_t = [b for a, b in tr_v if a == t]
            if len(tg_t) < 3: continue
            cand = np.array([pos[g] for g in genes if g != t and (t, g) not in anyreg and (g, t) not in anyreg])
            take = rng.choice(cand, min(NEG_PER_POS * len(tg_t), cand.size), replace=False)
            leg += [(pos[t], pos[b], 1) for b in tg_t] + [(pos[t], int(k), 0) for k in take]
        leg = np.array(leg, dtype=np.int64)
        np.savez(d / "h3_legacy_pairs.npz", tf=leg[:, 0], tg=leg[:, 1], y=leg[:, 2].astype(np.int8))
        # independent ChIP-only positives
        ch_v = sorted(p for p in chip_indep if p[0] in V and p[1] in V)
        indeg_c = {}
        for a, b in ch_v: indeg_c[b] = indeg_c.get(b, 0) + 1
        other_c = {tf: {b for a, b in ch_v if a != tf} | {b for a, b in tr_v if a != tf} for tf in {a for a, _ in ch_v}}
        rc = build(ch_v, genes, pos, anyreg - set(ch_v), indeg_c, det_dec, rng, other_c)
        if rc is not None:
            np.savez(d / "h3_chip_pairs.npz", tf=rc[0], tg=rc[1], y=rc[2], match_level=rc[3])
        m = dict(context=ctx, n_V=len(genes), n_lambert_tf=int(sum(g in lambert for g in genes)),
                 trrust_within_V=len(tr_v), h3_tfs=int(np.unique(tf).size), h3_pos=int(y.sum()),
                 h3_neg=int((1 - y).sum()), h3_neg_fully_matched=float(np.mean(lev[y == 0] == 0)),
                 legacy_pos=int(leg[:, 2].sum()), legacy_neg=int((1 - leg[:, 2]).sum()),
                 chip_pos=int(rc[2].sum()) if rc is not None else 0,
                 chip_tfs=int(np.unique(rc[0]).size) if rc is not None else 0)
        meta[ctx] = m; print(json.dumps(m), flush=True)
    (C / "pairs_meta.json").write_text(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
