"""E4b shared definitions: regulator / target sets, direction-test pairs, cell tiers.

Specification (fixed before any GRN result was seen):
  contexts      lung, lung_external, kidney_B, immune, kidney_A
  regulators R  = TFs appearing in cache/<ctx>/h3_pairs.npz  U  Lambert TFs in V (input features)
  targets T     = every gene appearing as a target in h3_pairs (positives and candidate non-edges)
                  U the endpoints of the direction-test pairs (needed so that importance(j -> i) exists;
                  these extra targets are used ONLY by the direction test)
  tiers         T1 = cache/<ctx>/expr_extraction.npy (the exact extraction cells), seeds 0..4
                T2 = cache/<ctx>/expr.npy (3,000 cells), seeds 0..1
  GENIE3        ExtraTreesRegressor(n_estimators=300, max_features='sqrt', random_state=seed)
  GRNBoost2     GradientBoostingRegressor(learning_rate=0.01, n_estimators=1000, max_depth=3, subsample=0.9,
                max_features=0.1, n_iter_no_change=25, validation_fraction=0.1, random_state=seed)
  importance    feature_importances_ normalised to sum 1 per target; a regulator is never its own feature
                (self column = 0, never scored because no h3 pair has TF == target)
  pair score    importance(TF -> target); importances that are exactly 0 stay 0 (ties -> mid-rank)
  headline      per (method, tier): seed-ensemble score = mean of the per-seed normalised importances;
                per-seed AUROCs reported with mean and SD (run-to-run variability)
  direction     TRRUST i->j with i, j in R, (j, i) in no regulatory set (TRRUST U DoRothEA any level, both
                files = cache any_regulatory.npy); compare importance(i->j) with importance(j->i)
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
R3 = HERE.parent
CACHE = R3 / "cache"
EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"
CONTEXTS = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
SEEDS = {"T1": [0, 1, 2, 3, 4], "T2": [0, 1]}
METHODS = ["genie3", "grnboost2"]
FROZEN = R3 / "FROZEN.txt"


def convention_hash() -> str:
    for line in FROZEN.read_text().splitlines():
        if "sha256" in line:
            return line.split()[-1]
    raise RuntimeError("hash not found")


def load_sets(ctx: str) -> dict:
    d = CACHE / ctx
    genes = json.loads((d / "obs_genes.json").read_text())
    h = np.load(d / "h3_pairs.npz")
    nl = np.load(d / "node_labels.npz")
    regs = np.array(sorted(set(np.unique(h["tf"]).tolist()) | set(np.flatnonzero(nl["is_tf_lambert"]).tolist())),
                    dtype=np.int64)
    rset = set(regs.tolist())
    # direction-test pairs
    tr = np.load(d / "edges_trrust.npy")
    anyreg = np.load(d / "any_regulatory.npy")
    anyset = {(int(a), int(b)) for a, b in anyreg}
    trset = {(int(a), int(b)) for a, b in tr}
    dirp = sorted((a, b) for a, b in trset if a in rset and b in rset and a != b and (b, a) not in anyset)
    dirp = np.array(dirp, dtype=np.int64).reshape(-1, 2)
    h3_targets = np.unique(h["tg"])
    targets = np.unique(np.concatenate([h3_targets, dirp.ravel()])).astype(np.int64)
    return dict(genes=genes, tf=h["tf"].astype(np.int64), tg=h["tg"].astype(np.int64), y=h["y"].astype(bool),
                regs=regs, targets=targets, h3_targets=h3_targets, dir_pairs=dirp)


def load_expr(ctx: str, tier: str) -> np.ndarray:
    f = "expr_extraction.npy" if tier == "T1" else "expr.npy"
    return np.asarray(np.load(CACHE / ctx / f), dtype=np.float32)


def importance_lookup(imp: np.ndarray, targets: np.ndarray, regs: np.ndarray, src: np.ndarray, dst: np.ndarray):
    """importance(src -> dst) from an (n_targets x n_regs) matrix; 0 if dst not a target or src not a regulator."""
    tpos = np.full(int(max(targets.max(), dst.max())) + 1, -1); tpos[targets] = np.arange(targets.size)
    rpos = np.full(int(max(regs.max(), src.max())) + 1, -1); rpos[regs] = np.arange(regs.size)
    ti = tpos[dst]; ri = rpos[src]; ok = (ti >= 0) & (ri >= 0)
    out = np.zeros(src.size)
    out[ok] = imp[ti[ok], ri[ok]]
    return out, ok
