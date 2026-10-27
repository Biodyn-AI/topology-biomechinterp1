"""E8 step 2 - train the co-occurrence-only instrument (gensim Word2Vec) on the step-1 corpora.

Task specification: sg=1, vector_size=256, window=512, min_count=1, epochs=10, negative=10, workers=2, 5 seeds;
token order shuffled per cell per epoch. Every unspecified parameter is the gensim 4.4.0 default (recorded in
each model's JSON).

COMPUTE DEVIATION (documented in SUMMARY.txt). With window=512 every token is trained against ~170 partners
(x 11 dot products of length 256). Measured on this machine at load ~200: 189 s per epoch for 256 lung cells,
i.e. ~14 days for the specified 5 seeds x 10 epochs x 26 M population tokens. Because token order is shuffled
independently per cell per pass, a window of w samples, per epoch, a uniformly random subset of each token's
within-cell co-occurrence partners; the expected SGNS objective is the same within-cell co-occurrence objective
(up to a per-cell weight ~ 1/L, L = tokens in the cell, that only matters where L varies). Production models
therefore use window=5 (the gensim default) with per-cell-per-pass shuffling; the specified window=512 is run
on the extraction-cell corpus as a calibration (W512 vs W5 with the number of pair updates matched,
`epochs_matched`), and the comparison is reported in SUMMARY.md.

Corpora: full = every cell of the extraction population; extr = the A2-reproduced extraction cells only
(lower bound). Cell order is also shuffled per pass (h5ad rows are grouped by donor / cell type; a fixed
order would couple the linear learning-rate decay to cell identity).
Seeds: gensim seed = 1000 + s; stream rng = default_rng([20261026, s, pass]). gensim with workers > 1 is not
bit-reproducible (thread scheduling); seeds fix initialisation and the corpus stream only.

Writes models/<ctx>__<corpus>__w<window>__e<epochs>__s<seed>.npy (n_V x 256 float32, rows in
cache/<ctx>/obs_genes.json order; NaN rows for V genes absent from the Word2Vec vocabulary) and .json.
Resumable. Usage: python e8_train.py <process-id>   (job lists below)
"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import numpy as np
import gensim
from gensim.models import Word2Vec

HERE = Path(__file__).resolve().parent
R3 = HERE.parent
CTXS = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
SEEDS = [0, 1, 2, 3, 4]
BASE = dict(sg=1, vector_size=256, min_count=1, negative=10, workers=2)
W_PROD, E_PROD, W_SPEC = 5, 10, 512
OUT = HERE / "models"


def pairs_per_token(lengths, window):
    """Expected number of (centre, context) updates per token under gensim's shrink_windows (effective
    window b ~ U{1..window}) for sentences of the given lengths (order-free: depends on L only)."""
    tot = 0.0; ntok = 0
    b = np.arange(1, window + 1)
    for L in np.unique(lengths):
        n_L = int((lengths == L).sum()); i = np.arange(L)
        cnt = np.minimum(b[None, :], i[:, None]) + np.minimum(b[None, :], (L - 1 - i)[:, None])
        tot += n_L * cnt.mean(1).sum(); ntok += n_L * L
    return tot / ntok


def jobs_for(pid):
    J = []
    if pid == 1:
        for s in SEEDS: J += [("lung", "full", W_PROD, E_PROD, s)]
        for s in SEEDS: J += [(c, "extr", W_PROD, E_PROD, s) for c in CTXS]
    elif pid == 2:
        for s in SEEDS: J += [("lung_external", "full", W_PROD, E_PROD, s), ("kidney_B", "full", W_PROD, E_PROD, s)]
    elif pid == 3:
        for s in SEEDS: J += [("immune", "full", W_PROD, E_PROD, s), ("kidney_A", "full", W_PROD, E_PROD, s)]
    elif pid == 4:      # calibration: the specified window on the extraction corpus ...
        J += [(c, "extr", W_SPEC, E_PROD, 0) for c in ["lung", "immune"]]
    elif pid == 5:      # ... and window 5 with the number of pair updates matched to it
        J += [(c, "extr", W_PROD, "matched", 0) for c in ["lung", "immune"]]
    return J


class ShuffledCorpus:
    """Re-iterable corpus; each pass (vocab build = pass 0, then one per epoch) shuffles the cell order and
    every cell's tokens."""
    def __init__(self, tok, off, cells, names, seed):
        self.tok, self.off, self.cells, self.names, self.seed = tok, off, cells, names, seed
        self.passes = 0

    def __iter__(self):
        rng = np.random.default_rng([20261026, self.seed, self.passes]); self.passes += 1
        for c in rng.permutation(self.cells):
            t = self.tok[self.off[c]:self.off[c + 1]]
            yield self.names[rng.permutation(t)].tolist()


def run(job, data):
    ctx, corpus, window, epochs, s = job
    tok, off, obs_idx, ext, names = data[ctx]
    cells = np.sort(ext) if corpus == "extr" else np.arange(off.size - 1)
    lengths = np.diff(off)[cells]
    ppt_w = pairs_per_token(lengths, window)
    if epochs == "matched":
        epochs = int(round(E_PROD * pairs_per_token(lengths, W_SPEC) / ppt_w))
    stem = OUT / f"{ctx}__{corpus}__w{window}__e{epochs}__s{s}"
    if stem.with_suffix(".npy").exists() and stem.with_suffix(".json").exists():
        return
    t0 = time.time()
    cor = ShuffledCorpus(tok, off, cells, names, s)
    m = Word2Vec(seed=1000 + s, window=window, epochs=epochs, **BASE)
    m.build_vocab(cor)
    m.train(cor, total_examples=m.corpus_count, epochs=m.epochs)
    W = np.full((obs_idx.size, BASE["vector_size"]), np.nan, np.float32); inv = np.zeros(obs_idx.size, bool)
    for r, gi in enumerate(obs_idx):
        k = f"g{gi}"
        if k in m.wv.key_to_index:
            W[r] = m.wv[k]; inv[r] = True
    np.save(stem.with_suffix(".npy"), W)
    info = dict(context=ctx, corpus=corpus, seed=s, gensim_seed=1000 + s, n_cells=int(cells.size),
                n_tokens=int(lengths.sum()), vocab=len(m.wv), n_V=int(obs_idx.size), n_V_in_vocab=int(inv.sum()),
                passes=cor.passes, pairs_per_token_per_epoch=ppt_w,
                pair_updates_total=float(ppt_w * lengths.sum() * epochs), finite=bool(np.isfinite(W[inv]).all()),
                params={k: getattr(m, k) for k in ("sg", "vector_size", "window", "min_count", "epochs", "negative",
                                                     "workers", "alpha", "min_alpha", "sample", "hs", "ns_exponent",
                                                     "cbow_mean", "shrink_windows", "max_vocab_size")},
                gensim=gensim.__version__, seconds=round(time.time() - t0, 1))
    stem.with_suffix(".json").write_text(json.dumps(info, indent=1, default=str))
    print(time.strftime("%H:%M:%S"), json.dumps({k: info[k] for k in ("context", "corpus", "seed", "n_cells", "vocab",
          "n_V_in_vocab", "passes", "pair_updates_total", "seconds")}), "window", window, "epochs", epochs, flush=True)


def main(pid):
    OUT.mkdir(exist_ok=True)
    J = jobs_for(pid); data = {}
    for c in sorted({j[0] for j in J}):
        z = np.load(HERE / "corpus" / f"{c}.npz"); d = R3 / "cache" / c
        obs_idx = np.load(d / "obs_idx.npy"); ext = np.load(d / "extraction_cells.npy")
        G = max(int(z["tok"].max()), int(obs_idx.max())) + 1
        data[c] = (z["tok"], z["off"], obs_idx, ext, np.array([f"g{i}" for i in range(G)], dtype=object))
    print("process", pid, "jobs", len(J), flush=True)
    for j in J:
        run(j, data)
    print("process", pid, "DONE", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]))
