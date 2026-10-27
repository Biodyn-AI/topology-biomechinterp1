"""E8 step 1 - tokenised corpora for the co-occurrence-only instrument (Word2Vec).

For each context, every cell of the extraction population (the h5ad the extraction used; A2 found no obs
subsetting for lung, lung_external, kidney_B, immune, kidney_A) is tokenised EXACTLY like the scGPT
extraction, using A2_axis_verification/a2_lib.tokenised_genes (non-zero genes; ScGPTDataset raises if a var
gene is missing from the vocabulary, so every var gene with X > 0 is a candidate token; top max_genes = 512 by
expression with numpy's default argsort tie handling on the float32 dense row).

Validation (must hold, else abort): the per-gene token counts over the A2-reproduced extraction cells
(cache/<ctx>/extraction_cells.npy) equal cache/<ctx>/nobs_reproduced.npy on obs_idx and are zero elsewhere,
and equal layer_gene_embedding_count.npy[0] of the run.

Writes corpus/<ctx>.npz: tok (int32 var indices, concatenated), off (int64 offsets, n_cells+1), and
corpus/<ctx>_meta.json. Resumable (skips contexts already written).
"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
R3 = HERE.parent
sys.path.insert(0, str(R3)); sys.path.insert(0, str(R3 / "A2_axis_verification"))
import contexts as X   # noqa: E402
import a2_lib as L     # noqa: E402

CTXS = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
MAX_GENES = 512
OUT = HERE / "corpus"


def main():
    OUT.mkdir(exist_ok=True)
    for ctx in CTXS:
        if (OUT / f"{ctx}.npz").exists() and (OUT / f"{ctx}_meta.json").exists():
            print(ctx, "exists, skip", flush=True); continue
        t0 = time.time()
        h5 = X.CONTEXTS[ctx]["h5ad"]
        R = L.H5Rows(h5)
        R.data = R.data[:]; R.indices = R.indices[:]          # in-memory copies; row_dense code unchanged
        n, G = R.shape
        toks, off = [], [0]
        for r in range(n):
            gi = L.tokenised_genes(R.row_dense(r), MAX_GENES)
            toks.append(gi.astype(np.int32)); off.append(off[-1] + gi.size)
        R.close()
        tok = np.concatenate(toks); off = np.array(off, dtype=np.int64)
        # ---- validation against A2 and the extraction's own count vector
        d = R3 / "cache" / ctx
        cells = np.load(d / "extraction_cells.npy"); obs_idx = np.load(d / "obs_idx.npy")
        nobs_rep = np.load(d / "nobs_reproduced.npy")
        cnt = np.zeros(G, np.int64)
        for c in cells: np.add.at(cnt, tok[off[c]:off[c + 1]], 1)
        ref = np.load(X.CONTEXTS[ctx]["emb"].parent / "layer_gene_embedding_count.npy")[0].astype(np.int64)
        ok_rep = bool(np.array_equal(cnt[obs_idx], nobs_rep)) and int(np.delete(cnt, obs_idx).sum()) == 0
        ok_ref = bool(np.array_equal(cnt, ref))
        if not (ok_rep and ok_ref):
            raise SystemExit(f"{ctx}: tokenisation does not reproduce the extraction counts")
        vocab_full = np.unique(tok)
        meta = dict(context=ctx, h5ad=str(h5), n_cells=int(n), G=int(G), n_tokens=int(tok.size),
                    tokens_per_cell_median=float(np.median(np.diff(off))), tokens_per_cell_max=int(np.diff(off).max()),
                    n_cells_truncated=int((np.diff(off) == MAX_GENES).sum()),
                    n_vocab_full_population=int(vocab_full.size), n_extraction_cells=int(cells.size),
                    extraction_counts_equal_nobs_reproduced=ok_rep, extraction_counts_equal_saved_count=ok_ref,
                    V_subset_of_full_vocab=bool(np.isin(obs_idx, vocab_full).all()), max_genes=MAX_GENES,
                    seconds=round(time.time() - t0, 1))
        np.savez(OUT / f"{ctx}.npz", tok=tok, off=off)
        (OUT / f"{ctx}_meta.json").write_text(json.dumps(meta, indent=1))
        print(json.dumps(meta), flush=True)


if __name__ == "__main__":
    main()
