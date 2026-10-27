"""Verifier step 1 (independent of e8_*.py and a2_lib): corpus tokens, embedding copies, gensim downsampling.

(a) Re-tokenise cells straight from the h5ad with anndata (not h5py / a2_lib): extraction cells + 300 random
    other cells per context; compare the token SET per cell with corpus/<ctx>.npz.
(b) Per-gene extraction-cell counts from my tokens vs layer_gene_embedding_count.npy[0] of the run.
(c) Local embedding copy vs source file: byte size and SHA-256 of layers 0 and 11 for V rows.
(d) gensim sample=1e-3 frequent-token downsampling: keep prob = (sqrt(f/s)+1) s/f per gene; how many V genes
    have keep prob < 1 in each corpus (full / extraction).
Writes verify/v1_inputs.json
"""
import json, sys, hashlib
from pathlib import Path
import numpy as np, anndata as ad
HERE = Path(__file__).resolve().parent; E8 = HERE.parent; R3 = E8.parent
sys.path.insert(0, str(R3)); import contexts as X
CTXS = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
out = {}
for ctx in CTXS:
    A = ad.read_h5ad(X.CONTEXTS[ctx]["h5ad"], backed="r")
    z = np.load(E8 / "corpus" / f"{ctx}.npz"); tok, off = z["tok"], z["off"]
    d = R3 / "cache" / ctx; ext = np.load(d / "extraction_cells.npy"); obs_idx = np.load(d / "obs_idx.npy")
    rng = np.random.default_rng(777)
    others = rng.choice(np.setdiff1d(np.arange(A.n_obs), ext), 300, replace=False)
    cells = np.concatenate([ext, others]); mism = 0; cnt = np.zeros(A.n_vars, np.int64); ties = 0
    for i, c in enumerate(cells):
        row = A.X[int(c)]
        row = np.asarray(row.toarray()).ravel() if hasattr(row, "toarray") else np.asarray(row).ravel()
        nz = np.flatnonzero(row > 0)
        if nz.size > 512:
            v = np.sort(row[nz])[::-1]; ties += int(v[511] == v[512])
            mine = nz[np.argsort(-row[nz])[:512]]
        else:
            mine = nz
        if set(mine.tolist()) != set(tok[off[c]:off[c + 1]].tolist()): mism += 1
        if i < ext.size: np.add.at(cnt, mine, 1)
    ref = np.load(X.CONTEXTS[ctx]["emb"].parent / "layer_gene_embedding_count.npy")[0]
    # (c) embedding copy
    src = np.load(X.CONTEXTS[ctx]["emb"], mmap_mode="r"); loc = np.load(Path.home() / "biomechinterp-local/subproject41_r3/emb" / f"{ctx}.npy", mmap_mode="r")
    hs = [hashlib.sha256(np.ascontiguousarray(e[l][obs_idx]).tobytes()).hexdigest() for e in (src, loc) for l in (0, 11)]
    # (d) downsampling
    L = np.diff(off); ds = {}
    for name, cs in (("full", np.arange(L.size)), ("extr", ext)):
        T = np.concatenate([tok[off[c]:off[c + 1]] for c in cs]); f = np.bincount(T, minlength=A.n_vars)[obs_idx] / T.size
        s = 1e-3; keep = np.minimum(1.0, (np.sqrt(f / s) + 1) * s / np.maximum(f, 1e-300))
        ds[name] = dict(n_tokens=int(T.size), n_V_keep_lt_1=int((keep < 1).sum()), min_keep=float(keep.min()),
                        frac_V_tokens_dropped=float(((1 - keep) * f).sum() / f.sum()))
    out[ctx] = dict(n_cells_checked=int(cells.size), n_extraction=int(ext.size), token_set_mismatches=mism, ties_at_512=ties,
                    extraction_counts_equal_saved=bool(np.array_equal(cnt, ref)),
                    emb_shape_equal=src.shape == loc.shape, emb_V_layers_0_11_sha_equal=(hs[0] == hs[2] and hs[1] == hs[3]),
                    downsampling=ds)
    print(ctx, json.dumps(out[ctx]), flush=True)
(HERE / "v1_inputs.json").write_text(json.dumps(out, indent=1))
