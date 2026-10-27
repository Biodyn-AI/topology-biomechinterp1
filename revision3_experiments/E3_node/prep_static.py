"""Extract the static scGPT token embedding and ESM2-15B (UCE per-protein; 5,120-d) protein embedding rows for every context's
observed genes V (exact symbol match). Writes partial/static_<ctx>.npz: idx_scgpt (positions in V),
emb_scgpt, idx_esm2, emb_esm2."""
import json, numpy as np
import e3_lib as L

def load(p):
    z = np.load(p, allow_pickle=False); g = [str(x) for x in z["gene_names"]]
    first = {}
    for i, x in enumerate(g): first.setdefault(x, i)
    return first, z["embeddings"]

meta = {}
for name, p in (("scgpt", L.STATIC_SCGPT), ("esm2", L.STATIC_ESM2)):
    pos, E = load(p)
    for ctx in L.CONTEXTS:
        genes = json.loads((L.CACHE / ctx / "obs_genes.json").read_text())
        iv = np.array([k for k, g in enumerate(genes) if g in pos], dtype=np.int64)
        rows = np.array([pos[genes[k]] for k in iv], dtype=np.int64)
        f = L.HERE / "partial" / f"static_{name}_{ctx}.npz"
        np.savez(f, idx=iv, emb=np.asarray(E[rows], dtype=np.float32))
        meta.setdefault(ctx, {})[name] = dict(n_V=len(genes), n_matched=int(iv.size))
        print(name, ctx, len(genes), iv.size, flush=True)
    del E
(L.HERE / "partial" / "static_meta.json").write_text(json.dumps(meta, indent=1))
