"""Shared loaders for E11 (compartment re-test). Read-only use of contexts.py, stats3.py and cache/."""
from __future__ import annotations
import hashlib, json, pickle, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
R3 = HERE.parent
sys.path.insert(0, str(R3))
import contexts as X   # noqa: E402
import stats3 as S     # noqa: E402

LOCAL_EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"
GENE2GO = X.REPO / "single_cell_mechinterp/data/perturb/gene2go_all.pkl"
GO_OBO = X.REPO / "subproject_61_two_axes_revision/data/go-basic.obo"
CTXS = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
TERMS = {"GO:0005615": "extracellular_space", "GO:0005829": "cytosol",
         "GO:0005739": "mitochondrion", "GO:0005788": "ER_lumen"}
PRIMARY = ["GO:0005615", "GO:0005829"]
SECONDARY = ["GO:0005739", "GO:0005788"]
POLE_FRAC = 0.27
N_LAYERS, N_DIRS = 12, 12


def convention_hash():
    txt = (R3 / "FROZEN.txt").read_text().split()
    frozen = txt[txt.index("sha256") + 1]
    actual = hashlib.sha256((R3 / "CONVENTION3.md").read_bytes()).hexdigest()
    return frozen, actual, frozen == actual


def load_emb(ctx):
    p = LOCAL_EMB / f"{ctx}.npy"
    return np.load(p if p.exists() else X.CONTEXTS[ctx]["emb"], mmap_mode="r")


def load_cache(ctx):
    d = R3 / "cache" / ctx
    genes = json.loads((d / "obs_genes.json").read_text())
    obs_idx = np.load(d / "obs_idx.npy")
    nl = np.load(d / "node_labels.npz")
    return genes, obs_idx, nl, d


def load_gene2go():
    with open(GENE2GO, "rb") as f:
        return pickle.load(f)


def go_namespaces():
    ns, cur = {}, None
    for line in open(GO_OBO):
        line = line.strip()
        if line == "[Term]": cur = {}
        elif line.startswith("id: GO:") and cur is not None: cur["id"] = line[4:]
        elif line.startswith("namespace:") and cur is not None and "id" in cur: ns[cur["id"]] = line.split(": ", 1)[1]
        elif line.startswith("alt_id:") and cur is not None and "id" in cur:
            ns.setdefault("__alt__" + line[8:], cur["id"])
        elif line == "" and cur is not None: cur = None
    alt = {k[7:]: v for k, v in ns.items() if k.startswith("__alt__")}
    ns = {k: v for k, v in ns.items() if not k.startswith("__alt__")}
    for a, prim in alt.items():
        if prim in ns: ns[a] = ns[prim]
    return ns


def sv_coords(emb, obs_idx, n_dirs=N_DIRS):
    """(n_obs, 12 layers * n_dirs) observed-gene SVD coordinates U_k S_k; column = L*n_dirs + k.
    Also returns singular values per layer (12 x 512)."""
    n = obs_idx.size
    C = np.empty((n, N_LAYERS * n_dirs)); Ss = []
    for L in range(N_LAYERS):
        US, Sv, _ = S.svd_observed(emb, L, obs_idx)
        C[:, L * n_dirs:(L + 1) * n_dirs] = US[:, :n_dirs]; Ss.append(Sv)
    return C, np.array(Ss)


def null_summary(x):
    x = np.asarray(x, float)
    q = np.quantile(x, [0.025, 0.5, 0.95, 0.975, 0.99])
    return dict(N=int(x.size), mean=float(x.mean()), sd=float(x.std(ddof=1)),
                q025=float(q[0]), q50=float(q[1]), q95=float(q[2]), q975=float(q[3]), q99=float(q[4]))
