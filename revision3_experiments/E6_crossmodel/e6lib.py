"""E6 cross-model comparison: shared data layer and fast exact paths.

All primary statistics come from ../stats3.py (tf_scores_loo, auroc, tf_auroc_null, StratifiedPairs,
vertex_subsample_ci, maslov_sneppen_chain, max_t_adjust). Three vectorised paths are added here ONLY for
resampling/null loops, and each is checked against the stats3 function on live data every time it runs
(`check_*` functions; a mismatch > 1e-9 aborts):
  * node_boot      gene bootstrap of the LOO-cosine TF AUROC through a precomputed cosine block G
                   (leave-gene-out on the row-duplicated matrix; equals stats3.tf_scores_loo at unit counts);
  * FastStrat      weighted within-TF stratified AUROC (equals stats3.StratifiedPairs.auroc(s, w)) and the
                   within-TF label-permutation null through fixed within-TF average ranks;
  * auroc_sorted   tie-aware AUROC against a fixed, pre-sorted negative set (equals stats3.auroc).
"""
from __future__ import annotations

import json
import os
import pickle
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
R3 = HERE.parent
sys.path.insert(0, str(R3))
import contexts as X          # noqa: E402
import stats3 as S            # noqa: E402
import prep_pairs as PP       # noqa: E402

WORK = HERE / "work"; WORK.mkdir(exist_ok=True)
CACHE = R3 / "cache"
LOCAL_EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"
REPO = X.REPO
FROZEN = (R3 / "FROZEN.txt").read_text().split()[2]

SEED = 20261026
CORE_CTX = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
BANDS = [(a, b) for a in range(1, 13) for b in range(a, 13)]            # 78 contiguous bands SVa..SVb
BAND_NAMES = [f"SV{a}-SV{b}" if a != b else f"SV{a}" for a, b in BANDS]
FIXED_BANDS = {"SV2-SV7": (2, 7), "SV1-SV6": (1, 6)}

# key -> (npz path, display label, checkpoint / provenance note, kind)
STATIC = {
    "scgpt_static": ("subproject_53_scgpt_gpl_replication/embeddings/scgpt_gene_embeddings.npz",
                     "scGPT static token embedding",
                     "scGPT whole-human best_model.pt, encoder nn.Embedding weights (60,694 x 512)", "static token"),
    "gf_v2_316m": ("subproject_54_geneformer_gpl_replication/embeddings/geneformer_gene_embeddings.npz",
                   "Geneformer V2-316M",
                   "ctheodoris/Geneformer HF snapshot 05fcbeb8 root model.safetensors (Geneformer-V2-316M per "
                   "subproject_54 PLAN/phase0: BERT, 18 layers, d=1152, gc104M tokenizer); input embedding + embedding LayerNorm",
                   "static token"),
    "gf_v1_10m": ("subproject_55_geneformer_v1_10m_gpl_replication/embeddings/geneformer_v1_10m_gene_embeddings.npz",
                  "Geneformer V1-10M", "ctheodoris/Geneformer Geneformer-V1-10M/model.safetensors, gc30M tokenizer, d=256",
                  "static token"),
    "gf_v2_104m": ("subproject_56_geneformer_v2_104m_gpl_replication/embeddings/geneformer_v2_104m_gene_embeddings.npz",
                   "Geneformer V2-104M", "ctheodoris/Geneformer Geneformer-V2-104M/model.safetensors, gc104M tokenizer, d=768",
                   "static token"),
    "aido_100m": ("subproject_57_aido_cell_100m_gpl_replication/embeddings/aido_cell_100m_gene_embeddings.npz",
                  "AIDO.Cell-100M", "genbio-ai/AIDO.Cell-100M; NOT a token table: last hidden state of one forward pass "
                  "on the HPA mean-expression cell (d=640)", "forward pass, neutral cell"),
    "scfoundation": ("subproject_58_scfoundation_gpl_replication/embeddings/scfoundation_gene_embeddings.npz",
                     "scFoundation", "genbio-ai/scFoundation models.ckpt; NOT a token table: official encoder+decoder "
                     "'gene' output of one forward pass on the HPA mean-expression cell (d=512)", "forward pass, neutral cell"),
    "tgpt": ("subproject_59_tgpt_gpl_replication/embeddings/tgpt_gene_embeddings.npz",
             "tGPT", "lixiangchun/transcriptome-gpt-1024-8-16-64 transformer.wte.weight (d=1024)", "static token"),
    "esm2_uce": ("subproject_60_uce_esm2_control_gpl_replication/embeddings/uce_esm2_human_gene_embeddings.npz",
                 "ESM2 (UCE protein file; 'ESM2-3B' in submitted version)",
                 "minwoosun/uce-misc Homo_sapiens.GRCh38.gene_symbol_to_embedding_ESM2.pt, 5120-d per-gene protein "
                 "embedding. 5120 is the hidden width of ESM2-15B (ESM2-3B is 2560-d); subproject_61's registry labels "
                 "this file ESM2-15B. The submitted manuscript called it ESM2-3B.", "protein LM"),
    "esm2_650m": ("subproject_61_two_axes_revision/data/derived_embeddings/esm2_650m_own.npz",
                  "ESM2-650M", "facebook/esm2_t33_650M_UR50D, final layer mean-pooled over residues of the longest "
                  "Swiss-Prot canonical sequence per gene symbol (subproject_61 e3b), d=1280", "protein LM"),
    "seq_features": ("subproject_61_two_axes_revision/data/derived_embeddings/seq_features.npz",
                     "32 sequence descriptors", "subproject_61 hand-built protein-sequence descriptors (d=32); "
                     "181 all-zero rows treated as missing", "hand-built"),
    "hpa_coexp": ("subproject_61_two_axes_revision/data/derived_embeddings/hpa_coexpression.npz",
                  "HPA co-expression", "subproject_61 e13: log1p HPA single-cell nCPM over 154 cell types, centred, "
                  "PCA to 128 components", "empirical expression"),
}
STATIC_KEYS = list(STATIC)


def load_static(key):
    """(genes, E float32) with all-zero rows dropped (treated as not in the vocabulary)."""
    z = np.load(REPO / STATIC[key][0])
    g = np.array([str(x) for x in z["gene_names"]]); E = np.asarray(z["embeddings"], dtype=np.float32)
    nz = np.abs(E).sum(1) > 0
    return list(g[nz]), E[nz]


def vocab(key):
    p = WORK / f"vocab_{key}.json"
    if p.exists(): return json.loads(p.read_text())
    g, _ = load_static(key); p.write_text(json.dumps(g)); return g


def static_spectral(key):
    """Top-12 U*S coordinates of the model's FULL (zero-row-free) vocabulary, mean-centred; cached."""
    p = WORK / f"svd12_{key}.npz"
    if p.exists():
        z = np.load(p); return [str(x) for x in z["genes"]], z["US"], z["S"]
    g, E = load_static(key)
    Xc = E.astype(np.float64); Xc -= Xc.mean(0, keepdims=True)
    if Xc.shape[1] > 2000:                      # ESM2 5120-d: eigendecomposition of X^T X (exact, cheaper)
        C = Xc.T @ Xc; w, V = np.linalg.eigh(C); o = np.argsort(w)[::-1]; w, V = w[o], V[:, o]
        S_ = np.sqrt(np.maximum(w, 0)); US = Xc @ V[:, :12]
    else:
        U, S_, Vt = np.linalg.svd(Xc, full_matrices=False); US = U[:, :12] * S_[:12]
    sgn = np.sign(US[np.abs(US).argmax(0), np.arange(12)]); sgn[sgn == 0] = 1; US = US * sgn
    np.savez(p, genes=np.array(g), US=US, S=S_[:min(len(S_), 512)])
    return g, US, S_


def ctx_obs(ctx):
    genes = json.loads((CACHE / ctx / "obs_genes.json").read_text())
    return genes, np.load(CACHE / ctx / "obs_idx.npy")


def ctx_emb(ctx):
    p = LOCAL_EMB / f"{ctx}.npy"
    return np.load(p if p.exists() else X.CONTEXTS[ctx]["emb"], mmap_mode="r")


def string700():
    p = WORK / "string700.pkl"
    if p.exists(): return pickle.loads(p.read_bytes())
    d = X.load_string_genomewide(700); p.write_bytes(pickle.dumps(d)); return d


def refs():
    p = WORK / "refs.pkl"
    if p.exists(): return pickle.loads(p.read_bytes())
    trrust, anyreg, chip, lambert = PP.load_refs()
    t, _ = X.load_trrust()
    tfs_g = set(t.TF.astype(str)); tgt_only_g = set(t.target.astype(str)) - tfs_g   # submitted definition
    r = dict(trrust=trrust, anyreg=anyreg, lambert=lambert, trrust_tf=tfs_g, trrust_target_only=tgt_only_g)
    p.write_bytes(pickle.dumps(r)); return r


# =====================================================================================
# gene-set benchmark construction
# =====================================================================================
def build_geneset(name, genes, det_dec, rf, stg, seed=SEED):
    """genes: sorted list; det_dec: detection decile per gene (int array). Returns dict of benchmark arrays."""
    pos = {g: k for k, g in enumerate(genes)}; G = set(genes)
    lam = np.array([g in rf["lambert"] for g in genes])
    ttf = np.array([g in rf["trrust_tf"] for g in genes]); ttg = np.array([g in rf["trrust_target_only"] for g in genes])
    st = sorted((a, b) for (a, b) in stg if a in G and b in G)
    pi = np.array([pos[a] for a, b in st], dtype=np.int64); pj = np.array([pos[b] for a, b in st], dtype=np.int64)
    # one Maslov-Sneppen rewiring of the gene-set-induced graph -> non-pairs (true pairs removed)
    a, b, acc = next(S.maslov_sneppen_chain(pi, pj, len(genes), 1, seed=seed, burn_factor=10))
    true = set(zip(np.minimum(pi, pj).tolist(), np.maximum(pi, pj).tolist()))
    k = np.array([(int(x), int(y)) not in true for x, y in zip(a, b)])
    ni, nj = a[k], b[k]
    # regulatory edges: prep_pairs rules restricted to the gene set
    trrust = rf["trrust"]
    tr_v = sorted(p for p in trrust if p[0] in G and p[1] in G)
    tfs_tr = {x for x, _ in tr_v}
    indeg = {}
    for x, y in tr_v: indeg[y] = indeg.get(y, 0) + 1
    other = {tf: {y for x, y in tr_v if x != tf} for tf in tfs_tr}
    rng = np.random.default_rng(seed)
    res = PP.build(tr_v, genes, pos, rf["anyreg"], indeg, np.asarray(det_dec), rng, other)
    tf, tg, y, lev = res if res is not None else (np.zeros(0, int),) * 4
    tp = sorted(set(g for p in tr_v for g in p))
    return dict(name=name, genes=np.array(genes), lambert=lam, trrust_tf=ttf, trrust_tgt=ttg,
                pi=pi, pj=pj, ni=ni, nj=nj, ms_acc=acc, n_rewired_true_removed=int((~k).sum()),
                e_tf=np.asarray(tf), e_tg=np.asarray(tg), e_y=np.asarray(y), e_lev=np.asarray(lev),
                det_dec=np.asarray(det_dec), n_trrust_within=len(tr_v), n_trrust_participating=len(tp))


def save_geneset(gs):
    np.savez(WORK / f"geneset_{gs['name']}.npz", **{k: v for k, v in gs.items() if k != "name"})


def load_geneset(name):
    z = np.load(WORK / f"geneset_{name}.npz", allow_pickle=False)
    d = {k: z[k] for k in z.files}; d["name"] = name; d["genes"] = [str(x) for x in d["genes"]]; return d


# =====================================================================================
# node task: fast paired bootstrap through a cosine block
# =====================================================================================
def node_boot(E, tf_pos, neg_pos, B, seed):
    """Gene bootstrap: labelled genes resampled with replacement; the TF centroid is rebuilt from the drawn TF
    copies and each gene is scored leave-GENE-out (all copies of itself excluded, so a duplicated TF is never
    compared with its own copy -- otherwise duplicates inflate TF scores by cos(self)=1). Same draws for every
    embedding with the same (n_lab, seed), so draws are paired across embeddings."""
    lab = np.concatenate([tf_pos, neg_pos]); ntf = tf_pos.size; nlab = lab.size
    En = S.l2norm_rows(np.asarray(E, dtype=np.float64))
    Gm = En[lab] @ En[tf_pos].T; selfsim = (En[lab] * En[lab]).sum(1)
    rng = np.random.default_rng(seed); out = np.empty(B); counts0 = None
    for r in range(B):
        c = np.bincount(rng.integers(0, nlab, nlab), minlength=nlab)
        if r == 0: counts0 = c
        ctf = c[:ntf]; K = ctf.sum()
        tot = Gm @ ctf
        sc_tf = (tot[:ntf] - ctf * selfsim[:ntf]) / np.maximum(K - ctf, 1); sc_ng = tot[ntf:] / K
        out[r] = S.auroc(np.repeat(sc_tf, ctf), np.repeat(sc_ng, c[ntf:]))
    return out, counts0


def check_node_boot(E, tf_pos, neg_pos, counts, val):
    """Re-derive replicate 0 explicitly on the row-duplicated matrix (cosine to every TF copy of a DIFFERENT
    gene, averaged), and check that with all counts = 1 it equals stats3.tf_scores_loo."""
    lab = np.concatenate([tf_pos, neg_pos]); ntf = tf_pos.size
    En = S.l2norm_rows(np.asarray(E, dtype=np.float64))
    idx = np.repeat(np.arange(lab.size), counts); is_tf = idx < ntf
    Ed = En[lab][idx]; Cm = Ed @ Ed[is_tf].T
    same = idx[:, None] == idx[is_tf][None, :]
    Cm[same] = 0.0; den = (~same).sum(1)
    s = Cm.sum(1) / np.maximum(den, 1)
    ref = S.auroc(s[is_tf], s[~is_tf])
    assert abs(ref - val) < 1e-9, ("node_boot mismatch", ref, val)
    s1 = S.tf_scores_loo(np.asarray(E, float)[lab], np.arange(ntf))
    Cm1 = En[lab] @ En[tf_pos].T; Cm1[np.arange(ntf), np.arange(ntf)] = 0
    s2 = Cm1.sum(1) / np.where(np.arange(lab.size) < ntf, ntf - 1, ntf)
    assert np.allclose(s1, s2, atol=1e-10), "leave-gene-out != tf_scores_loo at unit counts"


def fast_tf_auroc_null(E, tf_pos, tgt_pos, n_perm, rng, batch=500):
    """stats3.tf_auroc_null with identical rng draws, scores and ranks; the second argsort (rank = inverse of the
    sort permutation) is replaced by a scatter, which yields the identical rank matrix at half the cost."""
    lab = np.concatenate([tf_pos, tgt_pos]); k = tf_pos.size; n = lab.size
    En = S.l2norm_rows(np.asarray(E, dtype=np.float64)[lab]); selfsim = (En * En).sum(1)
    out = np.empty(n_perm); done = 0; ar = np.arange(1, n + 1, dtype=np.float64)
    while done < n_perm:
        b = min(batch, n_perm - done)
        M = np.zeros((n, b)); idx = np.argsort(rng.random((n, b)), axis=0)[:k]
        M[idx, np.arange(b)[None, :]] = 1.0
        tot = En @ (En.T @ M)
        sc = np.where(M > 0, (tot - selfsim[:, None]) / max(k - 1, 1), tot / k)
        o = sc.argsort(axis=0); ranks = np.empty_like(sc); np.put_along_axis(ranks, o, np.broadcast_to(ar[:, None], sc.shape), axis=0)
        pm = M > 0; n1 = pm.sum(0); n0 = n - n1
        out[done:done + b] = ((ranks * pm).sum(0) - n1 * (n1 + 1) / 2) / (n1 * n0)
        done += b
    return out


def null_from_perms(E, tf_pos, tgt_pos, perms):
    """Same statistic as fast_tf_auroc_null for given permutations (perms: k x N positions within lab), so one set of
    permutations is shared by every band (joint null for max-T) and generated once."""
    lab = np.concatenate([tf_pos, tgt_pos]); k = tf_pos.size; n = lab.size
    En = S.l2norm_rows(np.asarray(E, dtype=np.float64)[lab]); selfsim = (En * En).sum(1)
    N = perms.shape[1]; out = np.empty(N); ar = np.arange(1, n + 1, dtype=np.float64)
    for c0 in range(0, N, 500):
        idx = perms[:, c0:c0 + 500]; b = idx.shape[1]
        M = np.zeros((n, b)); M[idx, np.arange(b)[None, :]] = 1.0
        tot = En @ (En.T @ M)
        sc = np.where(M > 0, (tot - selfsim[:, None]) / max(k - 1, 1), tot / k)
        o = sc.argsort(axis=0); ranks = np.empty_like(sc); np.put_along_axis(ranks, o, np.broadcast_to(ar[:, None], sc.shape), axis=0)
        pm = M > 0; n1 = pm.sum(0); n0 = n - n1
        out[c0:c0 + b] = ((ranks * pm).sum(0) - n1 * (n1 + 1) / 2) / (n1 * n0)
    return out


def check_fast_null(E, tf_pos, tgt_pos, seed):
    a = S.tf_auroc_null(E, tf_pos, tgt_pos, 500, np.random.default_rng(seed))
    b = fast_tf_auroc_null(E, tf_pos, tgt_pos, 500, np.random.default_rng(seed))
    assert np.max(np.abs(a - b)) < 1e-12, "fast_tf_auroc_null mismatch"


# =====================================================================================
# edge task: exact vectorised within-TF stratified AUROC (weights) and permutation null
# =====================================================================================
class FastStrat:
    def __init__(self, sp, s):
        s = np.asarray(s, float)
        o = np.lexsort((s, sp.tf_code)); self.o = o
        code = sp.tf_code[o]; ss = s[o]; self.y = sp.y[o]; n = s.size
        newg = np.r_[True, (code[1:] != code[:-1]) | (ss[1:] != ss[:-1])]
        gid = np.cumsum(newg) - 1; gs = np.flatnonzero(newg); ge = np.r_[gs[1:], n]
        newt = np.r_[True, code[1:] != code[:-1]]; ts = np.flatnonzero(newt)
        self.code = code; self.gs = gs[gid]; self.ge = ge[gid]; self.tstart = ts[code]; self.ts = ts
        # within-TF average ranks (1-based) for the unweighted permutation null
        pos_in_tf = np.arange(n) - self.tstart
        self.rank = (pos_in_tf[self.gs] + 1 + pos_in_tf[self.ge - 1] + 1) / 2.0     # avg of tie group ranks

    def auroc_w(self, W):
        """W: (n_pairs, B) weights in the ORIGINAL pair order -> (B,) weighted stratified AUROC."""
        W = np.asarray(W, float)[self.o]
        y = self.y[:, None]; Wn = W * (~y); Wp = W * y
        cs = np.vstack([np.zeros((1, W.shape[1])), np.cumsum(Wn, 0)])
        lo = cs[self.gs] - cs[self.tstart]; hi = cs[self.ge] - cs[self.tstart]
        num = (Wp * (lo + 0.5 * (hi - lo))).sum(0)
        den = (np.add.reduceat(Wp, self.ts, 0) * np.add.reduceat(Wn, self.ts, 0)).sum(0)
        return num / den

    def perm_null(self, Yperm):
        """Yperm: (n_pairs, N) bool labels in ORIGINAL order, each column a within-TF permutation."""
        Y = np.asarray(Yperm, float)[self.o]
        n1 = np.add.reduceat(Y, self.ts, 0); ntot = np.diff(np.r_[self.ts, Y.shape[0]])[:, None]
        n0 = ntot - n1
        num = self.rank @ Y - (n1 * (n1 + 1) / 2).sum(0)
        return num / (n1 * n0).sum(0)


def two_way_weights(sp, B, rng):
    """The exact weights stats3.StratifiedPairs.bootstrap(two_way=True) draws, in the same rng order
    (call repeatedly with the same rng to stream the draws in column chunks)."""
    W = np.empty((sp.tf.size, B))
    for r in range(B):
        wtf = np.bincount(rng.integers(0, sp.tfs.size, sp.tfs.size), minlength=sp.tfs.size)
        w = wtf[sp.tf_code].astype(float)
        wtg = np.bincount(rng.integers(0, sp.tgs.size, sp.tgs.size), minlength=sp.tgs.size)
        W[:, r] = w * wtg[sp.tg_code]
    return W


def within_tf_perms(sp, N, seed):
    rng = np.random.default_rng(seed); Y = np.empty((sp.y.size, N), bool)
    for r in range(N): Y[:, r] = sp.permute_within_tf(rng)
    return Y


def check_faststrat(sp, s, W, Y):
    fs = FastStrat(sp, s)
    a0 = sp.auroc(s); assert abs(fs.auroc_w(np.ones((s.size, 1)))[0] - a0) < 1e-9
    for r in range(min(3, W.shape[1])):
        assert abs(fs.auroc_w(W[:, [r]])[0] - sp.auroc(s, W[:, r])) < 1e-9, "weighted strat mismatch"
    for r in range(min(3, Y.shape[1])):
        sp2 = S.StratifiedPairs(sp.tf, sp.tg, Y[:, r])
        assert abs(fs.perm_null(Y[:, [r]])[0] - sp2.auroc(s)) < 1e-9, "perm strat mismatch"


# =====================================================================================
# PPI task
# =====================================================================================
def pair_cos(En, a, b, chunk=None):
    chunk = chunk or max(256, int(4e6 / max(En.shape[1], 1)))
    out = np.empty(a.size)
    for s0 in range(0, a.size, chunk):
        out[s0:s0 + chunk] = np.einsum("ij,ij->i", En[a[s0:s0 + chunk]], En[b[s0:s0 + chunk]])
    return out


def auroc_sorted(pos, neg_sorted):
    lo = np.searchsorted(neg_sorted, pos, side="left"); hi = np.searchsorted(neg_sorted, pos, side="right")
    return float((lo + 0.5 * (hi - lo)).sum() / (pos.size * neg_sorted.size))


def band_cos(C, a, b):
    """Cosines of pairs (a, b) on every contiguous band of the 12 columns of C: (m, 78)."""
    P = np.cumsum(np.c_[np.zeros(a.size), C[a] * C[b]], 1)
    Qa = np.cumsum(np.c_[np.zeros(a.size), C[a] ** 2], 1); Qb = np.cumsum(np.c_[np.zeros(a.size), C[b] ** 2], 1)
    out = np.empty((a.size, len(BANDS)))
    for k, (lo, hi) in enumerate(BANDS):
        na = np.sqrt(Qa[:, hi] - Qa[:, lo - 1]); nb = np.sqrt(Qb[:, hi] - Qb[:, lo - 1])
        out[:, k] = (P[:, hi] - P[:, lo - 1]) / ((na + 1e-12) * (nb + 1e-12))
    return out


def ms_draws_cached(setname, pi, pj, n, N, seed):
    """N chained Maslov-Sneppen rewirings (burn-in 10m, thinning m) with true pairs removed, cached per gene set.
    Returns list of (a, b) index arrays. Identical to iterating stats3.maslov_sneppen_chain(pi, pj, n, N, seed)."""
    f = WORK / f"msdraws_{setname}_{seed}_{N}.npz"
    if not f.exists():
        true = set(zip(np.minimum(pi, pj).tolist(), np.maximum(pi, pj).tolist()))
        A, Bv, off, accs = [], [], [0], []
        for a, b, acc in S.maslov_sneppen_chain(pi, pj, n, N, seed=seed, burn_factor=10):
            k = np.fromiter(((int(x), int(y)) not in true for x, y in zip(a, b)), bool, a.size)
            A.append(a[k].astype(np.uint16)); Bv.append(b[k].astype(np.uint16)); off.append(off[-1] + int(k.sum())); accs.append(acc)
        tmp = f.with_name(f.stem + f".{os.getpid()}.tmp.npz")
        np.savez(tmp, a=np.concatenate(A), b=np.concatenate(Bv), off=np.array(off), acc=np.array(accs))
        tmp.rename(f)
    z = np.load(f); a, b, off = z["a"].astype(np.int64), z["b"].astype(np.int64), z["off"]
    return [(a[off[r]:off[r + 1]], b[off[r]:off[r + 1]]) for r in range(off.size - 1)]


def edge_perms_cached(setname, sp, N, seed):
    """N within-TF label permutations (stats3.StratifiedPairs.permute_within_tf, one rng), cached per gene set."""
    f = WORK / f"edgeperms_{setname}_{seed}_{N}.npz"
    if not f.exists():
        Y = within_tf_perms(sp, N, seed)
        tmp = f.with_name(f.stem + f".{os.getpid()}.tmp.npz"); np.savez(tmp, bits=np.packbits(Y, axis=0), n=Y.shape[0]); tmp.rename(f)
    z = np.load(f)
    return np.unpackbits(z["bits"], axis=0, count=int(z["n"])).astype(bool)


def auroc_sorted_fast(pos, neg_sorted, ties):
    """Exact tie-aware AUROC vs a fixed sorted negative set; queries are sorted first (cache-friendly search).
    ties=False skips the right-side search (valid when pos and neg share no values, i.e. continuous scores)."""
    ps = np.sort(pos)
    lo = np.searchsorted(neg_sorted, ps, side="left")
    if ties:
        hi = np.searchsorted(neg_sorted, ps, side="right"); return float((lo + 0.5 * (hi - lo)).sum() / (pos.size * neg_sorted.size))
    return float(lo.sum() / (pos.size * neg_sorted.size))


def null_summary(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    q = np.quantile(x, [0.025, 0.05, 0.5, 0.95, 0.975]) if x.size else [np.nan] * 5
    return dict(N=int(x.size), mean=float(x.mean()) if x.size else np.nan, sd=float(x.std(ddof=1)) if x.size > 1 else np.nan,
                q025=float(q[0]), q05=float(q[1]), q50=float(q[2]), q95=float(q[3]), q975=float(q[4]))


def ci(draws, est=None, level=0.95):
    d = np.asarray(draws, float); d = d[np.isfinite(d)]; a = (1 - level) / 2
    return float(np.quantile(d, a)), float(np.quantile(d, 1 - a))
