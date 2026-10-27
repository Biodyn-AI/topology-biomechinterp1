"""Round-2 canonical rebuild of every headline statistic. See CONVENTION.md."""
import sys, json, time
import numpy as np, pandas as pd
sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent))
import common as C

N_COPOLE, N_AUROC, SEED = 100_000, 10_000, 20260819
rng = np.random.default_rng(SEED)
t0 = time.time()

genes = C.h5ad_genes(); emb = C.load_emb()
nz0 = C.nonzero_mask(emb, 0)
cur = [i for i in C.curated_genes(genes) if nz0[i]]
sub = np.array(cur); names = [genes[i] for i in sub]; loc = {g: k for k, g in enumerate(names)}
n = len(names)
print(f"[setup] gene set n={n} (curated 209 minus {209-n} zero-norm at L0); "
      f"K={int(0.27*n)}; N_copole={N_COPOLE}; N_auroc={N_AUROC}", flush=True)

t = pd.read_csv(C.TRRUST_PATH, sep="\t", header=None, names=["TF","target","mode","pmid"])
tfs_all = set(t.TF.astype(str)); tgts_all = set(t.target.astype(str)) - tfs_all
tr = sorted({(a,b) for a,b in zip(t.TF.astype(str), t.target.astype(str)) if a in loc and b in loc and a!=b})
obj = json.loads(C.STRING_CACHE.read_text())
st = sorted({(p["g1"],p["g2"]) for p in obj["pairs"] if p["score"]>=0.7 and p["g1"] in loc and p["g2"] in loc and p["g1"]!=p["g2"]})
gene2go = C.load_gene2go()
tf_pos = np.array([k for k,g in enumerate(names) if g in tfs_all])
tgt_pos = np.array([k for k,g in enumerate(names) if g in tgts_all])
print(f"[setup] TRRUST pairs (deduplicated)={len(tr)}  STRING>=0.7={len(st)}  TF={tf_pos.size}  target={tgt_pos.size}", flush=True)

tri = np.array([loc[a] for a,b in tr]); trj = np.array([loc[b] for a,b in tr])
sti = np.array([loc[a] for a,b in st]); stj = np.array([loc[b] for a,b in st])

US = {L: C.layer_us(emb, L)[sub] for L in range(12)}
rows = []

# --- A. spectral collapse ---------------------------------------------------
print("[A] spectral collapse", flush=True)
def eff_rank(X):
    s = np.linalg.svd(X - X.mean(0, keepdims=True), compute_uv=False)
    p = s / s.sum(); p = p[p > 0]
    return float(np.exp(-(p*np.log(p)).sum()))
full_er = [eff_rank(np.asarray(emb[L], dtype=np.float64)) for L in range(12)]
sub_er  = [eff_rank(np.asarray(emb[L], dtype=np.float64)[sub]) for L in range(12)]
from scipy.stats import spearmanr
rho, p_rho = spearmanr(range(12), full_er)
Xs = np.asarray(emb[11], dtype=np.float64).copy()
for j in range(Xs.shape[1]): rng.shuffle(Xs[:, j])
print(f"    effective rank L0={full_er[0]:.2f} -> L11={full_er[11]:.2f} "
      f"({full_er[0]/full_er[11]:.1f}x)  Spearman rho={rho:.3f} p={p_rho:.2e}", flush=True)
print(f"    feature-shuffled L11 rank={eff_rank(Xs):.2f} (vs {full_er[11]:.2f}) -> compression is learned", flush=True)
rows.append(dict(claim="spectral_collapse_depth_trend", layer=None, stat=rho, null=None,
                 p=p_rho, p_kind="spearman_analytic", n_perm=None))

# --- B. SV1 GO compartment enrichment --------------------------------------
print("[B] SV1 GO compartment enrichment (Fisher, analytic)", flush=True)
def fisher_pole(US_L, comp, term, pole_frac=0.27):
    proj = US_L[:, comp]; top, bot, K = C.poles(proj, pole_frac)
    best = (0.0, 1.0)
    for pole in (top, bot):
        a = sum(1 for k in range(n) if k in pole and term in gene2go.get(names[k], set()))
        c = sum(1 for k in range(n) if k not in pole and term in gene2go.get(names[k], set()))
        o, p = C.fisher_2x2(a, len(pole)-a, c, (n-len(pole))-c)
        if p < best[1]: best = (o, p)
    return best
for term, lbl in [(C.GO_EXTRACELLULAR,"extracellular"), (C.GO_CYTOSOL,"cytosol")]:
    per = [fisher_pole(US[L], 0, term) for L in range(12)]
    ors = [x[0] for x in per]; ps = [x[1] for x in per]
    b = int(np.argmin(ps))
    print(f"    SV1 {lbl:14s} best L{b}: OR={ors[b]:.2f} p={ps[b]:.2e} | L11: OR={ors[11]:.2f} p={ps[11]:.2e}", flush=True)
    rows.append(dict(claim=f"SV1_{lbl}_L11", layer=11, stat=ors[11], null=None,
                     p=ps[11], p_kind="fisher_exact", n_perm=None))
    rows.append(dict(claim=f"SV1_{lbl}_best_layer", layer=b, stat=ors[b], null=None,
                     p=ps[b], p_kind="fisher_exact", n_perm=None))

# --- C. co-pole tests, N=100,000 gene-label shuffle ------------------------
def copole_claim(tag, comp, pi, pj):
    out = []
    for L in range(12):
        lab = C.pole_labels(US[L][:, comp])
        obs = C.copole_rate(lab, pi, pj)
        null = C.copole_null(lab, pi, pj, N_COPOLE, rng)
        p, b, N = C.emp_p(obs, null)
        out.append((L, obs, float(null.mean()), p, b))
    sig = sum(1 for _,_,_,p,_ in out if p < 0.05)
    bonf = sum(1 for _,_,_,p,_ in out if p < C.BONFERRONI_ALPHA)
    Lb = min(out, key=lambda r: r[3])
    print(f"    {tag:26s} best L{Lb[0]:<2} obs={Lb[1]:.3f} null={Lb[2]:.3f} p={Lb[3]:.2e} "
          f"(b={Lb[4]}/{N_COPOLE}) | p<0.05 {sig}/12 | p<alpha_bonf {bonf}/12", flush=True)
    for L,obs,nu,p,b in out:
        rows.append(dict(claim=f"{tag}_L{L}", layer=L, stat=obs, null=nu, p=p,
                         p_kind="perm_label_shuffle", n_perm=N_COPOLE))
    return out
print(f"[C] co-pole permutation tests (N={N_COPOLE:,})", flush=True)
copole_claim("SV2_TRRUST_copole", 1, tri, trj)
for comp, lbl in [(1,"SV2"), (2,"SV3"), (3,"SV4")]:
    copole_claim(f"{lbl}_STRING_copole", comp, sti, stj)

# --- D. TF/target AUROC, leave-one-out, N=10,000 ---------------------------
print(f"[D] TF/target AUROC SV2-SV7, leave-one-out (N={N_AUROC:,})", flush=True)
def auroc_perm(US_L, comps, tf_pos, tgt_pos, n_perm):
    E = US_L[:, comps]; E = E / (np.linalg.norm(E, axis=1, keepdims=True) + 1e-12)
    S = E @ E.T
    labelled = np.concatenate([tf_pos, tgt_pos]); k = tf_pos.size
    obs = C.tf_target_auroc(US_L, comps, tf_pos, tgt_pos, exclude_self=True)
    null = np.empty(n_perm)
    for r in range(n_perm):
        perm_tf = rng.choice(labelled, size=k, replace=False)
        tot = S[:, perm_tf].sum(axis=1)
        selfs = np.zeros(S.shape[0]); selfs[perm_tf] = 1.0
        istf = np.zeros(S.shape[0], dtype=bool); istf[perm_tf] = True
        sc = (tot - selfs) / np.where(istf, k-1, k)
        rest = np.setdiff1d(labelled, perm_tf, assume_unique=False)
        null[r] = C.auroc(sc[perm_tf], sc[rest])
    p, b, N = C.emp_p(obs, null)
    return obs, float(null.mean()), p, b
res = []
for L in range(12):
    o, nu, p, b = auroc_perm(US[L], slice(1,7), tf_pos, tgt_pos, N_AUROC)
    res.append((L,o,nu,p,b))
    rows.append(dict(claim=f"TF_target_AUROC_SV2_7_L{L}", layer=L, stat=o, null=nu, p=p,
                     p_kind="perm_label_shuffle", n_perm=N_AUROC))
mean_a = np.mean([r[1] for r in res]); Lb = max(res, key=lambda r: r[1])
sig = sum(1 for r in res if r[3] < 0.05); bonf = sum(1 for r in res if r[3] < C.BONFERRONI_ALPHA)
print(f"    mean AUROC={mean_a:.3f}  max={Lb[1]:.3f}@L{Lb[0]}  null_mean={np.mean([r[2] for r in res]):.3f}", flush=True)
print(f"    p<0.05 {sig}/12 | p<alpha_bonf {bonf}/12 | best p={min(r[3] for r in res):.2e}", flush=True)
Z, pZ = C.stouffer([r[3] for r in res])
print(f"    Stouffer across 12 layers: Z={Z:.2f} p={pZ:.3e}", flush=True)
rows.append(dict(claim="TF_target_AUROC_SV2_7_stouffer", layer=None, stat=mean_a, null=None,
                 p=pZ, p_kind="stouffer", n_perm=N_AUROC))

df = pd.DataFrame(rows)
df.to_csv(C.OUT/"headline_rebuilt.csv", index=False)
print(f"\n[done] {len(df)} statistics written to headline_rebuilt.csv  ({time.time()-t0:.0f}s)", flush=True)
