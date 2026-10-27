"""Census of every context: tissue identity, vocabulary, and how much reference signal it can test."""
import sys, json
import numpy as np, pandas as pd
sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
import contexts as X

MARK = {"kidney": ["LRP2","CUBN","SLC34A1","UMOD","PAX8","PAX2","SLC12A1","SLC22A8"],
        "immune": ["PTPRC","CD3E","CD3D","MS4A1","CD79A","LYZ","CD14","NKG7"],
        "lung":   ["SFTPC","SFTPB","SFTPA1","AGER","NKX2-1","SCGB1A1","EPCAM","CAV1"]}
t, tr = X.load_trrust()
tfs = set(t.TF)
stg = X.load_string_genomewide(700)
dch = set(X.load_dorothea(chipseq=True))
dcu = set(X.load_dorothea(chipseq=False))
trs = set(tr)
cur209 = set(pd.read_csv(X.S38/"cycle1_main/cycle1_edge_dataset.tsv", sep="\t")[["source","target"]].values.ravel())
rows = []
for name, c in X.CONTEXTS.items():
    emb, genes = X.load_context(name)
    cnt_path = c["emb"].parent / "layer_gene_embedding_count.npy"
    cnt = np.load(cnt_path)[0] if cnt_path.exists() else None
    gi = {g: i for i, g in enumerate(genes)}
    nz = np.linalg.norm(np.asarray(emb[0], dtype=np.float32), axis=1) > 0
    V = {g for g in genes if nz[gi[g]]}
    mk = {k: int(np.mean([cnt[gi[m]] for m in v if m in gi])) if cnt is not None and any(m in gi for m in v) else -1
          for k, v in MARK.items()}
    tr_in = [(a, b) for a, b in trs if a in V and b in V]
    tf_in = {a for a, _ in tr_in}
    tg_in = {b for _, b in tr_in} - tfs
    st_in = sum(1 for (a, b) in stg if a in V and b in V)
    rows.append(dict(context=name, tissue=c["tissue"], group=c["group"], G=len(genes), nonzero=len(V),
                     curated209=len(cur209 & V), TRRUST_edges=len(tr_in), TFs=len(tf_in), targets=len(tg_in),
                     DoRothEA_chip_noTRRUST=sum(1 for e in dch - trs if e[0] in V and e[1] in V),
                     DoRothEA_ABC_noTRRUST=sum(1 for e in dcu - trs if e[0] in V and e[1] in V),
                     STRING700=st_in, mk_kidney=mk["kidney"], mk_immune=mk["immune"], mk_lung=mk["lung"]))
    print(rows[-1], flush=True)
df = pd.DataFrame(rows)
df.to_csv(X.Path(__file__).parent/"00_context_census.csv", index=False)
pd.set_option("display.width", 250)
print("\n", df.to_string(index=False))
print(f"\nSTRING>=700 genome-wide pairs: {len(stg)} | TRRUST distinct pairs: {len(trs)} | "
      f"DoRothEA chip ABC: {len(dch)} | DoRothEA ABC: {len(dcu)}")
