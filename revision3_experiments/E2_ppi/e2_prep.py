"""E2 prep: genome-wide STRING v12 physical degree (combined >= 700) and, per context, every within-V pair with
ANY physical STRING link (used to exclude linked pairs from the co-expression-matched non-edge pool), plus
ComplexPortal / hu.MAP 2.0 co-complex pairs within V (lung, lung_external).
Writes _ckpt/genome_degree700.json, _ckpt/<ctx>_phys_any.npy, _ckpt/<ctx>_complex_{cp,humap}.npy, prep_meta.json.
"""
import csv, gzip, json, re, itertools
import numpy as np
import e2_common as E

X = E.X
CP = X.REPO / "subproject_61_two_axes_revision/data/complexportal_9606.tsv"
HUMAP = X.REPO / "subproject_61_two_axes_revision/data/humap2_complexes.txt"
UNIPROT = X.REPO / "subproject_61_two_axes_revision/data/uniprot_human_reviewed_canonical.fasta.gz"
ACC = re.compile(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}")


def main():
    meta = {}
    allp = X.load_string_genomewide(min_score=0)              # {(g1,g2): max combined score}, physical file
    meta["string_physical_pairs_any_score"] = len(allp)
    deg = {}
    n700 = 0
    for (a, b), s in allp.items():
        if s >= 700:
            n700 += 1; deg[a] = deg.get(a, 0) + 1; deg[b] = deg.get(b, 0) + 1
    meta["string_physical_pairs_ge700"] = n700
    (E.CKPT / "genome_degree700.json").write_text(json.dumps(deg))
    # uniprot -> symbol
    u2s = {}
    with gzip.open(UNIPROT, "rt") as f:
        for line in f:
            if line.startswith(">"):
                acc = line.split("|")[1]; m = re.search(r" GN=(\S+)", line)
                if m: u2s[acc] = m.group(1)
    cp_sets = []
    with open(CP) as f:
        r = csv.reader(f, delimiter="\t"); h = next(r)
        col = h.index("Expanded participant list")
        for row in r:
            accs = {a for a in ACC.findall(row[col])}
            syms = sorted({u2s[a] for a in accs if a in u2s})
            if len(syms) >= 2: cp_sets.append(syms)
    hm_sets = []
    with open(HUMAP) as f:
        r = csv.DictReader(f)
        for row in r:
            syms = sorted(set(row["genenames"].split()))
            if len(syms) >= 2: hm_sets.append((int(row["Confidence"]), syms))
    meta["complexportal_complexes_ge2_mapped"] = len(cp_sets)
    meta["humap2_complexes_ge2"] = len(hm_sets)
    for ctx in E.PRIMARY_CTX + E.ROBUST_CTX:
        c = E.load_ctx(ctx); pos = {g: k for k, g in enumerate(c["genes"])}
        ph = np.array(sorted((pos[a], pos[b]) for (a, b) in allp if a in pos and b in pos), dtype=np.int64).reshape(-1, 2)
        np.save(E.CKPT / f"{ctx}_phys_any.npy", ph)
        m = dict(n_V=c["n"], phys_any_within_V=int(len(ph)), string700_within_V=int(len(c["edges"])))
        if ctx in ("lung", "lung_external"):
            for tag, sets in (("cp", cp_sets), ("humap", [s for _, s in hm_sets])):
                pr = set()
                for sset in sets:
                    v = sorted(pos[g] for g in sset if g in pos)
                    pr.update(itertools.combinations(v, 2))
                arr = np.array(sorted(pr), dtype=np.int64).reshape(-1, 2)
                np.save(E.CKPT / f"{ctx}_complex_{tag}.npy", arr)
                st = {(min(a, b), max(a, b)) for a, b in c["edges"]}
                m[f"complex_{tag}_pairs_within_V"] = int(len(arr))
                m[f"complex_{tag}_overlap_with_string700"] = int(sum((a, b) in st for a, b in arr))
                m[f"complex_{tag}_nodes"] = int(np.unique(arr).size) if len(arr) else 0
        meta[ctx] = m; print(ctx, m, flush=True)
    (E.HERE / "prep_meta.json").write_text(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
