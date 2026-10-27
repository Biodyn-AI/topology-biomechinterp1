"""E2 sensitivity (g): the 120-statistic assortativity family on co-complex pairs within V (lung, lung_external),
Maslov-Sneppen N = 2,000. Edge sets (from e2_prep.py): cp = ComplexPortal (expanded participant list, UniProt ->
symbol), hm = hu.MAP 2.0 (all confidence levels); cpx / hmx = the same minus STRING-physical >= 700 edges
(co-complex evidence not already in the primary edge set).
"""
import numpy as np
import e2_common as E
from e2_nulls import run

if __name__ == "__main__":
    for ctx in ("lung", "lung_external"):
        st = {(min(a, b), max(a, b)) for a, b in E.load_ctx(ctx)["edges"]}
        for tag in ("cp", "hm"):
            e = np.load(E.CKPT / f"{ctx}_complex_{'cp' if tag == 'cp' else 'humap'}.npy")
            run(ctx, ["ms"], 2000, edges_override=e, tag=f"_{tag}")
            ex = np.array([p for p in map(tuple, e) if p not in st], dtype=np.int64).reshape(-1, 2)
            run(ctx, ["ms"], 2000, edges_override=ex, tag=f"_{tag}x")
