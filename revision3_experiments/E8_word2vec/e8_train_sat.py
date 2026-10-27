"""E8 step 2b - training-budget saturation check for the production (full-corpus, window 5) models.

The calibration (e8_train.py processes 4/5) showed that on the extraction corpus the window-5 model at 10 epochs is
far from the window-512 model, while window 5 with the number of pair updates matched to window 512 is close
to it (CKA 0.80 lung, 0.83 immune). The production full-corpus models (window 5, 10 epochs) carry ~1/58 (lung)
of the pair updates the specified window-512 full-corpus model would. This script trains the same production
recipe with 4x the epochs (40) for seed 0 in lung and immune, to see whether the three task statistics are
still moving with training budget on the full corpus. It reuses e8_train.run unchanged (same corpus stream,
seeds, gensim parameters); only `epochs` differs.

Usage: python e8_train_sat.py <ctx>      (writes models/<ctx>__full__w5__e40__s0.{npy,json})
"""
import sys
from pathlib import Path
import numpy as np
import e8_train as T

ctx = sys.argv[1]
z = np.load(T.HERE / "corpus" / f"{ctx}.npz"); d = T.R3 / "cache" / ctx
obs_idx = np.load(d / "obs_idx.npy"); ext = np.load(d / "extraction_cells.npy")
G = max(int(z["tok"].max()), int(obs_idx.max())) + 1
data = {ctx: (z["tok"], z["off"], obs_idx, ext, np.array([f"g{i}" for i in range(G)], dtype=object))}
print("saturation job", ctx, flush=True)
T.run((ctx, "full", T.W_PROD, 40, 0), data)
print("DONE", ctx, flush=True)
