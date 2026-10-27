"""E6 (resume agent): complete estimator (iii) for the STATIC models in Tier B.

The previous agent skipped the 78-band max-T (and the PPI band null) for static models in Tier B to save compute
(e6_run.SKIP_MAXT_STATIC_TIERB). Machine time allowed completing it. This driver re-runs e6_run.run() for those
55 (context, static model) jobs with the skip disabled. Every seed depends only on (set, task), so estimators (i)
and (ii) are recomputed bit-identically; the previous outputs are moved to work/_backup_tierB_static_nomaxt/ and
compared after the run (check_identical). Workers: python e6_tierb_static_maxt.py worker  (lock files *.lock2).
"""
import json, os, shutil, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import e6lib as L
import e6_run as R
R.SKIP_MAXT_STATIC_TIERB = False
BK = L.WORK / "_backup_tierB_static_nomaxt"
CTX = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]


def jobs():
    return [(f"tierB_{c}", k) for c in CTX for k in L.STATIC_KEYS]


def worker():
    for s, e in jobs():
        d = L.WORK / s
        try:
            fd = os.open(d / f"{e}.lock2", os.O_CREAT | os.O_EXCL | os.O_WRONLY); os.close(fd)
        except FileExistsError:
            continue
        b = BK / s; b.mkdir(parents=True, exist_ok=True)
        for ext in (".json", ".npz", ".p1d.npz"):
            f = d / f"{e}{ext}"
            if f.exists():
                z = np.load(f) if ext == ".npz" else None
                if ext == ".npz" and any(k.startswith("tobs_") for k in z.files):
                    continue                                  # already has band statistics (finished re-run)
                shutil.move(str(f), str(b / f.name))
        R.run(s, e)


def check_identical():
    bad = 0; n = 0
    for s, e in jobs():
        fo = BK / s / f"{e}.json"; fn = L.WORK / s / f"{e}.json"
        if not fo.exists() or not fn.exists(): print("missing", s, e); bad += 1; continue
        ro = {(r["task"], r["estimator"]): r for r in json.loads(fo.read_text())["rows"]}
        rn = {(r["task"], r["estimator"]): r for r in json.loads(fn.read_text())["rows"]}
        for k, r in ro.items():
            n += 1
            for f in ("value", "ci_lo", "ci_hi", "p"):
                a, b = r.get(f), rn[k].get(f)
                if a is None and b is None: continue
                if a is None or b is None or abs(a - b) > 1e-12:
                    if not (f == "p" and a is None): print("DIFF", s, e, k, f, a, b); bad += 1
    print(f"checked {n} rows of estimators (i)/(ii); {bad} differences")


if __name__ == "__main__":
    worker() if sys.argv[1] == "worker" else check_identical()
