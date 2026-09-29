#!/usr/bin/env python3
"""
check_blimmp_overlap.py

For each organism: how many of the modules BLIMMP calls at 0% removal
are still called at 90% removal?

Checks all four BLIMMP variants.
"""

from pathlib import Path
import pandas as pd

BASE_DIR = Path("/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP")
MOD_DIR  = BASE_DIR / "MODULE_MATRICES"

meta = pd.read_csv(MOD_DIR / "sample_metadata.csv", index_col="sample")
organisms = sorted(o for o in meta["organism"].unique()
                   if o != "Acinetobacter_baumannii")

MIN_REPS = 2

BLIMMP_VARIANTS = [
    ("BLIMMP_prior",  "BLIMMP prior"),
    ("BLIMMP_nosub",  "BLIMMP nosub"),
    ("BLIMMP",        "BLIMMP+sub"),
    ("BLIMMP_raw",    "BLIMMP raw"),
]

def majority_set(mat, org, pct):
    rows = meta[(meta["organism"] == org) & (meta["removal_pct"] == pct)]
    samples = [s for s in rows.index if s in mat.index]
    if not samples:
        return set()
    sub = mat.loc[samples]
    return set(sub.columns[sub.sum(axis=0) >= MIN_REPS])

for stem, label in BLIMMP_VARIANTS:
    p = MOD_DIR / f"{stem}_module_binary.csv"
    if not p.exists():
        print(f"\n{label}: matrix not found — run parse_module_matrices.py first")
        continue
    mat = pd.read_csv(p, index_col="sample")
    print(f"\n{'='*60}")
    print(f"{label}")
    print(f"{'='*60}")
    for org in organisms:
        at_0  = majority_set(mat, org, 0)
        at_90 = majority_set(mat, org, 90)
        overlap = at_0 & at_90
        missing = at_0 - at_90   # present at 0% but gone by 90%
        new_at_90 = at_90 - at_0 # new modules appearing at 90% (FPs)
        print(f"\n  {org}")
        print(f"    Modules at   0% removal : {len(at_0)}")
        print(f"    Modules at  90% removal : {len(at_90)}")
        print(f"    Still present at 90%    : {len(overlap)}  ({len(overlap)/len(at_0)*100:.1f}% of baseline)" if at_0 else "    (no baseline)")
        print(f"    Lost by 90%             : {len(missing)}")
        print(f"    New at 90% (not in 0%)  : {len(new_at_90)}")
        if missing:
            print(f"    Modules LOST: {sorted(missing)}")

print("\nDone.")
