#!/usr/bin/env python3

import sys
import pandas as pd
from pathlib import Path

CSVS_DIR = Path(".")

CSV_FILES = {
    "anvio":     CSVS_DIR / "blimmp_vs_anvio_scatter_data.csv",
    "metabolic": CSVS_DIR / "blimmp_vs_metabolic_scatter_data.csv",
    "kemet":     CSVS_DIR / "blimmp_vs_kemet_scatter_data.csv",
    "mpp":       CSVS_DIR / "blimmp_vs_mpp_violin_data.csv",
}

ORGANISMS = ["MED4", "SS120", "MIT9313", "Pseudomonas_fluorescens_SBW25"]
TOOLS     = ["anvio", "metabolic", "kemet", "mpp"]

# Category values that mean one tool disagrees with the other
DISCORDANT_CATS = {"blimmp_only", "anvio_only", "metabolic_only", "kemet_only", "mpp_only"}


frames = []
for tool, path in CSV_FILES.items():
    if not path.exists():
        print(f"[MISSING] {path.name} — skipping {tool}")
        continue
    df = pd.read_csv(path, usecols=["organism", "module", "category"])
    df["tool"] = tool
    frames.append(df)
    print(f"  Loaded {path.name}: {len(df):,} rows")

if not frames:
    sys.exit("No CSV files found. Run the scatter/violin scripts first, or set CSVS_DIR.")

combined = pd.concat(frames, ignore_index=True)
print(f"\nTotal rows loaded: {len(combined):,}")
print(f"Unique modules:    {combined['module'].nunique():,}")


combined["col"] = combined["organism"] + "_" + combined["tool"]

wide = (combined
        .pivot_table(index="module", columns="col",
                     values="category", aggfunc="first")
        .reset_index())

cat_cols = [c for c in wide.columns if c != "module"]
mask = wide[cat_cols].isin(DISCORDANT_CATS).any(axis=1)
filtered = wide[mask].copy()
print(f"\nDiscordant modules (blimmp_only or tool_only in ≥1 cell): {len(filtered):,}")
print(f"Consensus modules dropped (all 'both' or 'neither'):       {len(wide) - len(filtered):,}")

ordered_cols = ["module"]
for org in ORGANISMS:
    for tool in TOOLS:
        col = f"{org}_{tool}"
        if col in filtered.columns:
            ordered_cols.append(col)

# Include any unexpected columns at the end
extra = [c for c in filtered.columns if c not in ordered_cols]
filtered = filtered[ordered_cols + extra]

# Sort rows by module ID
filtered = filtered.sort_values("module").reset_index(drop=True)

def summarise(row):
    vals = [v for v in row[cat_cols] if pd.notna(v)]
    cats = set(vals)
    if "blimmp_only" in cats and any(c.endswith("_only") and c != "blimmp_only" for c in cats):
        return "mixed"
    if "blimmp_only" in cats:
        return "blimmp_only"
    if any(c.endswith("_only") and c != "blimmp_only" for c in cats):
        return "tool_only"
    return "other"

filtered["overall_pattern"] = filtered.apply(summarise, axis=1)


out = CSVS_DIR / "modules_disagreement_comparison.csv"
filtered.to_csv(out, index=False)
print(f"\nSaved → {out}")

# Quick breakdown
print("\nOverall pattern breakdown:")
print(filtered["overall_pattern"].value_counts().to_string())
