#!/usr/bin/env python3
"""
plot_ko_quadrant.py

At 0% removal, for each KO in the bacterial module universe:
  X = BLIMMP's max ko_probability (posterior) across 3 replicates at 0% removal.
      0 if the KO never appears in BLIMMP output.
  Y = number of non-BLIMMP tools (HMMsearch, KofamScan, anvi'o, METABOLIC,
      DRAM, KEMET) that stably detect this KO (≥2 of 3 replicates at 0%).
      Integer 0–6; jittered slightly for legibility.

Four quadrants divided at X=0.5 (BLIMMP threshold) and Y=0.5 (any other tool):
  Top-right:    Both BLIMMP and ≥1 other tool detect it     → Consensus
  Top-left:     ≥1 other tool detects it, BLIMMP prob < 0.5 → BLIMMP misses
  Bottom-right: BLIMMP prob ≥ 0.5, no other tool detects it → BLIMMP novel
  Bottom-left:  Neither BLIMMP nor any other tool detects it → Absent

KOs absent from all sources (X=0, Y=0) are excluded.

Reads:  KO_MATRICES/ (binary matrices for other tools)
        BLIMMP_SPLICED/ (raw __BLIMMP_dk.csv and *_substituted_dk.csv)
        BLIMMP_8Sep2026/BLIMMP/module_ko_reaction.json
Writes: PLOTS/ko_quadrant.pdf/.png
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = Path("/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP")
MAT_DIR  = BASE_DIR / "KO_MATRICES"
OUT_DIR  = BASE_DIR / "PLOTS"
OUT_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Bacterial module KO universe
# ---------------------------------------------------------------------------
MODULE_KO_JSON = BASE_DIR / "BLIMMP_8Sep2026" / "BLIMMP" / "module_ko_reaction.json"
with open(MODULE_KO_JSON) as f:
    module_data = json.load(f)
BACTERIAL_MODULE_KOS = set()
for ko_dict in module_data.values():
    BACTERIAL_MODULE_KOS.update(ko_dict.keys())
print(f"Module KO universe: {len(BACTERIAL_MODULE_KOS):,} KOs")

# ---------------------------------------------------------------------------
# Load metadata and binary matrices for non-BLIMMP tools
# ---------------------------------------------------------------------------
meta      = pd.read_csv(MAT_DIR / "sample_metadata.csv", index_col="sample")
organisms = sorted(meta["organism"].unique())

OTHER_STEMS = ["HMMsearch", "KofamScan", "anvio", "METABOLIC", "DRAM", "KEMET"]
N_OTHER     = len(OTHER_STEMS)

matrices = {}
for stem in OTHER_STEMS:
    p = MAT_DIR / f"{stem}_binary.csv"
    if p.exists():
        matrices[stem] = pd.read_csv(p, index_col="sample")
        print(f"  Loaded {stem}: {matrices[stem].shape}")
    else:
        print(f"  WARNING: {stem}_binary.csv not found — skipping")

available_other = [s for s in OTHER_STEMS if s in matrices]
N_AVAIL = len(available_other)

# ---------------------------------------------------------------------------
# Y-axis: count of other tools that stably detect each KO at 0% removal
# ---------------------------------------------------------------------------
def stable_kos(mat: pd.DataFrame, org: str) -> set:
    ref     = meta[(meta["organism"] == org) & (meta["removal_pct"] == 0)]
    samples = [s for s in ref.index if s in mat.index]
    if not samples:
        return set()
    sub = mat.loc[samples]
    return set(sub.sum(axis=0)[sub.sum(axis=0) >= 2].index)

print("\nComputing stable KOs for other tools...")
other_tool_counts = {}   # org → {ko: n_tools}
for org in organisms:
    ko_counts: dict = {}
    for stem in available_other:
        for ko in stable_kos(matrices[stem], org):
            ko_counts[ko] = ko_counts.get(ko, 0) + 1
    other_tool_counts[org] = ko_counts
    n_any = sum(1 for v in ko_counts.values() if v >= 1)
    print(f"  {org}: {n_any} KOs found by ≥1 other tool")

# ---------------------------------------------------------------------------
# X-axis: BLIMMP max ko_probability across 3 replicates at 0% removal
# Reads raw __BLIMMP_dk.csv and *_substituted_dk.csv (not the binary matrix).
# ---------------------------------------------------------------------------
def read_blimmp_probs(org: str) -> dict:
    """Return {ko_id: max_ko_probability} for this organism at 0% removal."""
    ref = meta[(meta["organism"] == org) & (meta["removal_pct"] == 0)]
    ko_probs: dict = {}
    for sample in ref.index:
        # Pass-1 posterior (direct + neighbor)
        p = BASE_DIR / "BLIMMP_SPLICED" / sample / f"{sample}__BLIMMP_dk.csv"
        if p.exists():
            df = pd.read_csv(p, usecols=["ko_id", "ko_probability"])
            df["ko_id"] = df["ko_id"].str.strip()
            for _, row in df.iterrows():
                if pd.notna(row["ko_probability"]):
                    ko   = row["ko_id"]
                    prob = float(row["ko_probability"])
                    ko_probs[ko] = max(ko_probs.get(ko, 0.0), prob)
        # Pass-2 substituted
        d = BASE_DIR / "BLIMMP_SPLICED" / sample
        if d.exists():
            for sub_file in d.glob("*_substituted_dk.csv"):
                df2 = pd.read_csv(sub_file, usecols=["ko_id", "ko_probability"])
                df2["ko_id"] = df2["ko_id"].str.strip()
                for _, row in df2.iterrows():
                    if pd.notna(row["ko_probability"]):
                        ko   = row["ko_id"]
                        prob = float(row["ko_probability"])
                        ko_probs[ko] = max(ko_probs.get(ko, 0.0), prob)
    return ko_probs

print("\nReading raw BLIMMP probabilities at 0% removal...")
blimmp_probs = {}
for org in organisms:
    blimmp_probs[org] = read_blimmp_probs(org)
    n_above = sum(1 for v in blimmp_probs[org].values() if v > 0.5)
    print(f"  {org}: {len(blimmp_probs[org])} KOs with any probability, "
          f"{n_above} above 0.5")

# ---------------------------------------------------------------------------
# Colors and style
# ---------------------------------------------------------------------------
Q_COLOR = {
    "both":         "#2a78d6",   # blue    — consensus
    "blimmp_only":  "#4a3aa7",   # violet  — BLIMMP novel
    "others_only":  "#eb6834",   # orange  — BLIMMP misses
    "neither":      "#c8c7c2",   # gray    — absent
}
Q_LABEL = {
    "both":        "Both (consensus)",
    "blimmp_only": "BLIMMP only (novel)",
    "others_only": "Other tools only (BLIMMP misses)",
    "neither":     "Neither detected",
}

SURFACE   = "#fcfcfb"
TEXT_PRI  = "#0b0b0b"
TEXT_SEC  = "#52514e"
GRID_COL  = "#e5e4de"
BLIMMP_T  = 0.5
rng = np.random.default_rng(42)

# ---------------------------------------------------------------------------
# Plot: 1 row × 4 columns (one per organism)
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 4, figsize=(18, 5), facecolor=SURFACE)
fig.suptitle(
    "BLIMMP probability vs. tool consensus at 0% removal\n"
    "(each point = one KO in the bacterial module universe)",
    color=TEXT_PRI, fontsize=12, fontweight="bold", y=1.03,
)

for ax, org in zip(axes, organisms):
    ax.set_facecolor(SURFACE)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GRID_COL)
    ax.tick_params(colors=TEXT_SEC, labelsize=7.5)
    ax.xaxis.label.set_color(TEXT_SEC)
    ax.yaxis.label.set_color(TEXT_SEC)
    ax.grid(color=GRID_COL, linewidth=0.5, zorder=0)

    probs  = blimmp_probs[org]
    counts = other_tool_counts[org]

    # Build arrays — skip KOs absent from everything
    xs, ys, qs = [], [], []
    for ko in BACTERIAL_MODULE_KOS:
        x = probs.get(ko, 0.0)
        y = counts.get(ko, 0)
        if x == 0.0 and y == 0:
            continue
        xs.append(x)
        ys.append(y)
        if   x >= BLIMMP_T and y >= 1: qs.append("both")
        elif x >= BLIMMP_T and y == 0: qs.append("blimmp_only")
        elif x <  BLIMMP_T and y >= 1: qs.append("others_only")
        else:                           qs.append("neither")

    xs = np.array(xs)
    ys = np.array(ys, dtype=float)
    ys_j = ys + rng.uniform(-0.25, 0.25, size=len(ys))   # Y jitter

    colors = [Q_COLOR[q] for q in qs]
    ax.scatter(xs, ys_j, c=colors, s=6, alpha=0.55, linewidths=0, zorder=3)

    # Quadrant dividers
    ax.axvline(BLIMMP_T, color=TEXT_SEC, linewidth=0.9, linestyle="--", zorder=4)
    ax.axhline(0.5,      color=TEXT_SEC, linewidth=0.9, linestyle="--", zorder=4)

    # Quadrant KO counts
    n = {q: qs.count(q) for q in Q_COLOR}
    kw = dict(fontsize=7, fontweight="bold", transform=ax.transData)
    ax.text(BLIMMP_T + 0.02, N_AVAIL - 0.15, f"n={n['both']:,}",
            color=Q_COLOR["both"], ha="left", va="top", **kw)
    ax.text(BLIMMP_T - 0.02, N_AVAIL - 0.15, f"n={n['others_only']:,}",
            color=Q_COLOR["others_only"], ha="right", va="top", **kw)
    ax.text(BLIMMP_T + 0.02, 0.3, f"n={n['blimmp_only']:,}",
            color=Q_COLOR["blimmp_only"], ha="left", va="top", **kw)
    ax.text(BLIMMP_T - 0.02, 0.3, f"n={n['neither']:,}",
            color=Q_COLOR["neither"], ha="right", va="top", **kw)

    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.6, N_AVAIL + 0.6)
    ax.set_yticks(range(N_AVAIL + 1))
    ax.set_yticklabels([str(i) for i in range(N_AVAIL + 1)], fontsize=7)
    ax.set_xlabel("BLIMMP ko_probability\n(max across 3 replicates at 0% removal)",
                  fontsize=8, color=TEXT_SEC)
    if ax is axes[0]:
        ax.set_ylabel(
            f"Other tools detecting KO stably (of {N_AVAIL})\n"
            "HMMsearch · KofamScan · anvi'o · METABOLIC · DRAM · KEMET",
            fontsize=7.5, color=TEXT_SEC,
        )
    ax.set_title(org.replace("_", " "), color=TEXT_PRI, fontsize=9,
                 fontweight="semibold", style="italic")

legend_handles = [
    Line2D([0], [0], marker="o", color="none",
           markerfacecolor=Q_COLOR[q], markersize=7, label=Q_LABEL[q])
    for q in ["both", "blimmp_only", "others_only", "neither"]
]
fig.legend(handles=legend_handles, loc="lower center", ncol=4,
           frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, -0.08),
           labelcolor=TEXT_SEC)

fig.text(
    0.5, -0.03,
    "BLIMMP threshold = 0.5 (vertical dashed line) | "
    "Horizontal dashed line: below = no other tool detects it; above = ≥1 other tool detects it stably",
    ha="center", va="top", fontsize=7, color=TEXT_SEC, style="italic",
)

fig.tight_layout(rect=[0, 0.10, 1, 1.0])
fig.savefig(OUT_DIR / "ko_quadrant.png", dpi=300,
            bbox_inches="tight", facecolor=SURFACE)
fig.savefig(OUT_DIR / "ko_quadrant.pdf",
            bbox_inches="tight", facecolor=SURFACE)
print(f"\nSaved: {OUT_DIR}/ko_quadrant.{{png,pdf}}")
plt.close("all")
