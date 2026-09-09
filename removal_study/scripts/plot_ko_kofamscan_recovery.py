#!/usr/bin/env python3
"""
plot_ko_kofamscan_recovery.py

Combined recall & precision vs KofamScan stable KOs at 0% removal.

Ground truth: KofamScan stable KOs at 0% removal (present in ≥2 of 3 replicates).

  Recall    = |tool KOs ∩ kofamscan_stable| / |kofamscan_stable|
              (fraction of KofamScan's set that this tool detects)
  Precision = |tool KOs ∩ kofamscan_stable| / |tool KOs|
              (fraction of tool's calls that are in KofamScan's set)

BLIMMP thresholds (ko_probability): direct > 0.01, neighbor > 0.01, +sub > 0.01.

Layout: 2 rows × 4 columns
  Top row:    Recall    — one panel per organism
  Bottom row: Precision — one panel per organism

Reads: KO_MATRICES/ (from parse_ko_matrices.py)
Writes: PLOTS/ko_kofamscan_recovery.pdf/.png
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = Path("/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP")
MAT_DIR  = BASE_DIR / "KO_MATRICES"
OUT_DIR  = BASE_DIR / "PLOTS"
OUT_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Tool configuration
# ---------------------------------------------------------------------------
TOOLS_CFG = [
    # (csv_stem,       display_name,          color,       linestyle,                 marker)
    ("HMMsearch",     "HMMsearch (E<1e-5)",   "#777777",   (0, (1, 1)),               "x"),
    ("KofamScan",     "KofamScan",            "#eb6834",   "-.",                      "s"),
    ("anvio",         "anvi'o",               "#1baf7a",   ":",                       "^"),
    ("METABOLIC",     "METABOLIC",            "#eda100",   (0, (3, 1, 1, 1)),         "D"),
    ("DRAM",          "DRAM",                 "#e87ba4",   (0, (3, 1, 1, 1, 1, 1)),  "P"),
    ("KEMET",         "KEMET",                "#008300",   (0, (8, 2)),               "h"),
    ("BLIMMP_prior",  "BLIMMP (direct)",      "#93b8e8",   (0, (5, 1)),               "o"),
    ("BLIMMP_nosub",  "BLIMMP (neighbor)",    "#2a78d6",   "-",                       "o"),
    ("BLIMMP",        "BLIMMP+sub",           "#4a3aa7",   "--",                      "o"),
]
STEMS       = [t[0] for t in TOOLS_CFG]
DISPLAY     = {t[0]: t[1] for t in TOOLS_CFG}
TOOL_COLOR  = {t[0]: t[2] for t in TOOLS_CFG}
TOOL_LS     = {t[0]: t[3] for t in TOOLS_CFG}
TOOL_MARKER = {t[0]: t[4] for t in TOOLS_CFG}

GROUND_TRUTH_STEM = "KofamScan"

# ---------------------------------------------------------------------------
# Load metadata and binary matrices
# ---------------------------------------------------------------------------
meta = pd.read_csv(MAT_DIR / "sample_metadata.csv", index_col="sample")
organisms      = sorted(meta["organism"].unique())
removal_levels = sorted(meta["removal_pct"].unique())
print(f"Organisms: {organisms}")
print(f"Removal levels: {removal_levels}")

matrices = {}
for stem in STEMS:
    path = MAT_DIR / f"{stem}_binary.csv"
    if not path.exists():
        print(f"  WARNING: {path} not found — skipping {stem}")
        continue
    matrices[stem] = pd.read_csv(path, index_col="sample")
    print(f"  Loaded {stem}: {matrices[stem].shape}")

if GROUND_TRUTH_STEM not in matrices:
    raise FileNotFoundError(
        f"Ground truth matrix '{GROUND_TRUTH_STEM}_binary.csv' not found in {MAT_DIR}"
    )

available_stems = [s for s in STEMS if s in matrices]
print(f"\n{len(available_stems)} tools loaded: {available_stems}")

# ---------------------------------------------------------------------------
# Step 1: KofamScan stable KO set per organism at 0% removal
# ---------------------------------------------------------------------------
def stable_kos(mat: pd.DataFrame, org: str) -> set:
    ref = meta[(meta["organism"] == org) & (meta["removal_pct"] == 0)]
    samples = [s for s in ref.index if s in mat.index]
    if not samples:
        return set()
    sub = mat.loc[samples]
    return set(sub.sum(axis=0)[sub.sum(axis=0) >= 2].index)


print(f"\nGround truth: {DISPLAY[GROUND_TRUTH_STEM]} stable KOs at 0% removal")
gt_kos = {}
for org in organisms:
    gt_kos[org] = stable_kos(matrices[GROUND_TRUTH_STEM], org)
    print(f"  {org}: {len(gt_kos[org])} KOs")

# ---------------------------------------------------------------------------
# Step 2: compute recall and precision per (tool, organism, removal_pct, replicate)
# ---------------------------------------------------------------------------
rows = []
for stem in available_stems:
    mat = matrices[stem]
    for org in organisms:
        gt = gt_kos[org]
        if not gt:
            continue
        org_meta = meta[meta["organism"] == org]
        for _, row in org_meta.iterrows():
            sample = row.name
            if sample not in mat.index:
                continue
            tool_kos  = set(mat.columns[mat.loc[sample] == 1])
            n_overlap = len(tool_kos & gt)
            recall    = n_overlap / len(gt)
            precision = n_overlap / len(tool_kos) if tool_kos else float("nan")
            rows.append(dict(
                organism=org, tool=stem,
                removal_pct=row["removal_pct"], replicate=row["replicate"],
                recall=recall, precision=precision,
            ))

pr_df = pd.DataFrame(rows)
summary = (
    pr_df.groupby(["organism", "tool", "removal_pct"])
    .agg(
        rec_mean=("recall",    "mean"), rec_sd=("recall",    "std"),
        pre_mean=("precision", "mean"), pre_sd=("precision", "std"),
    )
    .reset_index()
)

print("\nRecall and Precision at 0% removal (KofamScan ground truth):")
ref = summary[summary["removal_pct"] == 0]
for metric, label in [("rec_mean", "Recall"), ("pre_mean", "Precision")]:
    tbl = ref.pivot_table(index="tool", columns="organism", values=metric)
    print(f"\n  {label}:")
    print(tbl.applymap(lambda x: f"{x:.1%}" if pd.notna(x) else "n/a").to_string())

# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------
SURFACE  = "#fcfcfb"
TEXT_PRI = "#0b0b0b"
TEXT_SEC = "#52514e"
GRID_COL = "#e5e4de"

def style_ax(ax):
    ax.set_facecolor(SURFACE)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GRID_COL)
    ax.tick_params(colors=TEXT_SEC, labelsize=7.5)
    ax.xaxis.label.set_color(TEXT_SEC)
    ax.yaxis.label.set_color(TEXT_SEC)
    ax.grid(axis="y", color=GRID_COL, linewidth=0.6, zorder=0)
    ax.set_xticks(removal_levels)
    ax.set_xticklabels([f"{p}%" for p in removal_levels],
                       rotation=45, ha="right", fontsize=7)
    ax.set_ylim(-0.02, 1.08)
    ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1, decimals=0))

legend_handles = [
    Line2D([0], [0],
           color=TOOL_COLOR[s], linestyle=TOOL_LS[s],
           marker=TOOL_MARKER[s], linewidth=1.8, markersize=6,
           label=DISPLAY[s])
    for s in STEMS if s in matrices
]

# ---------------------------------------------------------------------------
# Figure: 2 rows (recall / precision) × 4 columns (organisms)
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 4, figsize=(18, 9), facecolor=SURFACE)

METRIC_CFG = [
    ("rec_mean", "rec_sd",
     "Recall\n(KofamScan KOs recovered / KofamScan stable, mean ± SD)"),
    ("pre_mean", "pre_sd",
     "Precision\n(KofamScan KOs recovered / tool KOs, mean ± SD)"),
]
ROW_LABELS = ["Recall", "Precision"]

for row_i, (mean_col, sd_col, ylabel) in enumerate(METRIC_CFG):
    for col_i, org in enumerate(organisms):
        ax = axes[row_i, col_i]
        style_ax(ax)
        n_gt   = len(gt_kos[org])
        org_df = summary[summary["organism"] == org]

        for stem in available_stems:
            t = org_df[org_df["tool"] == stem].sort_values("removal_pct")
            if t.empty:
                continue
            means = t[mean_col].clip(0, 1)
            sds   = t[sd_col].fillna(0)
            ax.plot(t["removal_pct"], means,
                    color=TOOL_COLOR[stem], linestyle=TOOL_LS[stem],
                    marker=TOOL_MARKER[stem], linewidth=1.8, markersize=5, zorder=3)
            ax.fill_between(
                t["removal_pct"],
                (means - sds).clip(lower=0),
                (means + sds).clip(upper=1),
                color=TOOL_COLOR[stem], alpha=0.12, zorder=2,
            )

        # Organism title only on first row
        if row_i == 0:
            ax.set_title(org.replace("_", " "), color=TEXT_PRI, fontsize=9,
                         fontweight="semibold", style="italic")
            ax.text(0.97, 0.04, f"KofamScan n={n_gt:,}",
                    transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=7, color=TEXT_SEC)

        # Y-axis label only on left-most column
        if col_i == 0:
            ax.set_ylabel(ylabel, fontsize=8, color=TEXT_SEC)

        ax.set_xlabel("Removal level", fontsize=8, color=TEXT_SEC)

        # Row badge (top-left corner of each panel)
        ax.text(0.03, 0.97, ROW_LABELS[row_i],
                transform=ax.transAxes, ha="left", va="top",
                fontsize=8, fontweight="bold", color=TEXT_PRI,
                bbox=dict(boxstyle="round,pad=0.2", fc=SURFACE, ec=GRID_COL, lw=0.8))

fig.suptitle(
    "Recall and Precision vs KofamScan stable KOs at 0% removal",
    color=TEXT_PRI, fontsize=13, fontweight="bold", y=1.01,
)

# BLIMMP threshold footnote
fig.text(
    0.5, -0.02,
    "BLIMMP thresholds (ko_probability): direct > 0.01, neighbor > 0.01, +sub > 0.01   "
    "|   Ground truth = KofamScan stable KOs at 0% removal (present in ≥2 of 3 replicates)",
    ha="center", va="top", fontsize=7.5, color=TEXT_SEC, style="italic",
)

fig.legend(handles=legend_handles, loc="lower center",
           ncol=len(legend_handles), frameon=False, fontsize=8,
           bbox_to_anchor=(0.5, -0.07), labelcolor=TEXT_SEC)
fig.tight_layout(rect=[0, 0.08, 1, 0.99])

fig.savefig(OUT_DIR / "ko_kofamscan_recovery.png", dpi=300,
            bbox_inches="tight", facecolor=SURFACE)
fig.savefig(OUT_DIR / "ko_kofamscan_recovery.pdf",
            bbox_inches="tight", facecolor=SURFACE)
print(f"\nSaved: {OUT_DIR}/ko_kofamscan_recovery.{{png,pdf}}")
plt.close("all")
