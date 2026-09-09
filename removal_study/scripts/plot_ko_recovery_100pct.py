#!/usr/bin/env python3
"""
plot_ko_recovery.py

Self-recovery: for each tool, what fraction of its own stable KO set at 0%
removal does it still detect at higher removal levels?

Ground truth (per tool, per organism):
    stable KOs = KOs present in ≥2 of 3 replicates at 0% removal.

For every (tool, organism, removal_pct, replicate):
    recovery = |tool KOs for sample ∩ tool_stable| / |tool_stable|

Mean ± SD across replicates is plotted per organism.

This answers: "How well does each tool maintain its own annotation as gene
content is progressively removed?" — a tool-neutral measure of robustness
that does not depend on KofamScan as an external reference.

BLIMMP thresholds (ko_probability): direct > 0.01, neighbor > 0.01, +sub > 0.01.

Reads: KO_MATRICES/ (from parse_ko_matrices.py)
Writes: PLOTS/ko_recovery.pdf/.png
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

available_stems = [s for s in STEMS if s in matrices]
print(f"\n{len(available_stems)} tools loaded: {available_stems}")

# ---------------------------------------------------------------------------
# Step 1: compute each tool's own stable KO set per organism at 0% removal
#   stable = present in ≥2 of 3 replicates
# ---------------------------------------------------------------------------
def stable_kos(mat: pd.DataFrame, org: str) -> set:
    """KOs present in ≥2 of 3 replicates at 0% removal for this organism."""
    ref = meta[(meta["organism"] == org) & (meta["removal_pct"] == 0)]
    samples = [s for s in ref.index if s in mat.index]
    if not samples:
        return set()
    sub = mat.loc[samples]
    rep_sum = sub.sum(axis=0)
    return set(rep_sum[rep_sum >= 2].index)


print("\nStable KO sets at 0% removal:")
tool_stable = {}   # {(stem, organism): set_of_KOs}
for stem in available_stems:
    for org in organisms:
        tool_stable[(stem, org)] = stable_kos(matrices[stem], org)
    n_by_org = {org: len(tool_stable[(stem, org)]) for org in organisms}
    print(f"  {stem}: " + "  ".join(f"{org}={n}" for org, n in n_by_org.items()))

# BLIMMP variants share the BLIMMP+sub stable set as a common denominator so
# that BLIMMP_prior / BLIMMP_nosub / BLIMMP+sub are directly comparable.
# Without this, each variant's smaller own stable set makes it appear to
# "recover" more even though it detects fewer absolute KOs.
BLIMMP_STEMS = {"BLIMMP_prior", "BLIMMP_nosub", "BLIMMP"}
if "BLIMMP" in matrices:
    print("\n  Overriding BLIMMP_prior / BLIMMP_nosub denominators → BLIMMP+sub stable set")
    for org in organisms:
        blimmp_gt = tool_stable[("BLIMMP", org)]
        for stem in ("BLIMMP_prior", "BLIMMP_nosub"):
            if stem in matrices:
                tool_stable[(stem, org)] = blimmp_gt

# ---------------------------------------------------------------------------
# Step 2: compute self-recovery per (tool, organism, removal_pct, replicate)
# ---------------------------------------------------------------------------
rows = []
for stem in available_stems:
    mat = matrices[stem]
    for org in organisms:
        gt = tool_stable[(stem, org)]
        if not gt:
            continue
        org_meta = meta[meta["organism"] == org]
        for _, row in org_meta.iterrows():
            sample = row.name
            if sample not in mat.index:
                continue
            detected = set(mat.columns[mat.loc[sample] == 1])
            rows.append(dict(
                organism=org, tool=stem,
                removal_pct=row["removal_pct"], replicate=row["replicate"],
                recovery=len(detected & gt) / len(gt),
                n_stable=len(gt),
            ))

rec_df = pd.DataFrame(rows)
summary = (
    rec_df.groupby(["organism", "tool", "removal_pct"])
    .agg(rec_mean=("recovery", "mean"), rec_sd=("recovery", "std"))
    .reset_index()
)

print("\nSelf-recovery at 0% removal (should be ~100% for all tools):")
ref_tbl = summary[summary["removal_pct"] == 0][["organism", "tool", "rec_mean"]]
print(ref_tbl.pivot(index="tool", columns="organism", values="rec_mean")
      .applymap(lambda x: f"{x:.1%}").to_string())

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
# Figure: 2×2 panels, one per organism
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(12, 9), facecolor=SURFACE)
fig.suptitle(
    "Self-recovery: each tool vs its own stable KOs at 0% removal\n"
    "(BLIMMP variants share BLIMMP+sub stable set as common denominator)",
    color=TEXT_PRI, fontsize=12, fontweight="bold", y=1.02,
)

for ax, org in zip(axes.flatten(), organisms):
    style_ax(ax)
    org_df = summary[summary["organism"] == org]

    for stem in available_stems:
        t = org_df[org_df["tool"] == stem].sort_values("removal_pct")
        if t.empty:
            continue
        means = t["rec_mean"].clip(0, 1)
        sds   = t["rec_sd"].fillna(0)
        n_stable = tool_stable[(stem, org)]
        ax.plot(t["removal_pct"], means,
                color=TOOL_COLOR[stem], linestyle=TOOL_LS[stem],
                marker=TOOL_MARKER[stem], linewidth=1.8, markersize=5, zorder=3)
        ax.fill_between(
            t["removal_pct"],
            (means - sds).clip(lower=0),
            (means + sds).clip(upper=1),
            color=TOOL_COLOR[stem], alpha=0.12, zorder=2,
        )

    ax.set_title(org.replace("_", " "), color=TEXT_PRI, fontsize=10,
                 fontweight="semibold", style="italic")
    ax.set_xlabel("Removal level", fontsize=8, color=TEXT_SEC)
    ax.set_ylabel("Own stable KOs recovered\n(mean ± SD, n=3)", fontsize=8, color=TEXT_SEC)

# BLIMMP threshold footnote
fig.text(
    0.5, -0.01,
    "Each tool's ground truth = own stable KOs at 0% removal (≥2 of 3 replicates).  "
    "BLIMMP_prior & BLIMMP_nosub use BLIMMP+sub's stable set as denominator for direct comparison.   "
    "|   BLIMMP thresholds (ko_probability): direct > 0.01, neighbor > 0.01, +sub > 0.01",
    ha="center", va="top", fontsize=7, color=TEXT_SEC, style="italic",
)

fig.legend(handles=legend_handles, loc="lower center",
           ncol=len(legend_handles), frameon=False, fontsize=9,
           bbox_to_anchor=(0.5, -0.06), labelcolor=TEXT_SEC)
fig.tight_layout(rect=[0, 0.07, 1, 0.99])

fig.savefig(OUT_DIR / "ko_recovery.png", dpi=300,
            bbox_inches="tight", facecolor=SURFACE)
fig.savefig(OUT_DIR / "ko_recovery.pdf",
            bbox_inches="tight", facecolor=SURFACE)
print(f"\nSaved: {OUT_DIR}/ko_recovery.{{png,pdf}}")
plt.close("all")
