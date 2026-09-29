#!/usr/bin/env python3
"""
plot_module_matrix.py

Baseline module detection heatmap: which KEGG modules does each tool call
on the complete (0% removal) genome?

Layout (one figure per organism):
  Rows:    KEGG modules present in ≥1 tool at 0% removal (majority vote ≥2/3 reps)
           sorted by number of tools detecting them (most detected on top)
  Columns: tools (one column per tool)
  Color:   tool color = present (≥2/3 replicates detect it), white = absent

This is a static snapshot of tool agreement on the complete genome —
NOT a removal-level recovery plot.

Reads: MODULE_MATRICES/ (built by parse_module_matrices.py)
Writes:
  PLOTS/module_matrix_{organism}.{pdf,png}
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = Path("/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP")
MOD_DIR  = BASE_DIR / "MODULE_MATRICES"
OUT_DIR  = BASE_DIR / "PLOTS"
OUT_DIR.mkdir(exist_ok=True)

meta      = pd.read_csv(MOD_DIR / "sample_metadata.csv", index_col="sample")
organisms = sorted(o for o in meta["organism"].unique()
                   if o != "Acinetobacter_baumannii")

# ---------------------------------------------------------------------------
# Tool configuration — saturated palette, same order as recovery plot
# ---------------------------------------------------------------------------
TOOLS_CFG = [
    ("anvio",           "anvi'o",          "#00996b"),
    ("BLIMMP_prior",    "BLIMMP prior",    "#a07bff"),
    ("BLIMMP_nosub",    "BLIMMP nosub",    "#7c3aed"),
    ("BLIMMP",          "BLIMMP+sub",      "#4338ca"),
    ("MetaPathPredict", "MetaPathPredict", "#f0347a"),
    ("METABOLIC",       "METABOLIC",       "#d97706"),
    ("DRAM",            "DRAM",            "#1d6fdb"),
    ("KEMET",           "KEMET",           "#dd1515"),
]

ORG_LABELS = {
    "MED4":                          r"$P.\ marinus$ MED4",
    "SS120":                         r"$P.\ marinus$ SS120",
    "MIT9313":                       r"$P.\ marinus$ MIT9313",
    "Pseudomonas_fluorescens_SBW25": r"$P.\ fluorescens$ SBW25",
}
def org_label(org):
    return ORG_LABELS.get(org, org.replace("_", " "))

# Load matrices
matrices = {}
for stem, label, color in TOOLS_CFG:
    p = MOD_DIR / f"{stem}_module_binary.csv"
    if p.exists():
        matrices[stem] = pd.read_csv(p, index_col="sample")
    else:
        print(f"  WARNING: {p.name} not found — skipping {stem}")

available  = [k for k, _, _ in TOOLS_CFG if k in matrices]
TOOL_COLOR = {k: c for k, _, c in TOOLS_CFG}
TOOL_LABEL = {k: l for k, l, _ in TOOLS_CFG}

MIN_REPS = 2   # majority vote: ≥2 of 3 replicates

# ---------------------------------------------------------------------------
# One figure per organism
# ---------------------------------------------------------------------------
for org in organisms:
    # For each tool: which modules are detected at 0% removal (majority vote)?
    tool_detected = {}
    for key in available:
        mat = matrices[key]
        rows_0 = meta[(meta["organism"] == org) & (meta["removal_pct"] == 0)]
        samples_0 = [s for s in rows_0.index if s in mat.index]
        if not samples_0:
            tool_detected[key] = set()
            continue
        sub = mat.loc[samples_0]
        tool_detected[key] = set(sub.columns[sub.sum(axis=0) >= MIN_REPS])

    # Union of all detected modules across tools
    union_modules = set()
    for key in available:
        union_modules |= tool_detected[key]

    if not union_modules:
        print(f"  {org}: no modules detected by any tool at 0% — skipping")
        continue

    # Sort: most tools detecting → top; alphabetical tiebreak
    def sort_key(m):
        n = sum(1 for k in available if m in tool_detected.get(k, set()))
        return (-n, m)

    union_modules = sorted(union_modules, key=sort_key)
    n_mods  = len(union_modules)
    n_tools = len(available)

    print(f"\n{org}: {n_mods} modules detected by ≥1 tool at 0% removal")
    for key in available:
        print(f"  {TOOL_LABEL[key]:20s}: {len(tool_detected[key])}")

    # Build binary matrix: modules × tools
    data = np.zeros((n_mods, n_tools), dtype=np.float32)
    for j, key in enumerate(available):
        for i, mod in enumerate(union_modules):
            data[i, j] = 1.0 if mod in tool_detected[key] else 0.0

    # Build a custom colormap: each column gets its own color for "present"
    # We'll draw each column separately as an image strip with its own cmap.

    fig_h = max(4.0, n_mods * 0.09 + 2.5)
    fig_w = max(6.0, n_tools * 1.0 + 1.5)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor="white")
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # Draw each tool column as a vertical strip
    col_width = 1.0  # normalized data units
    for j, key in enumerate(available):
        col_data = data[:, j].reshape(-1, 1)  # shape (n_mods, 1)
        color = TOOL_COLOR[key]
        cmap = mcolors.LinearSegmentedColormap.from_list(
            key, ["#f5f5f5", color], N=2
        )
        ax.imshow(col_data,
                  aspect="auto",
                  cmap=cmap,
                  vmin=0, vmax=1,
                  interpolation="nearest",
                  origin="upper",
                  extent=(j - 0.5, j + 0.5, n_mods - 0.5, -0.5))

    # Draw thin white grid lines between columns
    for j in range(n_tools - 1):
        ax.axvline(j + 0.5, color="white", lw=1.2)

    # Axes: tool names on x, module names on y
    ax.set_xticks(range(n_tools))
    ax.set_xticklabels(
        [f"{TOOL_LABEL[k]}\n(n={len(tool_detected.get(k, set()))})"
         for k in available],
        fontsize=8, color="#111111"
    )
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")

    ax.set_yticks(range(n_mods))
    ax.set_yticklabels(union_modules, fontsize=4.5, color="#222222")
    ax.set_ylabel("KEGG module", fontsize=8.5, color="#333333")
    ax.set_xlabel("Tool  (n = modules detected at 0% removal, ≥2/3 replicates)",
                  fontsize=8, color="#333333", labelpad=8)
    ax.xaxis.set_label_position("top")

    ax.set_xlim(-0.5, n_tools - 0.5)
    ax.set_ylim(n_mods - 0.5, -0.5)

    ax.tick_params(axis="x", length=0, pad=4)
    ax.tick_params(axis="y", length=0)

    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
        spine.set_color("#bbbbbb")

    # Legend: colored square = present
    legend_handles = [
        Patch(facecolor=TOOL_COLOR[k], edgecolor="none", label=TOOL_LABEL[k])
        for k in available
    ]
    legend_handles.append(Patch(facecolor="#f5f5f5", edgecolor="#bbbbbb", label="absent"))
    fig.legend(handles=legend_handles, loc="lower center",
               ncol=len(available) + 1, frameon=False,
               fontsize=7.5, bbox_to_anchor=(0.5, -0.03))

    fig.suptitle(
        f"Module detection at 0% removal — {org_label(org)}",
        fontsize=11, fontweight="bold", color="#111111", y=1.04
    )

    out_stem = f"module_matrix_{org}"
    for ext in ("pdf", "png"):
        fig.savefig(OUT_DIR / f"{out_stem}.{ext}", dpi=200,
                    bbox_inches="tight", facecolor="white")
    print(f"  Saved: {OUT_DIR}/{out_stem}.{{pdf,png}}")
    plt.close("all")

print("\nDone.")
