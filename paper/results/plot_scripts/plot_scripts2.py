#!/usr/bin/env python3
"""
plot_blimmp_vs_anvio_scatter.py

BLIMMP module_confidence (x) vs anvi'o pathwise_module_completeness (y)
at 0% removal, majority-voted across 3 replicates.
One panel per genome, four panels in a row.

Highlighted modules are drawn on top of every panel with a labelled
annotation so their position relative to the agreement quadrants is
immediately visible. Marker shape and color for each highlighted module
are assigned automatically from fixed palettes below, in the order the
module IDs are listed in HIGHLIGHT_MODULE_IDS, so adding or removing a
module ID never requires manually picking a new color/marker pair.

Palette (validated against light surface #fcfcfb):
  Both present  : #1e7a40  dark green
  BLIMMP only   : #1a5f9e  dark blue
  anvi'o only   : #7c3aed  violet
  Neither       : #111111  black (de-emphasis)
  Highlighted modules: auto-assigned, see HIGHLIGHT_COLOR_PALETTE /
                       HIGHLIGHT_MARKER_PALETTE below
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

BASE_DIR = Path("/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP")

ORGANISMS = [
    "MED4",
    "SS120",
    "MIT9313",
    "Pseudomonas_fluorescens_SBW25",
]

LABELS = {
    "MED4":                          r"$P.\ marinus$ MED4",
    "SS120":                         r"$P.\ marinus$ SS120",
    "MIT9313":                       r"$P.\ marinus$ MIT9313",
    "Pseudomonas_fluorescens_SBW25": r"$P.\ fluorescens$ SBW25",
}

REMOVAL_PCT   = 0
REPLICATES    = [1, 2, 3]
MIN_REPS      = 2
BLIMMP_THRESH = 0.60
ANVIO_THRESH  = 0.75

# ---------------------------------------------------------------------------
# Highlighted modules
# ---------------------------------------------------------------------------
# List module IDs and their display names here. Color and marker are
# assigned automatically (see build_highlight_styles below) by walking
# HIGHLIGHT_COLOR_PALETTE / HIGHLIGHT_MARKER_PALETTE in order, so the
# palettes only need to be as long as the largest highlight set you expect
# to use across the whole project (extend either list if you add more
# modules than there are palette entries; colors will cycle with a
# warning rather than crash).

HIGHLIGHT_MODULE_NAMES = {

    "M00034": "M00034 (Methionine salvage)",
    "M00126": "M00126 (Tetrahydrofolate biosynthesis)",
    "M00009": "M00009 (Citrate cycle, TCA cycle)",
}
HIGHLIGHT_MODULE_IDS = list(HIGHLIGHT_MODULE_NAMES.keys())

# Distinct, colorblind-conscious palette, validated against light surface
# #fcfcfb. Extend this list (and/or HIGHLIGHT_MARKER_PALETTE) if you add
# more than 8 highlighted modules.
HIGHLIGHT_COLOR_PALETTE = [
    "#f5c800",  # gold
    "#e76f51",  # burnt orange
    "#457b9d",  # steel blue
    "#2a9d8f",  # teal
    "#d62828",  # brick red
    "#6a4c93",  # purple
    "#8ac926",  # yellow-green
    "#ff6f9c",  # pink
]
# Distinct marker shapes; star gets extra render size (see HIGHLIGHT_SIZES).
HIGHLIGHT_MARKER_PALETTE = ["*", "D", "^", "s", "P", "v", "X", "o"]
HIGHLIGHT_SIZES = {"*": 200, "D": 90, "^": 100, "s": 90, "P": 110, "v": 100, "X": 110, "o": 90}


def build_highlight_styles(module_ids, color_palette, marker_palette):
    """Assign a (color, marker) pair to each module ID by cycling the
    given palettes in order. Warns (does not crash) if the number of
    modules exceeds the palette length, since colors/markers would repeat.
    """
    n = len(module_ids)
    if n > len(color_palette) or n > len(marker_palette):
        print(f"  [WARN] {n} highlighted modules but palette has "
              f"{len(color_palette)} colors / {len(marker_palette)} markers "
              f"-- some modules will repeat a color and/or marker.")
    styles = {}
    for i, mod_id in enumerate(module_ids):
        color = color_palette[i % len(color_palette)]
        marker = marker_palette[i % len(marker_palette)]
        label = HIGHLIGHT_MODULE_NAMES[mod_id]
        styles[mod_id] = (color, marker, label)
    return styles


HIGHLIGHT_MODULES = build_highlight_styles(
    HIGHLIGHT_MODULE_IDS, HIGHLIGHT_COLOR_PALETTE, HIGHLIGHT_MARKER_PALETTE
)

C_BOTH    = "#1e7a40"
C_BLIMMP  = "#1a5f9e"
C_ANVIO   = "#7c3aed"
C_NEITHER = "#111111"

SURFACE = "#fcfcfb"
TEXT    = "#2e2e2e"
GRID    = "#e5e4de"


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def sample_name(organism, pct, rep):
    return f"{organism}_{pct}percentremoved_replicate{rep}"


def load_blimmp(organism):
    frames = []
    for rep in REPLICATES:
        sname = sample_name(organism, REMOVAL_PCT, rep)
        path = (BASE_DIR / "BLIMMP_SPLICED" / sname
                / f"{sname}__BLIMMP_substituted_module_probabilities.csv")
        if not path.exists():
            print(f"  [MISSING] BLIMMP: {path.name}")
            continue
        df = pd.read_csv(path, usecols=["module", "module_confidence"])
        df["rep"] = rep
        frames.append(df)

    if not frames:
        return pd.DataFrame(columns=["module", "mean_confidence", "n_present"])

    combined = pd.concat(frames, ignore_index=True)
    return (combined.groupby("module")
            .agg(
                mean_confidence=("module_confidence", "mean"),
                n_present=("module_confidence",
                           lambda x: (x >= BLIMMP_THRESH).sum()),
            )
            .reset_index())


def load_anvio(organism):
    frames = []
    for rep in REPLICATES:
        sname = sample_name(organism, REMOVAL_PCT, rep)
        path = BASE_DIR / "ANVIO_SPLICED" / sname / f"{sname}_modules.txt"
        if not path.exists():
            print(f"  [MISSING] anvi'o: {path.name}")
            continue
        df = pd.read_csv(path, sep="\t", low_memory=False)

        if "pathwise_module_completeness" in df.columns:
            comp_col     = "pathwise_module_completeness"
            complete_col = "pathwise_module_is_complete"
        elif "stepwise_module_completeness" in df.columns:
            comp_col     = "stepwise_module_completeness"
            complete_col = "stepwise_module_is_complete"
            print(f"  [WARN] Using stepwise columns for {sname}")
        else:
            raise KeyError(f"No completeness column found in {path}")

        df = df[["module", comp_col, complete_col]].copy()
        df.columns = ["module", "completeness", "is_complete"]
        df["is_complete"] = df["is_complete"].astype(str).str.strip() == "True"
        df["rep"] = rep
        frames.append(df)

    if not frames:
        return pd.DataFrame(columns=["module", "mean_completeness", "n_complete"])

    combined = pd.concat(frames, ignore_index=True)
    return (combined.groupby("module")
            .agg(
                mean_completeness=("completeness", "mean"),
                n_complete=("is_complete", "sum"),
            )
            .reset_index())


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def assign_category(row):
    blimmp_present = row["n_present"]  >= MIN_REPS
    anvio_present  = row["n_complete"] >= MIN_REPS
    if blimmp_present and anvio_present:   return "both"
    if blimmp_present:                     return "blimmp_only"
    if anvio_present:                      return "anvio_only"
    return "neither"


def plot_panel(ax, organism):
    blimmp = load_blimmp(organism)
    anvio  = load_anvio(organism)

    merged = pd.merge(blimmp, anvio, on="module", how="outer")
    merged["mean_confidence"]   = merged["mean_confidence"].fillna(0)
    merged["mean_completeness"] = merged["mean_completeness"].fillna(0)
    merged["n_present"]         = merged["n_present"].fillna(0)
    merged["n_complete"]        = merged["n_complete"].fillna(0)
    merged["category"]          = merged.apply(assign_category, axis=1)

    # Separate highlighted modules before drawing background points
    hl_ids = list(HIGHLIGHT_MODULES.keys())
    bg     = merged[~merged["module"].isin(hl_ids)]
    targets = {m: merged[merged["module"] == m] for m in hl_ids}

    layers = [
        ("neither",    C_NEITHER, 0.35, 18, 1),
        ("blimmp_only", C_BLIMMP, 0.70, 18, 2),
        ("anvio_only",  C_ANVIO,  0.70, 18, 2),
        ("both",        C_BOTH,   0.85, 18, 3),
    ]

    for cat, color, alpha, size, zorder in layers:
        sub = bg[bg["category"] == cat]
        ax.scatter(
            sub["mean_confidence"],
            sub["mean_completeness"],
            s=size, alpha=alpha, color=color,
            linewidths=0, zorder=zorder,
        )

    # ── Highlighted modules ───────────────────────────────────────────────
    # Stagger annotation offsets so labels don't collide. With more
    # highlighted modules than offsets, the offset list cycles (labels
    # may collide if many highlighted modules land close together in a
    # given panel -- inspect output and add offsets here if needed).
    offsets = [
        (-0.13,  0.08), (0.05,  0.08), (-0.13, -0.11),
        (0.10, -0.10), (0.14,  0.04),
    ]
    for i, mod_id in enumerate(HIGHLIGHT_MODULE_IDS):
        color, marker, _ = HIGHLIGHT_MODULES[mod_id]
        sub = targets[mod_id]
        if sub.empty:
            print(f"  {mod_id}: not found in merged data for {organism}")
            continue
        row = sub.iloc[0]
        x, y = row["mean_confidence"], row["mean_completeness"]
        sz = HIGHLIGHT_SIZES.get(marker, 100)

        ax.scatter(x, y, s=sz, marker=marker,
                   color=color, edgecolors="#222222", linewidths=0.8,
                   zorder=10 + i)

        dx, dy = offsets[i % len(offsets)]
        # Flip horizontal offset if near right or left edge
        if x + dx < 0.02:  dx = abs(dx)
        if x + dx > 0.98:  dx = -abs(dx)
        ax.annotate(
            mod_id,
            xy=(x, y), xytext=(x + dx, y + dy),
            fontsize=7.5, color="#222222", fontweight="bold",
            arrowprops=dict(arrowstyle="-", color="#555555",
                            lw=0.7, shrinkA=0, shrinkB=3),
            zorder=13,
        )
        print(f"  {mod_id}: BLIMMP={x:.3f}  anvi'o={y:.3f}"
              f"  cat={row['category']}")

    # Threshold guide lines
    ax.axvline(BLIMMP_THRESH, color=C_BLIMMP, lw=0.8, ls="--", alpha=0.45, zorder=0)
    ax.axhline(ANVIO_THRESH,  color=C_ANVIO,  lw=0.8, ls="--", alpha=0.45, zorder=0)

    # Counts annotation
    counts  = bg["category"].value_counts()
    n_both  = counts.get("both", 0)
    n_blimmp = counts.get("blimmp_only", 0) + n_both
    n_anvio  = counts.get("anvio_only",  0) + n_both
    ax.text(0.97, 0.03,
            f"BLIMMP: {n_blimmp}\nanvi'o:  {n_anvio}\nboth:    {n_both}",
            transform=ax.transAxes, fontsize=7.5,
            va="bottom", ha="right", color="#444444")

    ax.set_xlim(-0.04, 1.04)
    ax.set_ylim(-0.04, 1.04)
    ax.set_aspect("equal")
    ax.set_facecolor(SURFACE)
    ax.set_title(LABELS[organism], fontsize=9.5, pad=6, color=TEXT)
    ax.set_xlabel("BLIMMP module confidence", fontsize=9, color=TEXT)
    ax.tick_params(labelsize=8, colors=TEXT)
    ax.grid(True, color=GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    for spine in ["left", "bottom"]:
        ax.spines[spine].set_color("#cccccc")
        ax.spines[spine].set_linewidth(0.6)

    return merged


def main():
    matplotlib.rcParams.update({
        "font.family":      "STIXGeneral",
        "mathtext.fontset": "stix",
    })

    fig, axes = plt.subplots(1, 4, figsize=(15, 4),
                             sharey=True, sharex=True,
                             facecolor=SURFACE,
                             gridspec_kw={"wspace": 0.06})
    fig.patch.set_facecolor(SURFACE)

    all_data = []
    for ax, organism in zip(axes, ORGANISMS):
        print(f"\n{organism}")
        merged = plot_panel(ax, organism)
        merged.insert(0, "organism", organism)
        all_data.append(merged)

    axes[0].set_ylabel("anvi'o pathwise module completeness", fontsize=9, color=TEXT)
    axes[0].tick_params(axis="y", labelsize=8)

    # ── Legend ───────────────────────────────────────────────────────────
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=C_BOTH,
               markersize=7, label="Both present",  alpha=0.85),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=C_BLIMMP,
               markersize=6, label="BLIMMP only",   alpha=0.70),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=C_ANVIO,
               markersize=6, label="anvi'o only",   alpha=0.70),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=C_NEITHER,
               markersize=5, label="Neither",        alpha=0.45),
    ]
    for mod_id in HIGHLIGHT_MODULE_IDS:
        color, marker, label = HIGHLIGHT_MODULES[mod_id]
        ms = 11 if marker == "*" else 7
        legend_elements.append(
            Line2D([0], [0], marker=marker, color="w",
                   markerfacecolor=color,
                   markeredgecolor="#222222", markeredgewidth=0.6,
                   markersize=ms, label=label)
        )

    hl_ids_str = ", ".join(HIGHLIGHT_MODULE_IDS)
    fig.legend(handles=legend_elements, fontsize=8.5,
               loc="lower center", ncol=4,
               bbox_to_anchor=(0.5, -0.13),
               framealpha=0.85, edgecolor="#cccccc", handletextpad=0.4)

    fig.suptitle(
        f"BLIMMP vs. anvi'o module completeness at 0% gene removal"
        f"  ({hl_ids_str} highlighted)",
        fontsize=11, y=1.02, color=TEXT,
    )

    for ext in ("pdf", "png"):
        out = f"blimmp_vs_anvio_scatter.{ext}"
        plt.savefig(out, dpi=200, bbox_inches="tight", facecolor=SURFACE)
        print(f"Saved -> {out}")

    # Companion CSV
    out_df = pd.concat(all_data, ignore_index=True)
    out_df = out_df.rename(columns={
        "mean_confidence":   "blimmp_mean_confidence",
        "mean_completeness": "anvio_mean_completeness",
    })
    out_df["is_highlight"] = out_df["module"].isin(HIGHLIGHT_MODULE_IDS)
    out_df[["organism", "module", "blimmp_mean_confidence",
            "anvio_mean_completeness", "category", "is_highlight"]].to_csv(
        "blimmp_vs_anvio_scatter_data.csv", index=False)
    print("Saved -> blimmp_vs_anvio_scatter_data.csv")


if __name__ == "__main__":
    main()
