#!/usr/bin/env python3
"""
plot_blimmp_vs_anvio_scatter.py

BLIMMP module_confidence (x) vs anvi'o pathwise_module_completeness (y)
at 0% removal, majority-voted across 3 replicates.
One panel per genome, four panels in a row.

Palette (validated against light surface #fcfcfb):
  Both present  : #c97d3a  orange
  BLIMMP only   : #4a3aa7  purple
  anvi'o only   : #0a8055  green
  Neither       : #aaaaaa  gray (intentionally recessive)
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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

REMOVAL_PCT      = 0
REPLICATES       = [1, 2, 3]
MIN_REPS         = 2          # majority vote threshold
BLIMMP_THRESH    = 0.60       # module_confidence cutoff for "present"
ANVIO_THRESH     = 0.75       # pathwise_module_completeness cutoff for "present"

C_BOTH    = "#1e7a40"   # dark green — both present
C_BLIMMP  = "#1a5f9e"   # dark blue  — BLIMMP only
C_ANVIO   = "#7c3aed"   # violet     — anvi'o only
C_NEITHER = "#111111"   # black      — neither (de-emphasis)


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def sample_name(organism, pct, rep):
    return f"{organism}_{pct}percentremoved_replicate{rep}"


def load_blimmp(organism):
    """
    Returns DataFrame: module | mean_confidence | n_present
    mean_confidence: mean module_confidence across replicates (NaN if file missing)
    n_present:       number of replicates where module_confidence >= BLIMMP_THRESH
    """
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
    agg = (combined.groupby("module")
           .agg(
               mean_confidence=("module_confidence", "mean"),
               n_present=("module_confidence", lambda x: (x >= BLIMMP_THRESH).sum()),
           )
           .reset_index())
    return agg


def load_anvio(organism):
    """
    Returns DataFrame: module | mean_completeness | n_complete
    mean_completeness: mean pathwise_module_completeness across replicates
    n_complete:        number of replicates where pathwise_module_is_complete == True

    Falls back to stepwise_module_completeness / stepwise_module_is_complete if
    pathwise columns are absent (older anvi'o output).
    """
    frames = []
    for rep in REPLICATES:
        sname = sample_name(organism, REMOVAL_PCT, rep)
        path = BASE_DIR / "ANVIO_SPLICED" / sname / f"{sname}_modules.txt"
        if not path.exists():
            print(f"  [MISSING] anvi'o: {path.name}")
            continue
        df = pd.read_csv(path, sep="\t", low_memory=False)

        # Resolve column names
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
    agg = (combined.groupby("module")
           .agg(
               mean_completeness=("completeness", "mean"),
               n_complete=("is_complete", "sum"),
           )
           .reset_index())
    return agg


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def assign_category(row):
    blimmp_present = row["n_present"]   >= MIN_REPS
    anvio_present  = row["n_complete"]  >= MIN_REPS
    if blimmp_present and anvio_present:
        return "both"
    elif blimmp_present:
        return "blimmp_only"
    elif anvio_present:
        return "anvio_only"
    else:
        return "neither"


def plot_panel(ax, organism):
    blimmp = load_blimmp(organism)
    anvio  = load_anvio(organism)

    merged = pd.merge(blimmp, anvio, on="module", how="outer")
    merged["mean_confidence"]  = merged["mean_confidence"].fillna(0)
    merged["mean_completeness"] = merged["mean_completeness"].fillna(0)
    merged["n_present"]  = merged["n_present"].fillna(0)
    merged["n_complete"] = merged["n_complete"].fillna(0)

    merged["category"] = merged.apply(assign_category, axis=1)

    layers = [
        ("neither",    C_NEITHER, 0.35, 18, 1),
        ("blimmp_only", C_BLIMMP, 0.70, 18, 2),
        ("anvio_only",  C_ANVIO,  0.70, 18, 2),
        ("both",        C_BOTH,   0.85, 18, 3),
    ]

    for cat, color, alpha, size, zorder in layers:
        sub = merged[merged["category"] == cat]
        ax.scatter(
            sub["mean_confidence"],
            sub["mean_completeness"],
            s=size, alpha=alpha, color=color,
            linewidths=0, zorder=zorder,
        )

    # Presence threshold lines
    ax.axvline(BLIMMP_THRESH, color=C_BLIMMP, lw=0.8, ls="--", alpha=0.45, zorder=0)
    ax.axhline(ANVIO_THRESH,  color=C_ANVIO,  lw=0.8, ls="--", alpha=0.45, zorder=0)

    # Module counts
    counts = merged["category"].value_counts()
    n_both   = counts.get("both", 0)
    n_blimmp = counts.get("blimmp_only", 0) + n_both
    n_anvio  = counts.get("anvio_only",  0) + n_both

    ax.text(0.97, 0.03,
            f"BLIMMP: {n_blimmp}\nanvi'o:  {n_anvio}\nboth:    {n_both}",
            transform=ax.transAxes, fontsize=7.5,
            va="bottom", ha="right", color="#444444")

    ax.set_xlim(-0.04, 1.04)
    ax.set_ylim(-0.04, 1.04)
    ax.set_aspect("equal")
    ax.set_title(LABELS[organism], fontsize=9.5, pad=6)
    ax.set_xlabel("BLIMMP module confidence", fontsize=9)

    ax.tick_params(labelsize=8)
    ax.grid(True, color="#e8e8e8", linewidth=0.5, zorder=0)
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)
        spine.set_color("#cccccc")

    return merged


def main():
    matplotlib.rcParams.update({
        "font.family":      "STIXGeneral",
        "mathtext.fontset": "stix",
    })

    fig, axes = plt.subplots(1, 4, figsize=(15, 4),
                             sharey=True, sharex=True,
                             gridspec_kw={"wspace": 0.06})

    all_data = []
    for ax, organism in zip(axes, ORGANISMS):
        print(f"\n{organism}")
        merged = plot_panel(ax, organism)
        merged = merged.copy()
        merged.insert(0, "organism", organism)
        all_data.append(merged)

    axes[0].set_ylabel("anvi'o pathwise module completeness", fontsize=9)

    # Legend (below all panels)
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=C_BOTH,
               markersize=7, label="Both present", alpha=0.85),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=C_BLIMMP,
               markersize=6, label="BLIMMP only", alpha=0.70),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=C_ANVIO,
               markersize=6, label="anvi'o only",  alpha=0.70),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=C_NEITHER,
               markersize=5, label="Neither",      alpha=0.45),
    ]
    fig.legend(handles=legend_elements, fontsize=8.5,
               loc="lower center", ncol=4,
               bbox_to_anchor=(0.5, -0.08),
               framealpha=0.85, edgecolor="#cccccc", handletextpad=0.4)

    fig.suptitle(
        "BLIMMP vs. anvi'o module completeness at 0% gene removal",
        fontsize=11, y=1.02,
    )

    out = "blimmp_vs_anvio_scatter.pdf"
    plt.savefig(out, dpi=200, bbox_inches="tight")
    print(f"\nSaved → {out}")

    out_png = "blimmp_vs_anvio_scatter.png"
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    print(f"Saved → {out_png}")

    # Companion CSV
    out_df = pd.concat(all_data, ignore_index=True)
    out_df = out_df.rename(columns={
        "mean_confidence":   "blimmp_mean_confidence",
        "mean_completeness": "anvio_mean_completeness",
    })
    out_df[["organism", "module", "blimmp_mean_confidence",
            "anvio_mean_completeness", "category"]].to_csv(
        "blimmp_vs_anvio_scatter_data.csv", index=False)
    print("Saved → blimmp_vs_anvio_scatter_data.csv")


if __name__ == "__main__":
    main()
