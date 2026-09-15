#!/usr/bin/env python3
"""
plot_blimmp_vs_metabolic_scatter.py

BLIMMP module_confidence (x) vs METABOLIC module completeness (y)
at 0% removal, majority-voted across 3 replicates.
One panel per genome, four panels in a row.

METABOLIC output:
  METABOLIC_SPLICED/METABOLIC_{sample}/METABOLIC_result_each_spreadsheet/
  METABOLIC_result_worksheet4.tsv
  (step-level: Module step × total Module step presence = Present/Absent)
  → completeness = fraction of steps "Present" per module per replicate
"""

import re
import warnings
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

ORGANISMS = ["MED4", "SS120", "MIT9313", "Pseudomonas_fluorescens_SBW25"]
LABELS = {
    "MED4":                          r"$P.\ marinus$ MED4",
    "SS120":                         r"$P.\ marinus$ SS120",
    "MIT9313":                       r"$P.\ marinus$ MIT9313",
    "Pseudomonas_fluorescens_SBW25": r"$P.\ fluorescens$ SBW25",
}

REMOVAL_PCT      = 0
REPLICATES       = [1, 2, 3]
MIN_REPS         = 2
BLIMMP_THRESH    = 0.60
METABOLIC_THRESH = 0.75

C_BOTH    = "#1e7a40"   # dark green — both present
C_BLIMMP  = "#1a5f9e"   # dark blue  — BLIMMP only
C_TOOL    = "#7c3aed"   # violet     — METABOLIC only
C_NEITHER = "#111111"   # black      — neither (de-emphasis)


def sample_name(organism, pct, rep):
    return f"{organism}_{pct}percentremoved_replicate{rep}"


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------
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
            .agg(mean_confidence=("module_confidence", "mean"),
                 n_present=("module_confidence", lambda x: (x >= BLIMMP_THRESH).sum()))
            .reset_index())


def load_metabolic(organism):
    """
    Fraction of steps "Present" per module from worksheet4.tsv, averaged
    across replicates.

    File: METABOLIC_SPLICED/METABOLIC_{sample}/
          METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet4.tsv
    Columns detected by partial name:
      step  col: contains "module step" but not "presence"
      pres  col: contains "step presence"
    Module IDs are extracted with regex M\\d{5} from the step column.
    """
    frames = []
    for rep in REPLICATES:
        sname = sample_name(organism, REMOVAL_PCT, rep)
        path = (BASE_DIR / "METABOLIC_SPLICED" / f"METABOLIC_{sname}"
                / "METABOLIC_result_each_spreadsheet"
                / "METABOLIC_result_worksheet4.tsv")
        if not path.exists():
            print(f"  [MISSING] METABOLIC worksheet4: METABOLIC_{sname}")
            continue
        df = pd.read_csv(path, sep="\t")
        step_cols = [c for c in df.columns
                     if "module step" in c.lower() and "presence" not in c.lower()]
        pres_cols = [c for c in df.columns if "step presence" in c.lower()]
        if not step_cols or not pres_cols:
            warnings.warn(f"  [METABOLIC ws4] {sname}: cols not found → {list(df.columns)[:8]}")
            continue
        df["mod"] = df[step_cols[0]].astype(str).str.extract(r"(M\d{5})", expand=False)
        df["present"] = (df[pres_cols[0]].astype(str).str.strip().str.lower()
                         == "present").astype(int)
        df = df.dropna(subset=["mod"])
        grp = df.groupby("mod")["present"].agg(["sum", "count"])
        grp["completeness"] = grp["sum"] / grp["count"].clip(lower=1)
        grp = grp.reset_index().rename(columns={"mod": "module"})
        grp["rep"] = rep
        frames.append(grp[["module", "completeness", "rep"]])

    if not frames:
        return pd.DataFrame(columns=["module", "mean_completeness", "n_complete"])
    combined = pd.concat(frames, ignore_index=True)
    return (combined.groupby("module")
            .agg(mean_completeness=("completeness", "mean"),
                 n_complete=("completeness", lambda x: (x >= METABOLIC_THRESH).sum()))
            .reset_index())


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def assign_category(row):
    b = row["n_present"]  >= MIN_REPS
    t = row["n_complete"] >= MIN_REPS
    if b and t:  return "both"
    if b:        return "blimmp_only"
    if t:        return "tool_only"
    return "neither"


def plot_panel(ax, organism):
    blimmp = load_blimmp(organism)
    tool   = load_metabolic(organism)

    merged = pd.merge(blimmp, tool, on="module", how="outer")
    merged["mean_confidence"]   = merged["mean_confidence"].fillna(0)
    merged["mean_completeness"] = merged["mean_completeness"].fillna(0)
    merged["n_present"]  = merged["n_present"].fillna(0)
    merged["n_complete"] = merged["n_complete"].fillna(0)
    merged["category"] = merged.apply(assign_category, axis=1)

    layers = [
        ("neither",     C_NEITHER, 0.35, 18, 1),
        ("blimmp_only", C_BLIMMP,  0.70, 18, 2),
        ("tool_only",   C_TOOL,    0.70, 18, 2),
        ("both",        C_BOTH,    0.85, 18, 3),
    ]
    for cat, color, alpha, size, zorder in layers:
        sub = merged[merged["category"] == cat]
        ax.scatter(sub["mean_confidence"], sub["mean_completeness"],
                   s=size, alpha=alpha, color=color, linewidths=0, zorder=zorder)

    ax.axvline(BLIMMP_THRESH,    color=C_BLIMMP, lw=0.8, ls="--", alpha=0.45, zorder=0)
    ax.axhline(METABOLIC_THRESH, color=C_TOOL,   lw=0.8, ls="--", alpha=0.45, zorder=0)

    counts = merged["category"].value_counts()
    n_both     = counts.get("both", 0)
    n_blimmp   = counts.get("blimmp_only", 0) + n_both
    n_metabolic = counts.get("tool_only",  0) + n_both
    ax.text(0.97, 0.03,
            f"BLIMMP:    {n_blimmp}\nMETABOLIC: {n_metabolic}\nboth:      {n_both}",
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

    axes[0].set_ylabel("METABOLIC pathway completeness", fontsize=9)

    legend_elements = [
        Line2D([0],[0], marker="o", color="w", markerfacecolor=C_BOTH,
               markersize=7, label="Both present", alpha=0.85),
        Line2D([0],[0], marker="o", color="w", markerfacecolor=C_BLIMMP,
               markersize=6, label="BLIMMP only",   alpha=0.70),
        Line2D([0],[0], marker="o", color="w", markerfacecolor=C_TOOL,
               markersize=6, label="METABOLIC only", alpha=0.70),
        Line2D([0],[0], marker="o", color="w", markerfacecolor=C_NEITHER,
               markersize=5, label="Neither",         alpha=0.45),
    ]
    fig.legend(handles=legend_elements, fontsize=8.5,
               loc="lower center", ncol=4,
               bbox_to_anchor=(0.5, -0.08),
               framealpha=0.85, edgecolor="#cccccc", handletextpad=0.4)
    fig.suptitle("BLIMMP vs. METABOLIC pathway completeness at 0% gene removal",
                 fontsize=11, y=1.02)

    plt.savefig("blimmp_vs_metabolic_scatter.pdf", dpi=200, bbox_inches="tight")
    plt.savefig("blimmp_vs_metabolic_scatter.png", dpi=200, bbox_inches="tight")
    print("\nSaved → blimmp_vs_metabolic_scatter.pdf / .png")

    # Companion CSV
    out_df = pd.concat(all_data, ignore_index=True)
    out_df = out_df.rename(columns={
        "mean_confidence":   "blimmp_mean_confidence",
        "mean_completeness": "metabolic_mean_completeness",
    })
    out_df[["organism", "module", "blimmp_mean_confidence",
            "metabolic_mean_completeness", "category"]].to_csv(
        "blimmp_vs_metabolic_scatter_data.csv", index=False)
    print("Saved → blimmp_vs_metabolic_scatter_data.csv")


if __name__ == "__main__":
    main()
