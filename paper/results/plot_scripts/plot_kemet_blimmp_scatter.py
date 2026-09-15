#!/usr/bin/env python3
"""
plot_blimmp_vs_kemet_scatter.py

BLIMMP module_confidence (x) vs KEMET module completeness (y)
at 0% removal, majority-voted across 3 replicates.
One panel per genome, four panels in a row.
"""

import re
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

REMOVAL_PCT   = 0
REPLICATES    = [1, 2, 3]
MIN_REPS      = 2
BLIMMP_THRESH = 0.60
KEMET_THRESH  = 0.75

C_BOTH    = "#1e7a40"   # dark green — both present
C_BLIMMP  = "#1a5f9e"   # dark blue  — BLIMMP only
C_TOOL    = "#7c3aed"   # violet     — KEMET only
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


def load_kemet(organism):
    """
    File: KEMET_SPLICED/{sample}/reports_tsv/reportKMC_*.tsv  (no header)
    Col 0: module ID (M#####),  Col 3: found__total (double-underscore ratio)
    """
    frames = []
    for rep in REPLICATES:
        sname = sample_name(organism, REMOVAL_PCT, rep)
        d = BASE_DIR / "KEMET_SPLICED" / sname
        if not d.exists():
            print(f"  [MISSING] KEMET dir: {sname}")
            continue
        report_files = list(d.glob("reports_tsv/reportKMC_*.tsv"))
        if not report_files:
            report_files = list(d.glob("**/reportKMC_*.tsv"))
        if not report_files:
            print(f"  [MISSING] KEMET reportKMC for {sname}")
            continue
        mod_ratios = {}
        for rpt in report_files:
            try:
                df = pd.read_csv(rpt, sep="\t", header=None,
                                 usecols=[0, 3], dtype=str)
                for mod_id, ratio_str in zip(df[0], df[3]):
                    mod_id = str(mod_id).strip()
                    if not re.match(r"^M\d{5}$", mod_id):
                        continue
                    ratio_str = str(ratio_str).strip()
                    if "__" not in ratio_str:
                        continue
                    found_s, total_s = ratio_str.split("__", 1)
                    try:
                        found = float(found_s)
                        total = float(total_s)
                        mod_ratios[mod_id] = found / total if total > 0 else 0.0
                    except ValueError:
                        pass
            except Exception as e:
                import warnings
                warnings.warn(f"KEMET {sname} {rpt.name}: {e}")
        if mod_ratios:
            grp = pd.DataFrame.from_dict(mod_ratios, orient="index",
                                         columns=["completeness"])
            grp.index.name = "module"
            grp["rep"] = rep
            frames.append(grp.reset_index())

    if not frames:
        return pd.DataFrame(columns=["module", "mean_completeness", "n_complete"])
    combined = pd.concat(frames, ignore_index=True)
    return (combined.groupby("module")
            .agg(mean_completeness=("completeness", "mean"),
                 n_complete=("completeness", lambda x: (x >= KEMET_THRESH).sum()))
            .reset_index())


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def assign_category(row, min_reps=MIN_REPS):
    b = row["n_present"]  >= min_reps
    t = row["n_complete"] >= min_reps
    if b and t:     return "both"
    if b:           return "blimmp_only"
    if t:           return "tool_only"
    return "neither"


def plot_panel(ax, organism):
    blimmp = load_blimmp(organism)
    tool   = load_kemet(organism)

    merged = pd.merge(blimmp, tool, on="module", how="outer")
    merged["mean_confidence"]   = merged["mean_confidence"].fillna(0)
    merged["mean_completeness"] = merged["mean_completeness"].fillna(0)
    merged["n_present"]  = merged["n_present"].fillna(0)
    merged["n_complete"] = merged["n_complete"].fillna(0)
    merged["category"] = merged.apply(assign_category, axis=1)

    layers = [
        ("neither",    C_NEITHER, 0.35, 18, 1),
        ("blimmp_only", C_BLIMMP, 0.70, 18, 2),
        ("tool_only",   C_TOOL,   0.70, 18, 2),
        ("both",        C_BOTH,   0.85, 18, 3),
    ]
    for cat, color, alpha, size, zorder in layers:
        sub = merged[merged["category"] == cat]
        ax.scatter(sub["mean_confidence"], sub["mean_completeness"],
                   s=size, alpha=alpha, color=color, linewidths=0, zorder=zorder)

    ax.axvline(BLIMMP_THRESH, color=C_BLIMMP, lw=0.8, ls="--", alpha=0.45, zorder=0)
    ax.axhline(KEMET_THRESH,  color=C_TOOL,   lw=0.8, ls="--", alpha=0.45, zorder=0)

    counts = merged["category"].value_counts()
    n_both   = counts.get("both", 0)
    n_blimmp = counts.get("blimmp_only", 0) + n_both
    n_tool   = counts.get("tool_only",   0) + n_both
    ax.text(0.97, 0.03,
            f"BLIMMP: {n_blimmp}\nKEMET:  {n_tool}\nboth:   {n_both}",
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

    axes[0].set_ylabel("KEMET module completeness", fontsize=9)

    legend_elements = [
        Line2D([0],[0], marker="o", color="w", markerfacecolor=C_BOTH,
               markersize=7, label="Both present", alpha=0.85),
        Line2D([0],[0], marker="o", color="w", markerfacecolor=C_BLIMMP,
               markersize=6, label="BLIMMP only", alpha=0.70),
        Line2D([0],[0], marker="o", color="w", markerfacecolor=C_TOOL,
               markersize=6, label="KEMET only",  alpha=0.70),
        Line2D([0],[0], marker="o", color="w", markerfacecolor=C_NEITHER,
               markersize=5, label="Neither",      alpha=0.45),
    ]
    fig.legend(handles=legend_elements, fontsize=8.5,
               loc="lower center", ncol=4,
               bbox_to_anchor=(0.5, -0.08),
               framealpha=0.85, edgecolor="#cccccc", handletextpad=0.4)
    fig.suptitle("BLIMMP vs. KEMET module completeness at 0% gene removal",
                 fontsize=11, y=1.02)

    plt.savefig("blimmp_vs_kemet_scatter.pdf", dpi=200, bbox_inches="tight")
    plt.savefig("blimmp_vs_kemet_scatter.png", dpi=200, bbox_inches="tight")
    print("\nSaved → blimmp_vs_kemet_scatter.pdf / .png")

    # Companion CSV
    out_df = pd.concat(all_data, ignore_index=True)
    out_df = out_df.rename(columns={
        "mean_confidence":   "blimmp_mean_confidence",
        "mean_completeness": "kemet_mean_completeness",
    })
    out_df[["organism", "module", "blimmp_mean_confidence",
            "kemet_mean_completeness", "category"]].to_csv(
        "blimmp_vs_kemet_scatter_data.csv", index=False)
    print("Saved → blimmp_vs_kemet_scatter_data.csv")


if __name__ == "__main__":
    main()
