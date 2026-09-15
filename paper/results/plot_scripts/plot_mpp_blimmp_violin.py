#!/usr/bin/env python3
"""
plot_blimmp_vs_mpp_violin.py

For the 100%-complete genome (0% removal), show how BLIMMP module
confidence is distributed for modules MetaPathPredict (MPP) calls
ABSENT vs PRESENT (majority vote across 3 replicates).

Because MPP outputs a binary 0/1 per module, a scatter plot collapses
onto two horizontal lines; a violin + strip plot reveals whether the
continuous BLIMMP confidence distinguishes the two groups.

Layout: 1 row × 4 panels (one per organism)
  X: MPP call (absent / present) — majority vote (≥ 2/3 replicates)
  Y: BLIMMP mean module_confidence (mean across 3 replicates)
  Violin body: distribution shape, neutral fill
  Strip dots: colored by agreement category
    both        — BLIMMP high & MPP present   "#1e7a40" (dark green)
    blimmp_only — BLIMMP high, MPP absent     "#1a5f9e" (dark blue)
    tool_only   — BLIMMP low,  MPP present    "#7c3aed" (violet)
    neither     — both low                    "#111111" (black, de-emphasis)

Palette validated (3 meaningful hues: green/blue/violet): all CVD checks pass.
Black for "neither" is de-emphasis (secondary encoding: tiny size + low alpha).
Font: Times New Roman (serif) with STIX math for italic species names.
Legend placed below all panels.

Output: blimmp_vs_mpp_violin.{pdf,png}
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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
MPP_THRESH    = 0.50   # binary output — 0.5 means ≥ 1 replicate called present

# Validated categorical palette — all 3 meaningful hues pass full CVD check:
#   green ↔ blue ΔE 19.0 (normal), 20.0 (deutan)
#   violet ↔ blue ΔE 18.8 (normal), 12.1 (deutan)
#   tritan ΔE 4.0 (floor band) → secondary encoding via size difference
# Black for "neither" is de-emphasis (not in the validated set): small dots, low alpha
C_BOTH    = "#1e7a40"   # both present   — dark green
C_BLIMMP  = "#1a5f9e"   # BLIMMP only    — dark blue
C_MPP     = "#7c3aed"   # MPP only       — violet (only CVD-safe hue vs green)
C_NEITHER = "#111111"   # neither        — black, de-emphasis via small size + low alpha

SURFACE = "#fcfcfb"
TEXT    = "#2e2e2e"
GRID    = "#e5e4de"


def sample_name(organism, pct, rep):
    return f"{organism}_{pct}percentremoved_replicate{rep}"


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------
def load_blimmp(organism):
    frames = []
    for rep in REPLICATES:
        sname = sample_name(organism, REMOVAL_PCT, rep)
        d = BASE_DIR / "BLIMMP_SPLICED" / sname
        cands = (list(d.glob("*_substituted_module_probabilities.csv"))
                 if d.exists() else [])
        if not cands:
            print(f"  [MISSING] BLIMMP: {sname}")
            continue
        df = pd.read_csv(cands[0], usecols=["module", "module_confidence"])
        df["rep"] = rep
        frames.append(df)
    if not frames:
        return pd.DataFrame(columns=["module", "mean_confidence", "n_present"])
    combined = pd.concat(frames, ignore_index=True)
    return (combined.groupby("module")
            .agg(mean_confidence=("module_confidence", "mean"),
                 n_present=("module_confidence", lambda x: (x >= BLIMMP_THRESH).sum()))
            .reset_index())


_MPP_CACHE = None

def _load_mpp_file():
    global _MPP_CACHE
    if _MPP_CACHE is not None:
        return _MPP_CACHE
    path = BASE_DIR / "METAPATHPREDICT_SPLICED" / "metapathpredict_all.tsv"
    if not path.exists():
        print(f"  [MISSING] MetaPathPredict output: {path}")
        _MPP_CACHE = pd.DataFrame()
        return _MPP_CACHE
    df = pd.read_csv(path, sep="\t", index_col=0)
    # Normalise index: extract sample name from path-like index entries
    new_idx = []
    for idx in df.index:
        parts = str(idx).replace("\\", "/").split("/")
        new_idx.append(parts[-2] if len(parts) >= 2 else parts[-1])
    df.index = new_idx
    print(f"  MetaPathPredict: {df.shape[0]} samples × {df.shape[1]} modules")
    _MPP_CACHE = df
    return _MPP_CACHE


def load_mpp(organism):
    df = _load_mpp_file()
    if df.empty:
        return pd.DataFrame(columns=["module", "mean_prediction", "n_complete"])
    frames = []
    for rep in REPLICATES:
        sname = sample_name(organism, REMOVAL_PCT, rep)
        if sname not in df.index:
            print(f"  [MISSING] MPP row: '{sname}'")
            continue
        row = df.loc[sname].astype(float)
        row = row[row.index.str.match(r"^M\d{5}$")]
        frame = pd.DataFrame({
            "module":      row.index.astype(str),
            "prediction":  row.values,
            "rep":         rep,
        })
        frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=["module", "mean_prediction", "n_complete"])
    combined = pd.concat(frames, ignore_index=True)
    return (combined.groupby("module")
            .agg(mean_prediction=("prediction", "mean"),
                 n_complete=("prediction", lambda x: (x >= MPP_THRESH).sum()))
            .reset_index())


# ---------------------------------------------------------------------------
# Category
# ---------------------------------------------------------------------------
def assign_category(row):
    b = row["n_present"]  >= MIN_REPS
    t = row["n_complete"] >= MIN_REPS
    if b and t:  return "both"
    if b:        return "blimmp_only"
    if t:        return "tool_only"
    return "neither"


# ---------------------------------------------------------------------------
# Panel
# ---------------------------------------------------------------------------
def plot_panel(ax, organism):
    blimmp = load_blimmp(organism)
    mpp    = load_mpp(organism)

    merged = pd.merge(blimmp, mpp, on="module", how="outer")
    merged["mean_confidence"]  = merged["mean_confidence"].fillna(0)
    merged["mean_prediction"]  = merged["mean_prediction"].fillna(0)
    merged["n_present"]        = merged["n_present"].fillna(0)
    merged["n_complete"]       = merged["n_complete"].fillna(0)
    merged["category"]         = merged.apply(assign_category, axis=1)
    merged["mpp_group"]        = merged["n_complete"].apply(
        lambda x: "present" if x >= MIN_REPS else "absent")

    absent  = merged.loc[merged["mpp_group"] == "absent",  "mean_confidence"].values
    present = merged.loc[merged["mpp_group"] == "present", "mean_confidence"].values

    # ---- Violin bodies -------------------------------------------------------
    data_groups = [d for d in [absent, present] if len(d) >= 3]
    positions   = [i for i, d in enumerate([absent, present]) if len(d) >= 3]

    if data_groups:
        parts = ax.violinplot(data_groups, positions=positions,
                              widths=0.55, showmeans=False,
                              showmedians=True, showextrema=False)
        for pc in parts["bodies"]:
            pc.set_facecolor("#ebebeb")
            pc.set_edgecolor("#aaaaaa")
            pc.set_linewidth(0.7)
            pc.set_alpha(0.55)
            pc.set_zorder(1)
        parts["cmedians"].set_color("#555555")
        parts["cmedians"].set_linewidth(1.6)
        parts["cmedians"].set_zorder(5)

    # ---- Strip plot (jittered) -----------------------------------------------
    rng = np.random.default_rng(42)
    layers = [
        # category, x_group, color, size, alpha, zorder
        ("neither",     "absent",  C_NEITHER, 18, 0.30, 2),
        ("neither",     "present", C_NEITHER, 18, 0.30, 2),
        ("blimmp_only", "absent",  C_BLIMMP,  18, 0.72, 4),
        ("tool_only",   "present", C_MPP,     18, 0.72, 4),
        ("both",        "present", C_BOTH,    18, 0.85, 5),
        ("blimmp_only", "present", C_BLIMMP,  18, 0.55, 3),
        ("tool_only",   "absent",  C_MPP,     18, 0.55, 3),
        ("both",        "absent",  C_BOTH,    18, 0.55, 3),
    ]
    for cat, group, color, size, alpha, zorder in layers:
        pos = 0 if group == "absent" else 1
        sub = merged[(merged["category"] == cat) & (merged["mpp_group"] == group)]
        if sub.empty:
            continue
        vals    = sub["mean_confidence"].values
        jitter  = rng.normal(0, 0.07, len(vals)).clip(-0.22, 0.22)
        ax.scatter(pos + jitter, vals,
                   s=size, alpha=alpha, color=color,
                   linewidths=0, zorder=zorder)

    # ---- Threshold line ------------------------------------------------------
    ax.axhline(BLIMMP_THRESH, color=C_BLIMMP, lw=0.9, ls="--", alpha=0.5, zorder=0)

    # ---- Counts annotation ---------------------------------------------------
    for pos, grp, vals_arr in [(0, "absent", absent), (1, "present", present)]:
        ax.text(pos, -0.1, f"n={len(vals_arr)}",
                ha="center", va="top", fontsize=7.5,
                color="#666666", transform=ax.get_xaxis_transform())

    # ---- Style ---------------------------------------------------------------
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(-0.04, 1.04)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["MPP\nabsent", "MPP\npresent"], fontsize=8.5)
    ax.tick_params(labelsize=8, colors=TEXT)
    ax.set_title(LABELS[organism], fontsize=9.5, pad=6, color=TEXT)
    ax.grid(axis="y", color=GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    for spine in ["left", "bottom"]:
        ax.spines[spine].set_color("#cccccc")
        ax.spines[spine].set_linewidth(0.6)

    return merged


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    matplotlib.rcParams.update({
        "font.family":      "STIXGeneral",
        "mathtext.fontset": "stix",
    })

    fig, axes = plt.subplots(1, 4, figsize=(14, 4.2),
                             sharey=True, sharex=False,
                             facecolor=SURFACE,
                             gridspec_kw={"wspace": 0.12})
    fig.patch.set_facecolor(SURFACE)

    all_data = []
    for ax, organism in zip(axes, ORGANISMS):
        ax.set_facecolor(SURFACE)
        print(f"\n{organism}")
        merged = plot_panel(ax, organism)
        merged = merged.copy()
        merged.insert(0, "organism", organism)
        all_data.append(merged)

    axes[0].set_ylabel("BLIMMP module confidence", fontsize=9, color=TEXT)
    axes[0].tick_params(axis="y", labelsize=8)

    # Threshold label on first panel
    axes[0].text(0.02, BLIMMP_THRESH + 0.02,
                 "BLIMMP thresh 0.60",
                 transform=axes[0].get_yaxis_transform(),
                 fontsize=6.5, color=C_BLIMMP, alpha=0.75, ha="left")

    # ---- Legend (below all panels) -------------------------------------------
    legend_patches = [
        mpatches.Patch(color=C_BOTH,    label="Both present"),
        mpatches.Patch(color=C_BLIMMP,  label="BLIMMP only"),
        mpatches.Patch(color=C_MPP,     label="MPP only"),
        mpatches.Patch(color=C_NEITHER, label="Neither", alpha=0.5),
    ]
    fig.legend(handles=legend_patches, fontsize=8.5,
               loc="lower center", ncol=4,
               bbox_to_anchor=(0.5, -0.08),
               framealpha=0.85, edgecolor="#cccccc",
               handlelength=1.2, handletextpad=0.5)

    fig.suptitle(
        "BLIMMP module confidence vs MetaPathPredict call at 0% gene removal",
        fontsize=11, y=1.02, color=TEXT)

    plt.savefig("blimmp_vs_mpp_violin.pdf", dpi=200, bbox_inches="tight",
                facecolor=SURFACE)
    plt.savefig("blimmp_vs_mpp_violin.png", dpi=200, bbox_inches="tight",
                facecolor=SURFACE)
    print("\nSaved → blimmp_vs_mpp_violin.pdf / .png")

    # Companion CSV
    out_df = pd.concat(all_data, ignore_index=True)
    out_df = out_df.rename(columns={
        "mean_confidence":  "blimmp_mean_confidence",
        "mean_prediction":  "mpp_mean_prediction",
    })
    out_df[["organism", "module", "blimmp_mean_confidence",
            "mpp_mean_prediction", "mpp_group", "category"]].to_csv(
        "blimmp_vs_mpp_violin_data.csv", index=False)
    print("Saved → blimmp_vs_mpp_violin_data.csv")


if __name__ == "__main__":
    main()
