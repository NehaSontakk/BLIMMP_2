#!/usr/bin/env python3
"""
plot_module_recovery.py

Self-recovery of KEGG module calls as genes are progressively removed.
Reads pre-built binary module matrices from MODULE_MATRICES/
(run parse_module_matrices.py first).

Ground truth per tool: modules it calls at 0% removal (≥2 of 3 replicates).
Metric: fraction / count of those baseline modules still called at each level.

X-axis truncated at 90% (100% = fully empty genome, trivially no modules).

Writes:
  PLOTS/module_recovery_pct.{png,pdf}
  PLOTS/module_recovery_counts.{png,pdf}
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
MOD_DIR  = BASE_DIR / "MODULE_MATRICES"
OUT_DIR  = BASE_DIR / "PLOTS"
OUT_DIR.mkdir(exist_ok=True)

meta = pd.read_csv(MOD_DIR / "sample_metadata.csv", index_col="sample")

# Explicit panel order: the three P. marinus strains grouped together
# first, followed by P. fluorescens SBW25. A. baumannii is excluded (as
# before). Using an explicit order here (rather than sorted(), which
# alphabetizes to MED4, MIT9313, Pseudomonas..., SS120 and both breaks
# up the marinus strains and pushes SS120 to the end) so panel order
# always matches this list regardless of what's present in the metadata.
ORGANISM_ORDER = [
    "MED4",
    "SS120",
    "MIT9313",
    "Pseudomonas_fluorescens_SBW25",
]
available_in_meta = set(meta["organism"].unique())
organisms = [o for o in ORGANISM_ORDER if o in available_in_meta]
missing_from_order = available_in_meta - set(ORGANISM_ORDER) - {"Acinetobacter_baumannii"}
if missing_from_order:
    print(f"  WARNING: organism(s) in metadata but not in ORGANISM_ORDER, "
          f"excluded from plots: {sorted(missing_from_order)}")

removal_levels = sorted(meta["removal_pct"].unique())
removal_levels_plot = [r for r in removal_levels if r < 100]

print(f"Organisms:      {organisms}")
print(f"Removal levels: {removal_levels_plot}  (100% excluded)\n")

# ---------------------------------------------------------------------------
# Tool configuration
# (csv_stem, display_label, color, linestyle, marker)
# Saturated palette
# ---------------------------------------------------------------------------
TOOLS_CFG = [
    ("anvio",           "anvi'o",          "#00996b", "-",                       "^"),
    ("BLIMMP_prior",    "BLIMMP prior",    "#a07bff", (0, (1, 1)),               "o"),
    ("BLIMMP_nosub",    "BLIMMP nosub",    "#7c3aed", (0, (4, 1)),               "o"),
    ("BLIMMP",          "BLIMMP+sub",      "#4338ca", "--",                      "o"),
    ("MetaPathPredict", "MetaPathPredict", "#f0347a", (0, (3, 1, 1, 1, 1, 1)), "P"),
    ("METABOLIC",       "METABOLIC",       "#d97706", "-.",                      "D"),
    ("DRAM",            "DRAM",            "#1d6fdb", (0, (5, 2)),               "v"),
    ("KEMET",           "KEMET",           "#dd1515", ":",                       "s"),
]

TOOL_LABEL = {k: lbl for k, lbl, *_ in TOOLS_CFG}
TOOL_COLOR = {k: col for k, _, col, *_ in TOOLS_CFG}
TOOL_LS    = {k: ls  for k, _, _, ls, _ in TOOLS_CFG}
TOOL_MARK  = {k: mk  for k, _, _, _, mk in TOOLS_CFG}

# Organism display labels
ORG_LABELS = {
    "MED4":                          r"$P.\ marinus$ MED4",
    "SS120":                         r"$P.\ marinus$ SS120",
    "MIT9313":                       r"$P.\ marinus$ MIT9313",
    "Pseudomonas_fluorescens_SBW25": r"$P.\ fluorescens$ SBW25",
    "Acinetobacter_baumannii":       r"$A.\ baumannii$",
}

def org_label(org):
    return ORG_LABELS.get(org, org.replace("_", " "))

# ---------------------------------------------------------------------------
# Load binary module matrices
# ---------------------------------------------------------------------------
matrices = {}
for stem, *_ in TOOLS_CFG:
    p = MOD_DIR / f"{stem}_module_binary.csv"
    if not p.exists():
        print(f"  WARNING: {p.name} not found — skipping {stem}")
        continue
    matrices[stem] = pd.read_csv(p, index_col="sample")
    print(f"  Loaded {stem}: {matrices[stem].shape}")

available = [k for k, *_ in TOOLS_CFG if k in matrices]
print(f"\n{len(available)} tools loaded: {available}\n")

# ---------------------------------------------------------------------------
# Stable module set: present in ≥2 of 3 replicates at 0% removal
# ---------------------------------------------------------------------------
MIN_REPS = 2

def stable_set(mat, org):
    ref = meta[(meta["organism"] == org) & (meta["removal_pct"] == 0)]
    samples = [s for s in ref.index if s in mat.index]
    if not samples:
        return set()
    sub = mat.loc[samples]
    return set(sub.columns[sub.sum(axis=0) >= MIN_REPS])

stable = {}
print("Stable module counts at 0% removal (≥2 of 3 replicates):")
for key in available:
    for org in organisms:
        stable[(key, org)] = stable_set(matrices[key], org)
    row = f"  {TOOL_LABEL[key]:20s}: " + "  ".join(
        f"{org.split('_')[0]}={len(stable[(key,org)])}" for org in organisms
    )
    print(row)
print()

# ---------------------------------------------------------------------------
# Recovery metrics
# ---------------------------------------------------------------------------
rows = []
for key in available:
    mat = matrices[key]
    for org in organisms:
        gt = stable[(key, org)]
        if not gt:
            continue
        gt_cols = [c for c in gt if c in mat.columns]
        org_meta = meta[meta["organism"] == org]
        for _, row in org_meta.iterrows():
            sample = row.name
            if sample not in mat.index:
                continue
            called_gt = mat.loc[sample, gt_cols].sum() if gt_cols else 0
            rows.append(dict(
                tool=key, organism=org,
                removal_pct=row["removal_pct"], replicate=row["replicate"],
                rec_frac=called_gt / len(gt),
                rec_count=int(called_gt),
                n_stable=len(gt),
            ))

df_rec = pd.DataFrame(rows)
summary = (
    df_rec.groupby(["tool", "organism", "removal_pct"])
    .agg(
        frac_mean=("rec_frac",  "mean"), frac_sd=("rec_frac",  "std"),
        cnt_mean =("rec_count", "mean"), cnt_sd =("rec_count", "std"),
        n_stable =("n_stable",  "first"),
    )
    .reset_index()
)

print("Recovery at 0% removal (should all be 1.00):")
ref = summary[summary["removal_pct"] == 0]
tbl = ref.pivot_table(index="tool", columns="organism", values="frac_mean")
print(tbl.map(lambda x: f"{x:.2f}" if pd.notna(x) else "n/a").to_string())
print()

# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------
n_orgs = len(organisms)

def style_ax(ax, ylim=None):
    ax.set_facecolor("white")
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#aaaaaa")
    ax.tick_params(colors="#222222", labelsize=8)
    ax.grid(axis="y", color="#dddddd", linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.set_xlim(removal_levels_plot[0] - 2, removal_levels_plot[-1] + 2)
    ax.set_xticks(removal_levels_plot)
    ax.set_xticklabels([f"{p}%" for p in removal_levels_plot],
                       rotation=45, ha="right", fontsize=7.5)
    if ylim is not None:
        ax.set_ylim(*ylim)


def draw_lines(ax, org_df, mean_col, sd_col):
    for key in available:
        t = org_df[(org_df["tool"] == key) &
                   (org_df["removal_pct"] < 100)].sort_values("removal_pct")
        if t.empty:
            continue
        means = t[mean_col].clip(lower=0)
        sds   = t[sd_col].fillna(0)
        ax.plot(t["removal_pct"], means,
                color=TOOL_COLOR[key], linestyle=TOOL_LS[key],
                marker=TOOL_MARK[key], linewidth=2.0, markersize=5.5, zorder=3)
        ax.fill_between(t["removal_pct"],
                        (means - sds).clip(lower=0), means + sds,
                        color=TOOL_COLOR[key], alpha=0.15, zorder=2)


def add_n_caption(ax, org):
    """Write baseline n per tool below the subplot (not inside it)."""
    half = (len(available) + 1) // 2
    line1 = "  ".join(
        f"{TOOL_LABEL[k]}: {len(stable.get((k, org), set()))}"
        for k in available[:half]
    )
    line2 = "  ".join(
        f"{TOOL_LABEL[k]}: {len(stable.get((k, org), set()))}"
        for k in available[half:]
    )
    text = f"n (baseline) — {line1}\n{line2}" if line2 else f"n (baseline) — {line1}"
    ax.text(0.5, -0.34, text,
            transform=ax.transAxes, ha="center", va="top",
            fontsize=5.5, color="#444444", clip_on=False,
            linespacing=1.6)


legend_handles = [
    Line2D([0], [0], color=TOOL_COLOR[k], linestyle=TOOL_LS[k],
           marker=TOOL_MARK[k], linewidth=2.0, markersize=6,
           label=TOOL_LABEL[k])
    for k in available
]

# ---------------------------------------------------------------------------
# Figure 1 — Fraction recovered
# ---------------------------------------------------------------------------
fig1, axes1 = plt.subplots(1, n_orgs,
                            figsize=(4.2 * n_orgs, 5.5),
                            sharey=True, sharex=True,
                            facecolor="white",
                            gridspec_kw={"wspace": 0.08})
if n_orgs == 1:
    axes1 = [axes1]
fig1.patch.set_facecolor("white")

fig1.suptitle(
    "Module self-recovery: fraction of each tool's 0%-removal modules still detected",
    fontsize=12, fontweight="bold", y=1.03, color="#111111",
)
for ax, org in zip(axes1, organisms):
    style_ax(ax, ylim=(-0.02, 1.08))
    ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1, decimals=0))
    draw_lines(ax, summary[summary["organism"] == org], "frac_mean", "frac_sd")
    ax.set_title(org_label(org), fontsize=9, fontweight="semibold", color="#111111")
    ax.set_xlabel("Gene removal", fontsize=8.5, color="#333333")
    if ax is axes1[0]:
        ax.set_ylabel("Module recovery (fraction of baseline, mean ± SD)",
                      fontsize=8.5, color="#333333")
    add_n_caption(ax, org)

fig1.legend(handles=legend_handles, loc="lower center", ncol=len(available),
            frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, -0.12),
            labelcolor="#222222")
fig1.tight_layout(rect=[0, 0.16, 1, 1.0])
fig1.savefig(OUT_DIR / "module_recovery_pct.png", dpi=300,
             bbox_inches="tight", facecolor="white")
fig1.savefig(OUT_DIR / "module_recovery_pct.pdf",
             bbox_inches="tight", facecolor="white")
print(f"Saved: {OUT_DIR}/module_recovery_pct.{{png,pdf}}")
plt.close("all")

# ---------------------------------------------------------------------------
# Figure 2 — Count recovered
# ---------------------------------------------------------------------------
fig2, axes2 = plt.subplots(1, n_orgs,
                            figsize=(4.2 * n_orgs, 5.5),
                            sharey=False, sharex=True,
                            facecolor="white",
                            gridspec_kw={"wspace": 0.12})
if n_orgs == 1:
    axes2 = [axes2]
fig2.patch.set_facecolor("white")

fig2.suptitle(
    "Module self-recovery: count of each tool's 0%-removal modules still detected",
    fontsize=12, fontweight="bold", y=1.03, color="#111111",
)
for ax, org in zip(axes2, organisms):
    style_ax(ax)
    draw_lines(ax, summary[summary["organism"] == org], "cnt_mean", "cnt_sd")
    ax.set_title(org_label(org), fontsize=9, fontweight="semibold", color="#111111")
    ax.set_xlabel("Gene removal", fontsize=8.5, color="#333333")
    if ax is axes2[0]:
        ax.set_ylabel("Modules recovered (count, mean ± SD)",
                      fontsize=8.5, color="#333333")
    add_n_caption(ax, org)

fig2.legend(handles=legend_handles, loc="lower center", ncol=len(available),
            frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, -0.12),
            labelcolor="#222222")
fig2.tight_layout(rect=[0, 0.16, 1, 1.0])
fig2.savefig(OUT_DIR / "module_recovery_counts.png", dpi=300,
             bbox_inches="tight", facecolor="white")
fig2.savefig(OUT_DIR / "module_recovery_counts.pdf",
             bbox_inches="tight", facecolor="white")
print(f"Saved: {OUT_DIR}/module_recovery_counts.{{png,pdf}}")
plt.close("all")