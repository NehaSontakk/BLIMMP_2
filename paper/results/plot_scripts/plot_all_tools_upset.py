#!/usr/bin/env python3
"""
plot_all_tools_upset_grouped.py

A grouped-bar variant of an UpSet plot: one shared x-axis of tool
combinations (the union of combinations seen in any genome, ranked by
total module count across all genomes), with one bar per genome per
combination (side-by-side, not stacked), colored by genome. A black
dot-matrix below the bars shows which tools are in each combination,
identical in spirit to a standard UpSet plot's matrix but drawn directly
with matplotlib rather than the upsetplot package, since upsetplot only
supports *stacked* bars-by-category (UpSet.add_stacked_bars), not
side-by-side grouped bars.

Reuses the same loaders, thresholds, and majority-vote-across-replicates
logic as plot_all_tools_consensus.py / plot_all_tools_upset.py.

Output: upset_grouped.{pdf,png}, all_tools_membership_data.csv (if not
already produced by plot_all_tools_upset.py)
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

BASE_DIR = Path("/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP")

ORGANISMS = ["MED4", "SS120", "MIT9313", "Pseudomonas_fluorescens_SBW25"]
GENOME_LABELS = {
    "MED4": "MED4",
    "SS120": "SS120",
    "MIT9313": "MIT9313",
    "Pseudomonas_fluorescens_SBW25": "SBW25",
}
GENOME_COLORS = {
    "MED4": "#4E79A7",
    "SS120": "#E15759",
    "MIT9313": "#59A14F",
    "Pseudomonas_fluorescens_SBW25": "#F28E2B",
}

REMOVAL_PCT   = 0
REPLICATES    = [1, 2, 3]
MIN_REPS      = 2
BLIMMP_THRESH = 0.60
ANVIO_THRESH  = 0.75
METABOLIC_THRESH = 0.75
KEMET_THRESH  = 0.75
MPP_THRESH    = 0.50

TOOL_NAMES = ["BLIMMP", "anvio", "METABOLIC", "KEMET", "MetaPathPredict"]


def sample_name(organism, pct, rep):
    return f"{organism}_{pct}percentremoved_replicate{rep}"


# ---------------------------------------------------------------------------
# Loaders (identical to plot_all_tools_consensus.py / plot_all_tools_upset.py)
# ---------------------------------------------------------------------------
def load_blimmp(organism):
    frames = []
    for rep in REPLICATES:
        sname = sample_name(organism, REMOVAL_PCT, rep)
        path = (BASE_DIR / "BLIMMP_SPLICED" / sname
                / f"{sname}__BLIMMP_substituted_module_probabilities.csv")
        if not path.exists():
            continue
        df = pd.read_csv(path, usecols=["module", "module_confidence"])
        df["rep"] = rep
        frames.append(df)
    if not frames:
        return pd.DataFrame(columns=["module", "present"])
    combined = pd.concat(frames, ignore_index=True)
    agg = (combined.groupby("module")["module_confidence"]
           .agg(lambda x: (x >= BLIMMP_THRESH).sum())
           .reset_index(name="n_present"))
    agg["present"] = agg["n_present"] >= MIN_REPS
    return agg[["module", "present"]]


def load_anvio(organism):
    frames = []
    for rep in REPLICATES:
        sname = sample_name(organism, REMOVAL_PCT, rep)
        path = BASE_DIR / "ANVIO_SPLICED" / sname / f"{sname}_modules.txt"
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t", low_memory=False)
        if "pathwise_module_completeness" in df.columns:
            comp_col = "pathwise_module_completeness"
        else:
            comp_col = "stepwise_module_completeness"
        df = df[["module", comp_col]].copy()
        df.columns = ["module", "completeness"]
        df["rep"] = rep
        frames.append(df)
    if not frames:
        return pd.DataFrame(columns=["module", "present"])
    combined = pd.concat(frames, ignore_index=True)
    agg = (combined.groupby("module")["completeness"]
           .agg(lambda x: (x >= ANVIO_THRESH).sum())
           .reset_index(name="n_present"))
    agg["present"] = agg["n_present"] >= MIN_REPS
    return agg[["module", "present"]]


def load_metabolic(organism):
    frames = []
    for rep in REPLICATES:
        sname = sample_name(organism, REMOVAL_PCT, rep)
        path = (BASE_DIR / "METABOLIC_SPLICED" / f"METABOLIC_{sname}"
                / "METABOLIC_result_each_spreadsheet"
                / "METABOLIC_result_worksheet4.tsv")
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        step_cols = [c for c in df.columns
                     if "module step" in c.lower() and "presence" not in c.lower()]
        pres_cols = [c for c in df.columns if "step presence" in c.lower()]
        if not step_cols or not pres_cols:
            continue
        df["mod"] = df[step_cols[0]].astype(str).str.extract(r"(M\d{5})", expand=False)
        df["present_step"] = (df[pres_cols[0]].astype(str).str.strip().str.lower()
                              == "present").astype(int)
        df = df.dropna(subset=["mod"])
        grp = df.groupby("mod")["present_step"].agg(["sum", "count"])
        grp["completeness"] = grp["sum"] / grp["count"].clip(lower=1)
        grp = grp.reset_index().rename(columns={"mod": "module"})
        grp["rep"] = rep
        frames.append(grp[["module", "completeness", "rep"]])
    if not frames:
        return pd.DataFrame(columns=["module", "present"])
    combined = pd.concat(frames, ignore_index=True)
    agg = (combined.groupby("module")["completeness"]
           .agg(lambda x: (x >= METABOLIC_THRESH).sum())
           .reset_index(name="n_present"))
    agg["present"] = agg["n_present"] >= MIN_REPS
    return agg[["module", "present"]]


def load_kemet(organism):
    import re
    frames = []
    for rep in REPLICATES:
        sname = sample_name(organism, REMOVAL_PCT, rep)
        d = BASE_DIR / "KEMET_SPLICED" / sname
        if not d.exists():
            continue
        report_files = list(d.glob("reports_tsv/reportKMC_*.tsv")) or list(d.glob("**/reportKMC_*.tsv"))
        if not report_files:
            continue
        mod_ratios = {}
        for rpt in report_files:
            df = pd.read_csv(rpt, sep="\t", header=None, usecols=[0, 3], dtype=str)
            for mod_id, ratio_str in zip(df[0], df[3]):
                mod_id = str(mod_id).strip()
                if not re.match(r"^M\d{5}$", mod_id):
                    continue
                ratio_str = str(ratio_str).strip()
                if "__" not in ratio_str:
                    continue
                found_s, total_s = ratio_str.split("__", 1)
                try:
                    found, total = float(found_s), float(total_s)
                    mod_ratios[mod_id] = found / total if total > 0 else 0.0
                except ValueError:
                    pass
        if mod_ratios:
            grp = pd.DataFrame.from_dict(mod_ratios, orient="index", columns=["completeness"])
            grp.index.name = "module"
            grp["rep"] = rep
            frames.append(grp.reset_index())
    if not frames:
        return pd.DataFrame(columns=["module", "present"])
    combined = pd.concat(frames, ignore_index=True)
    agg = (combined.groupby("module")["completeness"]
           .agg(lambda x: (x >= KEMET_THRESH).sum())
           .reset_index(name="n_present"))
    agg["present"] = agg["n_present"] >= MIN_REPS
    return agg[["module", "present"]]


_MPP_CACHE = None

def _load_mpp_file():
    global _MPP_CACHE
    if _MPP_CACHE is not None:
        return _MPP_CACHE
    path = BASE_DIR / "METAPATHPREDICT_SPLICED" / "metapathpredict_all.tsv"
    if not path.exists():
        _MPP_CACHE = pd.DataFrame()
        return _MPP_CACHE
    df = pd.read_csv(path, sep="\t", index_col=0)
    new_idx = []
    for idx in df.index:
        parts = str(idx).replace("\\", "/").split("/")
        new_idx.append(parts[-2] if len(parts) >= 2 else parts[-1])
    df.index = new_idx
    _MPP_CACHE = df
    return _MPP_CACHE


def load_mpp(organism):
    df = _load_mpp_file()
    if df.empty:
        return pd.DataFrame(columns=["module", "present"])
    frames = []
    for rep in REPLICATES:
        sname = sample_name(organism, REMOVAL_PCT, rep)
        if sname not in df.index:
            continue
        row = df.loc[sname].astype(float)
        row = row[row.index.str.match(r"^M\d{5}$")]
        frame = pd.DataFrame({"module": row.index.astype(str), "prediction": row.values, "rep": rep})
        frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=["module", "present"])
    combined = pd.concat(frames, ignore_index=True)
    agg = (combined.groupby("module")["prediction"]
           .agg(lambda x: (x >= MPP_THRESH).sum())
           .reset_index(name="n_present"))
    agg["present"] = agg["n_present"] >= MIN_REPS
    return agg[["module", "present"]]


def build_membership(organism):
    """One row per module, one boolean column per tool."""
    tools = {
        "BLIMMP": load_blimmp(organism),
        "anvio": load_anvio(organism),
        "METABOLIC": load_metabolic(organism),
        "KEMET": load_kemet(organism),
        "MetaPathPredict": load_mpp(organism),
    }
    merged = None
    for name, df in tools.items():
        df = df.rename(columns={"present": name})
        merged = df if merged is None else pd.merge(merged, df, on="module", how="outer")
    for name in tools:
        merged[name] = merged[name].fillna(False)
    return merged


def main():
    matplotlib.rcParams.update({"font.family": "STIXGeneral", "mathtext.fontset": "stix"})

    # ---- Step 1: per-genome combination counts, excluding all-absent modules ----
    per_genome_counts = {}
    all_membership = []
    for organism in ORGANISMS:
        merged = build_membership(organism)
        any_present = merged[TOOL_NAMES].any(axis=1)
        df = merged.loc[any_present, TOOL_NAMES]
        combo = df.apply(lambda row: tuple(t for t in TOOL_NAMES if row[t]), axis=1)
        per_genome_counts[organism] = combo.value_counts()

        merged_out = merged.copy()
        merged_out.insert(0, "organism", organism)
        all_membership.append(merged_out)

    # ---- Step 2: shared combination axis ----
    # Every combination observed in any genome is shown (no top-N cutoff),
    # except a combination is dropped entirely if it has zero modules in
    # every genome (i.e. it was never actually observed).
    all_combos = set()
    for counts in per_genome_counts.values():
        all_combos.update(counts.index)
    total_by_combo = {c: sum(per_genome_counts[o].get(c, 0) for o in ORGANISMS) for c in all_combos}
    all_combos = {c for c in all_combos if total_by_combo[c] > 0}

    # Order: combinations that include BLIMMP come first, ranked by the
    # number of tools in the combination (most tools first, i.e. broadest
    # agreement first), with ties broken by total module count (descending).
    # Combinations that do NOT include BLIMMP follow, ranked the same way.
    def sort_key(combo):
        has_blimmp = "BLIMMP" in combo
        return (0 if has_blimmp else 1, -len(combo), -total_by_combo[combo])

    ordered_combos = sorted(all_combos, key=sort_key)

    n_combos = len(ordered_combos)
    n_genomes = len(ORGANISMS)

    # ---- Step 3: plot ----
    fig, (ax_bar, ax_matrix) = plt.subplots(
        2, 1, figsize=(max(12, n_combos * 0.9), 7.5),
        gridspec_kw={"height_ratios": [3, 2], "hspace": 0.05},
        sharex=True,
    )

    bar_width = 0.8 / n_genomes
    x = np.arange(n_combos)

    # Alternating vertical shading behind every other combination column,
    # spanning both the bar chart and the matrix below it, so adjacent
    # combinations are visually separated even when bar heights are similar.
    # This mirrors the shading convention used by standard UpSet plots.
    for ci in range(n_combos):
        if ci % 2 == 1:
            ax_bar.axvspan(ci - 0.5, ci + 0.5, color="#f0f0f0", zorder=0)
            ax_matrix.axvspan(ci - 0.5, ci + 0.5, color="#f0f0f0", zorder=0)

    for gi, organism in enumerate(ORGANISMS):
        vals = [per_genome_counts[organism].get(c, 0) for c in ordered_combos]
        offset = (gi - (n_genomes - 1) / 2) * bar_width
        ax_bar.bar(x + offset, vals, bar_width,
                   color=GENOME_COLORS[organism], label=GENOME_LABELS[organism],
                   zorder=3)

    ax_bar.set_ylabel("Number of modules", fontsize=10)
    ax_bar.set_title(
        "Tool agreement on module presence, by genome (0% gene removal)\n"
        "BLIMMP-containing combinations first (most tools to fewest), then others; "
        "combinations with zero modules in every genome excluded",
        fontsize=11,
    )
    ax_bar.legend(title="Genome", frameon=False, fontsize=9, loc="upper right")
    for spine in ["top", "right"]:
        ax_bar.spines[spine].set_visible(False)
    ax_bar.grid(axis="y", color="#e0e0e0", linewidth=0.5, zorder=1)
    ax_bar.set_axisbelow(True)
    ax_bar.set_xlim(-0.6, n_combos - 0.4)

    # Black dot matrix, matching standard UpSet plot conventions
    for ci, combo in enumerate(ordered_combos):
        in_idx = [ti for ti, t in enumerate(TOOL_NAMES) if t in combo]
        for ti in range(len(TOOL_NAMES)):
            ax_matrix.scatter(ci, ti, s=140,
                              color="black" if ti in in_idx else "#d9d9d9",
                              zorder=3)
        if len(in_idx) > 1:
            ax_matrix.plot([ci, ci], [min(in_idx), max(in_idx)],
                           color="black", lw=1.5, zorder=2)

    ax_matrix.set_yticks(range(len(TOOL_NAMES)))
    ax_matrix.set_yticklabels(TOOL_NAMES, fontsize=9.5)
    ax_matrix.set_xticks(x)
    ax_matrix.set_xticklabels([])
    ax_matrix.set_ylim(-0.6, len(TOOL_NAMES) - 0.4)
    ax_matrix.invert_yaxis()
    ax_matrix.set_xlim(-0.6, n_combos - 0.4)
    for spine in ax_matrix.spines.values():
        spine.set_visible(False)

    plt.tight_layout()
    plt.savefig("upset_grouped.pdf", dpi=200, bbox_inches="tight")
    plt.savefig("upset_grouped.png", dpi=200, bbox_inches="tight")
    print("Saved -> upset_grouped.pdf / .png")

    out_df = pd.concat(all_membership, ignore_index=True)
    out_df.to_csv("all_tools_membership_data.csv", index=False)
    print("Saved -> all_tools_membership_data.csv")


if __name__ == "__main__":
    main()
