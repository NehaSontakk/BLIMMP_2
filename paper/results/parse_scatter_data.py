#!/usr/bin/env python3
"""
parse_scatter_data.py

Extract per-module continuous completeness scores and per-replicate
presence counts at 0% gene removal for BLIMMP and all comparison tools.
Values are averaged (or counted) across 3 replicates.

Output:
  MODULE_MATRICES/scatter_data_0pct.xlsx
    One sheet per organism.  Columns:
      module
      blimmp_mean_conf,      blimmp_n_present       (# reps with conf  >= 0.60)
      anvio_mean_comp,       anvio_n_complete        (# reps with comp  >= 0.75)
      kemet_mean_comp,       kemet_n_complete        (# reps with ratio >= 0.75)
      metabolic_mean_comp,   metabolic_n_complete    (# reps with frac  >= 0.75)
      dram_mean_comp,        dram_n_complete         (# reps with comp  >= 0.75)
      mpp_mean_pred,         mpp_n_complete          (# reps with pred  >= 0.50)
      cat_vs_anvio, cat_vs_kemet, cat_vs_metabolic,
      cat_vs_dram,  cat_vs_mpp
        values: both | blimmp_only | tool_only | neither
        (majority vote: both n values >= MIN_REPS = 2)

Data sources mirror what the scatter plot scripts read:
  BLIMMP    : BLIMMP_SPLICED/{s}/*_substituted_module_probabilities.csv
  anvi'o    : ANVIO_SPLICED/{s}/{s}_modules.txt  → stepwise_module_completeness
  KEMET     : KEMET_SPLICED/{s}/reports_tsv/reportKMC_*.tsv  → col 3 found__total
  METABOLIC : METABOLIC_SPLICED/METABOLIC_{s}/METABOLIC_result_each_spreadsheet/
              METABOLIC_result_worksheet4.tsv  → fraction of steps "Present"
  DRAM      : DRAM_SPLICED/{s}/distill/metabolism_summary.xlsx → module completeness
  MPP       : METAPATHPREDICT_SPLICED/metapathpredict_all.tsv  → binary 0/1

Usage:
    conda activate blimmp-work
    python parse_scatter_data.py
"""

import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
BASE_DIR = Path("/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP")
OUT_DIR  = BASE_DIR / "MODULE_MATRICES"
OUT_DIR.mkdir(exist_ok=True)

ORGANISMS   = ["MED4", "SS120", "MIT9313", "Pseudomonas_fluorescens_SBW25"]
REPLICATES  = [1, 2, 3]
REMOVAL_PCT = 0
MIN_REPS    = 2   # majority-vote threshold

THRESH = {
    "blimmp"   : 0.60,
    "anvio"    : 0.75,
    "kemet"    : 0.75,
    "metabolic": 0.75,
    "dram"     : 0.75,
    "mpp"      : 0.50,
}

def sname(org, rep):
    return f"{org}_{REMOVAL_PCT}percentremoved_replicate{rep}"

# ---------------------------------------------------------------------------
# Bacterial module universe
# ---------------------------------------------------------------------------
MODULE_KO_JSON = BASE_DIR / "BLIMMP_8Sep2026" / "BLIMMP" / "module_ko_reaction.json"
with open(MODULE_KO_JSON) as _f:
    BACTERIAL_MODULES = set(json.load(_f).keys())
print(f"Bacterial module universe: {len(BACTERIAL_MODULES):,} modules\n")


# ---------------------------------------------------------------------------
# Per-replicate loaders
# Each returns a pd.Series(float) indexed by module_id for one sample.
# Missing file → empty Series.
# ---------------------------------------------------------------------------

def load_blimmp_rep(org, rep) -> pd.Series:
    s = sname(org, rep)
    d = BASE_DIR / "BLIMMP_SPLICED" / s
    if not d.exists():
        return pd.Series(dtype=float)
    cands = (list(d.glob("*_substituted_module_probabilities.csv"))
             or list(d.glob("*_module_probabilities.csv")))
    if not cands:
        print(f"    [MISSING] BLIMMP csv for {s}")
        return pd.Series(dtype=float)
    df = pd.read_csv(cands[0], usecols=["module", "module_confidence"])
    return df.set_index("module")["module_confidence"]


def load_anvio_rep(org, rep) -> pd.Series:
    s = sname(org, rep)
    p = BASE_DIR / "ANVIO_SPLICED" / s / f"{s}_modules.txt"
    if not p.exists():
        print(f"    [MISSING] anvi'o modules.txt for {s}")
        return pd.Series(dtype=float)
    df = pd.read_csv(p, sep="\t", usecols=["module", "stepwise_module_completeness"])
    # Multiple completeness paths per module → take max
    return df.groupby("module")["stepwise_module_completeness"].max()


def load_kemet_rep(org, rep) -> pd.Series:
    """Col 3 of reportKMC_*.tsv: 'found__total' → ratio."""
    s = sname(org, rep)
    d = BASE_DIR / "KEMET_SPLICED" / s
    if not d.exists():
        print(f"    [MISSING] KEMET dir for {s}")
        return pd.Series(dtype=float)
    rpts = list(d.glob("reports_tsv/reportKMC_*.tsv"))
    if not rpts:
        rpts = list(d.glob("**/reportKMC_*.tsv"))
    if not rpts:
        print(f"    [MISSING] KEMET reportKMC for {s}")
        return pd.Series(dtype=float)
    ratios = {}
    for rpt in rpts:
        try:
            df = pd.read_csv(rpt, sep="\t", header=None, usecols=[0, 3], dtype=str)
            for mid, ratio_str in zip(df[0], df[3]):
                mid = str(mid).strip()
                if not re.match(r"^M\d{5}$", mid):
                    continue
                ratio_str = str(ratio_str).strip()
                if "__" not in ratio_str:
                    continue
                a, b = ratio_str.split("__", 1)
                try:
                    ratios[mid] = float(a) / float(b) if float(b) > 0 else 0.0
                except ValueError:
                    pass
        except Exception as e:
            warnings.warn(f"KEMET {s} {rpt.name}: {e}")
    return pd.Series(ratios, dtype=float)


def load_metabolic_rep(org, rep) -> pd.Series:
    """Fraction of steps 'Present' per module from worksheet4.tsv."""
    s = sname(org, rep)
    p = (BASE_DIR / "METABOLIC_SPLICED" / f"METABOLIC_{s}"
         / "METABOLIC_result_each_spreadsheet"
         / "METABOLIC_result_worksheet4.tsv")
    if not p.exists():
        print(f"    [MISSING] METABOLIC worksheet4 for {s}")
        return pd.Series(dtype=float)
    df = pd.read_csv(p, sep="\t")
    step_cols = [c for c in df.columns
                 if "module step" in c.lower() and "presence" not in c.lower()]
    pres_cols = [c for c in df.columns if "step presence" in c.lower()]
    if not step_cols or not pres_cols:
        warnings.warn(f"  [METABOLIC ws4] {s}: columns not found → {list(df.columns)[:8]}")
        return pd.Series(dtype=float)
    df["mod"] = df[step_cols[0]].astype(str).str.extract(r"(M\d{5})", expand=False)
    df["present"] = (df[pres_cols[0]].astype(str).str.strip().str.lower()
                     == "present").astype(int)
    df = df.dropna(subset=["mod"])
    grp = df.groupby("mod")["present"].agg(["sum", "count"])
    return (grp["sum"] / grp["count"].clip(lower=1)).rename_axis("module")


def load_dram_rep(org, rep) -> pd.Series:
    """Module completeness from metabolism_summary.xlsx."""
    s = sname(org, rep)
    path = BASE_DIR / "DRAM_SPLICED" / s / "distill" / "metabolism_summary.xlsx"
    if not path.exists():
        print(f"    [MISSING] DRAM metabolism_summary for {s}")
        return pd.Series(dtype=float)
    try:
        xl = pd.ExcelFile(path)
        # Prefer sheet with "module" in its name
        module_sheets = [sh for sh in xl.sheet_names if "module" in sh.lower()]
        sheet = module_sheets[0] if module_sheets else xl.sheet_names[0]
        df = xl.parse(sheet)
        # Auto-detect module ID column
        mod_col = next(
            (c for c in df.columns if c.lower() in
             ("module", "module_id", "kegg_module", "id")),
            next((c for c in df.columns if "module" in c.lower()), None)
        )
        # Auto-detect completeness column
        comp_col = next(
            (c for c in df.columns if c.lower() in
             ("completeness", "completion", "percent_steps_complete",
              "steps_complete", "fraction_complete")),
            next((c for c in df.columns
                  if any(k in c.lower() for k in
                         ["complet", "steps", "fraction", "score"])), None)
        )
        if mod_col and comp_col:
            sub = df[[mod_col, comp_col]].dropna()
            sub.columns = ["module", "completeness"]
            sub["completeness"] = pd.to_numeric(
                sub["completeness"], errors="coerce").fillna(0)
            return sub.set_index("module")["completeness"]
        # Fallback: try using sample name as column (wide format)
        if mod_col and s in df.columns:
            sub = df[[mod_col, s]].dropna()
            sub.columns = ["module", "completeness"]
            sub["completeness"] = pd.to_numeric(
                sub["completeness"], errors="coerce").fillna(0)
            return sub.set_index("module")["completeness"]
        warnings.warn(f"  [DRAM] {s}: sheet='{sheet}', "
                      f"cols={list(df.columns)[:8]} — cannot parse")
    except Exception as e:
        warnings.warn(f"  [DRAM] {s}: {e}")
    return pd.Series(dtype=float)


# MetaPathPredict: one TSV for all samples, loaded once
_MPP_DF = None

def _get_mpp() -> pd.DataFrame:
    global _MPP_DF
    if _MPP_DF is not None:
        return _MPP_DF
    p = BASE_DIR / "METAPATHPREDICT_SPLICED" / "metapathpredict_all.tsv"
    if not p.exists():
        print("  [MISSING] metapathpredict_all.tsv — MPP will be empty")
        _MPP_DF = pd.DataFrame()
        return _MPP_DF
    _MPP_DF = pd.read_csv(p, sep="\t", index_col=0)
    # Index may be a file path like KOFAM_SPLICED_FIXED/{sample}/{sample}_kofam.tsv
    # → extract second-to-last path component as sample name
    new_idx = []
    for idx in _MPP_DF.index:
        parts = str(idx).replace("\\", "/").split("/")
        new_idx.append(parts[-2] if len(parts) >= 2 else parts[-1])
    _MPP_DF.index = new_idx
    print(f"  MetaPathPredict: {_MPP_DF.shape[0]} samples × {_MPP_DF.shape[1]} modules")
    return _MPP_DF


def load_mpp_rep(org, rep) -> pd.Series:
    """Binary 0/1 prediction from MetaPathPredict all-sample TSV."""
    s = sname(org, rep)
    df = _get_mpp()
    if df.empty or s not in df.index:
        if not df.empty:
            print(f"    [MISSING] MPP row: '{s}'")
        return pd.Series(dtype=float)
    row = df.loc[s].astype(float)
    # Keep only columns that look like KEGG module IDs
    return row[row.index.str.match(r"^M\d{5}$")]


# ---------------------------------------------------------------------------
# Aggregate across replicates
# Returns DataFrame with columns: module, mean_col, n_col
# ---------------------------------------------------------------------------
def agg_tool(org, loader_fn, thresh, mean_col, n_col) -> pd.DataFrame:
    rep_series = []
    for rep in REPLICATES:
        ser = loader_fn(org, rep)
        # Restrict to bacterial module universe
        ser = ser[ser.index.isin(BACTERIAL_MODULES)] if not ser.empty else ser
        rep_series.append(ser.rename(rep))

    if all(s.empty for s in rep_series):
        return pd.DataFrame(columns=["module", mean_col, n_col])

    combined = pd.concat(rep_series, axis=1)   # index=modules, cols=rep numbers
    mean_vals = combined.mean(axis=1)
    n_vals    = (combined >= thresh).sum(axis=1)

    return pd.DataFrame({
        "module":  mean_vals.index,
        mean_col:  mean_vals.round(6).values,
        n_col:     n_vals.astype(int).values,
    })


# ---------------------------------------------------------------------------
# Category assignment (majority vote)
# ---------------------------------------------------------------------------
def assign_cat(blimmp_n, tool_n) -> str:
    b = int(blimmp_n) >= MIN_REPS
    t = int(tool_n)   >= MIN_REPS
    if b and t:  return "both"
    if b:        return "blimmp_only"
    if t:        return "tool_only"
    return "neither"


# ---------------------------------------------------------------------------
# Tool registry
# (key, loader_fn, threshold, mean_col_name, n_col_name)
# ---------------------------------------------------------------------------
TOOLS = [
    ("blimmp",    load_blimmp_rep,    THRESH["blimmp"],    "blimmp_mean_conf",    "blimmp_n_present"),
    ("anvio",     load_anvio_rep,     THRESH["anvio"],     "anvio_mean_comp",     "anvio_n_complete"),
    ("kemet",     load_kemet_rep,     THRESH["kemet"],     "kemet_mean_comp",     "kemet_n_complete"),
    ("metabolic", load_metabolic_rep, THRESH["metabolic"], "metabolic_mean_comp", "metabolic_n_complete"),
    ("dram",      load_dram_rep,      THRESH["dram"],      "dram_mean_comp",      "dram_n_complete"),
    ("mpp",       load_mpp_rep,       THRESH["mpp"],       "mpp_mean_pred",       "mpp_n_complete"),
]

TOOL_KEYS     = [t[0] for t in TOOLS]
TOOL_N_COLS   = {t[0]: t[4] for t in TOOLS}   # key → n_col name in merged DataFrame
COMPARISON_TOOLS = [k for k in TOOL_KEYS if k != "blimmp"]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
out_path = OUT_DIR / "scatter_data_0pct.xlsx"

with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
    for org in ORGANISMS:
        print(f"\n{'='*55}")
        print(f"Organism: {org}")

        # Build per-tool DataFrames and outer-join on module
        merged = None
        for key, loader_fn, thresh, mean_col, n_col in TOOLS:
            df = agg_tool(org, loader_fn, thresh, mean_col, n_col)
            n_nonzero = (df[mean_col] > 0).sum() if not df.empty and mean_col in df.columns else 0
            print(f"  {key:12s}: {len(df):4d} modules  ({n_nonzero} with non-zero value)")
            if merged is None:
                merged = df
            else:
                merged = pd.merge(merged, df, on="module", how="outer")

        if merged is None or merged.empty:
            print(f"  [SKIP] No data for {org}")
            continue

        # Fill numeric NAs with 0
        num_cols = [c for c in merged.columns if c != "module"]
        merged[num_cols] = merged[num_cols].fillna(0)

        # Category columns (BLIMMP vs each comparison tool)
        for tool_key in COMPARISON_TOOLS:
            n_col_t = TOOL_N_COLS[tool_key]
            merged[f"cat_vs_{tool_key}"] = merged.apply(
                lambda row, c=n_col_t: assign_cat(row["blimmp_n_present"], row[c]),
                axis=1,
            )

        # Restrict to bacterial module universe and sort
        merged = (merged[merged["module"].isin(BACTERIAL_MODULES)]
                  .sort_values("module")
                  .reset_index(drop=True))

        # Category summary
        print(f"\n  Modules in final table: {len(merged)}")
        for tool_key in COMPARISON_TOOLS:
            cat_col = f"cat_vs_{tool_key}"
            counts  = merged[cat_col].value_counts()
            print(f"  vs {tool_key:12s}: "
                  f"both={counts.get('both', 0):4d}  "
                  f"blimmp_only={counts.get('blimmp_only', 0):4d}  "
                  f"tool_only={counts.get('tool_only', 0):4d}  "
                  f"neither={counts.get('neither', 0):4d}")

        sheet = org[:31]   # Excel sheet name limit
        merged.to_excel(writer, sheet_name=sheet, index=False)
        print(f"\n  → Written to sheet '{sheet}'")

print(f"\nSaved: {out_path}")
