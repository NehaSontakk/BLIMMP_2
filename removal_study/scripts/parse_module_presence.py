#!/usr/bin/env python3
"""
parse_module_presence.py
Extracts module-level presence/absence from anvi'o, METABOLIC, and DRAM
across all organisms, removal levels, and replicates.

Module identifiers per tool:
  anvi'o    — M00XXX IDs  (stepwise_module_is_complete == True)
  METABOLIC — M00XXX IDs  (total Module presence == Present, worksheet3)
  DRAM      — product.tsv column names (coverage >= DRAM_THRESHOLD)
              DRAM's own names; no M00XXX mapping attempted here.
              This is self-consistent for the recovery plot.

Output: module_presence.csv (long format)
  organism, removal_pct, replicate, tool, module_id
"""

import os
import csv
import pandas as pd

BASE_DIR = "/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP"

ORGANISMS = [
    "Acinetobacter_baumannii",
    "MED4",
    "MIT9313",
    "Pseudomonas_fluorescens_SBW25",
]
REMOVAL_LEVELS = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
REPLICATES     = [1, 2, 3]
DRAM_THRESHOLD = 0.75


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def sample_name(organism, removal_pct, replicate):
    return f"{organism}_{removal_pct}percentremoved_replicate{replicate}"


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_anvio(organism, removal_pct, replicate):
    """M00XXX modules where stepwise_module_is_complete == True."""
    sname = sample_name(organism, removal_pct, replicate)
    path  = os.path.join(BASE_DIR, "ANVIO_SPLICED", sname, f"{sname}_modules.txt")
    if not os.path.exists(path):
        return None
    present = set()
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("stepwise_module_is_complete", "").strip() == "True":
                present.add(row["module"].strip())
    return present


def parse_metabolic(organism, removal_pct, replicate):
    """M00XXX modules where total Module presence == Present (worksheet3)."""
    sname = sample_name(organism, removal_pct, replicate)
    path  = os.path.join(
        BASE_DIR, "METABOLIC_SPLICED",
        f"METABOLIC_{sname}",
        "METABOLIC_result_each_spreadsheet",
        "METABOLIC_result_worksheet3.tsv",
    )
    if not os.path.exists(path):
        return None
    present = set()
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("total Module presence", "").strip() == "Present":
                present.add(row["Module ID"].strip())
    return present


def parse_dram(organism, removal_pct, replicate):
    """
    Modules from DRAM product.tsv with coverage >= DRAM_THRESHOLD.
    Uses DRAM's own column names as module identifiers (no M00XXX mapping).
    Columns with non-numeric values (CAZy, etc.) are skipped.
    """
    sname = sample_name(organism, removal_pct, replicate)
    path  = os.path.join(BASE_DIR, "DRAM_SPLICED", sname, "distill", "product.tsv")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:   # one row per genome
            present = set()
            for col, val in row.items():
                if col == "genome":
                    continue
                try:
                    if float(val) >= DRAM_THRESHOLD:
                        present.add(col.strip())
                except (ValueError, TypeError):
                    pass
            return present   # only one genome row per file
    return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    records = []
    missing = []

    for organism in ORGANISMS:
        print(f"Processing {organism}...")
        for removal_pct in REMOVAL_LEVELS:
            for replicate in REPLICATES:
                sname = sample_name(organism, removal_pct, replicate)
                results = {
                    "anvio":     parse_anvio(organism, removal_pct, replicate),
                    "metabolic": parse_metabolic(organism, removal_pct, replicate),
                    "dram":      parse_dram(organism, removal_pct, replicate),
                }
                for tool, mods in results.items():
                    if mods is None:
                        missing.append(f"{tool}:{sname}")
                        continue
                    for mod_id in mods:
                        records.append({
                            "organism":    organism,
                            "removal_pct": removal_pct,
                            "replicate":   replicate,
                            "tool":        tool,
                            "module_id":   mod_id,
                        })

    df  = pd.DataFrame(records)
    out = "module_presence.csv"
    df.to_csv(out, index=False)
    print(f"\nWrote {len(df)} records -> {out}")

    if missing:
        print(f"\n{len(missing)} missing input files:")
        for m in sorted(missing):
            print(f"  {m}")

    # Summary: modules at 0% removal per tool per organism
    df0      = df[df["removal_pct"] == 0]
    counts   = df0.groupby(["organism", "tool", "module_id"]).size().reset_index(name="n")
    stable   = counts[counts["n"] >= 2]
    baseline = stable.groupby(["organism", "tool"])["module_id"].nunique()
    print("\nBaseline modules (majority-vote, 0% removal):")
    print(baseline.to_string())


if __name__ == "__main__":
    main()
