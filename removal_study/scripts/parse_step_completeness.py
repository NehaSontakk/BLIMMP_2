#!/usr/bin/env python3
"""
parse_step_completeness.py
Extracts continuous step-level completeness scores from anvi'o and METABOLIC
for modules that are present at 0% removal (majority vote: >= 2/3 replicates
above COMPLETENESS_THRESHOLD).

anvi'o  : reads stepwise_module_completeness (float 0-1) from _modules.txt.
METABOLIC: reads worksheet4 (long format, one row per step).
          Columns: 'Module step' (e.g. M00001+01), 'total Module step presence'
          (Present/Absent).  Computes n_present / n_total per module.

Output: step_completeness.csv  (long format)
  organism, removal_pct, replicate, tool, module_id, completeness

Requires: pandas
  pip install pandas --break-system-packages
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
REMOVAL_LEVELS  = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
REPLICATES      = [1, 2, 3]
COMPLETENESS_THRESHOLD = 0.75   # fraction of steps for a module to count as "present"


def sample_name(organism, removal_pct, replicate):
    return f"{organism}_{removal_pct}percentremoved_replicate{replicate}"


# ---------------------------------------------------------------------------
# Parsers — each returns {module_id: completeness_fraction} or None
# ---------------------------------------------------------------------------

def parse_anvio(organism, removal_pct, replicate):
    """
    Reads stepwise_module_completeness (preferred) or module_completeness
    from anvi'o _modules.txt.  Returns the maximum completeness seen for
    each module (a module can appear on multiple rows if it has several paths).
    """
    sname = sample_name(organism, removal_pct, replicate)
    path  = os.path.join(BASE_DIR, "ANVIO_SPLICED", sname, f"{sname}_modules.txt")
    if not os.path.exists(path):
        return None

    result = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            module_id = row.get("module", "").strip()
            if not module_id:
                continue
            # Prefer stepwise completeness; fall back to module_completeness
            raw = (row.get("stepwise_module_completeness")
                   or row.get("module_completeness", "")).strip()
            try:
                completeness = float(raw)
            except (ValueError, TypeError):
                continue
            if module_id not in result or completeness > result[module_id]:
                result[module_id] = completeness
    return result


def parse_metabolic(organism, removal_pct, replicate):
    """
    Reads METABOLIC worksheet4 (long format, one row per module step).
    'Module step' column: e.g. M00001+01
    'total Module step presence' column: Present / Absent
    Returns {module_id: n_present_steps / n_total_steps}.
    """
    sname = sample_name(organism, removal_pct, replicate)
    path  = os.path.join(
        BASE_DIR, "METABOLIC_SPLICED",
        f"METABOLIC_{sname}",
        "METABOLIC_result_each_spreadsheet",
        "METABOLIC_result_worksheet4.tsv",
    )
    if not os.path.exists(path):
        return None

    import collections
    present_counts = collections.Counter()
    total_counts   = collections.Counter()

    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            module_step = row.get("Module step", "").strip()
            if not module_step or "+" not in module_step:
                continue
            module_id = module_step.split("+")[0]
            presence  = row.get("total Module step presence", "").strip().lower()
            total_counts[module_id] += 1
            if presence == "present":
                present_counts[module_id] += 1

    return {
        mod: present_counts[mod] / total
        for mod, total in total_counts.items()
        if total > 0
    }


# ---------------------------------------------------------------------------
# Baseline: modules whose completeness >= threshold in >= 2/3 replicates at 0%
# ---------------------------------------------------------------------------

def find_baseline(df0):
    """
    df0 : rows with removal_pct == 0.
    Returns {(organism, tool): set of baseline module_ids}.
    """
    import collections
    counts = collections.Counter()
    for _, row in df0.iterrows():
        if row["completeness"] >= COMPLETENESS_THRESHOLD:
            counts[(row["organism"], row["tool"], row["module_id"])] += 1

    baseline = {}
    for (org, tool, mod), cnt in counts.items():
        if cnt >= 2:
            baseline.setdefault((org, tool), set()).add(mod)
    return baseline


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    records = []
    missing = []

    parsers = {
        "anvio":     parse_anvio,
        "metabolic": parse_metabolic,
    }

    for organism in ORGANISMS:
        print(f"Processing {organism}...")
        for removal_pct in REMOVAL_LEVELS:
            for replicate in REPLICATES:
                sname = sample_name(organism, removal_pct, replicate)
                for tool, fn in parsers.items():
                    res = fn(organism, removal_pct, replicate)
                    if res is None:
                        missing.append(f"{tool}:{sname}")
                        continue
                    for module_id, completeness in res.items():
                        records.append(dict(
                            organism    = organism,
                            removal_pct = removal_pct,
                            replicate   = replicate,
                            tool        = tool,
                            module_id   = module_id,
                            completeness= completeness,
                        ))

    df = pd.DataFrame(records)

    # Identify baseline modules
    baseline = find_baseline(df[df["removal_pct"] == 0])

    df["in_baseline"] = df.apply(
        lambda r: r["module_id"] in baseline.get((r["organism"], r["tool"]), set()),
        axis=1,
    )
    df_base = df[df["in_baseline"]].copy()

    out = "step_completeness.csv"
    df_base.to_csv(out, index=False)
    print(f"\nWrote {len(df_base)} records -> {out}")

    if missing:
        print(f"\n{len(missing)} missing input files:")
        for m in sorted(missing):
            print(f"  {m}")

    # Summary
    summary = (
        df_base[df_base["removal_pct"] == 0]
        .groupby(["organism", "tool"])["module_id"]
        .nunique()
    )
    print("\nBaseline module counts (majority-vote, 0% removal):")
    print(summary.to_string())



if __name__ == "__main__":
    main()
