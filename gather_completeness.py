#!/usr/bin/env python3
"""
Gather BUSCO completeness percentages from all short_summary files found
anywhere under the current directory, and write them to a CSV.

Usage:
    python gather_completeness.py

No arguments needed — it recursively searches the current working directory
for files matching "short_summary*.txt" and writes "completeness_stats.csv"
in the current directory.

Output CSV columns: sample_name, completeness
  - sample_name: derived from the BUSCO run folder name (the "busco_<NAME>"
    directory containing the short_summary file), with the path relative to
    the current directory as a fallback prefix for uniqueness.
  - completeness: e.g. "98.4%"
"""

import os
import glob
import re
import csv
import sys

OUTPUT_CSV = "completeness_stats.csv"


def main():
    input_dir = os.getcwd()

    # Recursively find every BUSCO short summary file under the current dir
    pattern = os.path.join(input_dir, "**", "short_summary*.txt")
    short_summary_paths = sorted(glob.glob(pattern, recursive=True))

    if not short_summary_paths:
        print(f"ERROR: No BUSCO summary files found under: {input_dir}", file=sys.stderr)
        sys.exit(1)

    with open(OUTPUT_CSV, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["sample_name", "completeness"])

        for filepath in short_summary_paths:
            # e.g. .../Acinetobacter_baumannii/GCF_022459415_results/BUSCO_GCF_022459415/
            #        busco_GCF_022459415/short_summary_busco_GCF_022459415.txt

            # 1) Derive FA_BASE from the parent directory name "busco_<FA_BASE>"
            busco_dir = os.path.basename(os.path.dirname(filepath))
            if busco_dir.startswith("busco_"):
                fa_base = busco_dir[len("busco_"):]
            else:
                # Fallback: extract from filename "short_summary_busco_<FA_BASE>.txt"
                filename = os.path.basename(filepath)
                fa_base = filename.replace("short_summary_busco_", "").replace(
                    "short_summary_", ""
                ).replace(".txt", "")

            # 2) Build a unique sample name using the path relative to cwd,
            #    so genomes in different subfolders don't collide.
            rel_path = os.path.relpath(filepath, input_dir)
            rel_parts = rel_path.split(os.sep)
            # Use everything above the short_summary filename and busco_<FA_BASE>
            # dir as context, joined with fa_base for uniqueness.
            context_parts = rel_parts[:-2] if len(rel_parts) > 2 else rel_parts[:-1]
            if context_parts:
                sample_name = "/".join(context_parts + [fa_base])
            else:
                sample_name = fa_base

            # 3) Parse the "C:xx.x%" completeness line
            completeness = ""
            with open(filepath, "r") as fh:
                for line in fh:
                    m = re.match(r"^\s*C:(\d+\.\d+)%", line)
                    if m:
                        completeness = m.group(1) + "%"
                        break

            if not completeness:
                print(f"Warning: no completeness found in {filepath}", file=sys.stderr)

            writer.writerow([sample_name, completeness])

    print(f"Found {len(short_summary_paths)} BUSCO summary file(s).")
    print(f"Wrote completeness stats to: {os.path.join(input_dir, OUTPUT_CSV)}")


if __name__ == "__main__":
    main()
