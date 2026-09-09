#!/usr/bin/env python3
"""
Fix KofamScan detail-tsv column names for MetaPathPredict compatibility.

KofamScan (exec_annotation) outputs a tab-separated file where the header is:
    #\tgene name\tKO\tthrshld\tscore\tE-value\t"KO definition"
    (# is its own first field, followed by a separator line starting with #)

MetaPathPredict expects columns named:
    gene_name, KO, adaptive_threshold, score, evalue, KO definition

Creates fixed copies in KOFAM_SPLICED_FIXED/ leaving originals untouched.
"""
import glob
import os
from pathlib import Path

BASE_DIR = "/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP"
KOFAM_ROOT = f"{BASE_DIR}/KOFAM_SPLICED"
FIXED_ROOT = f"{BASE_DIR}/KOFAM_SPLICED_FIXED"

files = sorted(glob.glob(f"{KOFAM_ROOT}/**/*_kofam.tsv", recursive=True))
print(f"Found {len(files)} KofamScan TSV files to fix")

fixed = 0
skipped = 0

for src in files:
    rel = os.path.relpath(src, KOFAM_ROOT)
    dst = Path(FIXED_ROOT) / rel
    dst.parent.mkdir(parents=True, exist_ok=True)

    if dst.exists() and dst.stat().st_size > 0:
        skipped += 1
        continue

    with open(src) as f:
        lines = f.readlines()

    out_lines = []
    for line in lines:
        if line.startswith("#"):
            parts = line.rstrip("\n").split("\t")
            # Header line: ['#', 'gene name', 'KO', 'thrshld', 'score', 'E-value', '"KO definition"']
            if len(parts) > 1 and "gene" in parts[1]:
                # Rewrite header without leading '#', rename columns
                new_parts = parts[1:]  # drop the '#' field
                rename = {
                    "gene name": "gene_name",
                    "thrshld":   "adaptive_threshold",
                    "E-value":   "evalue",
                    '"KO definition"': "KO definition",
                }
                new_parts = [rename.get(p.strip(), p.strip()) for p in new_parts]
                out_lines.append("\t".join(new_parts) + "\n")
            # else: separator line (# ------- ...) — skip it
        else:
            out_lines.append(line)

    with open(dst, "w") as f:
        f.writelines(out_lines)

    fixed += 1

print(f"Fixed: {fixed}  Skipped (already done): {skipped}")
print(f"\nVerify one file:")
import subprocess
sample = Path(files[0])
rel = os.path.relpath(str(sample), KOFAM_ROOT)
fixed_file = Path(FIXED_ROOT) / rel
result = subprocess.run(
    ["head", "-2", str(fixed_file)],
    capture_output=True, text=True
)
print(result.stdout)
