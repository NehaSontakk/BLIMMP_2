#!/usr/bin/env python3
# make_lineage_specific_ko_matrices.py
import csv, os, glob
from pathlib import Path

BASE = Path("/xdisk/cgoubert/nsontakke/ATB_KO_Frequencies")
KO_MATRIX = BASE / "ko_matrix.tsv"
SAMPLES_DIR = BASE / "Lineage_Specific_Data" / "Lineage_Specific_SampleNames"
OUT_DIR = BASE / "Lineage_Specific_Data" / "KO_Matrices"

OUT_DIR.mkdir(parents=True, exist_ok=True)

# 1) Gather sample-id files
sample_files = sorted(glob.glob(str(SAMPLES_DIR / "sampleids_*_priors.csv")))
if not sample_files:
    raise SystemExit(f"No sample-id files found in {SAMPLES_DIR}")

# 2) Read header once, map sample name -> column index
with KO_MATRIX.open(newline="") as fin:
    reader = csv.reader(fin, delimiter="\t")
    header = next(reader)
    name_to_idx = {name: i for i, name in enumerate(header)}
    ko_idx = 0  # first column is KO id

# 3) Build writers (one per sample-id list), writing headers up front
writers = {}          # file_path -> (writer, keep_idxs)
out_handles = {}      # to close later
stats = {}            # out_path -> (n_requested, n_written)

for sf in sample_files:
    # load desired sample IDs for this group
    with open(sf) as f:
        want = {line.strip() for line in f if line.strip()}

    # columns to keep: KO id + any wanted samples present in header
    keep_idxs = [ko_idx] + [name_to_idx[s] for s in want if s in name_to_idx and name_to_idx[s] != ko_idx]

    out_path = OUT_DIR / f"ko_matrix_{Path(sf).name}"
    fout = out_path.open("w", newline="")
    out_handles[out_path] = fout
    w = csv.writer(fout, delimiter="\t")

    # write header
    w.writerow([header[i] for i in keep_idxs])
    writers[out_path] = (w, keep_idxs)

    stats[out_path] = (len(want), len(keep_idxs))  # requested vs written (with KO column)

# 4) Stream the big matrix ONCE and write selected columns for each output
with KO_MATRIX.open(newline="") as fin:
    reader = csv.reader(fin, delimiter="\t")
    next(reader)  # skip header already read above
    for row in reader:
        for out_path, (w, keep_idxs) in writers.items():
            w.writerow([row[i] for i in keep_idxs])

# 5) Close files
for fout in out_handles.values():
    fout.close()

# 6) Print verification summary
print("\nSummary of written files:")
for out_path, (n_req, n_written) in stats.items():
    print(f"{out_path.name}: requested {n_req} sample IDs -> wrote {n_written} columns (incl. KO)")
