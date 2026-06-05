#!/usr/bin/env python3
import csv, glob
from pathlib import Path

BASE = Path("/xdisk/cgoubert/nsontakke/ATB_KO_Frequencies")

COUNTS_DIR   = BASE / "Lineage_Specific_Data" / "KO_Counts_Lineage_Specific"
SAMPLES_DIR  = BASE / "Lineage_Specific_Data" / "Lineage_Specific_SampleNames"
OUT_FREQ_DIR = BASE / "Lineage_Specific_Data" / "KO_Frequencies_Lineage_Specific"

OUT_FREQ_DIR.mkdir(parents=True, exist_ok=True)

# find all count files
count_files = sorted(glob.glob(str(COUNTS_DIR / "ko_counts_ko_matrix_sampleids_*.tsv")))
if not count_files:
    raise SystemExit(f"No KO count files found in {COUNTS_DIR}")

for f_count in count_files:
    f_count = Path(f_count)
    base = f_count.name[len("ko_counts_"):-len(".tsv")]  # ko_matrix_sampleids_xxx
    # corresponding sampleids file
    f_samples = SAMPLES_DIR / (base.replace("ko_matrix_", "") + ".csv")

    if not f_samples.exists():
        print(f"[WARN] no sample file for {base} ({f_samples.name})")
        continue

    # count samples = number of non-empty lines
    with f_samples.open() as fs:
        n_samples = sum(1 for line in fs if line.strip())

    if n_samples == 0:
        print(f"[WARN] {f_samples.name} has 0 samples, skipping")
        continue

    out_path = OUT_FREQ_DIR / f"ko_freq_{base}.tsv"
    with f_count.open() as fin, out_path.open("w", newline="") as fout:
        rdr = csv.reader(fin, delimiter="\t")
        wtr = csv.writer(fout, delimiter="\t")
        header = next(rdr)
        # find columns
        try:
            idx_ko = header.index("KOs")
        except ValueError:
            idx_ko = 0
        try:
            idx_cnt = header.index("KO_count")
        except ValueError:
            idx_cnt = 1

        # write new header
        wtr.writerow(["KOs","KO_count","KO_frequency"])
        for row in rdr:
            if not row: continue
            ko = row[idx_ko]
            cnt = int(row[idx_cnt])
            freq = cnt / n_samples
            wtr.writerow([ko, cnt, f"{freq:.10f}"])

    print(f"[OK] {f_count.name} ({n_samples} samples) -> {out_path.name}")
