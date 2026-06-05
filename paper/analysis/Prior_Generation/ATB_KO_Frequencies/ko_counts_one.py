#!/usr/bin/env python3
import sys
from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path("/xdisk/cgoubert/nsontakke/ATB_KO_Frequencies")
OUT_CNT_DIR = BASE / "Lineage_Specific_Data" / "KO_Counts_Lineage_Specific"
OUT_CNT_DIR.mkdir(parents=True, exist_ok=True)

if len(sys.argv) != 2:
    sys.exit("usage: ko_counts_one.py /path/to/ko_matrix_*.csv")

fpath = Path(sys.argv[1])
base = fpath.name.rsplit(".", 1)[0]  # ko_matrix_*.csv -> ko_matrix_*

df = pd.read_csv(fpath, sep="\t", index_col=0)
n_kos, n_samples = df.shape

out_cnt = OUT_CNT_DIR / f"ko_counts_{base}.tsv"

if n_samples == 0:
    pd.DataFrame({"KO_count": np.zeros(n_kos, dtype=int)}, index=df.index).to_csv(out_cnt, sep="\t")
    print(f"[SKIP DATA] {fpath.name}: 0 samples -> wrote zero counts.")
    sys.exit(0)

bin_df = (df > 0).astype(np.uint8, copy=False)
ko_counts = bin_df.sum(axis=1).astype(int)
pd.DataFrame({"KO_count": ko_counts}, index=bin_df.index).to_csv(out_cnt, sep="\t")
print(f"[OK] {fpath.name}: {n_kos} KOs × {n_samples} samples -> {out_cnt.name}")
