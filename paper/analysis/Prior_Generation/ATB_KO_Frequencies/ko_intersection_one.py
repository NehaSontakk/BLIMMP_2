#!/usr/bin/env python3
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix

BASE = Path("/xdisk/cgoubert/nsontakke/ATB_KO_Frequencies")
OUT_INT_DIR = BASE / "Lineage_Specific_Data" / "KO_Intersections_Lineage_Specific"
OUT_INT_DIR.mkdir(parents=True, exist_ok=True)

if len(sys.argv) != 2:
    sys.exit("usage: ko_intersection_one.py /path/to/ko_matrix_*.csv")

fpath = Path(sys.argv[1])
base = fpath.name.rsplit(".", 1)[0]  # ko_matrix_*.csv -> ko_matrix_*

df = pd.read_csv(fpath, sep="\t", index_col=0)
n_kos, n_samples = df.shape

out_int = OUT_INT_DIR / f"ko_intersection_{base}.tsv"

if n_samples == 0:
    pd.DataFrame(index=df.index, columns=df.index, dtype=np.int32).to_csv(out_int, sep="\t")
    print(f"[SKIP DATA] {fpath.name}: 0 samples -> wrote empty intersection.")
    sys.exit(0)

bin_df = (df > 0).astype(np.uint8, copy=False)
X = csr_matrix(bin_df.to_numpy(dtype=np.uint8, copy=False))
M = (X @ X.T).astype(np.int32)

pd.DataFrame.sparse.from_spmatrix(M, index=bin_df.index, columns=bin_df.index).to_csv(out_int, sep="\t")
print(f"[OK] {fpath.name}: {n_kos} KOs × {n_samples} samples -> {out_int.name}")
