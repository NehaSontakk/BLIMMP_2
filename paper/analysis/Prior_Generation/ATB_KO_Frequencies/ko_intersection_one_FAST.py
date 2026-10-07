#!/usr/bin/env python3
"""
Faster variant of ko_intersection_one.py.

The computation (read -> binarize -> sparse matrix -> X @ X.T) is IDENTICAL
to the original, unchanged. Only the output-writing step is different:
instead of pandas.DataFrame.sparse.from_spmatrix(M, ...).to_csv(...) (slow,
and the suspected source of a confirmed correctness bug for at least one
KO), this writes each row directly from the CSR matrix's own arrays using
numpy, row by row, with a plain buffered file write.

If this produces correct values where the original didn't, that confirms
the bug was in pandas' sparse-to-csv path, not the matrix multiplication.

Usage: identical to the original
    python ko_intersection_one_FAST.py /path/to/ko_matrix_*.csv
"""
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

BASE = Path("/xdisk/twheeler/nsontakke/BLIMMP_2/paper/analysis/Prior_Generation/ATB_KO_Frequencies")
OUT_INT_DIR = BASE / "Lineage_Specific_Data" / "KO_Intersections_Lineage_Specific"
OUT_INT_DIR.mkdir(parents=True, exist_ok=True)

if len(sys.argv) != 2:
    sys.exit("usage: ko_intersection_one_FAST.py /path/to/ko_matrix_*.csv")

fpath = Path(sys.argv[1])
base = fpath.name.rsplit(".", 1)[0]

t0 = time.time()
print(f"Reading {fpath.name} ...", file=sys.stderr)

# Same dtype-safety fix we used in the debug script -- apply uint8 only to
# value columns, not the index column, avoiding the earlier cast error.
with open(fpath) as f:
    header_cols = f.readline().rstrip("\n").split("\t")
value_cols = header_cols[1:]
dtype_map = {col: np.uint8 for col in value_cols}
df = pd.read_csv(fpath, sep="\t", index_col=0, dtype=dtype_map)
n_kos, n_samples = df.shape
print(f"  shape: {n_kos} KOs x {n_samples} samples ({time.time()-t0:.1f}s)", file=sys.stderr)

out_int = OUT_INT_DIR / f"ko_intersection_{base}.tsv"

if n_samples == 0:
    with open(out_int, "w") as f:
        f.write("\t" + "\t".join(df.index) + "\n")
        for ko in df.index:
            f.write(ko + "\t" + "\t".join(["0"] * n_kos) + "\n")
    print(f"[SKIP DATA] {fpath.name}: 0 samples -> wrote empty intersection.")
    sys.exit(0)

print("Binarizing...", file=sys.stderr)
bin_df = (df > 0).astype(np.uint8, copy=False)

print("Building sparse matrix (int32, NOT uint8 -- uint8 silently overflows "
      "above 255 during X @ X.T's accumulation, which is the confirmed root "
      "cause of the K99001/K00234/K15052 bugs: wrong values were exactly "
      "true_value mod 256)...", file=sys.stderr)
X = csr_matrix(bin_df.to_numpy(dtype=np.int32, copy=False))

print("Computing X @ X.T (unchanged from the original)...", file=sys.stderr)
t1 = time.time()
M = (X @ X.T).astype(np.int32)
M = M.tocsr()  # ensure CSR for fast row slicing below
print(f"  done ({time.time()-t1:.1f}s)", file=sys.stderr)

ko_names = list(bin_df.index)
n = len(ko_names)

print(f"Writing {n} rows directly (numpy, no pandas sparse-to-csv)...", file=sys.stderr)
t2 = time.time()
with open(out_int, "w") as f:
    f.write("\t" + "\t".join(ko_names) + "\n")
    row_buf = np.zeros(n, dtype=np.int32)
    for i in range(n):
        row_buf[:] = 0
        start, end = M.indptr[i], M.indptr[i + 1]
        row_buf[M.indices[start:end]] = M.data[start:end]
        f.write(ko_names[i] + "\t" + "\t".join(map(str, row_buf)) + "\n")
        if (i + 1) % 2000 == 0:
            print(f"  ...{i+1}/{n} rows written ({time.time()-t2:.1f}s so far)",
                  file=sys.stderr)

print(f"[OK] {fpath.name}: {n_kos} KOs x {n_samples} samples -> {out_int.name} "
      f"(write took {time.time()-t2:.1f}s, total {time.time()-t0:.1f}s)")
