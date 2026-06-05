#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import pandas as pd
import numpy as np

INPUT_GLOB = "/xdisk/cgoubert/nsontakke/PROCESS_HMMER_OUTPUT/*_HMMER_grouped_hits_dedup.csv"
EVALUE_THRESHOLD = 1e-5
OUTPUT_TSV = "/xdisk/cgoubert/nsontakke/ATB_KO_Frequencies/ko_matrix.tsv"


def sample_from_filename(path: str) -> str:
    """Column name = filename without the trailing suffix '_HMMER_grouped_hits_dedup.csv'."""
    base = os.path.basename(path)
    suffix = "_HMMER_grouped_hits_dedup.csv"
    return base[:-len(suffix)] if base.endswith(suffix) else os.path.splitext(base)[0]


def load_ko_evalues(csv_path: str) -> pd.DataFrame:
    """
    Load a CSV with at least 'KO id' and 'E-value' columns.
    Returns a DataFrame with those two columns cleaned.
    """
    df = pd.read_csv(csv_path)
    # normalize column names just in case of stray spaces/casing
    cols = {c: c.strip() for c in df.columns}
    df = df.rename(columns=cols)

    if "KO id" not in df.columns or "E-value" not in df.columns:
        raise ValueError(f"Missing required columns in {csv_path}: found {list(df.columns)}")

    # Clean up
    df = df[["KO id", "E-value"]].copy()
    df["KO id"] = df["KO id"].astype(str).str.strip()
    df["E-value"] = pd.to_numeric(df["E-value"], errors="coerce")
    df = df.dropna(subset=["KO id"])
    return df


def main():
    files = sorted(glob.glob(INPUT_GLOB))
    if not files:
        raise SystemExit(f"No files matched: {INPUT_GLOB}")

    master_kos = set()
    sample_names = []
    for f in files:
        sample = sample_from_filename(f)
        sample_names.append(sample)
        df = load_ko_evalues(f)
        master_kos.update(df["KO id"].unique())

    # Build the output matrix initialized to 0
    index = sorted(master_kos)
    mat = pd.DataFrame(0, index=index, columns=sample_names, dtype=np.uint8)

    # -------- Pass 2: mark 1 where any hit for that KO has E-value < threshold --------
    for f, sample in zip(files, sample_names):
        df = load_ko_evalues(f)
        kos_pass = df.loc[df["E-value"] < EVALUE_THRESHOLD, "KO id"].dropna().astype(str).unique()
        if len(kos_pass) > 0:
            # Only set for those present in index (they should be)
            mat.loc[kos_pass, sample] = 1

    # Write result (row index is KO id)
    mat.index.name = "KO id"
    mat.to_csv(OUTPUT_TSV, sep="\t")
    print(f"Wrote {OUTPUT_TSV} with shape {mat.shape} "
          f"(KOs={mat.shape[0]}, files={mat.shape[1]}).")


if __name__ == "__main__":
    main()
