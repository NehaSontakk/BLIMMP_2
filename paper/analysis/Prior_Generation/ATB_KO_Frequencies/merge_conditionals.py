#!/usr/bin/env python3
import csv, glob, re
from pathlib import Path

BASE = Path("/xdisk/cgoubert/nsontakke/ATB_KO_Frequencies/Lineage_Specific_Data/KO_Conditionals_Pi_given_j")
OUT  = BASE / "merged_Pi_given_j.tsv"

# regex patterns to classify level
def detect_level(fname: str) -> str:
    if "phylum_level" in fname:
        return "phylum"
    elif "kingdom_level" in fname:
        return "kingdom"
    elif "domain_level" in fname:
        return "domain"
    else:
        return "other"

with OUT.open("w", newline="") as fout:
    w = csv.writer(fout, delimiter="\t")
    w.writerow(["KO_j","KO_i","Pi_given_j","Level","SourceFile"])
    for f in sorted(glob.glob(str(BASE / "neighbors_Pi_given_j_ko_matrix_sampleids_*_priors.tsv"))):
        level = detect_level(f)
        with open(f) as fin:
            rdr = csv.DictReader(fin, delimiter="\t")
            for row in rdr:
                w.writerow([row["KO_j"], row["KO_i"], row["P_i_given_j"], level, Path(f).name])
print(f"[OK] wrote merged file {OUT}")
