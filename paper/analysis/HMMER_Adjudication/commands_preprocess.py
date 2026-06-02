#!/usr/bin/env python3
from pathlib import Path

# Hard-coded paths (edit if needed)
INDIR = Path("/xdisk/cgoubert/nsontakke/ATB_Dechunked_HMMER")
OUTDIR = Path("/xdisk/cgoubert/nsontakke/PROCESS_HMMER_OUTPUT")
SCRIPT = "preprocess_domtblout.py"
OUTFILE = "commands.txt"

# Make sure output directory exists
OUTDIR.mkdir(parents=True, exist_ok=True)

commands = []
for f in sorted(INDIR.glob("*.domtblout")):
    base = f.stem
    out = OUTDIR / f"{base}_grouped_hits_dedup.csv"
    cmd = f"python {SCRIPT} {f} --out {out}"
    commands.append(cmd)

with open(OUTFILE, "w") as fh:
    fh.write("\n".join(commands) + "\n")

print(f"Wrote {len(commands)} commands to {OUTFILE}")
