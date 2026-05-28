#!/usr/bin/env python3
import os
import glob

INDIR = "ATB_HMMER_Run2"
OUTDIR = "ATB_Dechunked_HMMER_hmmout"

os.makedirs(OUTDIR, exist_ok=True)

files = sorted(glob.glob(os.path.join(INDIR, "*_chunk*.hmmout")))
if not files:
    raise SystemExit("No *_chunk*.hmmout files found")

groups = {}
for f in files:
    name = os.path.basename(f)
    prefix = name.split("_chunk")[0]
    groups.setdefault(prefix, []).append(f)

for prefix, flist in groups.items():
    out_path = os.path.join(OUTDIR, f"{prefix}.hmmout")
    print(f"Combining {len(flist)} files → {out_path}")

    if os.path.exists(out_path):
        print(f"Skipping {prefix} (already exists)")
        continue


    with open(out_path, "w") as out:
        for i, fn in enumerate(sorted(flist)):
            with open(fn) as f:
                out.write(f.read())
            if i != len(flist) - 1:
                out.write("\n\n")  # blank line between chunks

print("Done.")
