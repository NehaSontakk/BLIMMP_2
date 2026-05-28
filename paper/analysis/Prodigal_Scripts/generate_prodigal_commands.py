#!/usr/bin/env python3
import os
import glob

def main():
    INPUT_DIR  = "/xdisk/twheeler/nsontakke/ATB_Analysis_0725/ATB_TANTAN"
    OUTPUT_DIR = "/xdisk/twheeler/nsontakke/ATB_Analysis_0725/ATB_PRODIGAL"
    CMD_FILE   = "prodigal_cmds.txt"

    # 1) List all masked.fa files
    masked_paths = sorted(glob.glob(os.path.join(INPUT_DIR, "*_masked.fa")))
    print(f"Found masked-fa files: {len(masked_paths)}")

    # 2) Build a set of base names already processed (any of .faa, .fna or .gff exists)
    processed = set()
    for ext in ("faa", "fna", "gff"):
        pattern = os.path.join(OUTPUT_DIR, f"*__ORFs.{ext}")
        # note: we’ll strip off the extension including “_ORFs”
        for path in glob.glob(pattern.replace("__", "")):
            base = os.path.basename(path)
            # e.g. “SAMD00000344_ORFs.fna” → “SAMD00000344”
            sample = base.split("_ORFs.")[0]
            processed.add(sample)

    # 3) Emit one prodigal command per unprocessed masked-fa
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    with open(CMD_FILE, "w") as out:
        for mp in masked_paths:
            base = os.path.basename(mp).rsplit("_masked.fa", 1)[0]
            if base in processed:
                continue
            # build command
            out_base = os.path.join(OUTPUT_DIR, f"{base}_ORFs")
            cmd = (
                f"prodigal -i '{mp}' "
                f"-a '{out_base}.faa' "
                f"-d '{out_base}.fna' "
                f"-f gff "
                f"-o '{out_base}.gff' "
                f"-p single"
            )
            out.write(cmd + "\n")

    # 4) Report
    total_cmds = sum(1 for _ in open(CMD_FILE))
    print(f"Commands generated: {total_cmds} (≤ {len(masked_paths)})")

if __name__ == "__main__":
    main()
