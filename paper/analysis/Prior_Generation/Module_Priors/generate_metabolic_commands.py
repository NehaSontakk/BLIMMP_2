#!/usr/bin/env python3
from pathlib import Path

BASE_DIR  = Path("/xdisk/twheeler/nsontakke/ATB_Analysis_0725")
PRODIGAL_DIR  = BASE_DIR / "ATB_PRODIGAL"
METABOLIC_DIR = BASE_DIR / "ATB_METABOLIC"
COMMANDS_FILE = BASE_DIR / "metabolic_commands.txt"

METABOLIC_DIR.mkdir(parents=True, exist_ok=True)

faa_files = sorted(PRODIGAL_DIR.glob("*.faa"))
print(f"Found {len(faa_files)} .faa files")

with open(COMMANDS_FILE, "w") as cmd_file:
    for faa_file in faa_files:
        subdir_path = METABOLIC_DIR / f"{faa_file.stem}_METABOLIC"
        subdir_path.mkdir(parents=True, exist_ok=True)
        # Copy faa only if not already there
        dest = subdir_path / faa_file.name
        if not dest.exists():
            import shutil
            shutil.copy2(faa_file, dest)
        cmd = (
            f"perl /xdisk/twheeler/nsontakke/Software/METABOLIC_running_folder/"
            f"METABOLIC/METABOLIC-G.pl "
            f"-in {subdir_path} -o {subdir_path}\n"
        )
        cmd_file.write(cmd)

print(f"Wrote {len(faa_files)} commands to {COMMANDS_FILE}")
