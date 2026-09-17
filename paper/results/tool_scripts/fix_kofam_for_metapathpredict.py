#!/usr/bin/env python3
"""
Fix KofamScan detail-tsv for MetaPathPredict compatibility.

MetaPathPredict (utils.py read_kofamscan_detailed_tsv) reads files in binary
mode and keeps ONLY rows where split("\\t")[0] == "*".  It assigns column names
BY POSITION — never from the file header:
    0 = surpassed_threshold (the "*")
    1 = gene_identifier
    2 = k_number
    3 = adaptive_threshold
    4 = score
    5 = evalue
    6 = definition

So FIXED files must preserve the "*\\t" prefix on every significant-hit row.
Space-separated KofamScan files are converted to tab-separated format while
keeping the "*" indicator.  The header comment is passed through unchanged
(MetaPathPredict ignores it).

Reads from KOFAM_SPLICED, always rewrites KOFAM_SPLICED_FIXED (overwrites).
"""
import glob
import os
import subprocess
from pathlib import Path

BASE_DIR   = "/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP"
KOFAM_ROOT = f"{BASE_DIR}/KOFAM_SPLICED"
FIXED_ROOT = f"{BASE_DIR}/KOFAM_SPLICED_FIXED"

# Number of data columns AFTER the "*" indicator
# gene_identifier, k_number, adaptive_threshold, score, evalue, definition
N_DATA_COLS = 6

files = sorted(glob.glob(f"{KOFAM_ROOT}/**/*_kofam.tsv", recursive=True))
print(f"Found {len(files)} KofamScan TSV files\n")

fixed = errors = 0
warned = []

for src in files:
    rel      = os.path.relpath(src, KOFAM_ROOT)
    dst      = Path(FIXED_ROOT) / rel
    dst.parent.mkdir(parents=True, exist_ok=True)

    try:
        with open(src) as f:
            lines = f.readlines()
    except Exception as e:
        print(f"  [ERROR] Cannot read {rel}: {e}")
        errors += 1
        continue

    out_lines  = []
    star_count = 0

    for line in lines:
        stripped = line.rstrip("\n")
        if not stripped.strip():
            continue

        if stripped.startswith("#"):
            # Header / comment — pass through unchanged.
            # MetaPathPredict ignores these (first tab-field != "*").
            out_lines.append(stripped + "\n")

        elif "\t" in stripped:
            # ---------- Tab-separated row ----------
            # Pass through EXACTLY as-is so MetaPathPredict can filter on [0]=="*".
            out_lines.append(stripped + "\n")
            first_field = stripped.split("\t")[0]
            if first_field == "*":
                star_count += 1

        else:
            # ---------- Space-separated row (wrapper-script output) ----------
            # Detect and preserve the "*" indicator, then convert to tab-sep.
            has_star = stripped.startswith("*")
            rest     = stripped[1:].strip() if stripped[0] in ("*", " ") else stripped.strip()
            parts    = rest.split(None, N_DATA_COLS - 1)

            # Pad short rows (missing definition, etc.)
            while len(parts) < N_DATA_COLS:
                parts.append("")
            parts = parts[:N_DATA_COLS]

            indicator = "*" if has_star else " "
            out_lines.append(indicator + "\t" + "\t".join(parts) + "\n")
            if has_star:
                star_count += 1

    if star_count == 0:
        msg = f"[WARN] No '*' rows in {rel} — MetaPathPredict will crash (empty DataFrame)"
        print(f"  {msg}")
        warned.append(rel)

    try:
        with open(dst, "w") as f:
            f.writelines(out_lines)
        fixed += 1
    except Exception as e:
        print(f"  [ERROR] Cannot write {dst}: {e}")
        errors += 1

print(f"\nDone — Fixed: {fixed}  Errors: {errors}  Zero-star-row warnings: {len(warned)}")
if warned:
    print("Files with no '*' rows (these samples had no KO hits above threshold):")
    for w in warned:
        print(f"  {w}")

# ── Verification: simulate what MetaPathPredict actually does ────────────────
if files:
    sample  = files[0]
    rel     = os.path.relpath(sample, KOFAM_ROOT)
    fixed_f = Path(FIXED_ROOT) / rel
    print(f"\nVerify ({fixed_f.name}):")

    result = subprocess.run(["head", "-4", str(fixed_f)], capture_output=True, text=True)
    print(result.stdout)

    result2 = subprocess.run(
        ["grep", "-c", r"^\*", str(fixed_f)], capture_output=True, text=True
    )
    print(f"  '*' rows in fixed file: {result2.stdout.strip()}")

    # Replay the exact MetaPathPredict read loop
    try:
        import pandas as pd

        mp_lines = []
        with open(fixed_f, "rb") as fh:
            for row in fh:
                if row.decode(errors="replace").split("\t")[0] == "*":
                    mp_lines.append(row.decode(errors="replace").split("\t"))

        if not mp_lines:
            print("  [FAIL] MetaPathPredict simulation: 0 '*' rows found → KeyError will occur")
        else:
            data = pd.DataFrame(mp_lines)
            data.rename(columns={
                0: "surpassed_threshold",
                1: "gene_identifier",
                2: "k_number",
                3: "adaptive_threshold",
                4: "score",
                5: "evalue",
                6: "definition",
            }, inplace=True)
            print(f"  [OK] MetaPathPredict simulation: {len(data)} rows")
            print(f"  Columns: {data.columns.tolist()}")
            required = ["adaptive_threshold", "score", "evalue"]
            missing  = [c for c in required if c not in data.columns]
            if missing:
                print(f"  [FAIL] Missing required columns: {missing}")
            else:
                print("  [OK] All required MetaPathPredict columns present")
    except Exception as e:
        print(f"  [ERROR] Simulation failed: {e}")
