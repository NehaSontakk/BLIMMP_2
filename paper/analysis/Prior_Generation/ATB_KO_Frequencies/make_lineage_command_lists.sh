#!/usr/bin/env bash
set -euo pipefail

BASE="/xdisk/cgoubert/nsontakke/ATB_KO_Frequencies"
KO_DIR="${BASE}/Lineage_Specific_Data/KO_Matrices"

INT_PY="${BASE}/ko_intersection_one.py"
CNT_PY="${BASE}/ko_counts_one.py"

INT_CMDS="${BASE}/ko_intersection_lineage_commands.txt"
CNT_CMDS="${BASE}/ko_count_lineage_commands.txt"

: > "$INT_CMDS"
: > "$CNT_CMDS"

while IFS= read -r f; do
  echo "python \"$INT_PY\" \"$f\"" >> "$INT_CMDS"
  echo "python \"$CNT_PY\" \"$f\"" >> "$CNT_CMDS"
done < <(ls "${KO_DIR}"/ko_matrix_*.csv)

echo "Wrote:"
wc -l "$INT_CMDS" "$CNT_CMDS"
