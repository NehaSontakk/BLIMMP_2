#!/bin/bash
set -euo pipefail

BASE="/xdisk/twheeler/nsontakke/ATB_Analysis_0725/ATB_METABOLIC"
OUT="metabolic_sample_module_lastcol.tsv"

echo -e "Sample_ID\tmodule_id\tpresence" > "$OUT"

for sample_dir in "$BASE"/*_ORFs_METABOLIC; do
    [ -d "$sample_dir" ] || continue

    f="$sample_dir/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet3.tsv"
    [ -f "$f" ] || { echo "WARNING: missing $f"; continue; }

    sample_id="$(basename "$sample_dir" _ORFs_METABOLIC)"

    awk -v sid="$sample_id" 'BEGIN{FS=OFS="\t"}
        NR==1 { next }                          # skip header
        NF != 4 { next }                        # must have exactly 4 columns
        $1 !~ /^M[0-9]{5}$/ { next }           # col1 must be valid module ID
        $4 !~ /^(Present|Absent)$/ { next }     # col4 must be exactly Present or Absent
        { print sid, $1, $4 }
    ' "$f" >> "$OUT"
done

echo "Done. Output: $OUT"
echo "Unique samples: $(awk -F'\t' 'NR>1{print $1}' "$OUT" | sort -u | wc -l)"
echo "Total rows: $(wc -l < "$OUT")"
