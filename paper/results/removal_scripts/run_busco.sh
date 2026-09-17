#!/usr/bin/env bash
#SBATCH --job-name=busco_spliced
#SBATCH --output=logs/busco_%A_%a.out
#SBATCH --error=logs/busco_%A_%a.err
#SBATCH --mem-per-cpu=8GB
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=02:00:00
#SBATCH --account=twheeler
#SBATCH --partition=standard
#SBATCH --array=0-164

set -euo pipefail

### 1) Activate BUSCO env
source /home/u13/nsontakke/miniconda3/etc/profile.d/conda.sh
conda activate busco_env

### 2) Organisms and genome discovery
BASE_DIR="/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP"
mkdir -p "$BASE_DIR/logs"

BUSCO_OUT_ROOT="$BASE_DIR/BUSCO_SPLICED"
mkdir -p "$BUSCO_OUT_ROOT"

CENTRAL_TSV="$BUSCO_OUT_ROOT/busco_completeness.tsv"

ORGANISMS=(
    "Acinetobacter_baumannii"
    "MED4"
    "MIT9313"
    "Pseudomonas_fluorescens_SBW25"
    "SS120"
)

# Build sorted list of all spliced FASTAs across organisms
mapfile -t all_files < <(
    for org in "${ORGANISMS[@]}"; do
        find "$BASE_DIR/$org/${org}_SPLICED_GENOMES" \
             -maxdepth 1 -name "*_spliced.fna" 2>/dev/null
    done | sort
)

TOTAL=${#all_files[@]}
echo "Total spliced genomes found: $TOTAL"

if (( SLURM_ARRAY_TASK_ID >= TOTAL )); then
    echo "[$SLURM_ARRAY_TASK_ID] No file for this task index (TOTAL=$TOTAL)."
    exit 0
fi

input_fa="${all_files[$SLURM_ARRAY_TASK_ID]}"
base=$(basename "$input_fa" _spliced.fna)

echo "[$SLURM_ARRAY_TASK_ID] Sample : $base"
echo "[$SLURM_ARRAY_TASK_ID] Input  : $input_fa"

### 3) Per-genome BUSCO output dir
OUT_DIR="$BUSCO_OUT_ROOT/${base}"
mkdir -p "$OUT_DIR"
RUN_NAME="busco_${base}"

# Skip if already completed
if [[ -d "$OUT_DIR/$RUN_NAME" ]] && \
   find "$OUT_DIR/$RUN_NAME" -maxdepth 1 -name "short_summary*.txt" | grep -q .; then
    echo "[$SLURM_ARRAY_TASK_ID] SKIP: BUSCO already completed for $base"
else
    ### 4) Lineage selection
    case "$base" in
        MED4*|MIT9313*|SS120*)
            LINEAGE="cyanobacteria_odb10"
            ;;
        *)
            LINEAGE="bacteria_odb10"
            ;;
    esac

    echo "[$SLURM_ARRAY_TASK_ID] Lineage: $LINEAGE"
    echo "[$(date +%Y-%m-%dT%H:%M:%S)] Running BUSCO"

    cd "$OUT_DIR"

    busco \
        -i "$input_fa" \
        -o "$RUN_NAME" \
        -l "$LINEAGE" \
        -m genome \
        -c "$SLURM_NTASKS" \
        -f

    echo "[$(date +%Y-%m-%dT%H:%M:%S)] [$SLURM_ARRAY_TASK_ID] BUSCO done: $base"
fi

### 5) Extract completeness → decimal, append to central TSV
SUMMARY_FILE=$(find "$OUT_DIR/$RUN_NAME" -maxdepth 1 -name "short_summary*.txt" 2>/dev/null | head -n 1)

if [[ -n "${SUMMARY_FILE:-}" ]]; then
    completeness_pct=$(grep -oP 'C:\K[\d.]+(?=%)' "$SUMMARY_FILE" | head -1)
    completeness_dec=$(python3 -c "print(round(${completeness_pct}/100, 4))")

    echo "[$SLURM_ARRAY_TASK_ID] Completeness: ${completeness_pct}% → ${completeness_dec}"

    # flock prevents concurrent-write corruption across array tasks
    (
        flock -x 200
        echo "${base},${completeness_dec}" >> "$CENTRAL_TSV"
    ) 200>"${CENTRAL_TSV}.lock"
else
    echo "[$SLURM_ARRAY_TASK_ID] WARNING: could not locate short_summary for $base" >&2
fi
