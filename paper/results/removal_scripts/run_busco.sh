#!/usr/bin/env bash
#SBATCH --job-name=busco_array
#SBATCH --output=logs/busco_%A_%a.out
#SBATCH --error=logs/busco_%A_%a.err
#SBATCH --mem-per-cpu=8GB
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=02:00:00
#SBATCH --account=twheeler
#SBATCH --partition=standard
#SBATCH --array=0-3

set -euo pipefail

### 1) Activate BUSCO env
eval "$(conda shell.bash hook)"
conda activate busco_env

### 2) Base dir + genome list (same discovery logic as the tantan/prodigal script)
BASE_DIR="/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP"
mkdir -p "$BASE_DIR/logs"

mapfile -t all_files < <(find "$BASE_DIR" -mindepth 2 -maxdepth 2 -type f -name '*.fna' | sort)
TOTAL=${#all_files[@]}

if (( SLURM_ARRAY_TASK_ID >= TOTAL )); then
  echo "[$SLURM_ARRAY_TASK_ID] No file for this task index (TOTAL=$TOTAL)."
  exit 0
fi

input_fa="${all_files[$SLURM_ARRAY_TASK_ID]}"
genome_dir=$(dirname "$input_fa")
base=$(basename "$input_fa" .fna)

# Prefer the tantan-masked fasta if it exists (from the prior step), else fall back to raw input
masked_fa="$genome_dir/${base}_results/${base}_masked.fa"
if [[ -s "$masked_fa" ]]; then
  FASTA="$masked_fa"
else
  FASTA="$input_fa"
fi

### 3) Per-genome output dir + lineage selection
OUT_DIR="$genome_dir/${base}_results/BUSCO_${base}"
mkdir -p "$OUT_DIR"

# Set lineage per organism folder; defaults to bacteria_odb10.
# Prochlorococcus (MED4, MIT9313) and Synechococcus (SS120) are cyanobacteria,
# so switch those to cyanobacteria_odb10 for a more informative run.
case "$genome_dir" in
  *MED4*|*MIT9313*|*SS120*)
    LINEAGE="cyanobacteria_odb10"
    ;;
  *)
    LINEAGE="bacteria_odb10"
    ;;
esac

RUN_NAME="busco_${base}"

echo "[$SLURM_ARRAY_TASK_ID] Genome  : $base"
echo "[$SLURM_ARRAY_TASK_ID] Input   : $FASTA"
echo "[$SLURM_ARRAY_TASK_ID] Lineage : $LINEAGE"
echo "[$SLURM_ARRAY_TASK_ID] Out dir : $OUT_DIR"
echo "[$(date +%Y-%m-%dT%H:%M:%S)] Running BUSCO on ${base}"

cd "${OUT_DIR}"

### 4) Run BUSCO
busco \
  -i "${FASTA}" \
  -o "${RUN_NAME}" \
  -l "${LINEAGE}" \
  -m genome \
  -c ${SLURM_NTASKS} \
  -f

if [ $? -ne 0 ]; then
  echo "[$SLURM_ARRAY_TASK_ID] ERROR: BUSCO failed for $base" >&2
  exit 1
fi

### 5) Extract completeness
SUMMARY_FILE=$(find "${RUN_NAME}" -maxdepth 1 -name "short_summary*.txt" | head -n 1)
if [[ -n "$SUMMARY_FILE" ]]; then
  grep "^\s*C:" "$SUMMARY_FILE" > "${RUN_NAME}_completeness.txt"
else
  echo "[$SLURM_ARRAY_TASK_ID] WARNING: could not locate short_summary file for $base" >&2
fi

echo "[$(date +%Y-%m-%dT%H:%M:%S)] [$SLURM_ARRAY_TASK_ID] DONE: $base"
