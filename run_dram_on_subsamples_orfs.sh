#!/usr/bin/env bash
#SBATCH --job-name=dram_subsamples
#SBATCH --output=logs/dram_%A_%a.out
#SBATCH --error=logs/dram_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=6GB
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=twheeler
#SBATCH --partition=standard
#SBATCH --array=0-164

set -euo pipefail

# ---- Env ----
source /home/u13/nsontakke/miniconda3/etc/profile.d/conda.sh
conda activate my_dram_env

# ---- Config ----
BASE_DIR="/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP"
DRAM_SETUP_SCRIPT="/xdisk/twheeler/nsontakke/Software/DRAM/scripts/DRAM-setup.py"
DRAM_SCRIPT="/xdisk/twheeler/nsontakke/Software/DRAM/scripts/DRAM.py"
OUT_ROOT="$BASE_DIR/DRAM_SUBSAMPLES"

mkdir -p "$BASE_DIR/logs" "$OUT_ROOT"

# ---- Discover every subsampled ORF faa file across all genomes ----
# Matches the *_SUBSAMPLES output folders created by subsample_orfs.py, e.g.:
#   MED4/MED4_SUBSAMPLES/MED4_0percentremoved_replicate1.faa
#   Acinetobacter_baumannii/Acinetobacter_baumannii_SUBSAMPLES/Acinetobacter_baumannii_50percentremoved_replicate2.faa
mapfile -t all_faa < <(find "$BASE_DIR" -type d -name '*_SUBSAMPLES' -exec find {} -maxdepth 1 -type f -name '*.faa' \; | sort)
TOTAL=${#all_faa[@]}

echo "[$SLURM_ARRAY_TASK_ID] Total subsampled .faa files found: $TOTAL"

if (( SLURM_ARRAY_TASK_ID >= TOTAL )); then
  echo "[$SLURM_ARRAY_TASK_ID] No file for this task index (TOTAL=$TOTAL)."
  exit 0
fi

faa_path="${all_faa[$SLURM_ARRAY_TASK_ID]}"
base=$(basename "$faa_path" .faa)   # e.g. MED4_50percentremoved_replicate2

# ---- Per-sample output dir, flat under DRAM_SUBSAMPLES/ ----
OUT_DIR="$OUT_ROOT/$base"

echo "[$SLURM_ARRAY_TASK_ID] Input   : $faa_path"
echo "[$SLURM_ARRAY_TASK_ID] Out dir : $OUT_DIR"

if [[ -s "$OUT_DIR/annotations.tsv" ]]; then
  echo "[$SLURM_ARRAY_TASK_ID] SKIP: $base (annotations.tsv already exists)"
else
  # DRAM.py refuses to run if the output dir already exists, so remove any
  # partial/failed prior attempt before retrying.
  rm -rf "$OUT_DIR"

  python "$DRAM_SCRIPT" annotate_genes \
    -i "$faa_path" \
    -o "$OUT_DIR" \
    --threads "${SLURM_CPUS_PER_TASK}"
fi

echo "[$SLURM_ARRAY_TASK_ID] DONE: $base"
