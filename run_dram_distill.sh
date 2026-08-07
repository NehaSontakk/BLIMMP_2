#!/usr/bin/env bash
#SBATCH --job-name=dram_distill
#SBATCH --output=logs/dram_distill_%A_%a.out
#SBATCH --error=logs/dram_distill_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --cpus-per-task=2
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
DRAM_SCRIPT="/xdisk/twheeler/nsontakke/Software/DRAM/scripts/DRAM.py"
ANNOTATE_ROOT="$BASE_DIR/DRAM_SUBSAMPLES"

mkdir -p "$BASE_DIR/logs"

# ---- Same discovery as run_dram_array.sh, so indices line up ----
mapfile -t all_faa < <(find "$BASE_DIR" -type d -name '*_SUBSAMPLES' -exec find {} -maxdepth 1 -type f -name '*.faa' \; | sort)
TOTAL=${#all_faa[@]}

if (( SLURM_ARRAY_TASK_ID >= TOTAL )); then
  echo "[$SLURM_ARRAY_TASK_ID] No file for this task index (TOTAL=$TOTAL)."
  exit 0
fi

faa_path="${all_faa[$SLURM_ARRAY_TASK_ID]}"
base=$(basename "$faa_path" .faa)   # e.g. MED4_50percentremoved_replicate2

# ---- Skip 100%-removed samples: empty ORF set, nothing for DRAM to annotate ----
if [[ "$base" == *"_100percentremoved_"* ]]; then
  echo "[$SLURM_ARRAY_TASK_ID] SKIP: $base is a 100% removal sample (empty ORF set, no annotation to distill)."
  exit 0
fi

ANNOTATIONS_TSV="$ANNOTATE_ROOT/$base/annotations.tsv"
DISTILL_OUT="$ANNOTATE_ROOT/$base/distill"

echo "[$SLURM_ARRAY_TASK_ID] Sample         : $base"
echo "[$SLURM_ARRAY_TASK_ID] Annotations in : $ANNOTATIONS_TSV"
echo "[$SLURM_ARRAY_TASK_ID] Distill out    : $DISTILL_OUT"

if [[ ! -s "$ANNOTATIONS_TSV" ]]; then
  echo "[$SLURM_ARRAY_TASK_ID] ERROR: annotations.tsv not found or empty for $base -- run run_dram_array.sh first." >&2
  exit 1
fi

if [[ -s "$DISTILL_OUT/metabolism_summary.xlsx" ]]; then
  echo "[$SLURM_ARRAY_TASK_ID] SKIP: $base already distilled"
else
  rm -rf "$DISTILL_OUT"

  TRNAS_TSV="$ANNOTATE_ROOT/$base/trnas.tsv"
  RRNAS_TSV="$ANNOTATE_ROOT/$base/rrnas.tsv"

  EXTRA_ARGS=()
  [[ -s "$TRNAS_TSV" ]] && EXTRA_ARGS+=(--trna_path "$TRNAS_TSV")
  [[ -s "$RRNAS_TSV" ]] && EXTRA_ARGS+=(--rrna_path "$RRNAS_TSV")

  python "$DRAM_SCRIPT" distill \
    -i "$ANNOTATIONS_TSV" \
    -o "$DISTILL_OUT" \
    "${EXTRA_ARGS[@]}"
fi

echo "[$SLURM_ARRAY_TASK_ID] DONE: $base"
