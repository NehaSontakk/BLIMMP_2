#!/usr/bin/env bash
#SBATCH --job-name=metabolic_spliced
#SBATCH --output=logs/metabolic_spliced_%A_%a.out
#SBATCH --error=logs/metabolic_spliced_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=twheeler
#SBATCH --partition=standard
#SBATCH --array=0-164
set -euo pipefail

# ---- 1) Initialize conda ----
if [[ -f "${HOME}/mambaforge/etc/profile.d/conda.sh" ]]; then
  source "${HOME}/mambaforge/etc/profile.d/conda.sh"
elif [[ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]]; then
  source "${HOME}/miniconda3/etc/profile.d/conda.sh"
else
  echo "ERROR: cannot find conda.sh" >&2
  exit 1
fi

# The METABOLIC_v4.0 env's conda activation hook
# (activate-binutils_linux-64.sh) references an unset variable
# (ADDR2LINE), which is harmless under normal bash but fatal under
# `set -u`. Temporarily disable -u around the activate call, then
# restore it for the rest of the script.
set +u
conda activate METABOLIC_v4.0
set -u

# ---- 2) Config ----
BASE_DIR="/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP"
METABOLIC_SCRIPT="/xdisk/twheeler/nsontakke/Software/METABOLIC_running_folder/METABOLIC/METABOLIC-G.pl"
METABOLIC_ROOT="$BASE_DIR/METABOLIC_SPLICED"
mkdir -p "$BASE_DIR/logs" "$METABOLIC_ROOT"

# ---- 3) Discover every spliced genome fasta (mirrors the anvi'o/DRAM discovery pattern) ----
mapfile -t all_fna < <(find "$BASE_DIR" -type d -name '*_SPLICED_GENOMES' -exec find {} -maxdepth 1 -type f -name '*_spliced.fna' \; | sort)
TOTAL=${#all_fna[@]}
echo "[$SLURM_ARRAY_TASK_ID] Total spliced .fna files found: $TOTAL"
if (( SLURM_ARRAY_TASK_ID >= TOTAL )); then
  echo "[$SLURM_ARRAY_TASK_ID] No file for this task index (TOTAL=$TOTAL)."
  exit 0
fi

fna_path="${all_fna[$SLURM_ARRAY_TASK_ID]}"
base=$(basename "$fna_path" _spliced.fna)   # e.g. MED4_50percentremoved_replicate2

OUT_DIR="$METABOLIC_ROOT/METABOLIC_${base}"
echo "[$SLURM_ARRAY_TASK_ID] Sample     : $base"
echo "[$SLURM_ARRAY_TASK_ID] Input fna  : $fna_path"
echo "[$SLURM_ARRAY_TASK_ID] Out dir    : $OUT_DIR"

if [[ -d "$OUT_DIR" && "$(ls -A "$OUT_DIR" 2>/dev/null)" ]]; then
  echo "[$SLURM_ARRAY_TASK_ID] SKIP: $base (output already present in $OUT_DIR)"
  exit 0
fi
mkdir -p "$OUT_DIR"

# ---- METABOLIC-G.pl expects a folder of genome files, not a single file path ----
# -in-gn is genome mode: it wants real nucleotide fasta and calls prodigal
# itself. The file inside the folder needs a .fasta extension for
# METABOLIC to recognize it as genome-mode input.
TMPDIR="$METABOLIC_ROOT/temp_${base}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"
cp "$fna_path" "$TMPDIR/total.fasta"

set +e
perl "$METABOLIC_SCRIPT" \
  -in-gn "$TMPDIR" \
  -o "$OUT_DIR" \
  -cpu "${SLURM_CPUS_PER_TASK}"
METABOLIC_EXIT=$?
set -e

rm -rf "$TMPDIR"

if [[ $METABOLIC_EXIT -ne 0 ]]; then
  echo "[$SLURM_ARRAY_TASK_ID] ERROR: METABOLIC failed on $base (exit code $METABOLIC_EXIT)" >&2
  # Remove the (incomplete) output dir so a later rerun doesn't get skipped
  # by the "already present" check above.
  rm -rf "$OUT_DIR"
  exit 1
fi

echo "[$SLURM_ARRAY_TASK_ID] DONE: $base"
