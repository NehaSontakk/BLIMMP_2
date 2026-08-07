#!/usr/bin/env bash
#SBATCH --job-name=tantan_prodigal
#SBATCH --output=logs/tantan_prodigal_%A_%a.out
#SBATCH --error=logs/tantan_prodigal_%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=twheeler
#SBATCH --partition=standard
#SBATCH --array=0-3

set -euo pipefail

# ---- Tools ----
tantan_exec="/xdisk/twheeler/nsontakke/Software/tantan-49/bin/tantan"
prodigal_exec="$HOME/bin/prodigal"   # confirmed via `which prodigal` -> ~/bin/prodigal

# ---- Base dir (run this script from /xdisk/twheeler/nsontakke/Removal_Study_BLIMMP) ----
BASE_DIR="/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP"

mkdir -p "$BASE_DIR/logs"

# ---- List of input genomes (one per array task) ----
# Array index -> fasta file
mapfile -t all_files < <(find "$BASE_DIR" -mindepth 2 -maxdepth 2 -type f -name '*.fna' | sort)
TOTAL=${#all_files[@]}

if (( SLURM_ARRAY_TASK_ID >= TOTAL )); then
  echo "[$SLURM_ARRAY_TASK_ID] No file for this task index (TOTAL=$TOTAL)."
  exit 0
fi

input_fa="${all_files[$SLURM_ARRAY_TASK_ID]}"
genome_dir=$(dirname "$input_fa")      # e.g. .../Acinetobacter_baumannii
base=$(basename "$input_fa" .fna)      # e.g. GCF_022459415

# ---- Per-genome output folder ----
OUT_DIR="$genome_dir/${base}_results"
mkdir -p "$OUT_DIR"

masked_fa="$OUT_DIR/${base}_masked.fa"
orfs_faa="$OUT_DIR/${base}_ORFs.faa"

echo "[$SLURM_ARRAY_TASK_ID] Input     : $input_fa"
echo "[$SLURM_ARRAY_TASK_ID] Output dir: $OUT_DIR"

# ---- Step 1: tantan masking ----
if [[ -s "$masked_fa" ]]; then
  echo "[$SLURM_ARRAY_TASK_ID] SKIP tantan: $base (already masked)"
else
  echo "[$SLURM_ARRAY_TASK_ID] RUN tantan : $base"
  "$tantan_exec" -x N "$input_fa" > "$masked_fa"
fi

# ---- Step 2: prodigal gene prediction ----
if [[ -s "$orfs_faa" ]]; then
  echo "[$SLURM_ARRAY_TASK_ID] SKIP prodigal: $base (already predicted)"
else
  echo "[$SLURM_ARRAY_TASK_ID] RUN prodigal : $base"
  "$prodigal_exec" -i "$masked_fa" -a "$orfs_faa"
fi

echo "[$SLURM_ARRAY_TASK_ID] Done: $base"
