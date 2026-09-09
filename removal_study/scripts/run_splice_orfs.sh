#!/usr/bin/env bash
#SBATCH --job-name=splice_orfs
#SBATCH --output=logs/splice_orfs_%A_%a.out
#SBATCH --error=logs/splice_orfs_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=twheeler
#SBATCH --partition=standard
#SBATCH --array=0-3
set -euo pipefail

# ---- Base dir (run this script from /xdisk/twheeler/nsontakke/Removal_Study_BLIMMP) ----
BASE_DIR="/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP"
SCRIPT="$BASE_DIR/splice_orfs_from_genome.py"
mkdir -p "$BASE_DIR/logs"

# ---- Discover genomes the same way as the tantan/prodigal + subsample scripts ----
mapfile -t all_files < <(find "$BASE_DIR" -mindepth 2 -maxdepth 2 -type f -name '*.fna' | sort)
TOTAL=${#all_files[@]}
if (( SLURM_ARRAY_TASK_ID >= TOTAL )); then
  echo "[$SLURM_ARRAY_TASK_ID] No file for this task index (TOTAL=$TOTAL)."
  exit 0
fi

input_fa="${all_files[$SLURM_ARRAY_TASK_ID]}"
genome_dir=$(dirname "$input_fa")           # e.g. .../MED4
genome_name=$(basename "$genome_dir")       # e.g. MED4
base=$(basename "$input_fa" .fna)           # e.g. GCF_000011465.1_ASM1146v1_genomic

# ---- Locate the prodigal ORFs and the subsamples produced by earlier steps ----
ORFS_FAA="$genome_dir/${base}_results/${base}_ORFs.faa"
SUBSAMPLES_DIR="$genome_dir/${genome_name}_SUBSAMPLES"

if [[ ! -s "$ORFS_FAA" ]]; then
  echo "[$SLURM_ARRAY_TASK_ID] ERROR: ORF file not found or empty: $ORFS_FAA" >&2
  echo "[$SLURM_ARRAY_TASK_ID] (Run the tantan/prodigal step first.)" >&2
  exit 1
fi
if [[ ! -d "$SUBSAMPLES_DIR" ]]; then
  echo "[$SLURM_ARRAY_TASK_ID] ERROR: subsamples dir not found: $SUBSAMPLES_DIR" >&2
  echo "[$SLURM_ARRAY_TASK_ID] (Run subsample_orfs.py first.)" >&2
  exit 1
fi

# ---- Flat per-genome output folder, named after the genome directory ----
OUT_DIR="$genome_dir/${genome_name}_SPLICED_GENOMES"

echo "[$SLURM_ARRAY_TASK_ID] Genome        : $genome_name"
echo "[$SLURM_ARRAY_TASK_ID] Original FNA  : $input_fa"
echo "[$SLURM_ARRAY_TASK_ID] ORFs file     : $ORFS_FAA"
echo "[$SLURM_ARRAY_TASK_ID] Subsamples dir: $SUBSAMPLES_DIR"
echo "[$SLURM_ARRAY_TASK_ID] Out dir       : $OUT_DIR"

python3 "$SCRIPT" \
  -f "$input_fa" \
  -a "$ORFS_FAA" \
  -s "$SUBSAMPLES_DIR" \
  -n "$genome_name" \
  -o "$OUT_DIR"

echo "[$SLURM_ARRAY_TASK_ID] Done: $genome_name"
