#!/usr/bin/env bash
#SBATCH --job-name=dram_spliced
#SBATCH --output=logs/dram_spliced_%A_%a.out
#SBATCH --error=logs/dram_spliced_%A_%a.err
#SBATCH --time=05:00:00
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
DRAM_SCRIPT="/xdisk/twheeler/nsontakke/Software/DRAM/scripts/DRAM.py"
FILTER_SCRIPT="$BASE_DIR/filter_masked_contigs.py"
OUT_ROOT="$BASE_DIR/DRAM_SPLICED"
mkdir -p "$BASE_DIR/logs" "$OUT_ROOT"

# ---- Discover every spliced genome fasta (mirrors the anvi'o discovery pattern) ----
mapfile -t all_fna < <(find "$BASE_DIR" -type d -name '*_SPLICED_GENOMES' -exec find {} -maxdepth 1 -type f -name '*_spliced.fna' \; | sort)
TOTAL=${#all_fna[@]}
echo "[$SLURM_ARRAY_TASK_ID] Total spliced .fna files found: $TOTAL"
if (( SLURM_ARRAY_TASK_ID >= TOTAL )); then
  echo "[$SLURM_ARRAY_TASK_ID] No file for this task index (TOTAL=$TOTAL)."
  exit 0
fi

fna_path="${all_fna[$SLURM_ARRAY_TASK_ID]}"
base=$(basename "$fna_path" _spliced.fna)   # e.g. MED4_50percentremoved_replicate2

# ---- Skip 100%-removed samples: fully N-masked genome, nothing for DRAM to call or annotate ----
if [[ "$base" == *"_100percentremoved_"* ]]; then
  echo "[$SLURM_ARRAY_TASK_ID] SKIP: $base is a 100% removal sample (fully N-masked, nothing to annotate)."
  exit 0
fi

# ---- Per-sample output dir, flat under DRAM_SPLICED/ ----
OUT_DIR="$OUT_ROOT/$base"
ANNOTATIONS_TSV="$OUT_DIR/annotations.tsv"
DISTILL_OUT="$OUT_DIR/distill"

# Filtered input lives as a sibling of OUT_DIR (not inside it), since
# DRAM.py refuses to run if its output dir already exists.
FILTERED_FASTA="$OUT_ROOT/${base}_filtered.fna"
FILTERED_REPORT="$OUT_ROOT/${base}_filtered_report.tsv"

echo "[$SLURM_ARRAY_TASK_ID] Sample     : $base"
echo "[$SLURM_ARRAY_TASK_ID] Input fna  : $fna_path"
echo "[$SLURM_ARRAY_TASK_ID] Out dir    : $OUT_DIR"

# ---- Step A: drop fully/near-fully N-masked contigs ----
# DRAM's rRNA step (barrnap -> nhmmer) fails outright on a contig with no
# real nucleotide signal to detect an alphabet from ("Invalid alphabet type
# in target for nhmmer"). This is a composition problem, not a length
# problem, so --min_contig_size doesn't help. See filter_masked_contigs.py
# for the threshold; default is 0.99 (essentially all-N only).
if [[ -s "$FILTERED_FASTA" && -s "$FILTERED_REPORT" ]]; then
  echo "[$SLURM_ARRAY_TASK_ID] SKIP: $base (filtered FASTA already exists)"
else
  python "$FILTER_SCRIPT" \
    -i "$fna_path" \
    -o "$FILTERED_FASTA" \
    -r "$FILTERED_REPORT"
fi

# ---- Step B: annotate, letting DRAM call its own genes from the nucleotide fasta ----
# --min_contig_size defaults to 2500bp in DRAM, which would silently drop a
# lot of real content here since many spliced contigs are well under that
# (e.g. ~933bp average in the Acinetobacter test run). Set to 1 to keep
# every real contig regardless of length; the N-fraction filter above
# already handled the contigs that were actually a problem.
if [[ -s "$ANNOTATIONS_TSV" ]]; then
  echo "[$SLURM_ARRAY_TASK_ID] SKIP: $base (annotations.tsv already exists)"
else
  # DRAM.py refuses to run if the output dir already exists, so remove any
  # partial/failed prior attempt before retrying.
  rm -rf "$OUT_DIR"
  python "$DRAM_SCRIPT" annotate \
    -i "$FILTERED_FASTA" \
    -o "$OUT_DIR" \
    --min_contig_size 1 \
    --threads "${SLURM_CPUS_PER_TASK}"
fi

# ---- Step B: distill ----
if [[ ! -s "$ANNOTATIONS_TSV" ]]; then
  echo "[$SLURM_ARRAY_TASK_ID] ERROR: annotations.tsv not found or empty for $base after annotate step." >&2
  exit 1
fi

if [[ -s "$DISTILL_OUT/metabolism_summary.xlsx" ]]; then
  echo "[$SLURM_ARRAY_TASK_ID] SKIP: $base already distilled"
else
  rm -rf "$DISTILL_OUT"
  TRNAS_TSV="$OUT_DIR/trnas.tsv"
  RRNAS_TSV="$OUT_DIR/rrnas.tsv"
  EXTRA_ARGS=()
  [[ -s "$TRNAS_TSV" ]] && EXTRA_ARGS+=(--trna_path "$TRNAS_TSV")
  [[ -s "$RRNAS_TSV" ]] && EXTRA_ARGS+=(--rrna_path "$RRNAS_TSV")
  python "$DRAM_SCRIPT" distill \
    -i "$ANNOTATIONS_TSV" \
    -o "$DISTILL_OUT" \
    "${EXTRA_ARGS[@]}"
fi

echo "[$SLURM_ARRAY_TASK_ID] DONE: $base"
