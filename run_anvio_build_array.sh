#!/usr/bin/env bash
#SBATCH --job-name=anvio_subsamples
#SBATCH --output=logs/anvio_%A_%a.out
#SBATCH --error=logs/anvio_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=twheeler
#SBATCH --partition=standard
#SBATCH --array=0-164

set -euo pipefail

# ---- Env ----
source /home/u13/nsontakke/miniconda3/etc/profile.d/conda.sh
conda activate anvio-9

# ---- Config ----
BASE_DIR="/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP"
KEGG_DATA_DIR="/xdisk/twheeler/nsontakke/kegg-data"
GENE_CALLS_SCRIPT="$BASE_DIR/make_external_gene_calls.py"
OUT_ROOT="$BASE_DIR/ANVIO_SUBSAMPLES"

mkdir -p "$BASE_DIR/logs" "$OUT_ROOT"

# ---- Discover every subsampled ORF faa file (same pattern as DRAM array script) ----
mapfile -t all_faa < <(find "$BASE_DIR" -type d -name '*_SUBSAMPLES' -exec find {} -maxdepth 1 -type f -name '*.faa' \; | sort)
TOTAL=${#all_faa[@]}

echo "[$SLURM_ARRAY_TASK_ID] Total subsampled .faa files found: $TOTAL"

if (( SLURM_ARRAY_TASK_ID >= TOTAL )); then
  echo "[$SLURM_ARRAY_TASK_ID] No file for this task index (TOTAL=$TOTAL)."
  exit 0
fi

faa_path="${all_faa[$SLURM_ARRAY_TASK_ID]}"
base=$(basename "$faa_path" .faa)   # e.g. MED4_50percentremoved_replicate2

# ---- Per-sample output dir, flat under ANVIO_SUBSAMPLES/ ----
SAMPLE_DIR="$OUT_ROOT/$base"
mkdir -p "$SAMPLE_DIR"

GENE_CALLS_TSV="$SAMPLE_DIR/${base}_external_gene_calls.tsv"
PLACEHOLDER_FASTA="$SAMPLE_DIR/${base}_placeholder_contigs.fna"
CONTIGS_DB="$SAMPLE_DIR/${base}-contigs.db"

echo "[$SLURM_ARRAY_TASK_ID] Sample     : $base"
echo "[$SLURM_ARRAY_TASK_ID] Input faa  : $faa_path"
echo "[$SLURM_ARRAY_TASK_ID] Out dir    : $SAMPLE_DIR"

# ---- Step A: build external gene calls TSV + placeholder nucleotide FASTA from the faa ----
# anvi-gen-contigs-database always validates --contigs-fasta as nucleotide
# sequence, even when aa_sequence is supplied in the external-gene-calls TSV
# (it uses aa_sequence instead of translating, but the fasta itself still
# must pass the ACTGN character check). So we generate a placeholder
# nucleotide "contig" per protein rather than passing the raw .faa.
if [[ -s "$GENE_CALLS_TSV" && -s "$PLACEHOLDER_FASTA" ]]; then
  echo "[$SLURM_ARRAY_TASK_ID] SKIP: gene calls TSV + placeholder FASTA already exist"
else
  python "$GENE_CALLS_SCRIPT" -i "$faa_path" -o "$GENE_CALLS_TSV" -f "$PLACEHOLDER_FASTA"
fi

# ---- Step B: build contigs.db from the placeholder FASTA + external gene calls ----
if [[ -s "$CONTIGS_DB" ]]; then
  echo "[$SLURM_ARRAY_TASK_ID] SKIP: contigs.db already exists"
else
  anvi-gen-contigs-database \
    --contigs-fasta "$PLACEHOLDER_FASTA" \
    --project-name "$base" \
    --output-db-path "$CONTIGS_DB" \
    --external-gene-calls "$GENE_CALLS_TSV" \
    --ignore-internal-stop-codons \
    --num-threads "${SLURM_CPUS_PER_TASK}"
fi

# ---- Step C: annotate KOs (also stores module/BRITE membership) ----
if anvi-db-info "$CONTIGS_DB" 2>/dev/null | grep -q "KOfam"; then
  echo "[$SLURM_ARRAY_TASK_ID] SKIP: KOfam annotations already present"
else
  anvi-run-kegg-kofams \
    -c "$CONTIGS_DB" \
    --kegg-data-dir "$KEGG_DATA_DIR" \
    --num-threads "${SLURM_CPUS_PER_TASK}"
fi

echo "[$SLURM_ARRAY_TASK_ID] DONE: $base"
