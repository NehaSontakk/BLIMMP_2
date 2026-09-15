#!/usr/bin/env bash
#SBATCH --job-name=anvio_spliced
#SBATCH --output=logs/anvio_spliced_%A_%a.out
#SBATCH --error=logs/anvio_spliced_%A_%a.err
#SBATCH --time=02:00:00
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
OUT_ROOT="$BASE_DIR/ANVIO_SPLICED"
mkdir -p "$BASE_DIR/logs" "$OUT_ROOT"

# ---- Discover every spliced genome fasta (mirrors the earlier .faa discovery pattern) ----
mapfile -t all_fna < <(find "$BASE_DIR" -type d -name '*_SPLICED_GENOMES' -exec find {} -maxdepth 1 -type f -name '*_spliced.fna' \; | sort)
TOTAL=${#all_fna[@]}
echo "[$SLURM_ARRAY_TASK_ID] Total spliced .fna files found: $TOTAL"
if (( SLURM_ARRAY_TASK_ID >= TOTAL )); then
  echo "[$SLURM_ARRAY_TASK_ID] No file for this task index (TOTAL=$TOTAL)."
  exit 0
fi

fna_path="${all_fna[$SLURM_ARRAY_TASK_ID]}"
base=$(basename "$fna_path" _spliced.fna)   # e.g. MED4_50percentremoved_replicate2

# ---- Per-sample output dir, flat under ANVIO_SPLICED/ ----
SAMPLE_DIR="$OUT_ROOT/$base"
mkdir -p "$SAMPLE_DIR"

REFORMATTED_FASTA="$SAMPLE_DIR/${base}_reformatted.fna"
NAME_REPORT="$SAMPLE_DIR/${base}_name_conversion.txt"
CONTIGS_DB="$SAMPLE_DIR/${base}-contigs.db"

echo "[$SLURM_ARRAY_TASK_ID] Sample     : $base"
echo "[$SLURM_ARRAY_TASK_ID] Input fna  : $fna_path"
echo "[$SLURM_ARRAY_TASK_ID] Out dir    : $SAMPLE_DIR"

# ---- Step A: simplify contig names ----
# anvi-gen-contigs-database rejects contig names with characters outside
# letters, digits, and underscore. NCBI accessions like NZ_CP023029.1
# contain a period, which anvi'o will not accept. Simplify names up front
# and keep the conversion report so contig identities can be traced back.
if [[ -s "$REFORMATTED_FASTA" && -s "$NAME_REPORT" ]]; then
  echo "[$SLURM_ARRAY_TASK_ID] SKIP: reformatted FASTA already exists"
else
  anvi-script-reformat-fasta \
    "$fna_path" \
    -o "$REFORMATTED_FASTA" \
    --simplify-names \
    --report-file "$NAME_REPORT"
fi

# ---- Step B: build contigs.db, letting anvi'o call its own genes ----
if [[ -s "$CONTIGS_DB" ]]; then
  echo "[$SLURM_ARRAY_TASK_ID] SKIP: contigs.db already exists"
else
  anvi-gen-contigs-database \
    --contigs-fasta "$REFORMATTED_FASTA" \
    --project-name "$base" \
    --output-db-path "$CONTIGS_DB" \
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

# ---- Sanity check output mode names against this anvi'o version, once per array (task 0 only) ----
# Flag names and available modes have changed across anvi'o versions; if any
# of these three aren't recognized, fail with a clear message pointing at
# the fix rather than a cryptic downstream argparse error on all tasks.
if [[ "$SLURM_ARRAY_TASK_ID" == "0" ]]; then
  AVAILABLE_MODES=$(anvi-estimate-metabolism -c "$CONTIGS_DB" --kegg-data-dir "$KEGG_DATA_DIR" --list-available-modes 2>&1 || true)
  for mode in hits modules module_paths; do
    if ! echo "$AVAILABLE_MODES" | grep -q "$mode"; then
      echo "[$SLURM_ARRAY_TASK_ID] ERROR: output mode '$mode' not found in this anvi'o version's --list-available-modes output:" >&2
      echo "$AVAILABLE_MODES" >&2
      exit 1
    fi
  done
  echo "[$SLURM_ARRAY_TASK_ID] Confirmed output modes hits, modules, module_paths are available."
fi

# ---- Step D: estimate metabolism ----
# hits           -> every enzyme/KO annotation anvi'o found (gene_caller_id, KO accession, e-value, ...)
# modules        -> per-module completeness / presence-absence
# module_paths   -> the specific path(s) of KOs anvi'o considers this module to have
OUT_PREFIX="$SAMPLE_DIR/${base}"
if [[ -s "${OUT_PREFIX}_modules.txt" ]]; then
  echo "[$SLURM_ARRAY_TASK_ID] SKIP: metabolism estimation output already exists"
else
  anvi-estimate-metabolism \
    -c "$CONTIGS_DB" \
    --kegg-data-dir "$KEGG_DATA_DIR" \
    --output-modes hits,modules,module_paths \
    --output-file-prefix "$OUT_PREFIX" \
    --include-zeros
fi

echo "[$SLURM_ARRAY_TASK_ID] DONE: $base"
echo "[$SLURM_ARRAY_TASK_ID] Outputs:"
echo "[$SLURM_ARRAY_TASK_ID]   ${OUT_PREFIX}_hits.txt         (all KO/enzyme IDs anvi'o found)"
echo "[$SLURM_ARRAY_TASK_ID]   ${OUT_PREFIX}_modules.txt      (modules + completeness)"
echo "[$SLURM_ARRAY_TASK_ID]   ${OUT_PREFIX}_module_paths.txt (KO paths per module)"
