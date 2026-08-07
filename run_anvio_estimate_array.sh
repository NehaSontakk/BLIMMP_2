#!/usr/bin/env bash
#SBATCH --job-name=anvio_estimate
#SBATCH --output=logs/anvio_estimate_%A_%a.out
#SBATCH --error=logs/anvio_estimate_%A_%a.err
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
conda activate anvio-9

# ---- Config ----
BASE_DIR="/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP"
KEGG_DATA_DIR="/xdisk/twheeler/nsontakke/kegg-data"
OUT_ROOT="$BASE_DIR/ANVIO_SUBSAMPLES"

mkdir -p "$BASE_DIR/logs"

# ---- Same discovery as the build-array script, so indices line up ----
mapfile -t all_faa < <(find "$BASE_DIR" -type d -name '*_SUBSAMPLES' -exec find {} -maxdepth 1 -type f -name '*.faa' \; | sort)
TOTAL=${#all_faa[@]}

if (( SLURM_ARRAY_TASK_ID >= TOTAL )); then
  echo "[$SLURM_ARRAY_TASK_ID] No file for this task index (TOTAL=$TOTAL)."
  exit 0
fi

faa_path="${all_faa[$SLURM_ARRAY_TASK_ID]}"
base=$(basename "$faa_path" .faa)

SAMPLE_DIR="$OUT_ROOT/$base"
CONTIGS_DB="$SAMPLE_DIR/${base}-contigs.db"

if [[ ! -s "$CONTIGS_DB" ]]; then
  echo "[$SLURM_ARRAY_TASK_ID] ERROR: contigs.db not found for $base -- run run_anvio_build_array.sh first." >&2
  exit 1
fi

echo "[$SLURM_ARRAY_TASK_ID] Sample     : $base"
echo "[$SLURM_ARRAY_TASK_ID] Contigs DB : $CONTIGS_DB"

# ---- Sanity check output mode names against this anvi'o version, once per array (task 0 only) ----
# Flag names and available modes have changed across anvi'o versions; if any
# of these three aren't recognized, fail with a clear message pointing at
# the fix rather than a cryptic downstream argparse error on all 165 tasks.
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

# ---- Estimate metabolism: KO hits, module completeness, and module paths ----
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
