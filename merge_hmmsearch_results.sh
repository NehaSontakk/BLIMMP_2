#!/usr/bin/env bash
#SBATCH --job-name=merge_hmmsearch
#SBATCH --output=logs/merge_hmmsearch_%A_%a.out
#SBATCH --error=logs/merge_hmmsearch_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --mem=4GB
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=twheeler
#SBATCH --partition=standard
#SBATCH --array=0-4

set -euo pipefail

BASE_DIR="/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP"
OUT_ROOT="$BASE_DIR/HMMER_FULL_GENOME"
N_CHUNKS=50

mkdir -p "$BASE_DIR/logs"

# ---- Same genome discovery as the search-array script, so indices match ----
mapfile -t genome_faas < <(find "$BASE_DIR" -mindepth 3 -maxdepth 3 -type f -name '*_ORFs.faa' | sort)
N_GENOMES=${#genome_faas[@]}

if (( SLURM_ARRAY_TASK_ID >= N_GENOMES )); then
  echo "[$SLURM_ARRAY_TASK_ID] No genome for this task index (N_GENOMES=$N_GENOMES)."
  exit 0
fi

faa_path="${genome_faas[$SLURM_ARRAY_TASK_ID]}"
genome_name=$(basename "$(dirname "$(dirname "$faa_path")")")

CHUNK_DIR="$OUT_ROOT/$genome_name/chunks"
FINAL_TBLOUT="$OUT_ROOT/$genome_name/${genome_name}.tblout"
FINAL_DOMTBLOUT="$OUT_ROOT/$genome_name/${genome_name}.domtblout"

echo "[$SLURM_ARRAY_TASK_ID] Genome    : $genome_name"
echo "[$SLURM_ARRAY_TASK_ID] Chunk dir : $CHUNK_DIR"

mapfile -t chunk_tblouts < <(find "$CHUNK_DIR" -maxdepth 1 -type f -name "${genome_name}.chunk_*.tblout" | sort)
mapfile -t chunk_domtblouts < <(find "$CHUNK_DIR" -maxdepth 1 -type f -name "${genome_name}.chunk_*.domtblout" | sort)

if (( ${#chunk_domtblouts[@]} != N_CHUNKS )); then
  echo "[$SLURM_ARRAY_TASK_ID] ERROR: expected $N_CHUNKS domtblout chunks for $genome_name, found ${#chunk_domtblouts[@]}." >&2
  echo "[$SLURM_ARRAY_TASK_ID] (Some hmmsearch array tasks may still be running or failed.)" >&2
  exit 1
fi

# ---- Merge tblout: keep the 3-line HMMER comment header from the first chunk, then all data rows ----
head -3 "${chunk_tblouts[0]}" > "$FINAL_TBLOUT"
for f in "${chunk_tblouts[@]}"; do
  grep -v '^#' "$f" >> "$FINAL_TBLOUT" || true
done

# ---- Merge domtblout: same pattern ----
head -3 "${chunk_domtblouts[0]}" > "$FINAL_DOMTBLOUT"
for f in "${chunk_domtblouts[@]}"; do
  grep -v '^#' "$f" >> "$FINAL_DOMTBLOUT" || true
done

echo "[$SLURM_ARRAY_TASK_ID] Wrote: $FINAL_TBLOUT"
echo "[$SLURM_ARRAY_TASK_ID] Wrote: $FINAL_DOMTBLOUT"

# ---- Sanity check: total data-row count in merged file should equal sum across chunks ----
merged_rows=$(grep -vc '^#' "$FINAL_DOMTBLOUT" || true)
chunk_rows_sum=0
for f in "${chunk_domtblouts[@]}"; do
  rows=$(grep -vc '^#' "$f" || true)
  chunk_rows_sum=$(( chunk_rows_sum + rows ))
done

echo "[$SLURM_ARRAY_TASK_ID] Merged domtblout rows: $merged_rows (sum across chunks: $chunk_rows_sum)"
if (( merged_rows != chunk_rows_sum )); then
  echo "[$SLURM_ARRAY_TASK_ID] ERROR: row count mismatch after merge!" >&2
  exit 1
fi

echo "[$SLURM_ARRAY_TASK_ID] DONE: $genome_name"
