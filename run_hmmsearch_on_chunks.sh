#!/usr/bin/env bash
#SBATCH --job-name=hmmsearch_array
#SBATCH --output=logs/hmmsearch_%A_%a.out
#SBATCH --error=logs/hmmsearch_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=twheeler
#SBATCH --partition=standard
#SBATCH --array=0-249

set -euo pipefail

# ---- Config ----
BASE_DIR="/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP"
CHUNK_DIR="$BASE_DIR/HMM_CHUNKS"
OUT_ROOT="$BASE_DIR/HMMER_FULL_GENOME"
N_CHUNKS=50

mkdir -p "$BASE_DIR/logs" "$OUT_ROOT"

# ---- Genome list (5 genomes, matches your `*/*_results/*_ORFs.faa` listing) ----
mapfile -t genome_faas < <(find "$BASE_DIR" -mindepth 3 -maxdepth 3 -type f -name '*_ORFs.faa' | sort)
N_GENOMES=${#genome_faas[@]}

if (( N_GENOMES == 0 )); then
  echo "ERROR: no *_ORFs.faa files found under $BASE_DIR" >&2
  exit 1
fi

# ---- HMM chunk list ----
mapfile -t chunk_files < <(find "$CHUNK_DIR" -maxdepth 1 -type f -name 'chunk_*.hmm' | sort)
if (( ${#chunk_files[@]} != N_CHUNKS )); then
  echo "ERROR: expected $N_CHUNKS chunk files in $CHUNK_DIR, found ${#chunk_files[@]}." >&2
  echo "Run split_hmm_profiles.sh first." >&2
  exit 1
fi

# ---- Map SLURM_ARRAY_TASK_ID -> (genome_idx, chunk_idx) ----
# Layout: genome varies slowest, chunk fastest.
#   task 0..49   -> genome 0, chunks 0..49
#   task 50..99  -> genome 1, chunks 0..49
#   etc.
TASK_ID=$SLURM_ARRAY_TASK_ID
genome_idx=$(( TASK_ID / N_CHUNKS ))
chunk_idx=$(( TASK_ID % N_CHUNKS ))

if (( genome_idx >= N_GENOMES )); then
  echo "[$TASK_ID] No genome for this task index (genome_idx=$genome_idx >= N_GENOMES=$N_GENOMES)."
  exit 0
fi

faa_path="${genome_faas[$genome_idx]}"
chunk_path="${chunk_files[$chunk_idx]}"

# Genome name = the directory two levels up from the .faa
# e.g. .../Acinetobacter_baumannii/GCF_022459415_results/GCF_022459415_ORFs.faa -> Acinetobacter_baumannii
genome_name=$(basename "$(dirname "$(dirname "$faa_path")")")
chunk_tag=$(basename "$chunk_path" .hmm)   # e.g. chunk_007

OUT_DIR="$OUT_ROOT/$genome_name/chunks"
mkdir -p "$OUT_DIR"

tblout_out="$OUT_DIR/${genome_name}.${chunk_tag}.tblout"
domtblout_out="$OUT_DIR/${genome_name}.${chunk_tag}.domtblout"

echo "[$TASK_ID] Genome : $genome_name"
echo "[$TASK_ID] Chunk  : $chunk_tag ($chunk_path)"
echo "[$TASK_ID] FAA    : $faa_path"
echo "[$TASK_ID] Out    : $domtblout_out"

if [[ -s "$domtblout_out" ]]; then
  echo "[$TASK_ID] SKIP: output already exists"
else
  hmmsearch \
    --noali \
    --cpu "${SLURM_CPUS_PER_TASK}" \
    --tblout "$tblout_out" \
    --domtblout "$domtblout_out" \
    -o /dev/null \
    "$chunk_path" \
    "$faa_path"
fi

echo "[$TASK_ID] DONE"
