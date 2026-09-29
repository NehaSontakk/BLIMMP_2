#!/usr/bin/env bash
#SBATCH --job-name=hmmsearch
#SBATCH --output=logs/hmmsearch_%A_%a.out
#SBATCH --error=logs/hmmsearch_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=twheeler
#SBATCH --partition=standard
#SBATCH --array=0-131

set -eo pipefail
# ---- Env ----
source /home/u13/nsontakke/miniconda3/etc/profile.d/conda.sh


# ---- Config ----
BASE_DIR="/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP"
CHUNK_DIR="$BASE_DIR/HMM_CHUNKS"
KOFAM_ROOT="$BASE_DIR/KOFAM_SPLICED"
OUT_ROOT="$BASE_DIR/HMMSEARCH_SPLICED"
N_CHUNKS=50
mkdir -p "$BASE_DIR/logs"

# ---- Discover spliced FAA files (same order as kofamscan script) ----
mapfile -t all_faa < <(find "$KOFAM_ROOT" -name "*_spliced.faa" | sort)
TOTAL=${#all_faa[@]}
echo "[$SLURM_ARRAY_TASK_ID] Total spliced .faa files found: $TOTAL"

if (( SLURM_ARRAY_TASK_ID >= TOTAL )); then
    echo "[$SLURM_ARRAY_TASK_ID] No file for this task index."
    exit 0
fi

faa_path="${all_faa[$SLURM_ARRAY_TASK_ID]}"
base=$(basename "$faa_path" _spliced.faa)

SAMPLE_DIR="$OUT_ROOT/$base"
CHUNK_OUT="$SAMPLE_DIR/chunks"
MERGED="$SAMPLE_DIR/${base}_spliced.domtblout"
MERGED_ALN="$SAMPLE_DIR/${base}_spliced.aln"
mkdir -p "$CHUNK_OUT"

echo "[$SLURM_ARRAY_TASK_ID] Sample : $base"
echo "[$SLURM_ARRAY_TASK_ID] FAA    : $faa_path"
echo "[$SLURM_ARRAY_TASK_ID] Out    : $SAMPLE_DIR"

# ---- Skip if both merged outputs already exist ----
if [[ -s "$MERGED" && -s "$MERGED_ALN" ]]; then
    echo "[$SLURM_ARRAY_TASK_ID] SKIP: merged domtblout and alignment already exist"
    exit 0
fi

# ---- HMM chunk list ----
mapfile -t chunk_files < <(find "$CHUNK_DIR" -maxdepth 1 -name "chunk_*.hmm" | sort)
if (( ${#chunk_files[@]} != N_CHUNKS )); then
    echo "[$SLURM_ARRAY_TASK_ID] ERROR: expected $N_CHUNKS chunks, found ${#chunk_files[@]}" >&2
    exit 1
fi

# ---- Run hmmsearch on each chunk ----
for chunk_path in "${chunk_files[@]}"; do
    chunk_tag=$(basename "$chunk_path" .hmm)
    domtblout="$CHUNK_OUT/${base}.${chunk_tag}.domtblout"
    tblout="$CHUNK_OUT/${base}.${chunk_tag}.tblout"
    aln_out="$CHUNK_OUT/${base}.${chunk_tag}.aln"

    if [[ -s "$domtblout" && -s "$aln_out" ]]; then
        echo "[$SLURM_ARRAY_TASK_ID] SKIP chunk: $chunk_tag"
        continue
    fi

    echo "[$SLURM_ARRAY_TASK_ID] Running hmmsearch: $chunk_tag"
    hmmsearch \
        --cpu "${SLURM_CPUS_PER_TASK}" \
        --tblout "$tblout" \
        --domtblout "$domtblout" \
        -o "$aln_out" \
        "$chunk_path" \
        "$faa_path"
done

# ---- Merge all chunk domtblouts into one ----
echo "[$SLURM_ARRAY_TASK_ID] Merging $N_CHUNKS chunk domtblouts..."
first=1
for chunk_path in "${chunk_files[@]}"; do
    chunk_tag=$(basename "$chunk_path" .hmm)
    domtblout="$CHUNK_OUT/${base}.${chunk_tag}.domtblout"
    if (( first )); then
        cat "$domtblout" >> "$MERGED"
        first=0
    else
        grep -v "^#" "$domtblout" >> "$MERGED" || true
    fi
done
echo "[$SLURM_ARRAY_TASK_ID] Merged domtblout: $MERGED"
echo "[$SLURM_ARRAY_TASK_ID] Hits: $(grep -v '^#' "$MERGED" | wc -l)"

# ---- Merge all chunk alignment outputs into one ----
# Header comment lines (#) are kept from the first chunk only.
# Subsequent chunks: skip everything before the first "Query:" line,
# so query blocks concatenate cleanly.
echo "[$SLURM_ARRAY_TASK_ID] Merging $N_CHUNKS chunk alignments..."
first=1
for chunk_path in "${chunk_files[@]}"; do
    chunk_tag=$(basename "$chunk_path" .hmm)
    aln_out="$CHUNK_OUT/${base}.${chunk_tag}.aln"
    if (( first )); then
        # Strip the trailing [ok] from the first file so it doesn't appear mid-stream
        grep -v "^\[ok\]" "$aln_out" >> "$MERGED_ALN" || true
        first=0
    else
        # Skip header comment block; start from first Query: line
        awk '/^Query:/{found=1} found && !/^\[ok\]/{print}' "$aln_out" >> "$MERGED_ALN" || true
    fi
done
# Write a single [ok] footer
echo "[ok]" >> "$MERGED_ALN"
echo "[$SLURM_ARRAY_TASK_ID] Merged alignment: $MERGED_ALN"

# ---- Clean up chunk files to save space ----
rm -rf "$CHUNK_OUT"

echo "[$SLURM_ARRAY_TASK_ID] DONE: $base"
