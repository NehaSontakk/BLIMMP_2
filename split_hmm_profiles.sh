#!/usr/bin/env bash
#SBATCH --job-name=split_hmm
#SBATCH --output=logs/split_hmm_%j.out
#SBATCH --error=logs/split_hmm_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=twheeler
#SBATCH --partition=standard

set -euo pipefail

# ---- Config ----
HMM_DB="/xdisk/twheeler/nsontakke/kofamscan_data/profiles_hmmer/all_profiles.hmm"
CHUNK_DIR="/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP/HMM_CHUNKS"
N_CHUNKS=50

mkdir -p "$CHUNK_DIR" logs

echo "[split_hmm] Splitting $HMM_DB into $N_CHUNKS chunks -> $CHUNK_DIR"

# Each HMM record in an HMMER flatfile ends with a line containing only "//".
# We split on that boundary so no record is ever cut in half, then distribute
# whole records round-robin-free (in blocks) across N_CHUNKS output files.

TOTAL_PROFILES=$(grep -c '^NAME' "$HMM_DB")
echo "[split_hmm] Total profiles: $TOTAL_PROFILES"

awk -v chunk_dir="$CHUNK_DIR" -v n_chunks="$N_CHUNKS" -v total="$TOTAL_PROFILES" '
BEGIN {
    # ceil(total / n_chunks) profiles per chunk
    per_chunk = int((total + n_chunks - 1) / n_chunks)
    chunk_idx = 0
    profile_count_in_chunk = 0
    outfile = sprintf("%s/chunk_%03d.hmm", chunk_dir, chunk_idx)
}
{
    print $0 >> outfile
    if ($0 == "//") {
        profile_count_in_chunk++
        if (profile_count_in_chunk >= per_chunk && chunk_idx < n_chunks - 1) {
            close(outfile)
            chunk_idx++
            profile_count_in_chunk = 0
            outfile = sprintf("%s/chunk_%03d.hmm", chunk_dir, chunk_idx)
        }
    }
}
END {
    close(outfile)
}
' "$HMM_DB"

N_WRITTEN=$(ls "$CHUNK_DIR"/chunk_*.hmm | wc -l)
echo "[split_hmm] Wrote $N_WRITTEN chunk files"

# Sanity check: total NAME count across chunks should equal TOTAL_PROFILES
CHECK_TOTAL=$(grep -c '^NAME' "$CHUNK_DIR"/chunk_*.hmm | awk -F: '{sum+=$2} END{print sum}')
echo "[split_hmm] Profiles across all chunks: $CHECK_TOTAL (expected $TOTAL_PROFILES)"
if [[ "$CHECK_TOTAL" != "$TOTAL_PROFILES" ]]; then
  echo "[split_hmm] ERROR: profile count mismatch after splitting!" >&2
  exit 1
fi

echo "[split_hmm] DONE"
