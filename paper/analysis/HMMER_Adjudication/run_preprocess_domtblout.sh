#!/bin/bash
#SBATCH --job-name=domtblout_groups_serial
#SBATCH --output=logs/domtblout_groups_serial_%j.out
#SBATCH --error=logs/domtblout_groups_serial_%j.err
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00   # up to 7 days walltime (adjust as needed)
#SBATCH --account=twheeler
#SBATCH --partition=standard

set -euo pipefail

INDIR="/xdisk/cgoubert/nsontakke/ATB_Dechunked_HMMER"
OUTDIR="/xdisk/cgoubert/nsontakke/PROCESS_HMMER_OUTPUT"
SCRIPT="preprocess_domtblout.py"

mkdir -p logs "$OUTDIR"

# Load conda env
source /home/u13/nsontakke/miniconda3/etc/profile.d/conda.sh
conda activate test

# Iterate over every .domtblout file
for FILE in "$INDIR"/*.domtblout; do
    BASENAME=$(basename "$FILE" .domtblout)
    OUTFILE="$OUTDIR/${BASENAME}_grouped_hits_dedup.csv"

    if [[ -f "$OUTFILE" ]]; then
        echo "[$(date)] Skipping $FILE (already processed)"
        continue
    fi

    echo "[$(date)] Processing $FILE -> $OUTFILE"
    python "$SCRIPT" "$FILE" --out "$OUTFILE"
done

echo "[$(date)] All files processed."
