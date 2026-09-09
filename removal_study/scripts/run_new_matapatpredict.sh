#!/usr/bin/env bash
#SBATCH --job-name=metapathpredict
#SBATCH --output=logs/metapathpredict_%j.out
#SBATCH --error=logs/metapathpredict_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=twheeler
#SBATCH --partition=standard
set -eo pipefail

# ---- Env ----
source /home/u13/nsontakke/miniconda3/etc/profile.d/conda.sh
conda activate MetaPathPredict

# ---- Config ----
BASE_DIR="/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP"
KOFAM_ROOT="$BASE_DIR/KOFAM_SPLICED_FIXED"
OUT_DIR="$BASE_DIR/METAPATHPREDICT_SPLICED"
mkdir -p "$BASE_DIR/logs" "$OUT_DIR"

OUT_FILE="$OUT_DIR/metapathpredict_all.tsv"

# ---- Check KofamScan outputs exist ----
mapfile -t ko_files < <(find "$KOFAM_ROOT" -name "*_kofam.tsv" | sort)
echo "Found ${#ko_files[@]} KofamScan output files"

if (( ${#ko_files[@]} == 0 )); then
    echo "ERROR: No KofamScan outputs found in $KOFAM_ROOT"
    echo "Run run_kofamscan_spliced.sh first."
    exit 1
fi

if [[ -s "$OUT_FILE" ]]; then
    echo "SKIP: Output already exists: $OUT_FILE"
    exit 0
fi

echo "Running MetaPathPredict on all ${#ko_files[@]} samples..."
MetaPathPredict \
    -i "${ko_files[@]}" \
    -a kofamscan \
    -o "$OUT_FILE"

echo "Done. Output: $OUT_FILE"
echo "Rows: $(wc -l < "$OUT_FILE")"
echo "Cols (modules): $(head -1 "$OUT_FILE" | tr '\t' '\n' | wc -l)"
