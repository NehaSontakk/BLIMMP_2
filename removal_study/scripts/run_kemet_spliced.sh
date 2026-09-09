#!/usr/bin/env bash
#SBATCH --job-name=kemet_spliced
#SBATCH --output=logs/kemet_spliced_%A_%a.out
#SBATCH --error=logs/kemet_spliced_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=twheeler
#SBATCH --partition=standard
#SBATCH --array=0-131
set -eo pipefail

# ---- Env ----
source /home/u13/nsontakke/miniconda3/etc/profile.d/conda.sh
conda activate kemet

# ---- Config ----
BASE_DIR="/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP"
KEMET_DIR="/xdisk/twheeler/nsontakke/Software/KEMET"
KOFAM_ROOT="$BASE_DIR/KOFAM_SPLICED"
OUT_ROOT="$BASE_DIR/KEMET_SPLICED"
mkdir -p "$BASE_DIR/logs" "$OUT_ROOT"

# ---- Discover all spliced .fna files (same order as kofamscan script) ----
mapfile -t all_fna < <(find "$BASE_DIR" -type d -name '*_SPLICED_GENOMES' \
    -exec find {} -maxdepth 1 -type f -name '*_spliced.fna' \; | sort)
TOTAL=${#all_fna[@]}
echo "[$SLURM_ARRAY_TASK_ID] Total spliced .fna files found: $TOTAL"
if (( SLURM_ARRAY_TASK_ID >= TOTAL )); then
    echo "[$SLURM_ARRAY_TASK_ID] No file for this task index."
    exit 0
fi

fna_path="${all_fna[$SLURM_ARRAY_TASK_ID]}"
base=$(basename "$fna_path" _spliced.fna)

KO_TSV="$KOFAM_ROOT/$base/${base}_kofam.tsv"
SAMPLE_DIR="$OUT_ROOT/$base"
mkdir -p "$SAMPLE_DIR"

echo "[$SLURM_ARRAY_TASK_ID] Sample  : $base"
echo "[$SLURM_ARRAY_TASK_ID] KO file : $KO_TSV"
echo "[$SLURM_ARRAY_TASK_ID] Out dir : $SAMPLE_DIR"

# ---- Check KofamScan output exists ----
if [[ ! -s "$KO_TSV" ]]; then
    echo "[$SLURM_ARRAY_TASK_ID] ERROR: KofamScan output not found: $KO_TSV" >&2
    echo "[$SLURM_ARRAY_TASK_ID] Run run_kofamscan_spliced.sh first." >&2
    exit 1
fi

# ---- Skip if already done ----
if [[ -s "$SAMPLE_DIR/${base}_modules.tsv" ]]; then
    echo "[$SLURM_ARRAY_TASK_ID] SKIP: KEMET output already exists"
    exit 0
fi

# ---- KEMET needs the annotation file in a dedicated input folder
#      named to match the genome basename. Use a per-job temp folder
#      to avoid parallel jobs clobbering each other. ----
KEMET_INPUT="$SAMPLE_DIR/kemet_input"
mkdir -p "$KEMET_INPUT"
cp "$KO_TSV" "$KEMET_INPUT/${base}_spliced.tsv"

echo "[$SLURM_ARRAY_TASK_ID] Running KEMET..."
cd "$KEMET_DIR"
python kemet.py "$fna_path" \
    -a kofamkoala \
    -I "$KEMET_INPUT" \
    -O "$SAMPLE_DIR" \
    --skip_hmm \
    --skip_gsmm \
    --log

# ---- Move any reports KEMET wrote to internal dirs into the sample dir ----
for d in reports_tsv reports_txt; do
    if [[ -d "$KEMET_DIR/$d" ]]; then
        find "$KEMET_DIR/$d" -name "${base}*" -exec mv {} "$SAMPLE_DIR/" \; 2>/dev/null || true
    fi
done

echo "[$SLURM_ARRAY_TASK_ID] DONE: $base"
echo "[$SLURM_ARRAY_TASK_ID] Outputs in: $SAMPLE_DIR"
ls "$SAMPLE_DIR/"
