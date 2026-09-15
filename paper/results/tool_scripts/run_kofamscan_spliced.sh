#!/usr/bin/env bash
#SBATCH --job-name=kofamscan_spliced
#SBATCH --output=logs/kofamscan_spliced_%A_%a.out
#SBATCH --error=logs/kofamscan_spliced_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=twheeler
#SBATCH --partition=standard
#SBATCH --array=0-131
set -euo pipefail

# ---- Env ----
source /home/u13/nsontakke/miniconda3/etc/profile.d/conda.sh
conda activate kofam_microbe

# ---- Config ----
BASE_DIR="/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP"
KOFAM_DB="/xdisk/twheeler/nsontakke/kofamscan_data"
KOFAM_CFG="$KOFAM_DB/config.yml"
OUT_ROOT="$BASE_DIR/KOFAM_SPLICED"
mkdir -p "$BASE_DIR/logs" "$OUT_ROOT"

# ---- Write config if missing ----
if [[ ! -f "$KOFAM_CFG" ]]; then
    cat > "$KOFAM_CFG" <<EOF
profile: $KOFAM_DB/profiles
ko_list: $KOFAM_DB/ko_list
cpu: ${SLURM_CPUS_PER_TASK}
hmmer_options: "--cut_tc"
EOF
fi

# ---- Discover all spliced .fna files ----
mapfile -t all_fna < <(find "$BASE_DIR" -type d -name '*_SPLICED_GENOMES' \
    -exec find {} -maxdepth 1 -type f -name '*_spliced.fna' \; | sort)
TOTAL=${#all_fna[@]}
echo "[$SLURM_ARRAY_TASK_ID] Total spliced .fna files found: $TOTAL"
if (( SLURM_ARRAY_TASK_ID >= TOTAL )); then
    echo "[$SLURM_ARRAY_TASK_ID] No file for this task index."
    exit 0
fi

fna_path="${all_fna[$SLURM_ARRAY_TASK_ID]}"
base=$(basename "$fna_path" _spliced.fna)   # e.g. MED4_50percentremoved_replicate2

SAMPLE_DIR="$OUT_ROOT/$base"
mkdir -p "$SAMPLE_DIR"

FAA="$SAMPLE_DIR/${base}_spliced.faa"
KO_TSV="$SAMPLE_DIR/${base}_kofam.tsv"

echo "[$SLURM_ARRAY_TASK_ID] Sample  : $base"
echo "[$SLURM_ARRAY_TASK_ID] FNA     : $fna_path"
echo "[$SLURM_ARRAY_TASK_ID] Out dir : $SAMPLE_DIR"

# ---- Step 1: Prodigal → protein FASTA ----
if [[ -s "$FAA" ]]; then
    echo "[$SLURM_ARRAY_TASK_ID] SKIP prodigal: $FAA already exists"
else
    echo "[$SLURM_ARRAY_TASK_ID] Running Prodigal..."
    prodigal \
        -i "$fna_path" \
        -a "$FAA" \
        -p single \
        -f gff \
        -o /dev/null \
        -q
    echo "[$SLURM_ARRAY_TASK_ID] Prodigal done: $FAA"
fi

# ---- Step 2: KofamScan → KO annotations ----
if [[ -s "$KO_TSV" ]]; then
    echo "[$SLURM_ARRAY_TASK_ID] SKIP KofamScan: $KO_TSV already exists"
else
    echo "[$SLURM_ARRAY_TASK_ID] Running KofamScan..."
    TMP_DIR="$SAMPLE_DIR/tmp"
    mkdir -p "$TMP_DIR"
    exec_annotation \
        -c "$KOFAM_CFG" \
        --format detail-tsv \
        --cpu "$SLURM_CPUS_PER_TASK" \
        --tmp-dir "$TMP_DIR" \
        -o "$KO_TSV" \
        "$FAA"
    rm -rf "$TMP_DIR"
    echo "[$SLURM_ARRAY_TASK_ID] KofamScan done: $KO_TSV"
fi

echo "[$SLURM_ARRAY_TASK_ID] DONE: $base"
