#!/usr/bin/env bash
#SBATCH --job-name=tantan_array
#SBATCH --output=logs/tantan_%A_%a.out
#SBATCH --error=logs/tantan_%A_%a.err
#SBATCH --time=32:00:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=twheeler
#SBATCH --partition=standard
#SBATCH --array=0-43

set -euo pipefail

tantan_exec="/xdisk/twheeler/nsontakke/Software/tantan-49/bin/tantan"
prodigal_exec="/xdisk/twheeler/nsontakke/Software/prodigal/bin/prodigal"

INPUT_DIR="/xdisk/twheeler/nsontakke/ATB_Analysis_0725/ATB_FASTA"
MASKED_DIR="/xdisk/twheeler/nsontakke/ATB_Analysis_0725/ATB_TANTAN"
PRODIGAL_DIR="/xdisk/twheeler/nsontakke/ATB_Analysis_0725/PRODIGAL_RESULTS"

mkdir -p "$MASKED_DIR" "$PRODIGAL_DIR"

mapfile -t all_files < <(find "$INPUT_DIR" -maxdepth 1 -type f -name '*.fa' | sort)
TOTAL=${#all_files[@]}

CHUNK=500
START=$(( SLURM_ARRAY_TASK_ID * CHUNK ))
END=$(( START + CHUNK - 1 ))
if (( START >= TOTAL )); then
  echo "[$SLURM_ARRAY_TASK_ID] No files for this task (START=$START ≥ TOTAL=$TOTAL)."
  exit 0
fi
(( END = END < TOTAL-1 ? END : TOTAL-1 ))

for idx in $(seq "$START" "$END"); do
  input_fa="${all_files[$idx]}"
  base=$(basename "$input_fa" .fa)
  masked_fa="$MASKED_DIR/${base}_masked.fa"

  if [[ -s "$masked_fa" ]]; then
    echo "[$SLURM_ARRAY_TASK_ID] SKIP: $base (already masked)"
  else
    echo "[$SLURM_ARRAY_TASK_ID] MASK : $base"
    "$tantan_exec" -x N "$input_fa" > "$masked_fa"
  fi
done
echo "[$SLURM_ARRAY_TASK_ID] Done processing files $((START+1))–$((END+1)) of $TOTAL."

