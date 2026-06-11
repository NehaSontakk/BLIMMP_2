#!/bin/bash
#SBATCH --job-name=metabolic_batch
#SBATCH --output=logs2/metabolic_%A_%a.out
#SBATCH --error=logs2/metabolic_%A_%a.err
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --account=twheeler
#SBATCH --partition=standard
#SBATCH --array=1-100          # fixed: was 1-71

source /home/u13/nsontakke/miniconda3/etc/profile.d/conda.sh
conda activate METABOLIC_v4.0

BATCH_FILE="/xdisk/twheeler/nsontakke/ATB_Analysis_0725/metabolic_batches/batch_${SLURM_ARRAY_TASK_ID}.txt"

echo "Batch ID:    $SLURM_ARRAY_TASK_ID"
echo "Batch file:  $BATCH_FILE"
echo "Commands:    $(wc -l < $BATCH_FILE)"

COMPLETED=0
FAILED=0
LINE_NUM=0

while IFS= read -r COMMAND; do
    ((LINE_NUM++))
    [[ -z "$COMMAND" ]] && continue

    echo "[$(date +%H:%M:%S)] Command $LINE_NUM..."
    eval "$COMMAND"

    if [ $? -eq 0 ]; then
        ((COMPLETED++))
    else
        echo "ERROR: Command $LINE_NUM failed: $COMMAND"
        ((FAILED++))
    fi
done < "$BATCH_FILE"

echo "Batch $SLURM_ARRAY_TASK_ID Summary: $COMPLETED completed, $FAILED failed"
echo "========================================"
