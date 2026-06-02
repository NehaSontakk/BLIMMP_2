#!/usr/bin/env bash
#SBATCH --job-name=commands
#SBATCH --output=commands_%A_%a.out
#SBATCH --error=commands_%A_%a.err
#SBATCH --array=663-714%15  #1-715%10
#SBATCH --cpus-per-task=16
#SBATCH --account=twheeler
#SBATCH --partition=standard

module load parallel
source /home/u13/nsontakke/miniconda3/etc/profile.d/conda.sh
conda activate test


CMD_FILE="commands.txt"
TASK_ID=${SLURM_ARRAY_TASK_ID}
START_LINE=$(( (TASK_ID - 1) * 30 + 1 ))
END_LINE=$(( TASK_ID * 30 ))

sed -n "${START_LINE},${END_LINE}p" "$CMD_FILE" | parallel -j 16
