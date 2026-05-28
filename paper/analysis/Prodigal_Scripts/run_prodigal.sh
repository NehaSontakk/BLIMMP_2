#!/usr/bin/env bash
#SBATCH --job-name=prodigal_array
#SBATCH --output=/xdisk/twheeler/nsontakke/ATB_Analysis_0725/RUN_PRODIGAL/logs/prodigal_array_%A_%a.out
#SBATCH --error=/xdisk/twheeler/nsontakke/ATB_Analysis_0725/RUN_PRODIGAL/logs/prodigal_array_%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=twheeler
#SBATCH --partition=standard
#SBATCH --array=1-500%500
set -euo pipefail

CMD_FILE="prodigal_cmds.txt"
TOTAL=$(wc -l < "$CMD_FILE")

# Compute chunk size so 1000 tasks cover all TOTAL commands
CHUNK_SIZE=$(( (TOTAL + 500 - 1) / 500 ))  # = ceil(TOTAL/1000), here ≈13

# Determine this task’s slice
IDX=$SLURM_ARRAY_TASK_ID
START=$(( (IDX - 1) * CHUNK_SIZE + 1 ))
END=$(( START + CHUNK_SIZE - 1 ))
(( END > TOTAL )) && END=$TOTAL

if (( START > TOTAL )); then
  echo "[$IDX] No commands to run (start $START > total $TOTAL)."
  exit 0
fi

echo "[$IDX] Running commands $START–$END (of $TOTAL)"

# Extract and run each
sed -n "${START},${END}p" "$CMD_FILE" | while read -r cmd; do
  echo "[$IDX] $cmd"
  eval "$cmd"
done

echo "[$IDX] Done."
