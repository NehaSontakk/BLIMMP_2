#!/bin/bash
#SBATCH --job-name=ATB_untar
#SBATCH --output=logs/untar_%A_%a.out
#SBATCH --error=logs/untar_%A_%a.err
#SBATCH --array=1-56
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --partition=standard
#SBATCH --account=twheeler
mkdir -p logs

start=$(( (SLURM_ARRAY_TASK_ID - 1) * 500 + 1 ))
end=$(( SLURM_ARRAY_TASK_ID * 500 ))

# run the tar commands for that slice
sed -n "${start},${end}p" unzip_commands.txt | bash

