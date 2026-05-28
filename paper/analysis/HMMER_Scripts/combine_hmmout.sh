#!/bin/bash
#SBATCH --job-name=combine_hmmout
#SBATCH --output=combine_hmmout.%j.out
#SBATCH --error=combine_hmmout.%j.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=twheeler
#SBATCH --partition=standard

set -euo pipefail

echo "Job started at $(date)"
echo "Running on node: $(hostname)"

python combine_hmmout.py

echo "Job finished at $(date)"
