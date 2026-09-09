#!/bin/bash
#SBATCH --job-name=rm_metabolic_subsamples
#SBATCH --output=rm_metabolic_subsamples_%j.out
#SBATCH --error=rm_metabolic_subsamples_%j.err
#SBATCH --account=twheeler
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb
#SBATCH --time=02:00:00

set -euo pipefail

TARGET_DIR="/xdisk/twheeler/nsontakke/Removal_Study_BLIMMP/METABOLIC_SUBSAMPLES"

echo "Starting removal of: ${TARGET_DIR}"
echo "Job started at: $(date)"

if [ -d "${TARGET_DIR}" ]; then
    rm -rf "${TARGET_DIR}"
    echo "Successfully removed ${TARGET_DIR}"
else
    echo "Directory ${TARGET_DIR} does not exist. Nothing to do."
fi

echo "Job finished at: $(date)"
