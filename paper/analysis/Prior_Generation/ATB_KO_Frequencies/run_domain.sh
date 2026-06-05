#!/bin/bash
#SBATCH --job-name=ko_intersection_domain
#SBATCH --output=logs/ko_intersection_domain_%j.out
#SBATCH --error=logs/ko_intersection_domain_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G              # ← raise or lower as needed
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --account=twheeler

set -euo pipefail

source /home/u13/nsontakke/miniconda3/etc/profile.d/conda.sh
conda activate test

INPUT=/xdisk/cgoubert/nsontakke/ATB_KO_Frequencies/Lineage_Specific_Data/KO_Matrices/ko_matrix_sampleids_domain_level_priors.csv
/usr/bin/time -v python /xdisk/cgoubert/nsontakke/ATB_KO_Frequencies/ko_intersection_one.py "$INPUT"
