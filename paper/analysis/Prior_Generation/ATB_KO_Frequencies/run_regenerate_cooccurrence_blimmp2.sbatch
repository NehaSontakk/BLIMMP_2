#!/bin/bash
#SBATCH --job-name=regen_cooccurrence
#SBATCH --output=logs/regen_cooccurrence_%j.out
#SBATCH --error=logs/regen_cooccurrence_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --partition=standard
#SBATCH --account=twheeler

set -euo pipefail

cd /xdisk/twheeler/nsontakke/BLIMMP_2/paper/analysis/Prior_Generation/ATB_KO_Frequencies

python regenerate_cooccurrence_allvsall.py \
    --intersection-tsv Lineage_Specific_Data/KO_Intersections_Lineage_Specific/ko_intersection_ko_matrix_sampleids_domain_level_priors.tsv \
    --out CoOccurrence_AllVsAll_domain_level_priors.json
