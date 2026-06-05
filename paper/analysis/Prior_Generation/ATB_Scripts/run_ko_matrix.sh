#!/bin/bash
#SBATCH --job-name=ko_mat
#SBATCH --output=ko_mat.out
#SBATCH --error=ko_mat.err
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --account=twheeler
#SBATCH --partition=standard


source /home/u13/nsontakke/miniconda3/etc/profile.d/conda.sh
conda activate test
python create_ko_matrix.py

