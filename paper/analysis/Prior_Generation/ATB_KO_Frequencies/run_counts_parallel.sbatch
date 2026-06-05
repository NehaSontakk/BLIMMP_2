#!/bin/bash
#SBATCH --job-name=ko_counts
#SBATCH --output=logs/ko_counts_%j.out
#SBATCH --error=logs/ko_counts_%j.err
#SBATCH --partition=standard
#SBATCH --account=cgoubert
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G

module load parallel
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1

source /home/u13/nsontakke/miniconda3/etc/profile.d/conda.sh
conda activate test

BASE=/xdisk/cgoubert/nsontakke/ATB_KO_Frequencies
CMD_FILE=${BASE}/ko_count_lineage_commands.txt

parallel -a "$CMD_FILE" -j ${SLURM_CPUS_PER_TASK}
