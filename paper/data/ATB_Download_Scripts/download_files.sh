#!/bin/bash
#SBATCH --job-name=ATB_dl
#SBATCH --output=logs/download_%A_%a.out
#SBATCH --error=logs/download_%A_%a.err
#SBATCH --array=1-56              # 27599 ÷ 500 ≃ 55.2 → 56 tasks
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --partition=standard
#SBATCH --account=twheeler

mkdir -p ATB_Downloads logs

start=$(( (SLURM_ARRAY_TASK_ID - 1) * 500 + 1 ))
end=$(( SLURM_ARRAY_TASK_ID * 500 ))

sed -n "${start},${end}p" download_commands.txt | while read -r cmd; do
  # extract the output filename from the wget command
  # (assumes form: wget -O ATB_Downloads/<tar_xz> <url>)
  outfile=$(echo "$cmd" | awk '{print $3}')
  if [[ -f "$outfile" ]]; then
    echo "Skipping existing $outfile"
  else
    echo "Downloading $outfile"
    bash -c "$cmd"
  fi
done

