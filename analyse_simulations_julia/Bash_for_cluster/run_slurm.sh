#!/bin/bash
set -euo pipefail

for dir in 2025_07_24/h_*; do
  if [ -d "$dir" ]; then
    slurm_file=$(find "$dir" -maxdepth 1 -type f -name 'job_*.slurm' | head -n 1)
    if [ -n "$slurm_file" ]; then
      echo "Submitting $slurm_file"
      sbatch "$slurm_file"
    else
      echo "No job_*.slurm file found in $dir"
    fi
  fi
done


