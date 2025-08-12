#!/bin/bash
#SBATCH --job-name=rxr456_geneformer_TRM_revision_hyperopt_train
#SBATCH --output=geneformer_TRM_revision_hyperopt_train_output.txt
#SBATCH --error=geneformer_TRM_revision_hyperopt_train_error.txt
#SBATCH --time=24:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=40
#SBATCH --mem=128G
#SBATCH --mail-user=rxr456@case.edu
#SBATCH --mail-type=ALL

set -euo pipefail

# Conda
export PATH="/home/rxr456/miniconda3/bin:$PATH"
source /home/rxr456/miniconda3/etc/profile.d/conda.sh
conda activate geneformer

storage_dir=/mnt/vstor/SOM_PATH_DKB50/members/rxr456/Trapecar_geneformer_finetune

# run each subfolder’s first .py, serially, from its own dir
for subdir in "$storage_dir"/*/; do
  py_file=$(find "$subdir" -maxdepth 1 -name "*hyperopt.py" | head -n 1)
  if [[ -n "${py_file:-}" ]]; then
    echo "Running $py_file ..."
    ( cd "$subdir" && python -u "$(basename "$py_file")" )
  else
    echo "No hyperopt.py file found in $subdir"
  fi
done
