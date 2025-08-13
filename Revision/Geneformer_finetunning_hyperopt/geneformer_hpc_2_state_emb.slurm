#!/bin/bash
#SBATCH --job-name=rxr456_geneformer_TRM_revision_perturb5000
#SBATCH --output=geneformer_TRM_revision_perturb5000_output.txt
#SBATCH --error=geneformer_TRM_revision_perturb5000_error.txt
#SBATCH --time=3-00:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=128gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=65
#SBATCH --mail-user=rxr456@case.edu
#SBATCH --mail-type=ALL

set -euo pipefail

# Conda
export PATH="/home/rxr456/miniconda3/bin:$PATH"
source /home/rxr456/miniconda3/etc/profile.d/conda.sh
conda activate geneformer

storage_dir=/mnt/vstor/SOM_PATH_DKB50/members/rxr456/Trapecar_geneformer_finetune
repo_dir=/home/rxr456/Github/gut-liver-TRM/Revision/Geneformer_finetunning_hyperopt
repo_script="$repo_dir/geneformer_hpc_2_state_emb.py"   # single script in repo

# Run serially: same script, different working dir each time
for subdir in "$storage_dir"/*/; do
  # optional: only run if subdir contains a specific file
  # [[ -f "$subdir/config.yaml" ]] || { echo "Skip $subdir (no config)"; continue; }

  echo "Running $repo_script in $subdir ..."
  ( cd "$subdir" && python -u "$repo_script" )
done