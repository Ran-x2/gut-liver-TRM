#!/bin/bash
#SBATCH --job-name=rxr456_geneformer_gdT_perturb2000_debug
#SBATCH --output=geneformer_gdT_perturb2000_debug_output.txt
#SBATCH --error=geneformer_gdT_perturb2000_debug_error.txt
#SBATCH --time=15:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=128gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=65
#SBATCH --mail-user=rxr456@case.edu
#SBATCH --mail-type=ALL

# Use the Python from your Miniconda installation
export PATH="/home/rxr456/miniconda3/bin:$PATH"
source /home/rxr456/miniconda3/etc/profile.d/conda.sh
storage_dir=/mnt/vstor/SOM_PATH_DKB50/members/rxr456/gdT_geneformer_TRMTeff
conda activate geneformer
python $storage_dir/geneformer_hpc_2_state_emb.py