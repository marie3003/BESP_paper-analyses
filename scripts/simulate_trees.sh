#!/bin/bash
#SBATCH --job-name=simulate_trees
#SBATCH --output=simulate_trees.out
#SBATCH --error=simulate_trees.err
#SBATCH --time=01:00:00          # 1 hour
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=standard     # or your cluster's partition

# Load R module (if needed on your cluster)
module load stack/2024-06
module load gcc/12.2.0
module load r/4.4.0

# Run your R script
Rscript /cluster/work/stadler/beckermar/BESP_paper-analyses/scripts/simulate_trees.R
