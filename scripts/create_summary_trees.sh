#!/bin/bash
#SBATCH --job-name=create_summary_trees
#SBATCH --output=create_summary_trees.out
#SBATCH --error=create_summary_trees.err
#SBATCH --time=02:00:00          
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=standard 

module load jdk/8u141-b15
module load stack/2024-06
module load gcc/12.2.0
module load beast1/1.10.4
module load r/4.4.0

Rscript evaluate_mcmc.R

CSV_FILE="./successful_mcmc_runs.csv"

# Skip header and process each line
tail -n +2 "$CSV_FILE" | while IFS=',' read -r trees_path burnin; do
  # Build output path: same directory, same base name, .tree extension
  out_dir=$(dirname "$trees_path")
  base_name=$(basename "$trees_path" .trees)
  output_path="$out_dir/${base_name}.tree"

  # Run treeannotator
  treeannotator -burnin "$burnin" -heights median "$trees_path" "$output_path"
done