#!/bin/bash
#SBATCH --job-name=combine_runs
#SBATCH --output=logs/combine_runs_%A_%a.out
#SBATCH --error=logs/combine_runs_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=standard
#SBATCH --array=1-1800  # <-- Adjust to match number of lines - 1 in xml_list.txt

# Load modules
module load stack/2024-06
module load openjdk/21.0.3_9
module load gcc/12.2.0
module load beast1/1.10.4
module load libbeagle

XML_LIST="xml_list.txt"

# Get the line corresponding to this array task
xml_file=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$XML_LIST")

# Replace .xml with .log and .trees
base_log="${xml_file%.xml}.log"
base_trees="${xml_file%.xml}.trees"

# Input paths
log1="$base_log"
log2="${base_log/simulation_results/simulation_results_2}"
log3="${base_log/simulation_results/simulation_results_3}"

trees1="$base_trees"
trees2="${base_trees/simulation_results/simulation_results_2}"
trees3="${base_trees/simulation_results/simulation_results_3}"

# Output paths
out_log="${base_log/simulation_results/simulation_results_combined}"
out_trees="${base_trees/simulation_results/simulation_results_combined}"

# Ensure output directories exist
mkdir -p "$(dirname "$out_log")"
mkdir -p "$(dirname "$out_trees")"

echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Combining LOG files:"
echo "  $log1"
echo "  $log2"
echo "  $log3"
echo "→ $out_log"

logcombiner -burnin 1000 "$log1" "$log2" "$log3" "$out_log"

echo "Combining TREES files:"
echo "  $trees1"
echo "  $trees2"
echo "  $trees3"
echo "→ $out_trees"

logcombiner -trees -burnin 1000 "$trees1" "$trees2" "$trees3" "$out_trees"