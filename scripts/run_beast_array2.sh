#!/bin/bash
#SBATCH --job-name=beast_array
#SBATCH --output=logs/beast_%A_%a.out
#SBATCH --error=logs/beast_%A_%a.err
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-1800
#SBATCH --partition=standard

module load stack/2024-06
module load openjdk/21.0.3_9
module load gcc/12.2.0
module load beast1/1.10.4
module load libbeagle

# Path to the list of all XML files
FILE_LIST=/cluster/work/stadler/beckermar/BESP_paper-analyses/scripts/xml_list.txt

# Get the original XML file (still in simulation_results)
XML_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $FILE_LIST)

# Build relative path under simulation_results
REL_PATH=$(dirname "$(realpath --relative-to=/cluster/work/stadler/beckermar/BESP_paper-analyses/pop_size_simulations/simulation_results "$XML_FILE")")

# Determine new output folder and prefix path
OUT_XML_BASE="/cluster/work/stadler/beckermar/BESP_paper-analyses/pop_size_simulations/simulation_results_2/${REL_PATH%.xml}"

# Create output directory
mkdir -p "$(dirname "$OUT_XML_BASE")"

# Run BEAST using XML_FILE as input, but redirect output files to new location via -prefix
beast -seed 43 -prefix "$OUT_XML_BASE" "$XML_FILE"

