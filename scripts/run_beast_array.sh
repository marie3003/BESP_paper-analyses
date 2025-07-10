#!/bin/bash
#SBATCH --job-name=beast_array
#SBATCH --output=logs/dummy_%A_%a.out
#SBATCH --error=logs/dummy_%A_%a.err
#SBATCH --time=01:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-90
#SBATCH --partition=standard

# Define list of XML files
FILE_LIST=/cluster/work/stadler/beckermar/BESP_paper-analyses/scripts/xml_list.txt

# Get the XML file for this array task
XML_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $FILE_LIST)

# Extract just the filename (no path, no extension)
XML_NAME=$(basename "$XML_FILE" .xml)

# Load required modules
module load jdk/8u141-b15
module load stack/2024-06
module load gcc/12.2.0
module load beast1/1.10.4

# Define output paths
OUT_LOG=logs/${XML_NAME}.out
ERR_LOG=logs/${XML_NAME}.err

# Run BEAST with manual output redirection
beast -overwrite -seed 42 "$XML_FILE" > "$OUT_LOG" 2> "$ERR_LOG"