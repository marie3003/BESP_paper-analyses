#!/bin/bash
#SBATCH --job-name=beast_array
#SBATCH --output=logs/beast_%A_%a.out
#SBATCH --error=logs/beast_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-90
#SBATCH --partition=standard

module load jdk/8u141-b15
module load stack/2024-06
module load gcc/12.2.0
module load beast1/1.10.4

# Path to the list of all XML files
FILE_LIST=/cluster/work/stadler/beckermar/BESP_paper-analyses/scripts/xml_list.txt

# Get the XML file corresponding to the current task ID
XML_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $FILE_LIST)

# Run BEAST on that file
beast -overwrite -seed 42 -java -working "$XML_FILE"