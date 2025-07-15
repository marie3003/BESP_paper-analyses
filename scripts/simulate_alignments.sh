#!/bin/bash
#SBATCH --job-name=seqgen_snp
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err
#SBATCH --array=0-14
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:30:00

# Setup Conda in non-interactive shell
source ~/miniconda3/etc/profile.d/conda.sh  # Update path if needed
conda activate snp_sites

# Exit on error and print commands
set -euo pipefail

# Read task info from line in tree_jobs.txt
JOB_LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" tree_jobs.txt)
FOLDER=$(echo "$JOB_LINE" | cut -d ' ' -f 1)
INDEX=$(echo "$JOB_LINE" | cut -d ' ' -f 2)

BASE="../results/pop_size_simulations/independent_homochronous/${FOLDER}"
TREE_FILE="${BASE}/${FOLDER}.trees"

# Print debug info
echo "Running SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "Folder: $FOLDER"
echo "Tree index: $INDEX"
echo "Tree file: $TREE_FILE"

# Check if tree file exists
if [ ! -f "$TREE_FILE" ]; then
  echo "Error: Tree file not found: $TREE_FILE"
  exit 1
fi

# Extract the tree
TREE=$(sed -n "$((INDEX + 1))p" "$TREE_FILE" | tr -d '\r\n')

# Check if tree was extracted
if [ -z "$TREE" ]; then
  echo "Error: No tree found at index $INDEX in $TREE_FILE"
  exit 1
fi

# Define output paths
FASTA_FILE="${BASE}/${FOLDER}_${INDEX}.fasta"
SNP_FILE="${BASE}/${FOLDER}_${INDEX}_snps.fasta"

# Run seq-gen
echo "$TREE" | ./../../Seq-Gen-1.3.5/source/seq-gen \
  -mHKY -t0.5 -f0.25,0.25,0.25,0.25 \
  -l10000000 -s4.6e-8 -n1 > "$FASTA_FILE"

# Run snp-sites
snp-sites -o "$SNP_FILE" "$FASTA_FILE"