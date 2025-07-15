#!/bin/bash
#SBATCH --job-name=seqgen_snp
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err
#SBATCH --array=0-14
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:30:00

# conda activate snp_sites

set -e

# Read line corresponding to task ID
JOB_LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" tree_jobs.txt)
FOLDER=$(echo "$JOB_LINE" | cut -d ' ' -f 1)
INDEX=$(echo "$JOB_LINE" | cut -d ' ' -f 2)

BASE="../results/pop_size_simulations/independent_homochronous/${FOLDER}"
TREE_FILE="${BASE}/${FOLDER}.trees"
TREE=$(sed -n "$((INDEX + 1))p" "$TREE_FILE")
echo TREE_FILE: $TREE_FILE
echo TREE: $TREE

# Define output paths
FASTA_FILE="${BASE}/${FOLDER}_${INDEX}.fasta"
SNP_FILE="${BASE}/${FOLDER}_${INDEX}_snps.fasta"

# Run seq-gen
echo "$TREE" | ./../../Seq-Gen-1.3.5/source/seq-gen -mHKY -t0.5 -f0.25,0.25,0.25,0.25 -l10000000 -s4.6e-8 -n1 > "$FASTA_FILE"

# Run snp-sites
snp-sites -o "$SNP_FILE" "$FASTA_FILE"