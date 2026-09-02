#!/bin/bash
#SBATCH --job-name=check_subsampled_ess
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/check_subsampled_ess.out

module load stack/2024-06 r/4.4.0

BEAST_DIR="results/run1/beast_inference"
OUTPUT_TSV="${BEAST_DIR}/subsampled_ess_check.tsv"

Rscript scripts/check_subsampled_ess.R "${BEAST_DIR}" "${OUTPUT_TSV}"
