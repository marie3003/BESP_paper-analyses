#!/bin/bash
#SBATCH --job-name=test_compute_errors
#SBATCH --output=logs/compute_errors/test_compute_errors.log
#SBATCH --error=logs/compute_errors/test_compute_errors.log
#SBATCH --mem-per-cpu=8000
#SBATCH --cpus-per-task=1
#SBATCH --time=240
#SBATCH --partition=normal

mkdir -p logs/compute_errors

conda run -n beast_tools python scripts/compute_errors.py \
    --tsv results/run1/evaluation/linearconstant/uniform/med.tsv \
    --out_dir results/run1/evaluation/linearconstant/uniform/ \
    --max_reps 3
