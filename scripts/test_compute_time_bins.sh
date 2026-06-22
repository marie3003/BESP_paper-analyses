#!/bin/bash
#SBATCH --job-name=test_time_bins
#SBATCH --output=logs/compute_time_bins/test_compute_time_bins.log
#SBATCH --error=logs/compute_time_bins/test_compute_time_bins.log
#SBATCH --mem-per-cpu=8000
#SBATCH --cpus-per-task=1
#SBATCH --time=240
#SBATCH --partition=normal

mkdir -p logs/compute_time_bins

conda run -n beast_tools python -u scripts/compute_time_bins.py \
    --tsv results/run1/evaluation/linearconstant/uniform/med.tsv \
    --out_dir results/run1/evaluation/linearconstant/uniform/test/ \
    --config "10,400" \
    --max_reps 10
