#!/bin/bash
#SBATCH --job-name=snakemake_rerun_eval
#SBATCH --output=logs/snakemake_rerun_eval_%j.log
#SBATCH --error=logs/snakemake_rerun_eval_%j.log
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=2G

# Re-runs only step 10 (compute_errors) and step 11 (compute_time_bins).
# Run scripts/rerun_eval_dryrun.sh first and check its job list.
#
# Submit from the project root:
#   sbatch scripts/rerun_eval.sh

conda activate snakemake

snakemake all --profile profiles/euler \
    --forcerun compute_errors compute_time_bins \
    --allowed-rules compute_errors compute_time_bins \
    --rerun-triggers mtime \
    --keep-incomplete --keep-going --scheduler greedy
