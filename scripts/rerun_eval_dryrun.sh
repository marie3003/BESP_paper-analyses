#!/bin/bash
#SBATCH --job-name=snakemake_dryrun_eval
#SBATCH --output=logs/snakemake_dryrun_eval_%j.log
#SBATCH --error=logs/snakemake_dryrun_eval_%j.log
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4G

# Dry run: shows which jobs would run when re-running only
# step 10 (compute_errors) and step 11 (compute_time_bins).
# Any other rule that would need to run makes this fail instead.
#
# Submit from the project root:
#   sbatch scripts/rerun_eval_dryrun.sh

conda activate snakemake

snakemake all --profile profiles/euler -n \
    --forcerun compute_errors compute_time_bins \
    --allowed-rules compute_errors compute_time_bins \
    --rerun-triggers mtime
