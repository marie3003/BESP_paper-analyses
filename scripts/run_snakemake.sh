#!/bin/bash
#SBATCH --job-name=snakemake_main
#SBATCH --output=snakemake.log
#SBATCH --error=snakemake.log
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=2G

conda activate snakemake

snakemake all --keep-incomplete --scheduler greedy --profile profiles/euler
