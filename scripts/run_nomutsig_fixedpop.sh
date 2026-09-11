#!/bin/bash
#SBATCH --job-name=beast_nomutsig_fixedpop
#SBATCH --output=beast_nomutsig_fixedpopsize_%j.out
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1000

module load stack/2024-06
module load openjdk/21.0.3_9
module load gcc/12.2.0
module load beast1/1.10.4
module load libbeagle/3.1.2

beast -overwrite -seed 44 constcoal_linearconstant_uniform_nomutsig_fixedpopsize.T86.xml
