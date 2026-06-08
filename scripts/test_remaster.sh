#!/bin/bash
#SBATCH --job-name=test_remaster
#SBATCH --output=logs/test_remaster_%a.log
#SBATCH --error=logs/test_remaster_%a.log
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-7

module load stack/2024-06
module load gcc/12.2.0
module load openjdk/21.0.3_9
module load beast2/2.7.4

SAMPLING_TYPES=("independenthomochronous" "independenthomochronous" "independenthomochronous" "independenthomochronous"
                "linearconstant" "linearconstant" "linearconstant" "linearconstant")
POP_MODELS=("expgrowthfast" "expgrowthslow" "uniform" "bottleneck"
            "expgrowthfast" "expgrowthslow" "uniform" "bottleneck")

SAMPLING="${SAMPLING_TYPES[$SLURM_ARRAY_TASK_ID]}"
POPMODEL="${POP_MODELS[$SLURM_ARRAY_TASK_ID]}"
NREPLICATES=10
FILEBASE="results/run1/simulated_data"
SEED=42

TEMPLATE="results/run1/templates/remaster_${SAMPLING}_${POPMODEL}.xml"
TMPXML=$(mktemp --suffix=.xml)

sed "s|{\$nreplicates}|${NREPLICATES}|g; s|{\$filebase}|${FILEBASE}|g" "$TEMPLATE" > "$TMPXML"

mkdir -p "${FILEBASE}/${SAMPLING}/${POPMODEL}"

beast -overwrite -seed "$SEED" "$TMPXML"

rm "$TMPXML"
