#!/usr/bin/env bash
# Submits one SLURM job per (popmodel, mutsig) scenario (12 total, for
# sampling=linearconstant), running compute_errors.py on a single tree per
# scenario using the FULL (un-thinned) combined chain instead of the
# subsampled ~1000-sample files.
#
# Assumes the one-row scenario tsvs already exist under
# results/run1/single_tree_test_tsvs/ (built locally with
# make_single_tree_tsvs.py and uploaded via scp) - this script does NOT
# build them, it only submits the jobs.
#
# Usage (run on Euler, from the repo root):
#   ./scripts/submit_single_tree_full_chain.sh

set -euo pipefail

TSV_DIR="results/run1/single_tree_test_tsvs"
OUT_ROOT="results/run1/single_tree_full_chain_test"
SAMPLING="linearconstant"

POP_MODELS=(expgrowthfast expgrowthslow uniform bottleneck)
MUTSIGS=(low med high)

mkdir -p logs "$OUT_ROOT"

for popmodel in "${POP_MODELS[@]}"; do
  for mutsig in "${MUTSIGS[@]}"; do
    tsv="${TSV_DIR}/${SAMPLING}_${popmodel}_${mutsig}.tsv"
    if [ ! -s "$tsv" ]; then
        echo "SKIP (missing/empty tsv): $tsv"
        continue
    fi

    out_dir="${OUT_ROOT}/${popmodel}/${mutsig}"
    mkdir -p "$out_dir"

    echo "Submitting: ${SAMPLING}/${popmodel}/${mutsig}"
    sbatch --job-name="st_${popmodel}_${mutsig}" --time=02:00:00 --mem-per-cpu=48000 --cpus-per-task=1 \
        --output="logs/single_tree_full_chain_${popmodel}_${mutsig}_%j.out" \
        --wrap="conda run -n beast_tools python -u scripts/compute_errors.py \
            --tsv ${tsv} \
            --out_dir ${out_dir} \
            --traj_points 1000 \
            --log_suffix combined \
            --trees_suffix combined"
  done
done
