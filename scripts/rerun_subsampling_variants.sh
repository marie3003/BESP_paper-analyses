#!/usr/bin/env bash
# Redo subsampling 4 more times with different (additional) burnin values,
# each producing an independent evenly-spaced ~1000-sample subset, then
# rerun compute_errors.py on each variant into evaluation2..evaluation5.
#
# Does NOT touch results/run1/beast_inference/*.subsampled.{log,trees} or
# results/run1/evaluation/* (the original Snakemake outputs) - variant
# subsampled files get a distinct suffix in the same directory as the
# combined files, and error outputs go to sibling "evaluation{v}" dirs.
#
# Usage:
#   ./scripts/rerun_subsampling_variants.sh subsample   # step 1: logcombiner variants
#   ./scripts/rerun_subsampling_variants.sh errors       # step 2: compute_errors.py per variant
#   ./scripts/rerun_subsampling_variants.sh all          # both

set -euo pipefail

module load stack/2024-06 openjdk/21.0.3_9 gcc/12.2.0 beast1/1.10.4 libbeagle/3.1.2 2>/dev/null || true

BEAST_DIR="results/run1/beast_inference"
EVAL_DIR="results/run1/evaluation"
NREPLICATES=100

SAMPLING_TYPES=(independenthomochronous linearconstant)
POP_MODELS=(expgrowthfast expgrowthslow uniform bottleneck)
MUTSIGS=(lowmutsig medmutsig highmutsig)
declare -A MUTSIG_SHORT=( [lowmutsig]=low [medmutsig]=med [highmutsig]=high )

# Phase-shift the resampling grid: skip a handful of samples (not a big
# burnin chunk) before applying the same evenly-spaced resample interval,
# so each variant draws a different ~1000-sample subset of the SAME
# underlying ~27000 samples rather than discarding a meaningful chunk.
# Values are in units of SAMPLES (not states) - e.g. 5 means "start at the
# 5th sample" instead of the 0th.
VARIANTS=(2 3 4 5)
declare -A OFFSET_SAMPLES=( [2]=5 [3]=10 [4]=15 [5]=20 )

subsample_one() {
    local inlog=$1 intrees=$2 outlog=$3 outtrees=$4 offset=$5 logfile=$6

    STATE1=$(grep -Ev "^(#|Sample|State)" "$inlog" | awk 'NR==1{print $1}')
    STATE2=$(grep -Ev "^(#|Sample|State)" "$inlog" | awk 'NR==2{print $1}')
    N_SAMPLES=$(grep -Evc "^(#|Sample|State)" "$inlog" || true)
    LOG_FREQ=$(( STATE2 - STATE1 ))

    # burnin = just enough states to skip `offset` samples from the start
    BURNIN=$(( offset * LOG_FREQ ))

    REMAINING=$(( N_SAMPLES - offset ))
    if [ "$REMAINING" -lt 2 ]; then
        echo "SKIP (offset too large): $inlog offset=$offset" >> "$logfile"
        return
    fi

    RESAMPLE=$(( ((REMAINING * LOG_FREQ) / 1000 / LOG_FREQ) * LOG_FREQ ))
    if [ "$RESAMPLE" -lt "$LOG_FREQ" ]; then RESAMPLE=$LOG_FREQ; fi

    echo "$(basename "$inlog"): offset=$offset samples BURNIN=$BURNIN N_SAMPLES=$N_SAMPLES LOG_FREQ=$LOG_FREQ RESAMPLE=$RESAMPLE" >> "$logfile"
    logcombiner -burnin "$BURNIN" -resample "$RESAMPLE" "$inlog" "$outlog" >> "$logfile" 2>&1
    logcombiner -trees -burnin "$BURNIN" -resample "$RESAMPLE" "$intrees" "$outtrees" >> "$logfile" 2>&1
}

run_subsample() {
    mkdir -p logs/rerun_subsampling_variants
    for sampling in "${SAMPLING_TYPES[@]}"; do
      for popmodel in "${POP_MODELS[@]}"; do
        for mutsig in "${MUTSIGS[@]}"; do
          for i in $(seq 0 $((NREPLICATES-1))); do
            for model in constcoal skyline; do
              base="${BEAST_DIR}/${model}/${sampling}/${popmodel}/${mutsig}/${model}_${sampling}_${popmodel}_${mutsig}.T${i}"
              combined_log="${base}.combined.log"
              combined_trees="${base}.combined.trees"

              # Skip replicates that failed ESS originally (empty combined/subsampled outputs)
              [ -s "$combined_log" ] || continue
              [ -s "${base}.subsampled.log" ] || continue

              for v in "${VARIANTS[@]}"; do
                outlog="${base}.subsampled_v${v}.log"
                outtrees="${base}.subsampled_v${v}.trees"
                logfile="logs/rerun_subsampling_variants/v${v}_${sampling}_${popmodel}_${mutsig}_T${i}.log"
                subsample_one "$combined_log" "$combined_trees" "$outlog" "$outtrees" "${OFFSET_SAMPLES[$v]}" "$logfile"
              done
            done
          done
        done
      done
    done
}

run_errors() {
    for v in "${VARIANTS[@]}"; do
      for sampling in "${SAMPLING_TYPES[@]}"; do
        for popmodel in "${POP_MODELS[@]}"; do
          for mutsig in "${MUTSIGS[@]}"; do
            short=${MUTSIG_SHORT[$mutsig]}
            tsv="${EVAL_DIR}/${sampling}/${popmodel}/${short}.tsv"
            out_dir="results/run1/evaluation${v}/${sampling}/${popmodel}"
            mkdir -p "$out_dir" "logs/compute_errors_v${v}"
            echo "compute_errors variant $v: ${sampling}/${popmodel}/${mutsig}"
            conda run -n beast_tools python -u scripts/compute_errors.py \
                --tsv "$tsv" \
                --out_dir "$out_dir" \
                --traj_points 1000 \
                --log_suffix "subsampled_v${v}" \
                --trees_suffix "subsampled_v${v}" \
                > "logs/compute_errors_v${v}/${sampling}_${popmodel}_${mutsig}.log" 2>&1
          done
        done
      done
    done
}

case "${1:-all}" in
    subsample) run_subsample ;;
    errors)    run_errors ;;
    all)       run_subsample; run_errors ;;
    *) echo "Usage: $0 {subsample|errors|all}"; exit 1 ;;
esac
