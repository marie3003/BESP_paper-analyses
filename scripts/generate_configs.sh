#!/bin/bash

# Generates all BEAST config files for run1
# Run from the scripts/ directory

OUTDIR="../results/run1/config"
TEMPLATE_CONSTCOAL="../results/pop_size_simulations/templates/template_constantcoal.xml"
TEMPLATE_SKYLINE="../results/pop_size_simulations/templates/template_skyline.xml"

mkdir -p "$OUTDIR"

SAMPLING_TYPES=("independent_homochronous" "linear_constant")
POP_MODELS=("expgrowth_fast" "expgrowth_slow" "uniform" "bottleneck")

get_skyline_chainlength() {
    case $1 in
        lowmutsig)  echo 120000000 ;;
        medmutsig)  echo 80000000  ;;
        highmutsig) echo 80000000  ;;
    esac
}

get_skyline_logfreq() {
    case $1 in
        lowmutsig)  echo 12000 ;;
        medmutsig)  echo 8000  ;;
        highmutsig) echo 8000  ;;
    esac
}

get_nsites() {
    case $1 in
        lowmutsig)  echo 1000  ;;
        medmutsig)  echo 5000  ;;
        highmutsig) echo 10000 ;;
    esac
}

for sampling in "${SAMPLING_TYPES[@]}"; do
    for pop in "${POP_MODELS[@]}"; do
        for mutsig in lowmutsig medmutsig highmutsig; do

            NSITES=$(get_nsites $mutsig)
            TREES="../results/run1/simulated_data/${sampling}/${pop}/${pop}.trees"
            if [ "$sampling" = "independent_homochronous" ]; then
                HOMOCHRONOUS="True"
            else
                HOMOCHRONOUS="False"
            fi

            # --- constcoal ---
            sampling_short="${sampling//_/}"
            pop_short="${pop//_/}"
            OUTPATH="../results/run1/beast_inference/constcoal/${sampling}/${pop}/${mutsig}/"
            NAME="constcoal_${sampling_short}_${pop_short}_${mutsig}"
            FILE="${OUTDIR}/constcoal_${sampling_short}_${pop_short}_${mutsig}.cfg"

            cat > "$FILE" << EOF
# Input - paths relative to MakeBEASTXML.py
template : ${TEMPLATE_CONSTCOAL}
trees    : ${TREES}

# Scenario
homochronous: ${HOMOCHRONOUS}
nSites: ${NSITES}

# Parameters
clockRate      : 4.6e-8
initialPopSize : 1000

# Analysis
chainlength        : 10000000
filelog            : 1000
screenlog          : 1000
treelog            : 1000

# Output
name       : ${NAME}
outputpath : ${OUTPATH}
EOF

            # --- skyline ---
            SKYLINE_CHAIN=$(get_skyline_chainlength $mutsig)
            SKYLINE_LOG=$(get_skyline_logfreq $mutsig)
            OUTPATH="../results/run1/beast_inference/skyline/${sampling}/${pop}/${mutsig}/"
            NAME="skyline_${sampling_short}_${pop_short}_${mutsig}"
            FILE="${OUTDIR}/skyline_${sampling_short}_${pop_short}_${mutsig}.cfg"

            cat > "$FILE" << EOF
# Input - paths relative to MakeBEASTXML.py
template : ${TEMPLATE_SKYLINE}
trees    : ${TREES}

# Scenario
homochronous: ${HOMOCHRONOUS}
nSites: ${NSITES}

# Parameters
numPopulationBins  : 10
initialPopSize     : 1000
clockRate          : 4.6e-8

# Analysis
chainlength        : ${SKYLINE_CHAIN}
filelog            : ${SKYLINE_LOG}
screenlog          : ${SKYLINE_LOG}
treelog            : ${SKYLINE_LOG}

# Output
name       : ${NAME}
outputpath : ${OUTPATH}
EOF

        done
    done
done

echo "Generated $(ls $OUTDIR/*.cfg | wc -l) config files in $OUTDIR"
