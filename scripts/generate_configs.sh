#!/bin/bash

# Generates all BEAST config files for run1
# Run from the scripts/ directory

OUTDIR="../results/run1/config"
TEMPLATE_CONSTCOAL="results/run1/templates/template_constantcoal.xml"
TEMPLATE_SKYLINE="results/run1/templates/template_skyline.xml"

mkdir -p "$OUTDIR"

SAMPLING_TYPES=("independenthomochronous" "linearconstant")
POP_MODELS=("expgrowthfast" "expgrowthslow" "uniform" "bottleneck")

get_skyline_chainlength() {
    local sampling=$1 pop=$2 mutsig=$3
    if [ "$sampling" = "linearconstant" ]; then
        case "${pop}_${mutsig}" in
            expgrowthfast_lowmutsig)  echo 80000000  ;;
            expgrowthfast_medmutsig)  echo 80000000  ;;
            expgrowthfast_highmutsig) echo 80000000  ;;
            expgrowthslow_lowmutsig)  echo 120000000 ;;
            expgrowthslow_medmutsig)  echo 80000000  ;;
            expgrowthslow_highmutsig) echo 80000000  ;;
            uniform_lowmutsig)        echo 80000000  ;;
            uniform_medmutsig)        echo 80000000  ;;
            uniform_highmutsig)       echo 80000000  ;;
            bottleneck_lowmutsig)     echo 120000000 ;;
            bottleneck_medmutsig)     echo 120000000 ;;
            bottleneck_highmutsig)    echo 200000000 ;;
        esac
    else
        # independenthomochronous: per popmodel+mutsig
        case "${pop}_${mutsig}" in
            expgrowthfast_lowmutsig)  echo 200000000 ;;
            expgrowthfast_medmutsig)  echo 80000000  ;;
            expgrowthfast_highmutsig) echo 40000000  ;;
            expgrowthslow_lowmutsig)  echo 250000000 ;;
            expgrowthslow_medmutsig)  echo 80000000  ;;
            expgrowthslow_highmutsig) echo 40000000  ;;
            uniform_lowmutsig)        echo 200000000 ;;
            uniform_medmutsig)        echo 150000000 ;;
            uniform_highmutsig)       echo 40000000  ;;
            bottleneck_lowmutsig)     echo 200000000 ;;
            bottleneck_medmutsig)     echo 200000000 ;;
            bottleneck_highmutsig)    echo 200000000 ;;
        esac
    fi
}

get_skyline_logfreq() {
    local chain=$(get_skyline_chainlength $1 $2 $3)
    echo $((chain / 10000))
}

get_snpfraction() {
    case $1 in
        lowmutsig)  echo 0.1 ;;
        medmutsig)  echo 0.5 ;;
        highmutsig) echo 1.0 ;;
    esac
}

for sampling in "${SAMPLING_TYPES[@]}"; do
    for pop in "${POP_MODELS[@]}"; do
        for mutsig in lowmutsig medmutsig highmutsig; do

            SNPFRACTION=$(get_snpfraction $mutsig)
            TREES="results/run1/simulated_data/${sampling}/${pop}/${pop}.trees"
            if [ "$sampling" = "independenthomochronous" ]; then
                HOMOCHRONOUS="True"
            else
                HOMOCHRONOUS="False"
            fi

            # --- constcoal ---
            CONSTCOAL_CHAIN=30000000
            CONSTCOAL_LOG=3000
            OUTPATH="results/run1/beast_inference/constcoal/${sampling}/${pop}/${mutsig}/"
            NAME="constcoal_${sampling}_${pop}_${mutsig}"
            FILE="${OUTDIR}/constcoal_${sampling}_${pop}_${mutsig}.cfg"

            cat > "$FILE" << EOF
# Input - paths relative to MakeBEASTXML.py
template : ${TEMPLATE_CONSTCOAL}
trees    : ${TREES}

# Scenario
homochronous: ${HOMOCHRONOUS}
snpFraction: ${SNPFRACTION}

# Parameters
clockRate      : 4.6e-8
initialPopSize : 1000

# Analysis
chainlength        : ${CONSTCOAL_CHAIN}
filelog            : ${CONSTCOAL_LOG}
screenlog          : ${CONSTCOAL_LOG}
treelog            : ${CONSTCOAL_LOG}

# Output
name       : ${NAME}
outputpath : ${OUTPATH}
EOF

            # --- skyline ---
            SKYLINE_CHAIN=$(get_skyline_chainlength $sampling $pop $mutsig)
            SKYLINE_LOG=$(get_skyline_logfreq $sampling $pop $mutsig)
            OUTPATH="results/run1/beast_inference/skyline/${sampling}/${pop}/${mutsig}/"
            NAME="skyline_${sampling}_${pop}_${mutsig}"
            FILE="${OUTDIR}/skyline_${sampling}_${pop}_${mutsig}.cfg"

            cat > "$FILE" << EOF
# Input - paths relative to MakeBEASTXML.py
template : ${TEMPLATE_SKYLINE}
trees    : ${TREES}

# Scenario
homochronous: ${HOMOCHRONOUS}
snpFraction: ${SNPFRACTION}

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
