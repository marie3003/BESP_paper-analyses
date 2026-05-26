#!/bin/bash

# Generates BEAST config files for run_bottleneck
# Run from the scripts/ directory

OUTDIR="../results/run_bottleneck/config"
TEMPLATE_CONSTCOAL="results/run_bottleneck/templates/template_constantcoal.xml"
TEMPLATE_SKYLINE="results/run_bottleneck/templates/template_skyline.xml"

mkdir -p "$OUTDIR"

SAMPLING_TYPES=("independenthomochronous" "linearconstant")
POP_MODELS=("bottleneck20" "bottleneck50" "bottlenecklate" "bottlenecklatesampling")

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
            TREES="results/run_bottleneck/simulated_data/${sampling}/${pop}/${pop}.trees"
            if [ "$sampling" = "independenthomochronous" ]; then
                HOMOCHRONOUS="True"
            else
                HOMOCHRONOUS="False"
            fi

            # --- constcoal ---
            CONSTCOAL_CHAIN=30000000
            CONSTCOAL_LOG=3000
            OUTPATH="results/run_bottleneck/beast_inference/constcoal/${sampling}/${pop}/${mutsig}/"
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

            # --- skyline (using bottleneck chain lengths from run1) ---
            if [ "$sampling" = "independenthomochronous" ]; then
                SKYLINE_CHAIN=300000000
            else
                SKYLINE_CHAIN=200000000
            fi
            SKYLINE_LOG=$((SKYLINE_CHAIN / 10000))
            OUTPATH="results/run_bottleneck/beast_inference/skyline/${sampling}/${pop}/${mutsig}/"
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
