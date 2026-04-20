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
