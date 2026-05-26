#!/bin/bash
# Submit one BEAST job per XML file listed in the input txt file.
# Usage: bash scripts/run_skygrid_beast.sh <xml_list.txt>
#
# Each line in xml_list.txt should be an absolute or relative path to a BEAST XML file.
# Output logs and trees are written next to the XML file.

XML_LIST=${1:?"Usage: $0 <xml_list.txt>"}

if [ ! -f "$XML_LIST" ]; then
    echo "Error: file not found: $XML_LIST"
    exit 1
fi

while IFS= read -r xml || [ -n "$xml" ]; do
    # Skip empty lines and comments
    [[ -z "$xml" || "$xml" == \#* ]] && continue

    if [ ! -f "$xml" ]; then
        echo "Warning: XML not found, skipping: $xml"
        continue
    fi

    xml_abs=$(realpath "$xml")
    outdir=$(dirname "$xml_abs")
    name=$(basename "$xml_abs" .xml)
    logfile="${outdir}/${name}.log"
    treefile="${outdir}/${name}.trees"

    sbatch \
        --job-name="skygrid_${name}" \
        --output="${outdir}/${name}.slurm.log" \
        --error="${outdir}/${name}.slurm.log" \
        --time=16:00:00 \
        --mem-per-cpu=1G \
        --cpus-per-task=2 \
        --wrap="
            module load stack/2024-06 openjdk/21.0.3_9 gcc/12.2.0 beast1/1.10.4 libbeagle/3.1.2
            cd ${outdir}
            beast -overwrite -seed 45 ${xml_abs} > ${logfile} 2>&1
        "

    echo "Submitted: $name"
done < "$XML_LIST"
