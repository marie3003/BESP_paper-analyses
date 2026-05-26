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

    job_script="${outdir}/${name}.job.sh"
    cat > "$job_script" << EOF
#!/bin/bash
#SBATCH --job-name=skygrid_${name}
#SBATCH --output=${outdir}/${name}.slurm.log
#SBATCH --error=${outdir}/${name}.slurm.log
#SBATCH --time=16:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=2
#SBATCH --partition=normal

module load stack/2024-06
module load openjdk/21.0.3_9
module load gcc/12.2.0
module load beast1/1.10.4
module load libbeagle/3.1.2

set -euo pipefail

cd ${outdir}
beast -overwrite -seed 45 ${xml_abs}
EOF

    sbatch "$job_script"
    echo "Submitted: $name"
done < "$XML_LIST"
