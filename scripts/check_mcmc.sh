#!/bin/bash
#SBATCH --job-name=check_mcmc
#SBATCH --output=logs/check_mcmc.out
#SBATCH --error=logs/check_mcmc.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --partition=standard

module load stack/2024-06
module load r/4.4.0

# Step 1: Run the evaluation script to create CSV
Rscript evaluate_mcmc.R

CSV_FILE="./successful_mcmc_runs.csv"

# Step 2: Count number of successful runs (skip header)
N=$(tail -n +2 "$CSV_FILE" | wc -l)

echo "Detected $N successful MCMC runs."

# Step 3: Generate and submit the summary array script dynamically (overwriting any existing one)
cat << EOF > run_summary_array.sh
#!/bin/bash
#SBATCH --job-name=create_summary_trees
#SBATCH --output=logs/summary_%A_%a.out
#SBATCH --error=logs/summary_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=standard
#SBATCH --array=1-${N}

module load jdk/8u141-b15
module load stack/2024-06
module load gcc/12.2.0
module load beast1/1.10.4

CSV_FILE="./successful_mcmc_runs.csv"

# Get the line corresponding to the current array index
LINE=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" <(tail -n +2 "\$CSV_FILE"))
trees_path=\$(echo "\$LINE" | cut -d ',' -f 1)
burnin=\$(echo "\$LINE" | cut -d ',' -f 2)

# Build output path
out_dir=\$(dirname "\$trees_path")
base_name=\$(basename "\$trees_path" .trees)
output_path="\$out_dir/\${base_name}.tree"

echo "Annotating \$trees_path with burnin \$burnin"
treeannotator -burnin "\$burnin" -heights median "\$trees_path" "\$output_path"
EOF

chmod +x run_summary_array.sh
sbatch run_summary_array.sh
