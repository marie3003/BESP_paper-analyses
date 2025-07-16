# we are in scripts folder

# run r script to create simulation tree files


# Run python script to create simulation XML files
conda activate beast_tools
python MakeBEASTXML.py -i ../results/pop_size_simulations/config/


# run Beast on all XML files in the current directory
# adapt such that it runs in all subdirectories of simulation_results
# on mac I use beast command, on linux java -jar beast.jar
find /cluster/work/stadler/beckermar/BESP_paper-analyses/results/pop_size_simulations/simulation_results/ -name "*.xml" > xml_list.txt

module load jdk/8u141-b15
module load stack/2024-06
module load gcc/12.2.0
module load beast1/1.10.4

ls *.xml | parallel --delay 1 --jobs 75% --results outdir -I% --max-args 1 '"'/Volumes/BEAST\ v1.10.4\ 1/bin/beast'"' -overwrite -seed 42 % 

# run tracer on all log files in all subdirectories and delete trees that haven't converged

# combine .trees files into summary tree running TreeAnnotator

CSV_FILE="/Users/mariebecker/Documents/Uni/ETH/RotationStadler/BESP_paper-analyses/scripts/successful_mcmc_runs.csv"

# Skip header and process each line
tail -n +2 "$CSV_FILE" | while IFS=',' read -r trees_path burnin; do
  # Build output path: same directory, same base name, .tree extension
  out_dir=$(dirname "$trees_path")
  base_name=$(basename "$trees_path" .trees)
  output_path="$out_dir/${base_name}.tree"

  # Run treeannotator
  treeannotator -burnin "$burnin" -heights median "$trees_path" "$output_path"
done

# compare summary tree with true simulated tree


#SBATCH --output=logs/dummy_%A_%a.out
#SBATCH --error=logs/dummy_%A_%a.err