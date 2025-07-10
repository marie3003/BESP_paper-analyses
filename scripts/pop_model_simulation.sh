# we are in scripts folder

# run r script to create simulation tree files

# Run python script to create simulation XML files
python MakeBEASTXML.py -i ../results/pop_size_simulations/config/


# run Beast on all XML files in the current directory
# adapt such that it runs in all subdirectories of simulation_results
# on mac I use beast command, on linux java -jar beast.jar
module load jdk/8u141-b15
module load stack/2024-06
module load gcc/12.2.0
module load beast1/1.10.4

ls *.xml | parallel --delay 1 --jobs 75% --results outdir -I% --max-args 1 beast -overwrite -seed 42 %

# run tracer on all log files in all subdirectories and delete trees that haven't converged

# combine .trees files into summary tree running TreeAnnotator

# compare summary tree with true simulated tree