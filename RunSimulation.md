# Run Simulation
This repository contains a simulation study to estimate the impact that different population size models have on the estimation of the tMRCA in trees.


## Workflow
Here you can find an explanation on how to repoduce the simulation study.

### Simulate trees under a population model
The population models to choose from are specified in the `SimUtils.R` file. Specific models for simulation are selected in the `simulate_trees.R` file which can be run to simulate the trees.

### Make Beast xml files
For all simulated trees, sequences are generated based on a specified mutation rate and used as an input to Beast. Beast is then run with a constant coalescent population size prior and a Skyline Coalescent prior. To generate the xml files needed as input to Beast, `scripts/MakeBEASTXML.py` can be one. It uses templates specified in the `results/pop_size_simulations/templates` folder and adds parameters specified in the config files that can be found in the `results/pop_size_simulations/config` folder. `MakeBEASTXML.py` generates sequences using SeqGen. The path to the SeqGen installation needs to be specified in the script.

### Run Beast
To run beast, first a list of all beast xml files is needed. This list is generated with `find ../results/pop_size_simulations/ -name "*.xml" > xml_list.txt` run in the `scripts` folder.

Afterwards Beast is run on all xml files by running `run_beast_array.sh`.

### Continue with trees with sufficient ESS
The R-script `evaluate_mcmc.R` can be run which outputs a csv file specifying the paths to all successfully obtained `.trees` files with an ESS on all estimated values above $200$. 

### Create summary trees using treeannotator
Run the `create_summary_trees.sh` script to create summary trees using treeannotator.

