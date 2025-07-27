# Run Simulation
This repository contains a simulation study to estimate the impact that different population size models have on the estimation of the tMRCA in trees.


## Workflow
Here you can find an explanation on how to repoduce the simulation study.

### Simulate trees under a population model
The population models to choose from are specified in the `SimUtils.R` file. Specific models for simulation are selected in the `simulate_trees.R` file which can be run to simulate the trees. The simulation can be run using the `simulate_trees.sh` script (around 3h).

### Simulate alignment with SeqGen and get SNPs with `snp-sites`

Before the alignment can be generated, a list of all trees is needed. This is done by running `bash make_tree_jobs.sh`. Them, the snp-sites environment needs to be activated: `conda activate snp_sites` and the `simulate_alignments.sh` script from the `scripts` folder needs to be run (around 5min).

### Make Beast xml files
For all simulated trees, sequences are generated based on a specified mutation rate and used as an input to Beast. Beast is then run with a constant coalescent population size prior and a Skyline Coalescent prior. To generate the xml files needed as input to Beast, `scripts/MakeBEASTXML.py` can be one. It can be run running `scripts/run_make_beast_xml.sh`. It uses templates specified in the `results/pop_size_simulations/templates` folder and adds parameters specified in the config files that can be found in the `results/pop_size_simulations/config` folder. `MakeBEASTXML.py` subsamples sequences created before with SeqGen and snp sites.

### Run Beast
To run beast, first a list of all beast xml files is needed. This list is generated with `find ../results/pop_size_simulations/simulation_results_3 -name "*.xml" > xml_list_3.txt` run in the `scripts` folder.

Afterwards Beast is run on all xml files by running `run_beast_array.sh` (around 16h).

If you want to run Beast several times to combine the runs later, I copied the xml files to a new directory and reran Beast with a new seed (adapt `run_beast_array.sh`). A new `xml_list.txt` file needs to be specified as well. To copy run:
```
rsync -av --include='*/' --include='*.xml' --exclude='*' \
/cluster/work/stadler/beckermar/BESP_paper-analyses/results/pop_size_simulations/simulation_results/ \
/cluster/work/stadler/beckermar/BESP_paper-analyses/results/pop_size_simulations/simulation_results_3/
````


### Create summary trees of runs with sufficient ESS
The `check_mcmc.sh` script select all trees that have a sufficiently high ESS (above $200$) for all estimated parameters and runs treeannotator on all of them with a burnin of $10\%$. 

### Compare true simulated trees with summary trees


## Additional information:
- If you change the number of trees to simulate, also change the number of jobs to run in the array scripts