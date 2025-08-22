# Impact of molecular signal and population model assumptions on the estimation of node heights in time-resolved phylogenetic trees.

This repository contains the simulation study which is the second part of this project.

## Abstract

Time-calibrated phylogenetic trees provide critical insights into pathogen evolution, transmission, and epidemiology, but their accuracy depends on both the availability of mutational signal and the choice of coalescent priors. In this study, we investigated how assumptions about population size models affect estimates of node heights in dated phylogenies.
We derive analytic posteriors for the two-sequence tMRCA under constant size, exponential growth, exponential decline, and bottleneck models, and quantify prior influence using summary-based metrics such as mean shift, median shift, and mode shift, together with the Wasserstein-1 distance and the coalescent information ratio $\Omega$.
By mapping parameter regimes to exemplary pathogens, _Influenza_ with high mutational signal and _Mycobacterium tuberculosis_ with low signal, we show theoretically that prior impact declines with increasing mutational signal and effective population size. Yet, prior influence can be substantial when signal is weak or when population dynamics depart from constancy. Exponential growth can move the posterior to earlier or later times depending on the growth rate, whereas bottlenecks induce non-monotonic effects and secondary posterior peaks.
We then simulate trees under constant and exponential growth population models, generate sequences, and re-estimate times on fixed topologies in BEAST using either a constant coalescent prior or a piecewise-constant skyline prior. The skyline prior recovers temporal population size changes and slightly reduces root-height and deep-branch errors under low signal with fast growth, whereas both priors perform similarly when signal is medium to high or when population models are more similar to constant models.
Overall, demographic misspecification can bias time estimates in bacteria-like settings with limited signal. When the true demography is uncertain, flexible skyline priors are preferable, while constant priors remain adequate for high-signal viral analyses or large effective population sizes in near constant settings.

## Code template

Scripts to simulate trees and create Beast xml files were adapted from the BESP_paper-analysis repository by KV Parag, L du Plessis and OG Pybus with the corresponding paper _Jointly Inferring the Dynamics of Population Size and Sampling Intensity from Molecular Sequences_ ([https://doi.org/10.1093/molbev/msaa016](https://doi.org/10.1093/molbev/msaa016)).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3649734.svg)](https://doi.org/10.5281/zenodo.3649734)

## Dependencies
Used environments are summarized in the `environment_info` folder (for python environments beast_tools and snp_sites both used on the cluster, and beast-tools used for the evaluation notebooks) and the `renv` folder as well as the `renv.lock` file (for the used R environment).  

### BEAST1
BEAST version 1.10.4 was used and needs to be installed. 

### Other süecific packages
The `beastio` packages arise developed by Louis du Plessis and not available on CRAN yet. It can be installed from Github using `devtools::install_github()`.
- [beastio](https://github.com/laduplessis/beastio): 
	- _R-package with functions for pre- and post-processing of BEAST and BEAST2 files (good for automating the functionality in Tracer or Logcombiner e.g. checking convergence in hundreds of replicates from a simulation study or combining many chains)._

## Workflow

The following provides an explenation which scripts and notebooks need to be run to reproduce the simulation study. 

### Simulate trees under a population model
The population models to choose from are specified in the `SimUtils.R` file. Specific models for simulation are selected in the `simulate_trees.R` file which can be run to simulate the trees. The simulation can be run using the `simulate_trees.sh` script (around 3h).

### Simulate alignment with SeqGen and get SNPs with `snp-sites`

Before the alignment can be generated, a list of all trees is needed. This is done by running `bash make_tree_jobs.sh`. Them, the snp-sites environment needs to be activated: `conda activate snp_sites` and the `simulate_alignments.sh` script from the `scripts` folder needs to be run (around 5min).

### Make Beast xml files
For all simulated trees, sequences are generated based on a specified mutation rate and used as an input to Beast. Beast is then run with a constant coalescent population size prior and a Skyline Coalescent prior. To generate the xml files needed as input to Beast, `scripts/MakeBEASTXML.py` can be used. It can be run running `scripts/run_make_beast_xml.sh` after activating the `beast_tools` conda environment. It uses templates specified in the `results/pop_size_simulations/templates` folder and adds parameters specified in the config files that can be found in the `results/pop_size_simulations/config` folder. `MakeBEASTXML.py` subsamples sequences created before with SeqGen and snp sites.

### Run Beast
To run beast, first a list of all beast xml files is needed. This list is generated with `find ../results/pop_size_simulations/simulation_results_3 -name "*.xml" > xml_list_3.txt` run in the `scripts` folder.

Afterwards Beast is run on all xml files by running `run_beast_array.sh` (around 16h).

If you want to run Beast several times to combine the runs later, I copied the xml files to a new directory and reran Beast with a new seed (adapt `run_beast_array.sh`). A new `xml_list.txt` file needs to be specified as well. To copy run:
```
rsync -av --include='*/' --include='*.xml' --exclude='*' \
/cluster/work/stadler/beckermar/BESP_paper-analyses/results/pop_size_simulations/simulation_results/ \
/cluster/work/stadler/beckermar/BESP_paper-analyses/results/pop_size_simulations/simulation_results_3/
````
### Combine runs
To combine runs, the `combine_runs.sh` script needs to be run. The script assumes that all runs were repeated three times and saved in the `simulation_results`, `simulation_results_2` and `simulation_results_3` folders. A burnin of $10\%$ is applied.

### Create summary trees of runs with sufficient ESS
The `check_mcmc.sh` script select all combined runs that have a sufficiently high ESS (above $200$) considering the combined log and trees files. The successful combined runs are then saved in the `successful_mcmc_runs.csv` file that is created in the `scripts` folder. This file was renamed when downloading to `successful_combined_runs.csv`.  The `check_mcmc.sh`file then automatically runs treeannotator on the successful runs creating summary trees.

### Evaluate results
Results are evaluated in the `evaluate_populationsize.ipynb` (only considers population size) and the `evaluate_results.ipynb` (root times and branch lengths). Helper functions are loaded from the `evaluation_helpers.py` file.

### Additional information:
- If you change the number of trees to simulate, also change the number of jobs to run in the array scripts