from Bio import Nexus, Phylo, SeqIO
from collections import defaultdict
from io import StringIO

import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd
import numpy as np

import re
from pathlib import Path

from itertools import combinations


def determine_sim_tree_path(row, sim_tree_mapping):
    for key, path in sim_tree_mapping.items():
        if key in row["tree_path"]:
            return str(path)
    return None

def extract_tree_index(tree_path):
    match = re.search(r"\.T(\d+)\.", tree_path)
    return int(match.group(1)) if match else None

def extract_model_components(tree_name):
    model = "skyline" if "skyline" in tree_name else "constcoal"
    growth_model = "expgrowth_fast" if "expgrowth_fast" in tree_name else \
                   "expgrowth_slow" if "expgrowth_slow" in tree_name else "uniform"
    mutsig = "low" if "lowmutsig" in tree_name else \
             "med" if "medmutsig" in tree_name else "high"
    return pd.Series([model, growth_model, mutsig])

def assign_growth_rate(pop_model, growth_rates):
        if pop_model == "expgrowth_fast":
            return growth_rates[2]
        elif pop_model == "expgrowth_slow":
            return growth_rates[1]
        else:  # uniform
            return growth_rates[0]

def process_results(input_csv, present_pop_size = 2000, growth_rates = [0, 0.001, 0.002]):
    """
    Use all successful runs from the paths csv file.
    """

    df = pd.read_csv(input_csv)

    # Base repo path (local)
    repo_base = Path("/Users/mariebecker/Documents/Uni/ETH/RotationStadler/BESP_paper-analyses")

    # Adjust all paths to be relative to local repo and replace .trees with .tree
    df["tree_path"] = df["trees_path"].apply(
        lambda p: str((repo_base / Path(p).relative_to("/cluster/work/stadler/beckermar/BESP_paper-analyses")).with_suffix(".tree"))
    )

    df["log_path"] = df["trees_path"].apply(
        lambda p: str((repo_base / Path(p).relative_to("/cluster/work/stadler/beckermar/BESP_paper-analyses")).with_suffix(".log"))
    )

    # Extract model components
    tree_name = df["tree_path"].apply(lambda p: Path(p).stem)
    df[["model", "population_model", "mutation_signal"]] = tree_name.apply(extract_model_components)

    df = df.drop(columns=["trees_path", "burnin"])
    df["tree_index"] = df["tree_path"].apply(extract_tree_index)

    # Determine sim_tree_path BEFORE pivoting
    sim_tree_base_path = "/Users/mariebecker/Documents/Uni/ETH/RotationStadler/BESP_paper-analyses/results/pop_size_simulations/independent_homochronous"
    sim_tree_mapping = {
        "expgrowth_fast": f"{sim_tree_base_path}/expgrowth_fast/expgrowth_fast.trees",
        "expgrowth_slow": f"{sim_tree_base_path}/expgrowth_slow/expgrowth_slow.trees",
        "uniform": f"{sim_tree_base_path}/uniform/uniform.trees"
    }
    df["sim_tree_path"] = df.apply(determine_sim_tree_path, axis=1, args=(sim_tree_mapping,))

    # Add constant present population size
    df["present_pop_size"] = present_pop_size

    # Add growth rate based on population model
    df["growth_rate"] = df["population_model"].apply(lambda pop_model: assign_growth_rate(pop_model, growth_rates))

    # Define the shared columns to group by
    group_cols = ['population_model', 'mutation_signal', 'tree_index', 'sim_tree_path', 'present_pop_size', 'growth_rate']

    # Pivot the tree and log paths for each model (constcoal or skyline)
    pivoted = df.pivot_table(
        index=group_cols,
        columns='model',
        values=['tree_path', 'log_path'],
        aggfunc='first'  # in case there are multiple matches, take the first
    ).reset_index()

    # Flatten the column MultiIndex
    pivoted.columns = ['_'.join(col).strip('_') for col in pivoted.columns.values]

    # Rename for clarity
    pivoted = pivoted.rename(columns={
        'tree_path_constcoal': 'tree_path_constcoal',
        'log_path_constcoal': 'log_path_constcoal',
        'tree_path_skyline': 'tree_path_skyline',
        'log_path_skyline': 'log_path_skyline'
    })

    # Optional: sort the result
    pivoted = pivoted.sort_values(by=group_cols).reset_index(drop=True)

    return pivoted

def get_mean_population_size(log_path, burnin=0.1, mode = "constcoal"):
    """
    Parses a BEAST log file and computes the average population size 
    after discarding burn-in.

    Parameters:
        log_path (str): Path to the BEAST .log file.
        burnin (int or float): Burn-in can be given as:
            - a float between 0 and 1 (fraction of the samples to discard), or
            - an int (number of initial states to discard)
        mode (str): The model type, either "constcoal" or "skyline". Determines whether to look for one constant population size or multiple population sizes in the log file.

    Returns:
        float: Average of 'constant.popSize' after burn-in or array of averages for skyline model.
    """
    if pd.isna(log_path):
        return None

    df = pd.read_csv(log_path, comment='#', sep='\t')

    burnin_rows = int(len(df) * burnin)
    df_post_burnin = df.iloc[burnin_rows:]

    if mode == "constcoal":
        mean_pop_size = df_post_burnin['constant.popSize'].mean()
    elif mode == "skyline":
        # For skyline model, we assume 'skyline.popSize' is the column with population sizes
        mean_pop_size = np.array(df_post_burnin.filter(like='skyline.popSize').mean())
    return mean_pop_size

def get_skyline_group_boundaries(tree, num_groups=10):
    """
    Given a dated BEAST tree with internal node times, compute the time points 
    (backwards in time) where Bayesian Skyline population groups end.
    
    Parameters:
        tree (str): Tree structuref
        num_groups (int): Number of skyline population groups (default: 10).
        
    Returns:
        list: Time points (in units of tree height) where each population group ends,
              ordered from present (0) backwards in time.
    """

    root_height = max(tree.depths().values())
    node_times = []

    for clade in tree.get_nonterminals():
        height = root_height - tree.distance(clade)
        node_times.append(height)

    # Sort times ascending (from recent to ancient)
    node_times.sort()

    total_nodes = len(node_times)
    base_group_size = total_nodes // num_groups
    remainder = total_nodes % num_groups

    # Compute how many nodes each group gets
    group_sizes = [base_group_size + 1 if i < remainder else base_group_size for i in range(num_groups)]

    boundaries = []
    index = 0
    for i in range(num_groups):
        index += group_sizes[i]
        if i < num_groups - 1:
            boundaries.append(node_times[index - 1])
        else:
            # Last interval ends at the root
            boundaries.append(root_height)

    return boundaries