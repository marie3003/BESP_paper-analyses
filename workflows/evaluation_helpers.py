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

def process_results(input_csv):
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

    # Path to simulated trees file
    sim_tree_base_path = "/Users/mariebecker/Documents/Uni/ETH/RotationStadler/BESP_paper-analyses/results/pop_size_simulations/independent_homochronous"
    sim_tree_mapping = {
        "expgrowth_fast": f"{sim_tree_base_path}/expgrowth_fast/expgrowth_fast.trees",
        "expgrowth_slow": f"{sim_tree_base_path}/expgrowth_slow/expgrowth_slow.trees",
        "uniform": f"{sim_tree_base_path}/uniform/uniform.trees"
    }

    df["sim_tree_path"] = df.apply(determine_sim_tree_path, axis=1, args=(sim_tree_mapping,))
    df["sim_tree_index"] = df["tree_path"].apply(extract_tree_index)

    tree_name = df["tree_path"].apply(lambda p: Path(p).stem)
    df[["model", "population_model", "mutation_signal"]] = tree_name.apply(extract_model_components)

    return df