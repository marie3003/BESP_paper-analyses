from Bio import Nexus, Phylo, SeqIO
from math import comb, exp, inf
from collections import defaultdict
from io import StringIO

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from matplotlib.container import ErrorbarContainer
from matplotlib import gridspec

import seaborn as sns

import pandas as pd
import numpy as np

import re
from pathlib import Path

import itertools
from itertools import combinations

from scipy.stats import norm, gaussian_kde


def hpd(vals, credible_mass=0.95):
    vals = np.sort(vals)
    n = len(vals)
    width = int(np.ceil(credible_mass * n))
    spans = vals[width - 1:] - vals[:n - width + 1]
    i = np.argmin(spans)
    return float(vals[i]), float(vals[i + width - 1])


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
    pop_match = re.search(r"(expgrowthfast|expgrowthslow|bottleneck[a-z0-9]*|uniform)", tree_name)
    growth_model = pop_match.group(1) if pop_match else "unknown"
    mutsig = "low" if "lowmutsig" in tree_name else \
             "med" if "medmutsig" in tree_name else "high"
    sampling = "linearconstant" if "linearconstant" in tree_name else "independenthomochronous"
    return pd.Series([model, growth_model, mutsig, sampling])

def assign_model_params(pop_model):
    if pop_model == "expgrowthfast":
        return pd.Series({"present_pop_size": 2000, "growth_rate": 0.02,
                          "bottleneck_size": None, "bottleneck_start": None, "bottleneck_end": None})
    elif pop_model == "expgrowthslow":
        return pd.Series({"present_pop_size": 2000, "growth_rate": 0.01,
                          "bottleneck_size": None, "bottleneck_start": None, "bottleneck_end": None})
    elif pop_model == "bottleneck":
        return pd.Series({"present_pop_size": 1000, "growth_rate": None,
                          "bottleneck_size": 100, "bottleneck_start": 20, "bottleneck_end": 23})
    elif pop_model == "bottleneck20":
        return pd.Series({"present_pop_size": 1000, "growth_rate": None,
                          "bottleneck_size": 20, "bottleneck_start": 10, "bottleneck_end": 13})
    elif pop_model == "bottleneck50":
        return pd.Series({"present_pop_size": 1000, "growth_rate": None,
                          "bottleneck_size": 50, "bottleneck_start": 10, "bottleneck_end": 13})
    elif pop_model == "bottleneck100":
        return pd.Series({"present_pop_size": 1000, "growth_rate": None,
                          "bottleneck_size": 100, "bottleneck_start": 10, "bottleneck_end": 13})
    elif pop_model == "bottleneck200":
        return pd.Series({"present_pop_size": 1000, "growth_rate": None,
                          "bottleneck_size": 200, "bottleneck_start": 10, "bottleneck_end": 13})
    elif pop_model in ("bottlenecklate", "bottlenecklatesampling"):
        return pd.Series({"present_pop_size": 1000, "growth_rate": None,
                          "bottleneck_size": 10, "bottleneck_start": 50, "bottleneck_end": 53})
    elif pop_model == "bottleneckmid50":
        return pd.Series({"present_pop_size": 1000, "growth_rate": None,
                          "bottleneck_size": 50, "bottleneck_start": 20, "bottleneck_end": 23})
    elif pop_model == "bottleneckmid100":
        return pd.Series({"present_pop_size": 1000, "growth_rate": None,
                          "bottleneck_size": 100, "bottleneck_start": 20, "bottleneck_end": 23})
    else:  # uniform
        return pd.Series({"present_pop_size": 1000, "growth_rate": None,
                          "bottleneck_size": None, "bottleneck_start": None, "bottleneck_end": None})

def process_results(input_csv, population_models=None):
    """
    Use all successful runs from the paths csv file.
    """
    if population_models is None:
        population_models = ["expgrowthfast", "expgrowthslow", "uniform", "bottleneck"]

    repo_base = Path("/Users/mariebecker/Documents/Uni/ETH/RotationStadler/BESP_paper-analyses")

    # infer run folder from the csv path (e.g. results/run_bottleneck/... -> run_bottleneck)
    csv_parts = Path(input_csv).parts
    run_folder = next((p for p in csv_parts if p.startswith("run")), "run1")

    df = pd.read_csv(input_csv)

    def to_local_path(p, suffix):
        rel = p.replace("/cluster/work/stadler/beckermar/BESP_paper-analyses/", "")
        rel = rel.lstrip("/")
        return str(repo_base / rel.replace(".combined.trees", suffix))

    df["tree_path"] = df["trees_path"].apply(lambda p: to_local_path(p, ".combined_summary.tree"))
    df["log_path"]  = df["trees_path"].apply(lambda p: to_local_path(p, ".combined.log"))
    df = df.drop(columns=["trees_path"])

    tree_name = df["tree_path"].apply(lambda p: Path(p).stem)
    df[["model", "population_model", "mutation_signal", "sampling"]] = tree_name.apply(extract_model_components)
    df["tree_index"] = df["tree_path"].apply(extract_tree_index)

    # Split by model and merge on shared keys
    shared_cols = ['population_model', 'mutation_signal', 'sampling', 'tree_index']
    constcoal = df[df["model"] == "constcoal"][shared_cols + ["tree_path", "log_path"]].rename(
        columns={"tree_path": "tree_path_constcoal", "log_path": "log_path_constcoal"})
    skyline = df[df["model"] == "skyline"][shared_cols + ["tree_path", "log_path"]].rename(
        columns={"tree_path": "tree_path_skyline", "log_path": "log_path_skyline"})

    merged = constcoal.merge(skyline, on=shared_cols, how="outer")

    # Build complete expected index and left-join to find completely missing runs
    all_combos = pd.DataFrame(
        list(itertools.product(
            population_models,
            ["low", "med", "high"],
            ["independenthomochronous", "linearconstant"],
            range(100)
        )),
        columns=["population_model", "mutation_signal", "sampling", "tree_index"]
    )

    merged = all_combos.merge(merged, on=shared_cols, how="left")

    # Add sim_tree_path
    sim_base = repo_base / "results" / run_folder / "simulated_data"
    merged["sim_tree_path"] = merged.apply(
        lambda r: str(sim_base / r["sampling"] / r["population_model"] / f"{r['population_model']}.trees"),
        axis=1
    )

    # Add model-specific parameters
    params = merged["population_model"].apply(assign_model_params)
    merged = pd.concat([merged, params], axis=1)
    merged = merged.sort_values(by=shared_cols).reset_index(drop=True)

    # Split into successful and unsuccessful
    successful = merged[merged["tree_path_constcoal"].notna() & merged["tree_path_skyline"].notna()].reset_index(drop=True)

    unsuccessful = merged[merged["tree_path_constcoal"].isna() | merged["tree_path_skyline"].isna()].copy()
    unsuccessful["failed"] = unsuccessful.apply(
        lambda r: "both" if pd.isna(r["tree_path_constcoal"]) and pd.isna(r["tree_path_skyline"])
                  else "constcoal" if pd.isna(r["tree_path_constcoal"])
                  else "skyline", axis=1
    )
    unsuccessful = unsuccessful[shared_cols + ["failed"]].reset_index(drop=True)

    return successful, unsuccessful

def get_median_population_size(log_path, burnin=0, mode="constcoal"):
    """
    Parses a BEAST log file and computes the median population size 
    after discarding burn-in, along with 95% confidence intervals.

    Parameters:
        log_path (str): Path to the BEAST .log file.
        burnin (int or float): Burn-in can be given as:
            - a float between 0 and 1 (fraction of the samples to discard), or
            - an int (number of initial states to discard)
        mode (str): Either "constcoal" or "skyline".

    Returns:
        If mode == "constcoal":
            tuple: (median, lower_95, upper_95) of constant.popSize
        If mode == "skyline":
            np.ndarray: shape (n_intervals, 3), each row is (median, lower_95, upper_95)
    """
    df = pd.read_csv(log_path, comment='#', sep='\t')

    #burnin_rows = int(len(df) * burnin) if isinstance(burnin, float) else burnin
    #df_post_burnin = df.iloc[burnin_rows:]

    if mode == "constcoal":
        values = df['constant.popSize']
        median = np.median(values)
        lower = np.percentile(values, 2.5)
        upper = np.percentile(values, 97.5)
        samples = np.random.choice(values, size=100, replace=False)
        return median, lower, upper, samples

    elif mode == "skyline":
        skyline_df = df.filter(like='skyline.popSize')
        medians = skyline_df.median().values
        lowers = skyline_df.quantile(0.025).values
        uppers = skyline_df.quantile(0.975).values
        sampled_df = skyline_df.sample(n=100)
        return medians, lowers, uppers, sampled_df
        #return np.stack([medians, lowers, uppers], axis=1)  # shape: (n_intervals, 3)

    else:
        raise ValueError("mode must be either 'constcoal' or 'skyline'")


def process_skygrid_logs(skygrid_dir, burnin=0.1):
    """
    Scan a skygrid beast_inference directory, parse each .log file, and return
    a DataFrame with per-grid-point median and 95% HPD population sizes.

    Directory structure expected:
        skygrid_dir / sampling / pop_model / mutsig / *.log

    Parameters:
        skygrid_dir: path to the skygrid beast_inference folder
        burnin: fraction of samples to discard as burn-in (default 0.1 = 10%)

    Returns:
        DataFrame with columns:
            sampling, population_model, mutation_signal, tree_index,
            grid_point (0-indexed), time (years before present),
            median, hpd_lower, hpd_upper (population size, linear scale)
    """
    skygrid_dir = Path(skygrid_dir)
    records = []

    for log_path in sorted(skygrid_dir.rglob("*.log")):
        if "slurm" in log_path.name:
            continue

        # parse metadata from filename
        components = extract_model_components(log_path.stem)
        model, pop_model, mutsig, sampling = components

        match = re.search(r"\.T(\d+)$", log_path.stem)
        if not match:
            continue
        tree_index = int(match.group(1))

        df = pd.read_csv(log_path, comment="#", sep="\t")

        # apply burn-in
        n_burnin = int(len(df) * burnin) if isinstance(burnin, float) else burnin
        df = df.iloc[n_burnin:]

        logpop_cols = sorted(
            [c for c in df.columns if c.startswith("skygrid.logPopSize")],
            key=lambda c: int(re.search(r"\d+$", c).group())
        )
        if not logpop_cols:
            continue

        cut_off = df["skygrid.cutOff"].iloc[0]
        n_points = len(logpop_cols)
        # left boundaries: 0, 5, 10, ..., cutOff (evenly spaced)
        # last interval extends from cutOff to the root — store median root height separately
        times = np.linspace(0, cut_off, n_points)
        median_root_height = np.median(df["treeModel.rootHeight"].values)

        row = {
            "sampling": sampling,
            "population_model": pop_model,
            "mutation_signal": mutsig,
            "tree_index": tree_index,
            "cut_off": cut_off,
            "median_root_height": median_root_height,
        }

        for i, col in enumerate(logpop_cols):
            pop_samples = np.exp(df[col].values)
            median = np.median(pop_samples)
            sorted_samples = np.sort(pop_samples)
            n = len(sorted_samples)
            hpd_window = int(np.floor(0.95 * n))
            min_width = np.inf
            hpd_lower, hpd_upper = sorted_samples[0], sorted_samples[-1]
            for j in range(n - hpd_window):
                width = sorted_samples[j + hpd_window] - sorted_samples[j]
                if width < min_width:
                    min_width = width
                    hpd_lower = sorted_samples[j]
                    hpd_upper = sorted_samples[j + hpd_window]

            row[f"median_{i}"] = median
            row[f"hpd_lower_{i}"] = hpd_lower
            row[f"hpd_upper_{i}"] = hpd_upper
            row[f"time_{i}"] = times[i]

        records.append(row)

    return pd.DataFrame(records)

    




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


### CALCULATE NODE HEIGHTS AND BRANCH LENGTHS
def get_branch_info(tree, simulated=False):
    """
    Extract branch length and node height information from a tree.

    For simulated trees (simulated=True): plain newick, no annotations.
      - Branch lengths from clade.branch_length
      - Heights computed as root_height - distance(root, clade), where
        root_height = max root-to-tip distance (correct since times are
        shifted so most recent sample is at t=0)
      - No CI intervals available

    For estimated trees (simulated=False): BEAST summary trees with
      TreeAnnotator annotations.
      - Branch lengths and heights from length_median / height_median
      - 95% HPD intervals from length_95%_HPD / height_95%_HPD

    Returns a DataFrame with columns:
    node, bl, height, bl_lower_ci, bl_upper_ci, height_lower_ci, height_upper_ci, internal
    """
    rows = []
    i = 0

    if simulated:
        root_height = max(tree.distance(tree.root, leaf) for leaf in tree.get_terminals())

    for clade in tree.find_clades(order="preorder"):
        if clade.name:
            clade_id = clade.name
            internal = False
        else:
            clade_id = f"internal_{i}"
            i += 1
            internal = True

        if simulated:
            assert clade.branch_length is not None or clade == tree.root, \
                f"Missing branch length on non-root node {clade_id}"
            bl = clade.branch_length if clade.branch_length is not None else 0.0
            height = root_height - tree.distance(tree.root, clade)
            rows.append({
                "node": clade_id, "bl": bl, "height": height,
                "bl_lower_ci": None, "bl_upper_ci": None,
                "height_lower_ci": None, "height_upper_ci": None,
                "internal": internal
            })
        else:
            assert clade.comment is not None or clade == tree.root, \
                f"Missing annotations on non-root node {clade_id}"
            bl = 0.0
            height = 0.0
            bl_lower_ci, bl_upper_ci = None, None
            height_lower_ci, height_upper_ci = None, None

            if clade.comment:
                bl_median_match = re.search(r'length_median=([\d\.eE+-]+)', clade.comment)
                if bl_median_match:
                    bl = float(bl_median_match.group(1))

                h_median_match = re.search(r'height_median=([\d\.eE+-]+)', clade.comment)
                if h_median_match:
                    height = float(h_median_match.group(1))

                bl_match = re.search(r'length_95%_HPD=\{([\d\.eE+-]+),([\d\.eE+-]+)\}', clade.comment)
                if bl_match:
                    bl_lower_ci = float(bl_match.group(1))
                    bl_upper_ci = float(bl_match.group(2))

                h_match = re.search(r'height_95%_HPD=\{([\d\.eE+-]+),([\d\.eE+-]+)\}', clade.comment)
                if h_match:
                    height_lower_ci = float(h_match.group(1))
                    height_upper_ci = float(h_match.group(2))

            rows.append({
                "node": clade_id, "bl": bl, "height": height,
                "bl_lower_ci": bl_lower_ci, "bl_upper_ci": bl_upper_ci,
                "height_lower_ci": height_lower_ci, "height_upper_ci": height_upper_ci,
                "internal": internal
            })

    return pd.DataFrame(rows)


def compare_tree_metrics(tree_sim, tree_estimate):
    """
    Compare branch lengths and node heights between a simulated tree and a BEAST
    summary tree (constcoal or skyline) node by node.

    Assumes both trees have the same topology and node ordering (as guaranteed when
    BEAST is run with the true simulated tree as starting tree).

    Returns a DataFrame with one row per node containing:
    - bl_sim / bl_estimate: median branch lengths from simulated and estimated tree
    - bl_ci_lower/upper_estimate: 95% HPD for estimated branch length
    - bl_inside_ci: whether the true branch length falls within the 95% HPD
    - height_sim / height_estimate: median node heights
    - height_ci_lower/upper_estimate: 95% HPD for estimated node height
    - height_inside_ci: whether the true height falls within the 95% HPD
    - bl_diff, bl_relative_error, bl_abs_relative_error: branch length errors
    - height_diff, height_relative_error, height_abs_relative_error: height errors
    - internal: whether the node is an internal node (True) or tip (False)
    """
    sim_info = get_branch_info(tree_sim, simulated=True).rename(columns={"bl": "bl_sim", "height": "height_sim"})
    est_info = get_branch_info(tree_estimate, simulated=False).rename(columns={
        "bl": "bl_estimate", "height": "height_estimate",
        "bl_lower_ci": "bl_ci_lower_estimate", "bl_upper_ci": "bl_ci_upper_estimate",
        "height_lower_ci": "height_ci_lower_estimate", "height_upper_ci": "height_ci_upper_estimate",
    })

    branch_length_df = pd.concat([
        sim_info[["node", "bl_sim", "height_sim"]],
        est_info[["bl_estimate", "bl_ci_lower_estimate", "bl_ci_upper_estimate",
                  "height_estimate", "height_ci_lower_estimate", "height_ci_upper_estimate", "internal"]]
    ], axis=1)

    branch_length_df["bl_inside_ci"] = (
        (branch_length_df["bl_ci_lower_estimate"] <= branch_length_df["bl_sim"]) &
        (branch_length_df["bl_sim"] <= branch_length_df["bl_ci_upper_estimate"])
    ).where(branch_length_df["bl_ci_lower_estimate"].notna())

    branch_length_df["height_inside_ci"] = (
        (branch_length_df["height_ci_lower_estimate"] <= branch_length_df["height_sim"]) &
        (branch_length_df["height_sim"] <= branch_length_df["height_ci_upper_estimate"])
    ).where(branch_length_df["height_ci_lower_estimate"].notna())

    def _abs_ci(signed_lower, signed_upper):
        crosses_zero = (signed_lower <= 0) & (signed_upper >= 0)
        abs_lower = np.where(crosses_zero, 0.0, np.minimum(np.abs(signed_lower), np.abs(signed_upper)))
        abs_upper = np.maximum(np.abs(signed_lower), np.abs(signed_upper))
        return abs_lower, abs_upper

    # Errors for branch lengths
    bl_sim = branch_length_df["bl_sim"]
    branch_length_df["bl_diff"]                = branch_length_df["bl_estimate"] - bl_sim
    branch_length_df["bl_diff_lower"]          = branch_length_df["bl_ci_lower_estimate"] - bl_sim
    branch_length_df["bl_diff_upper"]          = branch_length_df["bl_ci_upper_estimate"] - bl_sim
    branch_length_df["bl_relative_error"]       = branch_length_df["bl_diff"] / bl_sim
    branch_length_df["bl_relative_error_lower"] = branch_length_df["bl_diff_lower"] / bl_sim
    branch_length_df["bl_relative_error_upper"] = branch_length_df["bl_diff_upper"] / bl_sim
    branch_length_df["bl_abs_relative_error"] = np.abs(branch_length_df["bl_relative_error"])
    bl_abs_lo, bl_abs_up = _abs_ci(branch_length_df["bl_relative_error_lower"],
                                    branch_length_df["bl_relative_error_upper"])
    branch_length_df["bl_abs_relative_error_lower"] = bl_abs_lo
    branch_length_df["bl_abs_relative_error_upper"] = bl_abs_up

    # Errors for node heights
    height_sim = branch_length_df["height_sim"]
    internal   = branch_length_df["internal"]
    branch_length_df["height_diff"]                = branch_length_df["height_estimate"] - height_sim
    branch_length_df["height_diff_lower"]          = branch_length_df["height_ci_lower_estimate"] - height_sim
    branch_length_df["height_diff_upper"]          = branch_length_df["height_ci_upper_estimate"] - height_sim
    safe_h = height_sim.replace(0, np.nan)
    branch_length_df["height_relative_error"]       = np.where(internal, branch_length_df["height_diff"] / safe_h, 0.0)
    branch_length_df["height_relative_error_lower"] = np.where(internal, branch_length_df["height_diff_lower"] / safe_h, 0.0)
    branch_length_df["height_relative_error_upper"] = np.where(internal, branch_length_df["height_diff_upper"] / safe_h, 0.0)
    branch_length_df["height_abs_relative_error"] = np.abs(branch_length_df["height_relative_error"])
    h_abs_lo, h_abs_up = _abs_ci(branch_length_df["height_relative_error_lower"],
                                   branch_length_df["height_relative_error_upper"])
    branch_length_df["height_abs_relative_error_lower"] = np.where(internal, h_abs_lo, np.nan)
    branch_length_df["height_abs_relative_error_upper"] = np.where(internal, h_abs_up, np.nan)

    return branch_length_df

def tree_metrics_all_trees(path_df):
    """
    For each replicate in path_df, load the true simulated tree and both BEAST
    summary trees (constcoal and skyline), then compute node-level branch length
    and height errors via compare_tree_metrics.

    Input columns required in path_df:
    - tree_path_constcoal, tree_path_skyline: local paths to BEAST summary trees
    - sim_tree_path: path to newick file containing all simulated trees
    - tree_index: which tree (0-indexed line) to extract from sim_tree_path

    Returns a long-format DataFrame with one row per node per replicate per model,
    with columns for errors, CIs, model metadata (model, growth_model, mutsig, sampling).
    """
    results = []
    for sim_tree_path, group in path_df.groupby("sim_tree_path"):
        sim_trees = list(Phylo.parse(sim_tree_path, "newick"))
        for _, row in group.iterrows():
            tree_sim      = sim_trees[row["tree_index"]]
            tree_constcoal = Phylo.read(row["tree_path_constcoal"], "nexus")
            tree_skyline   = Phylo.read(row["tree_path_skyline"], "nexus")

            eval_df_constcoal = compare_tree_metrics(tree_sim, tree_constcoal)
            eval_df_skyline   = compare_tree_metrics(tree_sim, tree_skyline)

            eval_df_constcoal["tree_name"]  = Path(row["tree_path_constcoal"]).stem
            eval_df_skyline["tree_name"]    = Path(row["tree_path_skyline"]).stem
            eval_df_constcoal["tree_index"] = row["tree_index"]
            eval_df_skyline["tree_index"]   = row["tree_index"]
            results.append(eval_df_constcoal)
            results.append(eval_df_skyline)

    # Save combined output
    combined_df = pd.concat(results, ignore_index=True)

    combined_df[["model", "population_model", "mutation_signal", "sampling"]] = combined_df["tree_name"].apply(extract_model_components)

    combined_df["height_abs_diff"] = np.abs(combined_df["height_diff"])
    combined_df["bl_abs_diff"] = np.abs(combined_df["bl_diff"])

    return combined_df

def root_height_all_trees(path_df):
    """
    For each replicate in path_df, compare the root height of the simulated tree
    against the constcoal and skyline BEAST summary trees.

    Returns a long-format DataFrame with one row per replicate per model, with
    columns: tree_name, tree_index, sim_root_height, estimated_root_height,
    model, growth_model, mutsig, and derived error metrics.
    """
    results = []
    for _, row in path_df.iterrows():
        tree_constcoal = Phylo.read(row["tree_path_constcoal"], "nexus")
        tree_skyline = Phylo.read(row["tree_path_skyline"], "nexus")
        sim_trees = list(Phylo.parse(row["sim_tree_path"], "newick"))
        tree_sim = sim_trees[row["tree_index"]]
        
        sim_root_height = max(tree_sim.distance(tree_sim.root, leaf) for leaf in tree_sim.get_terminals())
        constcoal_root_height = max(tree_constcoal.distance(tree_constcoal.root, leaf) for leaf in tree_constcoal.get_terminals())
        skyline_root_height = max(tree_skyline.distance(tree_skyline.root, leaf) for leaf in tree_skyline.get_terminals())

        
        # Tree name and index
        tree_name_constcoal = Path(row["tree_path_constcoal"]).stem
        tree_name_skyline = Path(row["tree_path_skyline"]).stem
        tree_index = row["tree_index"]

        # Add both model estimates as separate rows
        results.append({
            "tree_name": tree_name_constcoal,
            "tree_index": tree_index,
            "sim_root_height": sim_root_height,
            "estimated_root_height": constcoal_root_height
        })

        results.append({
            "tree_name": tree_name_skyline,
            "tree_index": tree_index,
            "sim_root_height": sim_root_height,
            "estimated_root_height": skyline_root_height
        })

    # Save combined output
    combined_df = pd.DataFrame(results)

    combined_df[["model", "population_model", "mutation_signal", "sampling"]] = combined_df["tree_name"].apply(extract_model_components)

    combined_df['diff_root_height'] = combined_df['estimated_root_height'] - combined_df['sim_root_height']
    combined_df['abs_diff_root_height'] = np.abs(combined_df['diff_root_height'])
    combined_df['rel_diff_root_height'] = combined_df['diff_root_height'] / combined_df['sim_root_height']
    combined_df['abs_rel_diff_root_height'] = np.abs(combined_df['rel_diff_root_height'])

    return combined_df

### POPULATION SIZE ERROR
def compute_cumulative_population(time_points, populations):
    """
    Compute the cumulative population size at each interval end time point,
    assuming population is constant within each interval.

    Parameters:
    - time_points: array of interval end time points (length N)
    - populations: array of population size within each interval (length N)

    Returns:
    - cumulative: array of cumulative population size at each end time point (length N)
    """
    time_points = np.asarray(time_points)
    
    # Compute interval widths
    start_points = np.concatenate(([0], time_points[:-1]))
    widths = time_points - start_points

    # Compute population contribution in each interval
    contributions = widths * populations

    # Cumulative sum, prepend 0 for the start
    cumulative = np.concatenate(([0], np.cumsum(contributions)))
    return cumulative

def cumulative_exp_pop_size(t, present_pop_size, growth_rate):
    if growth_rate != 0:
        return present_pop_size/growth_rate * (1 -  np.exp(-growth_rate * t))
    else:
        return present_pop_size * t

def extract_tree_info(row, num_groups=10, burnin=0):

    # Load tree
    tree = Phylo.read(row["tree_path_skyline"], "nexus")

    # Load simulated tree and compute branch length statistics
    sim_trees = list(Phylo.parse(row["sim_tree_path"], "newick"))
    tree_sim = sim_trees[row["tree_index"]]
    tip_names = {c.name for c in tree_sim.get_terminals()}
    total_bl_sim = sum(c.branch_length for c in tree_sim.find_clades() if c.branch_length is not None)
    tip_bl_sim = sum(c.branch_length for c in tree_sim.find_clades()
                     if c.branch_length is not None and c.name in tip_names)

    # Extract population size estimates
    skyline_times = get_skyline_group_boundaries(tree, num_groups=num_groups)
    skyline_medians, skyline_lowers, skyline_uppers, skyline_samples_df  = get_median_population_size(row["log_path_skyline"], burnin=burnin, mode="skyline")
    coalescent_median, coalescent_lower, coalescent_upper, coalescent_samples = get_median_population_size(row["log_path_constcoal"], burnin=burnin, mode="constcoal")

    skyline_cumulative_medians = compute_cumulative_population(skyline_times, skyline_medians)
    skyline_cumulative_lowers  = compute_cumulative_population(skyline_times, skyline_lowers)
    skyline_cumulative_uppers  = compute_cumulative_population(skyline_times, skyline_uppers)

    return pd.Series({
            "skyline_medians": skyline_medians,
            "skyline_times": skyline_times,
            "skyline_cumulative_medians": skyline_cumulative_medians,
            "skyline_cumulative_lowers": skyline_cumulative_lowers,
            "skyline_cumulative_uppers": skyline_cumulative_uppers,
            "coalescent_median": coalescent_median,
            "skyline_lowers": skyline_lowers,
            "skyline_uppers": skyline_uppers,
            "coalescent_lower": coalescent_lower,
            "coalescent_upper": coalescent_upper,
            "skyline_samples": skyline_samples_df,
            "coalescent_samples": coalescent_samples,
            "total_branch_length_sim": total_bl_sim,
            "tip_branch_length_sim": tip_bl_sim})


def add_tree_information(df, num_groups = 10, burnin = 0):
    estimates = df.apply(lambda row: extract_tree_info(row, num_groups=num_groups, burnin=burnin), axis=1)
    return pd.concat([df, estimates], axis=1)

def _abs_error_ci(signed_lower, signed_upper):
    """Compute CI bounds for absolute error given signed lower/upper bounds."""
    crosses_zero = signed_lower <= 0 <= signed_upper
    abs_lower = 0.0 if crosses_zero else min(abs(signed_lower), abs(signed_upper))
    abs_upper = max(abs(signed_lower), abs(signed_upper))
    return abs_lower, abs_upper


def calculate_population_size_error(
    t,
    true_traj,
    skyline_times=None,
    skyline_medians=None,
    skyline_lowers=None,
    skyline_uppers=None,
    coalescent_median=None,
    coalescent_lower=None,
    coalescent_upper=None,
):
    true_pop_size = true_traj(t)

    if coalescent_median is not None:
        est_pop_size = coalescent_median
        est_lower = coalescent_lower
        est_upper = coalescent_upper
    elif skyline_times is not None and skyline_medians is not None:
        idx = max(np.searchsorted(skyline_times, t, side="left") - 1, 0)
        est_pop_size = skyline_medians[idx]
        est_lower = skyline_lowers[idx] if skyline_lowers is not None else None
        est_upper = skyline_uppers[idx] if skyline_uppers is not None else None
    else:
        raise ValueError("Either coalescent_median or skyline_times + skyline_medians must be provided.")

    error = est_pop_size - true_pop_size
    rel_error = error / true_pop_size

    diff_lower = (est_lower - true_pop_size) if est_lower is not None else np.nan
    diff_upper = (est_upper - true_pop_size) if est_upper is not None else np.nan
    rel_lower  = diff_lower / true_pop_size if est_lower is not None else np.nan
    rel_upper  = diff_upper / true_pop_size if est_upper is not None else np.nan

    abs_diff_lo, abs_diff_up = _abs_error_ci(diff_lower, diff_upper) if est_lower is not None else (np.nan, np.nan)
    abs_rel_lo,  abs_rel_up  = _abs_error_ci(rel_lower,  rel_upper)  if est_lower is not None else (np.nan, np.nan)

    return {
        "true_pop_size": true_pop_size,
        "est_pop_size": est_pop_size,
        "diff_pop_size": error,
        "diff_pop_size_lower": diff_lower,
        "diff_pop_size_upper": diff_upper,
        "abs_diff_pop_size": abs(error),
        "abs_diff_pop_size_lower": abs_diff_lo,
        "abs_diff_pop_size_upper": abs_diff_up,
        "rel_diff_pop_size": rel_error,
        "rel_diff_pop_size_lower": rel_lower,
        "rel_diff_pop_size_upper": rel_upper,
        "abs_rel_diff_pop_size": abs(rel_error),
        "abs_rel_diff_pop_size_lower": abs_rel_lo,
        "abs_rel_diff_pop_size_upper": abs_rel_up,
    }

def cumulative_true_pop_size(t, pop_model, row):
    """Analytical cumulative integral of the true population size from 0 to t."""
    if pop_model in ("expgrowthfast", "expgrowthslow"):
        N0, r = row["present_pop_size"], row["growth_rate"]
        return N0 / r * (1 - np.exp(-r * t))
    elif pop_model == "uniform":
        return row["present_pop_size"] * t
    elif pop_model.startswith("bottleneck"):
        N0 = row["present_pop_size"]
        Nb = row["bottleneck_size"]
        bs = row["bottleneck_start"]
        be = row["bottleneck_end"]
        if t <= bs:
            return N0 * t
        elif t < be:
            return N0 * bs + Nb * (t - bs)
        else:
            return N0 * bs + Nb * (be - bs) + N0 * (t - be)
    else:
        raise ValueError(f"Unknown pop_model: {pop_model}")


def calculate_population_size_cumulative_error(t,
                                               skyline_times=None,
                                               skyline_cumulative_medians=None, skyline_medians=None,
                                               skyline_cumulative_lowers=None, skyline_lowers=None,
                                               skyline_cumulative_uppers=None, skyline_uppers=None,
                                               coalescent_median=None, coalescent_lower=None, coalescent_upper=None,
                                               pop_model=None, row=None):
    cum_pop_size_sim = cumulative_true_pop_size(t, pop_model, row)

    def _cum_est(cum_medians, medians):
        idx = np.searchsorted(skyline_times, t, side="left") - 1
        start_times = [0] + list(skyline_times)
        return cum_medians[idx] + medians[idx] * (t - start_times[idx])

    if coalescent_median is not None:
        cum_est     = coalescent_median * t
        cum_est_lo  = coalescent_lower * t if coalescent_lower is not None else np.nan
        cum_est_up  = coalescent_upper * t if coalescent_upper is not None else np.nan
    elif skyline_times is not None and skyline_cumulative_medians is not None:
        cum_est    = _cum_est(skyline_cumulative_medians, skyline_medians)
        cum_est_lo = _cum_est(skyline_cumulative_lowers, skyline_lowers) if skyline_cumulative_lowers is not None else np.nan
        cum_est_up = _cum_est(skyline_cumulative_uppers, skyline_uppers) if skyline_cumulative_uppers is not None else np.nan

    error    = cum_est    - cum_pop_size_sim
    error_lo = cum_est_lo - cum_pop_size_sim
    error_up = cum_est_up - cum_pop_size_sim
    abs_lo, abs_up = _abs_error_ci(error_lo, error_up) if not np.isnan(error_lo) else (np.nan, np.nan)

    return {
        "cum_pop_size_sim":          cum_pop_size_sim,
        "cum_pop_size_est":          cum_est,
        "cum_pop_size_error":        error,
        "cum_pop_size_error_lower":  error_lo,
        "cum_pop_size_error_upper":  error_up,
        "cum_pop_size_abs_error_lower": abs_lo,
        "cum_pop_size_abs_error_upper": abs_up,
    }

def add_population_size_errors(node_df, tree_df):
    """
    Merges node-level tree metrics with population size estimates from tree_df,
    then computes pointwise and cumulative population size errors at each node's
    time against the true trajectory for that pop model.
    """
    merged_df = node_df.merge(
        tree_df,
        on=["tree_index", "mutation_signal", "population_model", "sampling"],
        how="left",
        suffixes=("", "_tree")
    )

    def compute_errors(row):
        is_const = row["model"] == "constcoal"
        is_skyline = row["model"] == "skyline"
        true_traj = get_true_traj(row["population_model"], row)

        pointwise = calculate_population_size_error(
            t=row["height_sim"],
            true_traj=true_traj,
            coalescent_median=row["coalescent_median"] if is_const else None,
            coalescent_lower=row["coalescent_lower"] if is_const else None,
            coalescent_upper=row["coalescent_upper"] if is_const else None,
            skyline_times=row["skyline_times"] if is_skyline else None,
            skyline_medians=row["skyline_medians"] if is_skyline else None,
            skyline_lowers=row["skyline_lowers"] if is_skyline else None,
            skyline_uppers=row["skyline_uppers"] if is_skyline else None,
        )

        cumulative = calculate_population_size_cumulative_error(
            t=row["height_sim"],
            coalescent_median=row["coalescent_median"] if is_const else None,
            coalescent_lower=row["coalescent_lower"] if is_const else None,
            coalescent_upper=row["coalescent_upper"] if is_const else None,
            skyline_times=row["skyline_times"] if is_skyline else None,
            skyline_cumulative_medians=row["skyline_cumulative_medians"] if is_skyline else None,
            skyline_cumulative_lowers=row["skyline_cumulative_lowers"] if is_skyline else None,
            skyline_cumulative_uppers=row["skyline_cumulative_uppers"] if is_skyline else None,
            skyline_medians=row["skyline_medians"] if is_skyline else None,
            skyline_lowers=row["skyline_lowers"] if is_skyline else None,
            skyline_uppers=row["skyline_uppers"] if is_skyline else None,
            pop_model=row["population_model"],
            row=row
        )

        return pd.Series({**pointwise, **cumulative})

    error_df = merged_df.apply(compute_errors, axis=1)
    result = pd.concat([merged_df, error_df], axis=1)

    result['cum_pop_size_abs_error'] = np.abs(result['cum_pop_size_error'])
    result['cum_pop_size_error_relsim'] = result['cum_pop_size_error'] / result['cum_pop_size_sim']

    drop_cols = [
        # list/array columns
        'skyline_times', 'skyline_medians', 'skyline_cumulative_medians',
        'skyline_cumulative_lowers', 'skyline_cumulative_uppers',
        'skyline_lowers', 'skyline_uppers', 'skyline_samples', 'coalescent_samples',
        # path columns
        'tree_path_constcoal', 'log_path_constcoal', 'tree_path_skyline',
        'log_path_skyline', 'sim_tree_path',
        # scalar overview_df_ columns (kept there, not needed per-node)
        'root_height', 'total_branch_length_sim', 'tip_branch_length_sim',
        'present_pop_size', 'growth_rate', 'bottleneck_size', 'bottleneck_start', 'bottleneck_end',
        'coalescent_median', 'coalescent_lower', 'coalescent_upper',
    ]
    return result.drop(columns=[c for c in drop_cols if c in result.columns])

def get_root_height_df(tree_metrics_combined):
    """Extract root node rows and drop branch length columns, returning one row per replicate per model."""
    drop_cols = ['node', 'bl_sim', 'bl_estimate', 'bl_ci_lower_estimate', 'bl_ci_upper_estimate',
                 'bl_inside_ci', 'bl_diff', 'bl_relative_error', 'bl_abs_relative_error', 'bl_abs_diff']
    return (tree_metrics_combined[tree_metrics_combined.node == 'internal_0']
            .copy()
            .reset_index(drop=True)
            .drop(columns=[c for c in drop_cols if c in tree_metrics_combined.columns]))


def get_root_height_df_wide(root_height_df):
    """Pivot root height DataFrame so constcoal and skyline estimates are in separate columns."""
    # Columns that differ per model and should be pivoted
    pivot_cols = ['height_estimate', 'height_ci_lower_estimate', 'height_ci_upper_estimate',
                  'height_inside_ci', 'height_diff', 'height_relative_error',
                  'height_abs_relative_error', 'height_abs_diff',
                  'est_pop_size', 'diff_pop_size', 'abs_diff_pop_size',
                  'rel_diff_pop_size', 'abs_rel_diff_pop_size',
                  'cum_pop_size_est', 'cum_pop_size_error', 'cum_pop_size_abs_error', 'cum_pop_size_error_relsim', 'root_cum_pop_size_error',
                  'rel_cum_pop_size_error',
                  'ci_true', 'ci_estimated', 'ci_diff', 'ci_rel_error', 'ci_abs_diff', 'ci_abs_rel_error',
                  'rate_true', 'rate_estimated', 'rate_diff', 'rate_rel_error', 'rate_abs_diff', 'rate_abs_rel_error']

    # Columns identical across models — keep once as index
    shared_cols = ['tree_index', 'mutation_signal', 'population_model', 'sampling',
                   'height_sim', 'true_pop_size', 'cum_pop_size_sim']

    pivot_cols = [c for c in pivot_cols if c in root_height_df.columns]
    shared_cols = [c for c in shared_cols if c in root_height_df.columns]

    root_height_df_wide = root_height_df.pivot_table(
        index=shared_cols,
        columns='model',
        values=pivot_cols
    ).reset_index()

    root_height_df_wide.columns = ['_'.join(col).strip('_') for col in root_height_df_wide.columns.values]

    return root_height_df_wide

def evaluate_ci_overlap(root_height_df_wide, sampling):

    def estimate_sigma(lower, upper):
        return (upper - lower) / (2 * 1.96)

    def compute_overlap(mu1, sigma1, mu2, sigma2):
        d = abs(mu1 - mu2)
        denom = np.sqrt(2 * (sigma1**2 + sigma2**2))
        return 2 * norm.cdf(-d / denom)
    
    root_height_df_wide = root_height_df_wide[root_height_df_wide["sampling"] == sampling].copy()

    mutsigs = ['high', 'med', 'low']
    population_models = ['uniform', 'expgrowthslow', 'expgrowthfast', 'bottleneck']

    root_height_df_wide['mutation_signal'] = pd.Categorical(root_height_df_wide['mutation_signal'], categories=mutsigs, ordered=True)
    root_height_df_wide['population_model'] = pd.Categorical(root_height_df_wide['population_model'], categories=population_models, ordered=True)


    root_height_df_wide['sigma_skyline'] = estimate_sigma(root_height_df_wide['skyline_height_est_hpd_lower'], root_height_df_wide['skyline_height_est_hpd_upper'])
    root_height_df_wide['sigma_constcoal'] = estimate_sigma(root_height_df_wide['constcoal_height_est_hpd_lower'], root_height_df_wide['constcoal_height_est_hpd_upper'])

    root_height_df_wide['overlap'] = compute_overlap(
        root_height_df_wide['skyline_height_est_median'],
        root_height_df_wide['sigma_skyline'],
        root_height_df_wide['constcoal_height_est_median'],
        root_height_df_wide['sigma_constcoal']
    )

    overlaps = root_height_df_wide.groupby(['population_model', 'mutation_signal'], observed=True)['overlap'].mean()
    return overlaps


def plot_ci_coverage_heatmap(root_height_df_wide, sampling, title=None,
                             save_path=None, figsize=None,
                             cbar_label="Fraction of replicates"):
    """
    For each combination of mutation signal × population model, compute the fraction
    of replicates falling into each of 4 categories:
      - both inside CI
      - only constcoal inside CI
      - only skyline inside CI
      - neither inside CI
    """
    mutsig_order    = ["low", "med", "high"]
    mutsig_labels   = {"low": "Low mut. signal", "med": "Med. mut. signal", "high": "High mut. signal"}
    pop_model_order = ["uniform", "expgrowthslow", "expgrowthfast", "bottleneck"]
    title_map       = {"expgrowthfast": "Exp-growth (fast)", "expgrowthslow": "Exp-growth (slow)",
                       "uniform": "Uniform", "bottleneck": "Bottleneck"}

    _cmap_colors = ["#b5d2f2", "#7394c2", "#397398", "#80557e", "#d991b4", "#a6444f"]
    cmap   = plt.matplotlib.colors.LinearSegmentedColormap.from_list("event_cmap", _cmap_colors)
    cmap_r = cmap.reversed()

    df = root_height_df_wide[root_height_df_wide["sampling"] == sampling].copy()
    const_col = "constcoal_height_inside_ci" if "constcoal_height_inside_ci" in df.columns else "height_inside_ci_constcoal"
    sky_col   = "skyline_height_inside_ci"   if "skyline_height_inside_ci"   in df.columns else "height_inside_ci_skyline"
    ci_const = df[const_col].astype(bool)
    ci_sky   = df[sky_col].astype(bool)
    df["both"]           = ci_const & ci_sky
    df["only constcoal"] = ci_const & ~ci_sky
    df["only skyline"]   = ~ci_const & ci_sky
    df["neither"]        = ~ci_const & ~ci_sky

    nrows, ncols = len(mutsig_order), len(pop_model_order)
    default_figsize = (2.8 * ncols + 1.4, 2.6 * nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize or default_figsize, squeeze=False)

    for i, mutsig in enumerate(mutsig_order):
        for j, pop_model in enumerate(pop_model_order):
            ax = axes[i, j]
            subset = df[(df["mutation_signal"] == mutsig) & (df["population_model"] == pop_model)]
            n = len(subset)

            matrix = pd.DataFrame(
                [[subset["only skyline"].sum(), subset["both"].sum()],
                 [subset["neither"].sum(),      subset["only constcoal"].sum()]],
                index=["Skyline ✓", "Skyline ✗"],
                columns=["Constcoal ✗", "Constcoal ✓"],
                dtype=float
            ) / max(n, 1)

            annot = (matrix * 100).map(lambda v: f"{v:.0f}%")
            sns.heatmap(matrix, ax=ax, vmin=0, vmax=1, annot=annot, fmt="",
                        annot_kws={"fontsize": 11, "fontweight": "bold", "color": "white"},
                        cmap=cmap_r, linewidths=0.5, cbar=False)
            ax.set_aspect("equal")
            if i == 0:
                ax.set_title(title_map.get(pop_model, pop_model), fontsize=12)
            if j == 0:
                ax.set_ylabel(mutsig_labels[mutsig], fontsize=12)
            else:
                ax.set_ylabel("")
            ax.set_xlabel("")
            ax.tick_params(labelsize=10)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
            ax.set_yticklabels(ax.get_yticklabels(), rotation=90, va="center")

    suptitle = title if title is not None else f"Root height CI coverage — {sampling}"
    fig.suptitle(suptitle, fontsize=14)
    plt.tight_layout(rect=[0, 0, 0.91, 0.95])
    sm      = plt.cm.ScalarMappable(cmap=cmap_r, norm=plt.Normalize(vmin=0, vmax=100))
    cbar_ax = fig.add_axes([0.91, 0.15, 0.02, 0.65])
    cbar    = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_ticks([0, 25, 50, 75, 100])
    cbar.set_label(cbar_label, fontsize=13)
    cbar.ax.tick_params(labelsize=12)

    if save_path is not None:
        fig.savefig(save_path, bbox_inches="tight")
    plt.show()


def plot_ci_coverage_heatmap_nodes(tree_metrics_combined, mode="height", title=None):
    """
    CI coverage contingency heatmap from full tree_metrics_combined.

    mode="height": node height CI coverage, tips excluded (internal nodes only)
    mode="bl":     branch length CI coverage, root excluded (non-root nodes only)

    For each mutsig × popmodel condition, shows a 2×2 heatmap of fraction of
    nodes where constcoal and skyline CIs contain the true value.
    One plot per sampling type.
    """
    if mode == "height":
        ci_col = "height_inside_ci"
        df = tree_metrics_combined[tree_metrics_combined["internal"] == True].copy()
        default_title = "Node Height CI Coverage (internal nodes)"
    elif mode == "bl":
        ci_col = "bl_inside_ci"
        df = tree_metrics_combined[tree_metrics_combined["node"] != "internal_0"].copy()
        default_title = "Branch Length CI Coverage (non-root nodes)"
    else:
        raise ValueError("mode must be 'height' or 'bl'")

    title = title or default_title
    mutsig_order = ["low", "med", "high"]
    pop_model_order = ["uniform", "expgrowthslow", "expgrowthfast", "bottleneck"]
    samplings = ["independenthomochronous", "linearconstant"]

    _cmap_colors = ["#b5d2f2", "#7394c2", "#397398", "#80557e", "#d991b4", "#a6444f"]
    cmap_r = plt.matplotlib.colors.LinearSegmentedColormap.from_list("event_cmap", _cmap_colors).reversed()

    for sampling in samplings:
        df_s = df[df["sampling"] == sampling].copy()

        prefix_const = f"constcoal_{ci_col}"
        prefix_sky   = f"skyline_{ci_col}"

        if prefix_const in df_s.columns and prefix_sky in df_s.columns:
            # Wide format (new compute_errors.py output)
            merged = df_s[["population_model", "mutation_signal", "tree_index", "node_id",
                            prefix_const, prefix_sky]].copy()
            merged = merged.rename(columns={"node_id": "node",
                                            prefix_const: "ci_constcoal",
                                            prefix_sky:   "ci_skyline"})
        else:
            # Long format (old tree_metrics_combined)
            constcoal = df_s[df_s["model"] == "constcoal"][
                ["population_model", "mutation_signal", "tree_index", "node", ci_col]
            ].rename(columns={ci_col: "ci_constcoal"})
            skyline = df_s[df_s["model"] == "skyline"][
                ["population_model", "mutation_signal", "tree_index", "node", ci_col]
            ].rename(columns={ci_col: "ci_skyline"})
            merged = constcoal.merge(skyline, on=["population_model", "mutation_signal", "tree_index", "node"])

        merged["ci_constcoal"] = merged["ci_constcoal"].astype(bool)
        merged["ci_skyline"]   = merged["ci_skyline"].astype(bool)
        merged["both"]           = merged["ci_constcoal"] & merged["ci_skyline"]
        merged["only constcoal"] = merged["ci_constcoal"] & ~merged["ci_skyline"]
        merged["only skyline"]   = ~merged["ci_constcoal"] & merged["ci_skyline"]
        merged["neither"]        = ~merged["ci_constcoal"] & ~merged["ci_skyline"]

        fig, axes = plt.subplots(len(mutsig_order), len(pop_model_order),
                                 figsize=(2.5 * len(pop_model_order), 2.5 * len(mutsig_order)))

        for i, mutsig in enumerate(mutsig_order):
            for j, pop_model in enumerate(pop_model_order):
                ax = axes[i, j]
                subset = merged[(merged["mutation_signal"] == mutsig) & (merged["population_model"] == pop_model)]
                n = max(len(subset), 1)

                matrix = pd.DataFrame(
                    [[subset["only skyline"].sum(), subset["both"].sum()],
                     [subset["neither"].sum(), subset["only constcoal"].sum()]],
                    index=["Skyline ✓", "Skyline ✗"],
                    columns=["Constcoal ✗", "Constcoal ✓"],
                    dtype=float
                ) / n

                sns.heatmap(matrix, ax=ax, vmin=0, vmax=1, annot=True, fmt=".2f",
                            cmap=cmap_r, linewidths=0.5, cbar=False)
                if i == 0: ax.set_title(pop_model, fontsize=9)
                if j == 0: ax.set_ylabel(f"{mutsig}", fontsize=9)
                else: ax.set_ylabel("")
                ax.set_xlabel("")
                ax.tick_params(labelsize=8)
                ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

        fig.suptitle(f"{title} ({sampling})", fontsize=14)
        plt.tight_layout(rect=[0, 0, 0.92, 1])
        sm = plt.cm.ScalarMappable(cmap=cmap_r, norm=plt.Normalize(vmin=0, vmax=1))
        cbar_ax = fig.add_axes([0.94, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(sm, cax=cbar_ax)
        cbar.set_label("Fraction of nodes", fontsize=9)
        cbar.ax.tick_params(labelsize=8)
        plt.show()


### PLOTTING

def plot_population_summary_ax(ax, skyline_all_times, skyline_all_medians, constant_all_estimates,
                               true_traj,
                               color_exp="#a6444f", color_sky="#80557e", color_const="#397398",
                               time_horizon=0, y_range = None, mode = 'both'):
    # Compute root height
    root_height = max([max(times) for times in skyline_all_times])
    t_max = root_height if time_horizon == 0 else min(time_horizon, root_height)
    t_vals = np.linspace(0, t_max, 1000)
    N_true = np.array([true_traj(t) for t in t_vals])

    # Interpolate all skyline estimates onto a common grid
    interpolated = []
    for times, medians in zip(skyline_all_times, skyline_all_medians):
        step_times = [0.0] + times[:-1]
        step_vals = np.zeros_like(t_vals)
        for start, end, val in zip(step_times, times, medians):
            mask = (t_vals >= start) & (t_vals < end)
            step_vals[mask] = val
        step_vals[t_vals >= times[-1]] = medians[-1]
        interpolated.append(step_vals)

    interpolated = np.array(interpolated)
    skyline_median = np.median(interpolated, axis=0)
    skyline_lower = np.percentile(interpolated, 2.5, axis=0)
    skyline_upper = np.percentile(interpolated, 97.5, axis=0)

    # Plot
    ax.plot(t_vals, N_true, color=color_exp, linewidth=2, label="True Population Size", zorder = 6)

    const_median = np.median(constant_all_estimates)
    const_lower = np.percentile(constant_all_estimates, 2.5)
    const_upper = np.percentile(constant_all_estimates, 97.5)

    if mode == 'both':
        ax.hlines(const_median, 0, t_max, color=color_const, linestyle=':', label="Constant Pop. Estimate", linewidth = 2, zorder = 5)
        ax.fill_between(t_vals, const_lower, const_upper, color=color_const, alpha=0.3, zorder = 4)

        ax.plot(t_vals, skyline_median, color=color_sky, linestyle='--', label="Skyline Estimate", linewidth = 2, zorder = 5)
        ax.fill_between(t_vals, skyline_lower, skyline_upper, color=color_sky, alpha=0.3, zorder = 4)
    elif mode == 'skyline':
        ax.plot(t_vals, skyline_median, color=color_sky, linestyle='--', label="Skyline Estimate", linewidth = 2, zorder = 5)
        ax.fill_between(t_vals, skyline_lower, skyline_upper, color=color_sky, alpha=0.3, zorder = 4)
    elif mode == 'constcoal':
        ax.hlines(const_median, 0, t_max, color=color_const, linestyle=':', label="Constant Pop. Estimate", linewidth = 2, zorder = 5)
        ax.fill_between(t_vals, const_lower, const_upper, color=color_const, alpha=0.3, zorder = 4)
    
    if y_range:
        ax.set_ylim(y_range)
    
    ax.set_xlabel("Time before present")
    ax.set_ylabel("Population Size")
    ax.invert_xaxis()

def plot_population_summary(path_info_df, sampling, x_range=None, title="", y_range=None, add_samples=False, mode='both', pop_models=None):
    """
    Summary plot per condition (mutation signal × population model) for a given sampling type,
    showing median and 95% CI for skyline and constant-coalescent estimates.
    """
    if pop_models is None:
        pop_models = sorted(path_info_df["population_model"].unique())
    mut_signals = ["low", "med", "high"]
    ncols, nrows = len(pop_models), len(mut_signals)

    df = path_info_df[path_info_df["sampling"] == sampling]
    time_horizon = x_range[1] if x_range else 0

    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.5 * nrows))

    if nrows == 1:
        axes = np.expand_dims(axes, axis=0)
    if ncols == 1:
        axes = np.expand_dims(axes, axis=1)

    for i, mut_sig in enumerate(mut_signals):
        for j, pop_model in enumerate(pop_models):
            ax = axes[i][j]
            subset = df[
                (df["mutation_signal"] == mut_sig) &
                (df["population_model"] == pop_model)
            ]

            if subset.empty:
                continue

            true_traj = get_true_traj(pop_model, subset.iloc[0])
            skyline_all_times = subset["skyline_times"].tolist()
            skyline_all_medians = subset["skyline_medians"].tolist()
            constant_all_estimates = subset["coalescent_median"].tolist()

            if add_samples:
                skyline_sample_dfs   = subset["skyline_samples"]
                skyline_times_sample_dfs = subset["skyline_times_samples"] if "skyline_times_samples" in subset.columns else None
                constant_samples = subset["coalescent_samples"]

                for k, samples_df in enumerate(skyline_sample_dfs):
                    times_df = skyline_times_sample_dfs.iloc[k] if skyline_times_sample_dfs is not None else None
                    for row_index, (l, row) in enumerate(samples_df.iterrows()):
                        sky_times = list(times_df.iloc[row_index]) if times_df is not None else skyline_all_times[k]
                        plot_population_trajectories_ax(
                            ax,
                            skyline_times=sky_times,
                            skyline_medians=row,
                            constant_pop_estimate=constant_samples.iloc[k].iloc[row_index],
                            true_traj=true_traj,
                            alpha=0.01,
                            first_plot=False,
                            time_horizon=time_horizon,
                            color_const='lightgray',
                            color_sky='lightgray',
                            mode=mode,
                        )

            plot_population_summary_ax(
                ax,
                skyline_all_times=skyline_all_times,
                skyline_all_medians=skyline_all_medians,
                constant_all_estimates=constant_all_estimates,
                true_traj=true_traj,
                time_horizon=time_horizon,
                y_range=y_range,
                mode=mode,
            )

            if x_range:
                ax.set_xlim(x_range[1], x_range[0])  # inverted x-axis

            if j == 0:
                ax.set_ylabel(f"{mut_sig.capitalize()} mut. signal\nPopulation Size")
            if i == 0:
                ax.set_title(pop_model, fontsize=12)

    handles, labels = axes[0][0].get_legend_handles_labels()
    if add_samples:
        if mode == 'skyline':
            grey_line = Line2D([0], [0], color='lightgrey', linestyle='--', label='Sample Skyline Trajectories')
            handles.append(grey_line)
            labels.append('Sample Skyline Trajectories')
        elif mode == 'constcoal':
            grey_line = Line2D([0], [0], color='lightgrey', linestyle=':', label='Sample Constant Trajectories')
            handles.append(grey_line)
            labels.append('Sample Constant Trajectories')
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.02), ncol=3)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.suptitle(f"{title} ({sampling})")
    plt.show()


def plot_skygrid_summary(skygrid_df, sampling, true_traj_fn=None, x_range=None, y_range=None,
                         title="", pop_models=None):
    """
    Plot skygrid population size estimates for each scenario (mut signal × pop model),
    one line per tree replicate, plus the true trajectory.

    Parameters:
        skygrid_df: DataFrame from process_skygrid_logs()
        sampling: "independenthomochronous" or "linearconstant"
        true_traj_fn: optional dict {pop_model: callable(t)} overriding get_true_traj
        x_range: (min, max) time axis — x-axis is inverted (time before present)
        y_range: (min, max) y axis
        title: plot title
        pop_models: list of pop models to show; inferred from data if None
    """
    df = skygrid_df[skygrid_df["sampling"] == sampling]
    if pop_models is None:
        pop_models = sorted(df["population_model"].unique())
    mut_signals = ["low", "med", "high"]

    ncols, nrows = len(pop_models), len(mut_signals)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.5 * nrows))
    if nrows == 1:
        axes = np.expand_dims(axes, axis=0)
    if ncols == 1:
        axes = np.expand_dims(axes, axis=1)

    color_skygrid = "#4e9a8a"
    color_true    = "#a6444f"

    for i, mut_sig in enumerate(mut_signals):
        for j, pop_model in enumerate(pop_models):
            ax = axes[i][j]
            subset = df[(df["mutation_signal"] == mut_sig) & (df["population_model"] == pop_model)]

            if subset.empty:
                ax.set_visible(False)
                continue

            # grid times are the same for all rows; read from first row
            row0 = subset.iloc[0]
            time_cols = sorted([c for c in subset.columns if re.match(r"time_\d+$", c)],
                               key=lambda c: int(c.split("_")[-1]))
            times = row0[time_cols].values.astype(float)

            # one line per tree replicate
            for _, row in subset.iterrows():
                median_cols = sorted([c for c in subset.columns if re.match(r"median_\d+$", c)],
                                     key=lambda c: int(c.split("_")[-1]))
                lower_cols  = sorted([c for c in subset.columns if c.startswith("hpd_lower_")],
                                     key=lambda c: int(c.split("_")[-1]))
                upper_cols  = sorted([c for c in subset.columns if c.startswith("hpd_upper_")],
                                     key=lambda c: int(c.split("_")[-1]))

                medians = row[median_cols].values.astype(float)
                lowers  = row[lower_cols].values.astype(float)
                uppers  = row[upper_cols].values.astype(float)

                # append root height as the right boundary of the last interval
                root_height = row["median_root_height"]
                plot_times = np.append(times, root_height)
                plot_medians = np.append(medians, medians[-1])
                plot_lowers  = np.append(lowers,  lowers[-1])
                plot_uppers  = np.append(uppers,  uppers[-1])

                ax.step(plot_times, plot_medians, where="post", color=color_skygrid, alpha=0.6, linewidth=1)
                ax.fill_between(plot_times, plot_lowers, plot_uppers, step="post", color=color_skygrid, alpha=0.1)

            # true trajectory
            t_max = x_range[1] if x_range else times[-1]
            t_vals = np.linspace(0, t_max, 500)
            if true_traj_fn and pop_model in true_traj_fn:
                true_traj = true_traj_fn[pop_model]
            else:
                true_traj = get_true_traj(pop_model, assign_model_params(pop_model))
            ax.plot(t_vals, [true_traj(t) for t in t_vals],
                    color=color_true, linewidth=2, label="True", zorder=6)

            ax.invert_xaxis()
            if x_range:
                ax.set_xlim(x_range[1], x_range[0])
            if y_range:
                ax.set_ylim(y_range)
            if j == 0:
                ax.set_ylabel(f"{mut_sig.capitalize()} mut. signal\nPopulation Size")
            if i == 0:
                ax.set_title(pop_model, fontsize=12)
            if i == nrows - 1:
                ax.set_xlabel("Time before present")

    legend_handles = [
        Line2D([0], [0], color=color_true,    linewidth=2,            label="True population size"),
        Line2D([0], [0], color=color_skygrid, linewidth=1, alpha=0.6, label="Skygrid median (per tree)"),
    ]
    fig.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, 1.02), ncol=2)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.suptitle(f"{title} ({sampling})")
    plt.show()


def plot_population_trajectories_ax(ax, skyline_times, skyline_medians, constant_pop_estimate,
                                     true_traj,
                                     color_exp="#a6444f", color_sky="#80557e", color_const="#397398",
                                     alpha=0.4, label_prefix="", first_plot=True, time_horizon=0, y_range = (0, 6000), mode = 'both'):

    root_height = skyline_times[-1]
    t_max = root_height if time_horizon == 0 else min(time_horizon, root_height)
    t_vals = np.linspace(0, t_max, 1000)
    N_true = np.array([true_traj(t) for t in t_vals])

    # Stepwise skyline
    skyline_start_times = [0.0] + skyline_times[:-1]
    step_times, step_values = [], []
    for start, end, value in zip(skyline_start_times, skyline_times, skyline_medians):
        if start > t_max:
            break
        end = min(end, t_max)
        step_times.extend([start, end])
        step_values.extend([value, value])
    
    if first_plot:
        ax.plot(t_vals, N_true, color=color_exp, alpha=1, linewidth= 2, label="True Population Size", zorder = 4)
        ax.hlines(constant_pop_estimate, 0, t_max, color=color_const, linestyle=':', 
              alpha=alpha, label = "Constant Pop. Estimate")
        ax.plot(step_times, step_values, drawstyle='steps-post', linestyle='--', color=color_sky, 
            alpha=alpha, label="Skyline Estimate")
    else:
        if mode == 'both':
            ax.hlines(constant_pop_estimate, 0, t_max, color=color_const, linestyle=':', 
                alpha=alpha)
            ax.plot(step_times, step_values, drawstyle='steps-post', linestyle='--', color=color_sky, 
                alpha=alpha)
        elif mode == 'skyline':
            ax.plot(step_times, step_values, drawstyle='steps-post', linestyle='--', color='lightgrey', 
                alpha=alpha)
        elif mode == 'constcoal':
            ax.hlines(constant_pop_estimate, 0, t_max, color='lightgrey', linestyle=':', 
                alpha=alpha)

    if first_plot:
        if y_range: 
            ax.set_ylim(y_range)
        ax.set_xlabel("Time before present")
        ax.set_ylabel("Population Size")
        ax.invert_xaxis()



def plot_summary_population_grid(path_info_df, sampling, num_groups=10, burnin=0.1, time_horizon=0, title="", y_range=(0, 6000)):
    """
    Plots all population trajectories for all trees in a subplot grid by
    population model (columns) and mutation signal (rows).
    """
    pop_models = ["uniform", "expgrowthfast", "expgrowthslow", "bottleneck"]
    mut_signals = ["low", "med", "high"]
    ncols, nrows = len(pop_models), len(mut_signals)

    df = path_info_df[path_info_df["sampling"] == sampling]

    fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 3.5*nrows))

    if nrows == 1:
        axes = np.expand_dims(axes, axis=0)
    if ncols == 1:
        axes = np.expand_dims(axes, axis=1)

    first_plot_per_ax = {(i, j): True for i in range(nrows) for j in range(ncols)}

    for idx, row in df.iterrows():
        pop_model = row["population_model"]
        mut_sig = row["mutation_signal"]
        tree_index = row["tree_index"]

        if pop_model not in pop_models or mut_sig not in mut_signals:
            continue

        row_idx = mut_signals.index(mut_sig)
        col_idx = pop_models.index(pop_model)
        ax = axes[row_idx][col_idx]

        true_traj = get_true_traj(pop_model, row)

        plot_population_trajectories_ax(
            ax,
            skyline_times=row["skyline_times"],
            skyline_medians=row["skyline_medians"],
            constant_pop_estimate=row["coalescent_median"],
            true_traj=true_traj,
            alpha=0.5,
            first_plot=first_plot_per_ax[(row_idx, col_idx)],
            time_horizon=time_horizon,
            y_range=y_range,
        )
        first_plot_per_ax[(row_idx, col_idx)] = False

    for i, mut_sig in enumerate(mut_signals):
        axes[i][0].set_ylabel(f"{mut_sig.capitalize()} mut.\nsignal\nPopulation Size")

    for j, pop_model in enumerate(pop_models):
        axes[0][j].set_title(pop_model, fontsize=12)

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.02), ncol=3)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.suptitle(f"{title} ({sampling})")
    plt.show()

def boxplot_branch_length_errors(df_complete, sampling, metric="bl_abs_relative_error", title="", logscale=True, plot_type='box', y_range=None):

    df = df_complete[df_complete["sampling"] == sampling].copy()

    if metric.startswith("bl"):
        df = df.dropna(subset=[metric])
    elif metric.startswith("height"):
        df = df[df.internal == True]

    hue_order = ["constcoal", "skyline"]
    palette = ["#a6444f", "#397398"]

    mutsig_order = ["low", "med", "high"]
    growth_model_order = ["uniform", "expgrowthfast", "expgrowthslow", "bottleneck"]
    condition_order = [f"{g}/{m}" for g in growth_model_order for m in mutsig_order]

    df["condition"] = df["population_model"] + "/" + df["mutation_signal"]

    plt.figure(figsize=(16, 6))

    if plot_type == 'box':
        ax = sns.boxplot(
            data=df, x="condition", y=metric, hue="model",
            palette=palette, hue_order=hue_order,
            showfliers=False, order=condition_order,
        )
    elif plot_type == 'violin':
        ax = sns.violinplot(
            data=df, x="condition", y=metric, hue="model",
            palette=palette, hue_order=hue_order,
            order=condition_order, cut=2, inner="box", density_norm='width'
        )

    sns.stripplot(
        data=df, x="condition", y=metric, hue="model",
        palette=palette, hue_order=hue_order, order=condition_order,
        dodge=True, alpha=0.2, size=2, jitter=True, linewidth=0
    )

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:len(hue_order)], labels[:len(hue_order)],
              title="Model", loc='center left', bbox_to_anchor=(1.01, 0.5),
              frameon=False, fontsize=12, title_fontsize=13)

    if logscale:
        ax.set_yscale("log")
        ax.yaxis.set_major_locator(plt.LogLocator(base=10, numticks=15))
        ax.yaxis.set_minor_locator(plt.NullLocator())
    if y_range:
        ax.set_ylim(y_range)
    ax.set_ylabel(f"{metric}", fontsize=12)
    ax.set_xlabel("")
    ax.set_xticks(range(len(condition_order)))
    ax.set_xticklabels([cond.split("/")[1] for cond in condition_order], rotation=0, fontsize=12)

    for i, growth_model in enumerate(growth_model_order):
        mid = i * 3 + 1
        ax.text(mid, -0.08, growth_model, ha='center', va='top', fontsize=12,
                transform=ax.get_xaxis_transform())

    plt.title(f"{title} ({sampling})")
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout(rect=[0, 0.03, 0.95, 1])
    plt.show()

def plot_tree_comparison(branch_length_df):

    # Drop rows with missing or non-positive values (log scale cannot handle <= 0)
    df_log = branch_length_df.dropna(subset=["bl_ci_lower_constcoal", "bl_ci_upper_constcoal"])
    df_log = df_log[(df_log["bl_sim"] > 0) & (df_log["bl_constcoal"] > 0)]

    # Compute asymmetric error bars for y-axis (const-coal tree)
    yerr_lower = df_log["bl_constcoal"] - df_log["bl_ci_lower_constcoal"]
    yerr_upper = df_log["bl_ci_upper_constcoal"] - df_log["bl_constcoal"]
    yerr = np.array([yerr_lower, yerr_upper])

    # Determine unified log-space limits
    log_min = df_log[["bl_sim", "bl_constcoal"]].min().min()
    log_max = df_log[["bl_sim", "bl_constcoal"]].max().max()
    margin = 0.1  # for padding
    log_min = log_min * (1 - margin)
    log_max = log_max * (1 + margin)

    # Plot
    plt.figure(figsize=(8, 6))
    plt.errorbar(
    df_log["bl_sim"],
    df_log["bl_constcoal"],
    yerr=yerr,
    fmt='o',
    ecolor='gray',
    color='steelblue',
    markersize=4,
    elinewidth=0.6,
    capsize=3,
    alpha=1,
    label="Branch lengths"
    )

    # Add diagonal x = y line
    x_vals = np.logspace(np.log10(log_min), np.log10(log_max), 100)
    plt.plot(x_vals, x_vals, 'k--')

    # Set log scale
    plt.xscale("log")
    plt.yscale("log")

    # Set same limits and equal aspect ratio
    plt.xlim(log_min, log_max)
    plt.ylim(log_min, log_max)
    plt.gca().set_aspect('equal', adjustable='box')

    # Labels and formatting
    plt.xlabel("Branch length [years] (simulated tree, log scale)")
    plt.ylabel("Branch length [years] (const-coalescent tree, log scale)")
    plt.title("Log-Scaled Branch Length Comparison with Credible Intervals")
    plt.legend()
    plt.grid(True, which='both', linestyle=':', linewidth=0.5)
    plt.show()


def plot_error_by_mutsig(root_height_df_wide, sampling,
                         y_col="height_rel_error_median",
                         y_label=None, y_range=None,
                         pop_models=None, figsize=None,
                         title=None, save_path=None,
                         violin=False):
    """
    Four subplots side by side (one per population model). For each mutation
    signal (low, med, high) two jittered scatter columns are shown — constcoal
    (blue) and skyline (purple) — with individual replicate values in grey and
    the per-group median overlaid as a larger coloured dot.

    y_col is the base column name without the model prefix, e.g.
    'height_rel_error_median'. The function looks up
    '{model}_{y_col}' in root_height_df_wide.
    """
    colors      = {"constcoal": "#397398", "skyline": "#80557e"}
    labels      = {"constcoal": "Const. coalescent", "skyline": "Skyline"}
    mutsig_order    = ["low", "med", "high"]
    mutsig_labels   = {"low": "Low", "med": "Med.", "high": "High"}
    title_map       = {"expgrowthfast": "Exp-growth (fast)", "expgrowthslow": "Exp-growth (slow)",
                       "uniform": "Uniform", "bottleneck": "Bottleneck"}
    pop_model_order = pop_models or ["uniform", "expgrowthslow", "expgrowthfast", "bottleneck"]

    df = root_height_df_wide[root_height_df_wide["sampling"] == sampling].copy()

    # x positions: each mutsig gets a slot of width 1; constcoal left, skyline right
    offsets = {"constcoal": -0.15, "skyline": 0.15}
    x_base  = {ms: i for i, ms in enumerate(mutsig_order)}

    ncols = len(pop_model_order)
    fig, axes = plt.subplots(1, ncols, figsize=figsize or (4.5 * ncols, 4.5), sharey=True)
    if ncols == 1:
        axes = [axes]

    rng = np.random.default_rng(42)

    for ax, pop_model in zip(axes, pop_model_order):
        for model in ("constcoal", "skyline"):
            col = f"{model}_{y_col}"
            if col not in df.columns:
                continue
            color = colors[model]
            for ms in mutsig_order:
                sub = df[(df["population_model"] == pop_model) & (df["mutation_signal"] == ms)].dropna(subset=[col])
                x_center = x_base[ms] + offsets[model]
                jitter = rng.uniform(-0.06, 0.06, size=len(sub))
                xs = x_center + jitter
                ax.scatter(xs, sub[col].values, color="#aaaaaa", s=15,
                           alpha=0.8, zorder=2, linewidths=0)
                vals = sub[col].dropna().values
                if len(vals) >= 2:
                    hpd_lo, hpd_hi = hpd(np.sort(vals))
                    ax.plot([x_center, x_center], [hpd_lo, hpd_hi],
                            color=color, linewidth=1.5, zorder=3, solid_capstyle="round")
                    cap = 0.05
                    ax.plot([x_center - cap, x_center + cap], [hpd_lo, hpd_lo],
                            color=color, linewidth=1.5, zorder=3)
                    ax.plot([x_center - cap, x_center + cap], [hpd_hi, hpd_hi],
                            color=color, linewidth=1.5, zorder=3)
                    if violin:
                        kde = gaussian_kde(vals, bw_method=0.4)
                        ys = np.linspace(hpd_lo, hpd_hi, 200)
                        xs_kde = kde(ys)
                        xs_kde = xs_kde / xs_kde.max() * 0.12
                        ax.fill_betweenx(ys, x_center - xs_kde, x_center + xs_kde,
                                         color=color, alpha=0.35, zorder=3, linewidth=0)
                ax.scatter(x_center, sub[col].median(), color=color, s=75,
                           zorder=4, edgecolors="white", linewidths=0.5)

        ax.set_xticks(range(len(mutsig_order)))
        ax.set_xticklabels([mutsig_labels[ms] for ms in mutsig_order], fontsize=14)
        ax.set_xlabel("Mutation signal", fontsize=14)
        ax.set_title(title_map.get(pop_model, pop_model), fontsize=15)
        ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, zorder=1)
        ax.grid(True, color="lightgray", linewidth=0.4, alpha=0.5, zorder=0)
        ax.tick_params(labelsize=14, labelleft=True)
        ax.set_ylabel(y_label or y_col, fontsize=15, labelpad=2)
        if y_range:
            ax.set_ylim(y_range)

    handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=c,
                          markersize=11, label=labels[m])
               for m, c in colors.items()]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.99, 0.97),
               ncol=2, frameon=False, fontsize=14)

    suptitle = title if title is not None else sampling
    plt.suptitle(suptitle, fontsize=17)
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.subplots_adjust(wspace=0.35)

    if save_path is not None:
        fig.savefig(save_path, bbox_inches="tight")
    plt.show()


def plot_height_vs_popsize_error(root_height_df_wide, sampling,
                                  x_col="pop_diff_median", y_col="height_diff_median",
                                  x_range=None, y_range=None, alpha=0.6, s=10,
                                  title=None, x_label=None, y_label=None,
                                  pop_models=None, mutsigs=None,
                                  figsize=None, show_zero_lines=True,
                                  legend_inside=True,
                                  save_path=None):
    """
    Scatterplot of y_col vs x_col for constcoal and skyline, one subplot per
    population model (columns) × mutation signal (rows), filtered by sampling type.
    Draws per-replicate HPD intervals as error bars.
    Column names are expected with model prefix: {model}_{col}_median/hpd_lower/hpd_upper.

    Parameters
    ----------
    save_path : str or Path, optional
        If provided, save figure to this path (PDF/PNG inferred from extension)
        instead of displaying it.
    x_label, y_label : str, optional
        Axis labels; default to the column name.
    title : str, optional
        Suptitle; default to sampling type.
    show_zero_lines : bool
        Draw dashed zero reference lines on each subplot.
    figsize : tuple, optional
        Override default figure size.
    """
    df = root_height_df_wide[root_height_df_wide["sampling"] == sampling].copy()

    mutsig_order       = mutsigs    if mutsigs    is not None else ["low", "med", "high"]
    growth_model_order = pop_models if pop_models is not None else ["uniform", "expgrowthslow", "expgrowthfast", "bottleneck"]

    mutsig_labels = {"low": "Low mut. signal", "med": "Med. mut. signal", "high": "High mut. signal"}
    title_map     = {"expgrowthfast": "Exp-growth (fast)", "expgrowthslow": "Exp-growth (slow)",
                     "uniform": "Uniform", "bottleneck": "Bottleneck"}
    legend_labels = {"constcoal": "Const. coalescent", "skyline": "Skyline"}
    colors        = {"constcoal": "#397398", "skyline": "#80557e"}
    markers       = {"constcoal": "o",       "skyline": "^"}

    nrows, ncols = len(mutsig_order), len(growth_model_order)
    default_figsize = (4 * ncols, 3.2 * nrows)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                             figsize=figsize or default_figsize,
                             squeeze=False)

    def _col(model, base):
        if f"{model}_{base}" in df.columns:
            return f"{model}_{base}"
        if f"{base}_{model}" in df.columns:
            return f"{base}_{model}"
        return None

    def _lo(base): return base.replace("_median", "_hpd_lower")
    def _hi(base): return base.replace("_median", "_hpd_upper")

    for i, mutsig in enumerate(mutsig_order):
        for j, pop_model in enumerate(growth_model_order):
            ax = axes[i, j]
            subset = df[(df["mutation_signal"] == mutsig) & (df["population_model"] == pop_model)]

            if subset.empty:
                ax.axis("off")
                continue

            for model, color in colors.items():
                xm_key = _col(model, x_col)
                ym_key = _col(model, y_col)
                if xm_key is None or ym_key is None:
                    continue

                xm = subset[xm_key]
                ym = subset[ym_key]

                xl_key = _col(model, _lo(x_col))
                xh_key = _col(model, _hi(x_col))
                yl_key = _col(model, _lo(y_col))
                yh_key = _col(model, _hi(y_col))

                has_x_ci = xl_key is not None and xh_key is not None
                has_y_ci = yl_key is not None and yh_key is not None

                if has_x_ci or has_y_ci:
                    xerr = ([xm - subset[xl_key], subset[xh_key] - xm] if has_x_ci else None)
                    yerr = ([ym - subset[yl_key], subset[yh_key] - ym] if has_y_ci else None)
                    eb = ax.errorbar(xm, ym,
                                     xerr=xerr, yerr=yerr,
                                     fmt=markers[model], color=color, ecolor=color,
                                     elinewidth=0.8, capsize=2, markersize=4,
                                     label=legend_labels[model])
                    eb[0].set_alpha(alpha)
                    for line in eb[1]: line.set_alpha(0.15)
                    for line in eb[2]: line.set_alpha(0.15)
                else:
                    ax.scatter(xm, ym, color=color, marker=markers[model],
                               alpha=alpha, s=s, label=legend_labels[model])

            ax.grid(True, color="lightgray", linewidth=0.4, alpha=0.5, zorder=0)
            if show_zero_lines:
                ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, zorder=1)
                ax.axvline(0, color="gray", linestyle="--", linewidth=0.8, zorder=1)
            if x_range: ax.set_xlim(x_range)
            if y_range: ax.set_ylim(y_range)
            ax.tick_params(labelsize=13)

            ax.set_xlabel(x_label or x_col, fontsize=14)

            if i == 0:
                ax.set_title(title_map.get(pop_model, pop_model), fontsize=15)
            ax.set_ylabel(y_label or y_col, fontsize=14, labelpad=2)
            if legend_inside:
                leg = ax.legend(loc="upper right", frameon=True, fontsize=11,
                                markerscale=1.2, handletextpad=0.4, borderpad=0.4)
                leg.get_frame().set_linewidth(0.6)
                leg.get_frame().set_edgecolor("gray")

    if not legend_inside:
        handles = []
        for m, c in colors.items():
            marker_line = Line2D([0], [0], marker=markers[m], color=c, markersize=9,
                                 linestyle="none", alpha=1.0)
            caplines    = (Line2D([0], [0], color=c, linewidth=1.0),
                           Line2D([0], [0], color=c, linewidth=1.0))
            barlines    = (LineCollection([[[0, -0.3], [0, 0.3]]], colors=[c], linewidths=[0.8]),
                           LineCollection([[[-0.3, 0], [0.3, 0]]], colors=[c], linewidths=[0.8]))
            handle      = ErrorbarContainer((marker_line, caplines, barlines),
                                            has_xerr=True, has_yerr=True,
                                            label=legend_labels[m])
            handles.append(handle)
        fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.99, 0.97),
                   ncol=2, frameon=False, fontsize=15, borderpad=0.8, handletextpad=0.6)

    suptitle = title if title is not None else sampling
    plt.suptitle(suptitle, fontsize=17)
    plt.tight_layout(rect=[0, 0, 1, 0.94 if not legend_inside else 0.97])

    # Place mutsig labels just left of the ylabel, anchored to its rendered position.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    for j, mutsig in enumerate(mutsig_order):
        ax = axes[j][0]
        ylabel_bb = ax.yaxis.label.get_window_extent(renderer)
        x_fig = inv.transform((ylabel_bb.x0, ylabel_bb.y0 + ylabel_bb.height / 2))[0]
        pos = ax.get_position()
        y_fig = (pos.y0 + pos.y1) / 2
        fig.text(x_fig - 0.01, y_fig, mutsig_labels.get(mutsig, mutsig),
                 fontsize=15, ha="center", va="center", rotation=90,
                 transform=fig.transFigure)

    if save_path is not None:
        fig.savefig(save_path, bbox_inches="tight")
    plt.show()

def plot_height_errors_by_time_bin(ax, df_sub, t_max, bins=10, bin_width=None, error_col="height_abs_relative_error", x_col = 'height_sim', y_range = None, plot_type = 'box', add_legend = True, skyline_col = "#397398", constcoal_col =  "#a6444f"):
    df_sub = df_sub.copy()

    if x_col == 'height_sim':
        if bin_width is not None:
            # Fixed-width bins up to t_max, last bin extends to at least t_max
            bin_edges = list(range(0, t_max, bin_width))
            open_end = max(df_sub["height_sim"].max() + 1e-6, t_max + 1e-6)
            bin_edges.append(open_end)
        else:
            bin_edges = list(np.linspace(0, t_max, bins + 1))
            max_height = df_sub["height_sim"].max()
            if max_height > t_max:
                bin_edges.append(max_height + 1e-6)
    else:
        min_val = df_sub[x_col].min()
        max_val = df_sub[x_col].max()
        bin_edges = list(np.linspace(min_val, max_val, bins + 1))
    
    df_sub["bin"] = pd.cut(df_sub[x_col], bins=bin_edges, include_lowest=True)

    # Drop NaNs
    df_sub = df_sub.dropna(subset=["bin", error_col])

    # Plot
    if plot_type == 'box':
        sns.boxplot(
            data=df_sub,
            x="bin",
            y=error_col,
            hue="model",
            ax=ax,
            showfliers=False,
            palette=[constcoal_col, skyline_col],
            alpha = 0.7
        )
    elif plot_type == 'violin':
        sns.violinplot(
            data=df_sub,
            x="bin",
            y=error_col,
            hue="model",
            ax=ax,
            palette=[constcoal_col, skyline_col],
            cut=0,
            bw_adjust=0.5,
            density_norm='width', # area such that each box has the same area, width for same width
            inner="box",
            alpha = 0.7,
        )
        

    ax.axhline(0, color='gray', linestyle='--', linewidth=1)
    
    # Formatting
    #ax.set_ylim(-5, 5)
    if x_col == 'height_sim':
        ax.set_xlabel("Time before present")
    else:
        ax.set_xlabel(x_col)
    tick_labels = [f"{int(interval.left)}-{int(interval.right)}" for interval in df_sub["bin"].cat.categories]
    ax.set_xticks(range(len(tick_labels)))
    ax.set_xticklabels(tick_labels, rotation=45)
    if y_range:
        ax.set_ylim(y_range)
    if add_legend:
        ax.legend(title="Model", loc="upper right")
    else:
        legend = ax.get_legend()
        if legend is not None:
            legend.remove()



def plot_population_summary_ax_95cf(ax, skyline_all_times, skyline_all_medians, constant_all_estimates,
                               skyline_all_lowers, skyline_all_uppers, constant_all_lowers, constant_all_uppers,
                               true_traj,
                               color_exp="#a6444f", color_sky="#80557e", color_const="#397398",
                               time_horizon=0, mode = 'skyline', plot_true_size = True):
    # Compute root height
    root_height = max([max(times) for times in skyline_all_times])
    t_max = root_height if time_horizon == 0 else min(time_horizon, root_height)
    t_vals = np.linspace(0, t_max, 1000)
    N_true = np.array([true_traj(t) for t in t_vals])

    if plot_true_size:
        ax.plot(t_vals, N_true, color=color_exp, linewidth=2, label="True Population Size")

    # Interpolate all skyline estimates onto a common grid
    def plot_population_all_trees(skyline_all_values, constant_all_values, color_sky, color_const, skyline_label, constant_label):
        interpolated = []
        for times, medians in zip(skyline_all_times, skyline_all_values):
            step_times = [0.0] + times[:-1]
            step_vals = np.zeros_like(t_vals)
            for start, end, val in zip(step_times, times, medians):
                mask = (t_vals >= start) & (t_vals < end)
                step_vals[mask] = val
            step_vals[t_vals >= times[-1]] = medians[-1]
            interpolated.append(step_vals)

        interpolated = np.array(interpolated)
        skyline_median = np.median(interpolated, axis=0)
        skyline_lower = np.percentile(interpolated, 2.5, axis=0)
        skyline_upper = np.percentile(interpolated, 97.5, axis=0)

        # Plot
        const_median = np.median(constant_all_values)
        const_lower = np.percentile(constant_all_values, 2.5)
        const_upper = np.percentile(constant_all_values, 97.5)
        if mode == 'constcoal':
            ax.hlines(const_median, 0, t_max, color=color_const, linestyle=':', label=constant_label)
            ax.fill_between(t_vals, const_lower, const_upper, color=color_const, alpha=0.2)
        elif mode == 'skyline':
            ax.plot(t_vals, skyline_median, color=color_sky, linestyle='--', label=skyline_label)
            ax.fill_between(t_vals, skyline_lower, skyline_upper, color=color_sky, alpha=0.2)

    plot_population_all_trees(skyline_all_medians, constant_all_estimates, color_sky, color_const, "Skyline Estimate", "Constant Estimate")
    plot_population_all_trees(skyline_all_lowers, constant_all_lowers, "#d991b4", "#7394c2", "Skyline Estimate 95% CI", "Constant Estimate 95% CI")
    plot_population_all_trees(skyline_all_uppers, constant_all_uppers, "#d991b4", "#7394c2", None, None)

    # Repeat for upper confidence interval values
    
    ax.set_xlabel("Time before present")
    ax.set_ylabel("Population Size")
    ax.invert_xaxis()


def get_true_traj(pop_model, row):
    """Return true population trajectory function for a given pop model."""
    if pop_model == "expgrowthfast":
        return lambda t: row["present_pop_size"] * np.exp(-row["growth_rate"] * t)
    elif pop_model == "expgrowthslow":
        return lambda t: row["present_pop_size"] * np.exp(-row["growth_rate"] * t)
    elif pop_model.startswith("bottleneck"):
        ps, bs, bstart, bend = row["present_pop_size"], row["bottleneck_size"], row["bottleneck_start"], row["bottleneck_end"]
        def bottleneck_traj(t):
            if t <= bstart:
                return ps
            elif t < bend:
                return bs
            else:
                return ps
        return bottleneck_traj
    else:  # uniform
        return lambda t: row["present_pop_size"]


def plot_population_summary_95cf(path_info_df, sampling, time_horizon=0, title="", mode='skyline', plot_true_size=True):
    """
    Summary plot per condition (mutation signal × population model) for a given sampling type,
    showing median and 95% CI for skyline and constant-coalescent estimates.
    """
    pop_models = ["uniform", "expgrowthfast", "expgrowthslow", "bottleneck"]
    mut_signals = ["low", "med", "high"]
    ncols, nrows = len(pop_models), len(mut_signals)

    df = path_info_df[path_info_df["sampling"] == sampling]

    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.5 * nrows))

    if nrows == 1:
        axes = np.expand_dims(axes, axis=0)
    if ncols == 1:
        axes = np.expand_dims(axes, axis=1)

    for i, mut_sig in enumerate(mut_signals):
        for j, pop_model in enumerate(pop_models):
            ax = axes[i][j]
            subset = df[
                (df["mutation_signal"] == mut_sig) &
                (df["population_model"] == pop_model)
            ]

            if subset.empty:
                continue

            true_traj = get_true_traj(pop_model, subset.iloc[0])

            plot_population_summary_ax_95cf(
                ax,
                skyline_all_times=subset["skyline_times"].tolist(),
                skyline_all_medians=subset["skyline_medians"].tolist(),
                constant_all_estimates=subset["coalescent_median"].tolist(),
                skyline_all_lowers=subset["skyline_lowers"].tolist(),
                skyline_all_uppers=subset["skyline_uppers"].tolist(),
                constant_all_lowers=subset["coalescent_lower"].tolist(),
                constant_all_uppers=subset["coalescent_upper"].tolist(),
                true_traj=true_traj,
                time_horizon=time_horizon,
                mode=mode,
                plot_true_size=plot_true_size,
            )

            if j == 0:
                ax.set_ylabel(f"{mut_sig.capitalize()} mut.\nsignal\nPopulation Size")
            if i == 0:
                ax.set_title(pop_model, fontsize=12)

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.02), ncol=3)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.suptitle(f"{title} ({sampling})")
    plt.show()


def plot_bl_vs_popsize_error_by_time_bin(df, sampling, bin_width=50, time_col="height_sim",
                                         pop_error_col="rel_diff_pop_size", bl_error_col="bl_relative_error",
                                         x_range=None, y_range=None, show_errorbars=True, min_values_per_bin=1):
    df = df[df["sampling"] == sampling].copy()
    mutsig_order = ["low", "med", "high"]
    growth_model_order = ["uniform", "expgrowthslow", "expgrowthfast", "bottleneck"]
    hue_order = ["constcoal", "skyline"]
    colors = {"constcoal": "#a6444f", "skyline": "#397398"}

    fig, axes = plt.subplots(nrows=len(mutsig_order), ncols=len(growth_model_order), figsize=(5*len(growth_model_order), 3.5*len(mutsig_order)), sharex=False, sharey=False)

    for i, mutsig in enumerate(mutsig_order):
        for j, growth_model in enumerate(growth_model_order):
            ax = axes[i, j]
            df_sub = df[(df["mutation_signal"] == mutsig) & (df["population_model"] == growth_model)].copy()

            if df_sub.empty or df_sub[time_col].isna().all():
                continue

            df_sub["time_bin"] = pd.cut(
                df_sub[time_col],
                bins=np.arange(0, df_sub[time_col].max() + bin_width, bin_width),
                include_lowest=True
            )

            for model in hue_order:
                df_model = df_sub[df_sub["model"] == model]
                grouped = df_model.groupby("time_bin", observed=True)
                valid = {k: g for k, g in grouped
                         if g[pop_error_col].notna().sum() > min_values_per_bin
                         and g[bl_error_col].notna().sum() > min_values_per_bin}

                x_median = pd.Series({k: g[pop_error_col].median() for k, g in valid.items()})
                y_median = pd.Series({k: g[bl_error_col].median() for k, g in valid.items()})

                if show_errorbars:
                    x_lower = pd.Series({k: g[pop_error_col].quantile(0.025) for k, g in valid.items()})
                    x_upper = pd.Series({k: g[pop_error_col].quantile(0.975) for k, g in valid.items()})
                    y_lower = pd.Series({k: g[bl_error_col].quantile(0.025) for k, g in valid.items()})
                    y_upper = pd.Series({k: g[bl_error_col].quantile(0.975) for k, g in valid.items()})

                    ax.errorbar(
                        x_median.values,
                        y_median.values,
                        xerr=[(x_median - x_lower).values, (x_upper - x_median).values],
                        yerr=[(y_median - y_lower).values, (y_upper - y_median).values],
                        fmt='o',
                        label = model,
                        #label=model if i == 0 and j == 0 else None,
                        color=colors[model],
                        ecolor = 'lightgray',
                        alpha=0.6,
                        capsize=3,
                        markersize=5
                    )
                else:
                    ax.scatter(
                        x_median,
                        y_median,
                        label = model,
                        #label=model if i == 0 and j == 0 else None,
                        color=colors[model],
                        alpha=0.6,
                        s=40
                    )

            if x_range:
                ax.set_xlim(x_range)
            if y_range:
                ax.set_ylim(y_range)

            if i == 0:
                ax.set_title(growth_model.replace("expgrowthfast", "Exp-growth fast").replace("expgrowthslow", "Exp-growth slow").replace("uniform", "Uniform").replace("bottleneck", "Bottleneck"))
            if j == 0:
                ax.set_ylabel(f"{mutsig} mutation signal\n{bl_error_col}")
            else:
                ax.set_ylabel(bl_error_col)

            ax.axhline(0, color="gray", linestyle="--", linewidth=1)
            ax.axvline(0, color="gray", linestyle="--", linewidth=1)
            ax.set_xlabel(pop_error_col)
            ax.legend(loc = 'upper left')

    #handles, labels = axes[0, 0].get_legend_handles_labels()
    #fig.legend(handles, labels, loc='upper right', title="Model")
    fig.suptitle("Branch Length vs. Population Size Errors\n(Median ± 95% CI in Time Bins)", fontsize=16)
    plt.tight_layout(rect=[0, 0, 0.95, 0.95])
    plt.show()


def plot_bl_vs_popsize_error_timecoloring_dualmodel(df, sampling, bin_width=50, time_col="height_sim",
                                                    pop_error_col="rel_diff_pop_size", bl_error_col="bl_relative_error",
                                                    x_range=None, y_range=None, show_errorbars=True,
                                                    min_values_per_bin=1, time_limit=None, log_colorbar=True,
                                                    x_label=None, y_label=None,
                                                    pop_models=None, mutsigs=None):
    
    """
    Scatter plot of branch length error vs population size error for both constcoal
    and skyline, colored by time bin using separate colormaps per model.
    Each point is the median error in a time bin; error bars show 95% CI across trees.
    Bins with fewer than min_values_per_bin valid values are excluded.
    """
    df = df[df["sampling"] == sampling].copy()
    if time_limit:
        df = df[df[time_col] < time_limit]

    mutsig_order = mutsigs if mutsigs is not None else ["low", "med", "high"]
    growth_model_order = pop_models if pop_models is not None else ["uniform", "expgrowthslow", "expgrowthfast", "bottleneck"]
    if isinstance(mutsig_order, str): mutsig_order = [mutsig_order]
    if isinstance(growth_model_order, str): growth_model_order = [growth_model_order]

    title_map = {"expgrowthfast": "Exp-growth fast", "expgrowthslow": "Exp-growth slow",
                 "uniform": "Uniform", "bottleneck": "Bottleneck"}
    model_markers = {"constcoal": "o", "skyline": "^"}
    model_cmaps = {"constcoal": plt.cm.Greens, "skyline": plt.cm.Reds}

    single_pop_model = len(growth_model_order) == 1

    df["time_bin"] = pd.cut(df[time_col],
                            bins=np.arange(0, df[time_col].max() + bin_width, bin_width),
                            include_lowest=True)
    all_bins = df["time_bin"].cat.categories
    bin_centers = [(b.left + b.right) / 2 for b in all_bins]
    n_bins = len(all_bins)
    t_min_pos = max(bin_centers[0], 1)
    t_max_val = bin_centers[-1]
    norm = plt.matplotlib.colors.LogNorm(vmin=t_min_pos, vmax=t_max_val) if log_colorbar \
           else plt.Normalize(vmin=t_min_pos, vmax=t_max_val)
    bin_to_idx = {b: i for i, b in enumerate(all_bins)}

    if single_pop_model:
        nrows, ncols = 1, len(mutsig_order)
    else:
        nrows, ncols = len(mutsig_order), len(growth_model_order)

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                             figsize=(5 * ncols, 3.5 * nrows),
                             squeeze=False)

    for i, mutsig in enumerate(mutsig_order):
        for j, growth_model in enumerate(growth_model_order):
            ax = axes[0, i] if single_pop_model else axes[i, j]
            df_sub = df[(df["mutation_signal"] == mutsig) & (df["population_model"] == growth_model)]

            for model in ["constcoal", "skyline"]:
                df_model = df_sub[df_sub["model"] == model]
                grouped = df_model.groupby("time_bin", observed=True)

                valid_groups = {k: g for k, g in grouped
                                if g[pop_error_col].notna().sum() >= min_values_per_bin
                                and g[bl_error_col].notna().sum() >= min_values_per_bin}

                for b, g in valid_groups.items():
                    color = model_cmaps[model](norm(max(bin_centers[bin_to_idx[b]], t_min_pos)))
                    xm = g[pop_error_col].median()
                    ym = g[bl_error_col].median()
                    label = model if (i == 0 and j == 0 and b == next(iter(valid_groups))) else None

                    if show_errorbars:
                        ax.errorbar(xm, ym,
                                    xerr=[[xm - g[pop_error_col].quantile(0.025)],
                                          [g[pop_error_col].quantile(0.975) - xm]],
                                    yerr=[[ym - g[bl_error_col].quantile(0.025)],
                                          [g[bl_error_col].quantile(0.975) - ym]],
                                    fmt=model_markers[model], color=color, ecolor='lightgray',
                                    elinewidth=1, alpha=0.7, capsize=3, markersize=5, label=label)
                    else:
                        ax.scatter(xm, ym, color=color, marker=model_markers[model],
                                   s=40, alpha=0.6, label=label)

            if x_range: ax.set_xlim(x_range)
            if y_range: ax.set_ylim(y_range)
            ax.axhline(0, color="gray", linestyle="--", linewidth=1)
            ax.axvline(0, color="gray", linestyle="--", linewidth=1)
            ax.set_xlabel(x_label if x_label else pop_error_col, fontsize=12)
            ax.tick_params(axis="both", labelsize=11)
            if single_pop_model:
                ax.set_title(f"{mutsig} mut. signal", fontsize=13)
                if i == 0: ax.set_ylabel(y_label if y_label else bl_error_col, fontsize=12)
            else:
                if i == 0: ax.set_title(title_map.get(growth_model, growth_model), fontsize=14)
                if j == 0: ax.set_ylabel(f"{mutsig} mut. signal\n{y_label if y_label else bl_error_col}", fontsize=12)

    fig.suptitle(f"Branch Length vs. Population Size Error ({sampling})\nMedian ± 95% CI, colored by time bin", fontsize=14)
    fig.subplots_adjust(right=0.86)

    fig.subplots_adjust(right=0.86)
    cbar_const = fig.colorbar(plt.cm.ScalarMappable(cmap=model_cmaps["constcoal"], norm=norm),
                               cax=fig.add_axes([0.88, 0.15, 0.015, 0.7]))
    cbar_const.set_label("Time before present — constcoal", rotation=270, labelpad=15, fontsize=12)

    cbar_sky = fig.colorbar(plt.cm.ScalarMappable(cmap=model_cmaps["skyline"], norm=norm),
                             cax=fig.add_axes([0.95, 0.15, 0.015, 0.7]))
    cbar_sky.set_label("Time before present — skyline", rotation=270, labelpad=15, fontsize=12)

    plt.show()



def plot_hpd_width_vs_bl(df, sampling, tree_index,
                          pop_models=None, mutsigs=None,
                          x_range=None, y_range=None,
                          alpha=0.5, s=15, title=None):
    """
    Scatterplot of 95% HPD width (bl_ci_upper_estimate - bl_ci_lower_estimate)
    vs median estimated branch length (bl_estimate), for one specific tree.
    Columns = population models, rows = mutation signals.
    Constcoal = red, skyline = blue.
    """
    df = df[(df["sampling"] == sampling) & (df["tree_index"] == tree_index)].copy()
    df = df.dropna(subset=["bl_estimate", "bl_ci_lower_estimate", "bl_ci_upper_estimate"])
    df["hpd_width"] = df["bl_ci_upper_estimate"] - df["bl_ci_lower_estimate"]

    mutsig_order = mutsigs if mutsigs is not None else ["low", "med", "high"]
    growth_model_order = pop_models if pop_models is not None else ["uniform", "expgrowthslow", "expgrowthfast", "bottleneck"]
    if isinstance(mutsig_order, str): mutsig_order = [mutsig_order]
    if isinstance(growth_model_order, str): growth_model_order = [growth_model_order]

    title_map = {"expgrowthfast": "Exp-growth fast", "expgrowthslow": "Exp-growth slow",
                 "uniform": "Uniform", "bottleneck": "Bottleneck"}
    model_colors = {"constcoal": "#a6444f", "skyline": "#397398"}

    fig, axes = plt.subplots(len(mutsig_order), len(growth_model_order),
                             figsize=(4 * len(growth_model_order), 3.5 * len(mutsig_order)),
                             squeeze=False)

    for i, mutsig in enumerate(mutsig_order):
        for j, growth_model in enumerate(growth_model_order):
            ax = axes[i, j]
            df_sub = df[(df["mutation_signal"] == mutsig) & (df["population_model"] == growth_model)]

            for model, color in model_colors.items():
                df_m = df_sub[df_sub["model"] == model]
                if df_m.empty:
                    continue
                ax.scatter(df_m["bl_estimate"], df_m["hpd_width"],
                           color=color, s=s, alpha=alpha, label=model)

            if x_range: ax.set_xlim(x_range)
            if y_range: ax.set_ylim(y_range)
            ax.set_xlabel("Estimated branch length", fontsize=9)
            ax.tick_params(labelsize=8)
            if i == 0: ax.set_title(title_map.get(growth_model, growth_model), fontsize=11)
            if j == 0: ax.set_ylabel(f"{mutsig} mut. signal\nHPD width", fontsize=9)

    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=c,
                           markersize=7, label=m) for m, c in model_colors.items()]
    fig.legend(handles=handles, loc='upper right', bbox_to_anchor=(0.98, 1.02), ncol=2)
    plt.suptitle(title or f"Tree {tree_index} — HPD width vs bl_estimate ({sampling})", fontsize=13)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.show()


def plot_bl_error_vs_bl_sim(df, sampling, tree_index,
                             y_col="bl_relative_error",
                             pop_models=None, mutsigs=None,
                             x_range=None, y_range=None,
                             alpha=0.5, s=15, title=None):
    """
    Scatterplot of branch length error (y_col) vs simulated branch length (bl_sim)
    for one specific tree_index and sampling type.
    Columns = population models, rows = mutation signals.
    Constcoal = red, skyline = blue.
    Error bars from {y_col}_lower / {y_col}_upper if available.
    """
    df = df[(df["sampling"] == sampling) & (df["tree_index"] == tree_index)].copy()
    df = df.dropna(subset=["bl_sim", y_col])

    mutsig_order = mutsigs if mutsigs is not None else ["low", "med", "high"]
    growth_model_order = pop_models if pop_models is not None else ["uniform", "expgrowthslow", "expgrowthfast", "bottleneck"]
    if isinstance(mutsig_order, str): mutsig_order = [mutsig_order]
    if isinstance(growth_model_order, str): growth_model_order = [growth_model_order]

    title_map = {"expgrowthfast": "Exp-growth fast", "expgrowthslow": "Exp-growth slow",
                 "uniform": "Uniform", "bottleneck": "Bottleneck"}
    model_colors = {"constcoal": "#a6444f", "skyline": "#397398"}

    y_lower_col = f"{y_col}_lower" if f"{y_col}_lower" in df.columns else None
    y_upper_col = f"{y_col}_upper" if f"{y_col}_upper" in df.columns else None
    has_errorbars = y_lower_col is not None

    fig, axes = plt.subplots(len(mutsig_order), len(growth_model_order),
                             figsize=(4 * len(growth_model_order), 3.5 * len(mutsig_order)),
                             squeeze=False)

    for i, mutsig in enumerate(mutsig_order):
        for j, growth_model in enumerate(growth_model_order):
            ax = axes[i, j]
            df_sub = df[(df["mutation_signal"] == mutsig) & (df["population_model"] == growth_model)]

            for model, color in model_colors.items():
                df_m = df_sub[df_sub["model"] == model]
                if df_m.empty:
                    continue
                if has_errorbars:
                    yerr = [df_m[y_col] - df_m[y_lower_col],
                            df_m[y_upper_col] - df_m[y_col]]
                    ax.errorbar(df_m["bl_sim"], df_m[y_col],
                                yerr=yerr, fmt='o', color=color,
                                ecolor='lightgray', elinewidth=0.6,
                                alpha=alpha, capsize=2, markersize=3, label=model)
                else:
                    ax.scatter(df_m["bl_sim"], df_m[y_col],
                               color=color, s=s, alpha=alpha, label=model)

            ax.axhline(0, color="gray", linestyle="--", linewidth=0.8)
            if x_range: ax.set_xlim(x_range)
            if y_range: ax.set_ylim(y_range)
            ax.set_xlabel("Simulated branch length", fontsize=9)
            ax.tick_params(labelsize=8)
            if i == 0: ax.set_title(title_map.get(growth_model, growth_model), fontsize=11)
            if j == 0: ax.set_ylabel(f"{mutsig} mut. signal\n{y_col}", fontsize=9)

    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=c,
                           markersize=7, label=m) for m, c in model_colors.items()]
    fig.legend(handles=handles, loc='upper right', bbox_to_anchor=(0.98, 1.02), ncol=2)
    plt.suptitle(title or f"Tree {tree_index} — {y_col} vs bl_sim ({sampling})", fontsize=13)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.show()


def plot_single_tree_errors(df, sampling, tree_index,
                            x_col="diff_pop_size", y_col="bl_relative_error",
                            pop_models=None, mutsigs=None,
                            x_range=None, y_range=None,
                            time_col="height_sim", log_colorbar=True,
                            title=None):
    """
    Scatter plot of y_col vs x_col for a single tree_index, one point per node.
    Color = node height (time before present), colormapped via viridis.
    Error bars inferred from {col}_lower and {col}_upper columns if present.
    One subplot per population model × mutation signal.
    Separate panels per model (constcoal / skyline).
    """
    df = df[(df["sampling"] == sampling) & (df["tree_index"] == tree_index)].copy()

    mutsig_order = mutsigs if mutsigs is not None else ["low", "med", "high"]
    growth_model_order = pop_models if pop_models is not None else ["uniform", "expgrowthslow", "expgrowthfast", "bottleneck"]
    if isinstance(mutsig_order, str): mutsig_order = [mutsig_order]
    if isinstance(growth_model_order, str): growth_model_order = [growth_model_order]

    title_map = {"expgrowthfast": "Exp-growth fast", "expgrowthslow": "Exp-growth slow",
                 "uniform": "Uniform", "bottleneck": "Bottleneck"}
    model_markers = {"constcoal": "o", "skyline": "^"}
    model_cmaps   = {"constcoal": plt.cm.Greens, "skyline": plt.cm.Reds}

    t_vals = df[time_col].dropna()
    t_min_pos = max(t_vals.min(), 1)
    t_max_val = t_vals.max()
    norm = plt.matplotlib.colors.LogNorm(vmin=t_min_pos, vmax=t_max_val) if log_colorbar \
           else plt.Normalize(vmin=t_min_pos, vmax=t_max_val)

    # Check if lower/upper columns exist
    x_lower_col = f"{x_col}_lower" if f"{x_col}_lower" in df.columns else None
    x_upper_col = f"{x_col}_upper" if f"{x_col}_upper" in df.columns else None
    y_lower_col = f"{y_col}_lower" if f"{y_col}_lower" in df.columns else None
    y_upper_col = f"{y_col}_upper" if f"{y_col}_upper" in df.columns else None
    has_errorbars = any([x_lower_col, y_lower_col])

    nrows, ncols = len(mutsig_order), len(growth_model_order)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.5 * nrows), squeeze=False)

    for i, mutsig in enumerate(mutsig_order):
        for j, growth_model in enumerate(growth_model_order):
            ax = axes[i, j]
            df_sub = df[(df["mutation_signal"] == mutsig) & (df["population_model"] == growth_model)]

            for model in ["constcoal", "skyline"]:
                df_m = df_sub[df_sub["model"] == model].dropna(subset=[x_col, y_col, time_col])
                if df_m.empty:
                    continue
                cmap_m = model_cmaps[model]

                if has_errorbars:
                    for _, row in df_m.iterrows():
                        color = cmap_m(norm(max(row[time_col], t_min_pos)))
                        xerr = [[row[x_col] - row[x_lower_col]], [row[x_upper_col] - row[x_col]]] \
                               if x_lower_col else None
                        yerr = [[row[y_col] - row[y_lower_col]], [row[y_upper_col] - row[y_col]]] \
                               if y_lower_col else None
                        ax.errorbar(row[x_col], row[y_col],
                                    xerr=xerr, yerr=yerr,
                                    fmt=model_markers[model], color=color,
                                    ecolor='lightgray', elinewidth=0.8,
                                    alpha=0.7, capsize=2, markersize=4)
                else:
                    colors = [cmap_m(norm(max(t, t_min_pos))) for t in df_m[time_col]]
                    ax.scatter(df_m[x_col], df_m[y_col],
                               c=colors, marker=model_markers[model], s=20, alpha=0.7)

            ax.axhline(0, color="gray", linestyle="--", linewidth=0.8)
            ax.axvline(0, color="gray", linestyle="--", linewidth=0.8)
            if x_range: ax.set_xlim(x_range)
            if y_range: ax.set_ylim(y_range)
            ax.set_xlabel(x_col, fontsize=10)
            ax.tick_params(labelsize=9)
            if i == 0: ax.set_title(title_map.get(growth_model, growth_model), fontsize=11)
            if j == 0: ax.set_ylabel(f"{mutsig} mut. signal\n{y_col}", fontsize=10)

    fig.subplots_adjust(right=0.86)
    cbar_const = fig.colorbar(plt.cm.ScalarMappable(cmap=model_cmaps["constcoal"], norm=norm),
                               cax=fig.add_axes([0.88, 0.15, 0.015, 0.7]))
    cbar_const.set_label("Time before present — constcoal", rotation=270, labelpad=15, fontsize=10)

    cbar_sky = fig.colorbar(plt.cm.ScalarMappable(cmap=model_cmaps["skyline"], norm=norm),
                             cax=fig.add_axes([0.95, 0.15, 0.015, 0.7]))
    cbar_sky.set_label("Time before present — skyline", rotation=270, labelpad=15, fontsize=10)

    plt.suptitle(title or f"Tree {tree_index} — {x_col} vs {y_col} ({sampling})", fontsize=13)
    plt.show()


def plot_bl_vs_popsize_error_timecoloring(df, sampling, bin_width=50, time_col="height_sim",
                                         pop_error_col="rel_diff_pop_size", bl_error_col="bl_relative_error",
                                         x_range=None, y_range=None, show_errorbars=True,
                                         model='skyline', time_limit=None, log_colorbar=True):
    """
    Scatter plot of branch length error vs population size error, with points
    colored by time bin. One subplot per population model × mutation signal.
    """
    df = df[(df["sampling"] == sampling) & (df["model"] == model)].copy()
    if time_limit:
        df = df[df[time_col] < time_limit]

    mutsig_order = ["low", "med", "high"]
    growth_model_order = ["uniform", "expgrowthfast", "expgrowthslow", "bottleneck"]
    title_map = {"expgrowthfast": "Exp-growth fast", "expgrowthslow": "Exp-growth slow",
                 "uniform": "Uniform", "bottleneck": "Bottleneck"}
    cmap = plt.cm.viridis

    df["time_bin"] = pd.cut(df[time_col],
                            bins=np.arange(0, df[time_col].max() + bin_width, bin_width),
                            include_lowest=True)
    all_bins = df["time_bin"].cat.categories
    bin_centers = [(b.left + b.right) / 2 for b in all_bins]
    t_min_pos = max(bin_centers[0], 1)
    t_max_val = bin_centers[-1]
    norm = plt.matplotlib.colors.LogNorm(vmin=t_min_pos, vmax=t_max_val) if log_colorbar \
           else plt.Normalize(vmin=t_min_pos, vmax=t_max_val)
    bin_colors = {b: cmap(norm(max(c, t_min_pos))) for b, c in zip(all_bins, bin_centers)}

    fig, axes = plt.subplots(nrows=len(mutsig_order), ncols=len(growth_model_order),
                             figsize=(5 * len(growth_model_order), 3.5 * len(mutsig_order)))

    for i, mutsig in enumerate(mutsig_order):
        for j, growth_model in enumerate(growth_model_order):
            ax = axes[i, j]
            df_sub = df[(df["mutation_signal"] == mutsig) & (df["population_model"] == growth_model)]
            grouped = df_sub.groupby("time_bin", observed=True)

            x_med = grouped[pop_error_col].median()
            y_med = grouped[bl_error_col].median()

            for b in x_med.index:
                color = bin_colors[b]
                if show_errorbars:
                    x_lo, x_hi = grouped[pop_error_col].quantile(0.025)[b], grouped[pop_error_col].quantile(0.975)[b]
                    y_lo, y_hi = grouped[bl_error_col].quantile(0.025)[b], grouped[bl_error_col].quantile(0.975)[b]
                    ax.errorbar(x_med[b], y_med[b],
                                xerr=[[x_med[b] - x_lo], [x_hi - x_med[b]]],
                                yerr=[[y_med[b] - y_lo], [y_hi - y_med[b]]],
                                fmt='o', color=color, ecolor='lightgray',
                                elinewidth=1, alpha=0.7, capsize=3)
                else:
                    ax.scatter(x_med[b], y_med[b], color=color, alpha=0.6, s=40)

            ax.axhline(0, color="gray", linestyle="--", linewidth=1)
            ax.axvline(0, color="gray", linestyle="--", linewidth=1)
            if x_range: ax.set_xlim(x_range)
            if y_range: ax.set_ylim(y_range)
            ax.set_xlabel(pop_error_col)
            if i == 0: ax.set_title(title_map.get(growth_model, growth_model))
            if j == 0: ax.set_ylabel(f"{mutsig} mut. signal\n{bl_error_col}")

    # Logarithmic colorbar using actual time values
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.015, 0.7])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Time before present', rotation=270, labelpad=15)

    fig.suptitle(f"Branch Length vs. Population Size Error — {model} ({sampling})", fontsize=14)
    plt.show()


def plot_height_errors_over_time_bin(ax, df_sub, t_max, bins=10, bin_width=None, error_col="height_abs_relative_error", x_col='height_sim', y_range=None, add_legend=True):
    df_sub = df_sub.copy()

    # Only consider rows within the [0, t_max] time range
    df_sub = df_sub[(df_sub[x_col] >= 0) & (df_sub[x_col] <= t_max)]

    # Create bin edges
    if bin_width is not None:
        bin_edges = list(range(0, t_max, bin_width))
        bin_edges.append(max(df_sub[x_col].max() + 1e-6, t_max + 1e-6))
    else:
        bin_edges = np.linspace(0, t_max, bins + 1)
    df_sub["time_bin"] = pd.cut(df_sub[x_col], bins=bin_edges, include_lowest=True)
    df_sub = df_sub.dropna(subset=["time_bin", error_col])

    colors = {"skyline": "#80557e", "constcoal": "#397398"}
    models = df_sub["model"].unique()

    for model in models:
        df_model = df_sub[df_sub["model"] == model]
        grouped = df_model.groupby("time_bin", observed=True)

        # Get bin midpoints for x-axis
        x_medians = [interval.mid for interval in grouped.groups.keys()]
        y_medians = grouped[error_col].median()
        y_lower = grouped[error_col].quantile(0.025)
        y_upper = grouped[error_col].quantile(0.975)

        ax.errorbar(
            x=x_medians,
            y=y_medians,
            yerr=[y_medians - y_lower, y_upper - y_medians],
            fmt='o',
            label=model,
            color=colors.get(model, 'gray'),
            ecolor='lightgray',
            capsize=3,
            alpha=0.6,
            markersize=4
        )
    
    ax.axhline(0, color="gray", linestyle="--", linewidth=1)
    ax.set_xlim(t_max, 0)  # reverse x-axis

    ax.set_xlabel("Time before present")
    if y_range:
        ax.set_ylim(y_range)

    ax.tick_params(axis='x', rotation=45)



def plot_combined_population_and_error(df_population, df_error, sampling,
                                       time_horizon=800, pop_y_range=(0, 5000),
                                       error_y_range=(-1, 3), bins=10, bin_width=None,
                                       error_col='bl_relative_error', x_col='height_sim',
                                       mode='both', plot_type='violin',
                                       pop_models=None, mutsigs=None):
    df_population = df_population[df_population["sampling"] == sampling].copy()
    df_error = df_error[df_error["sampling"] == sampling].copy()

    mutsig_order = mutsigs if mutsigs is not None else ["low", "med", "high"]
    growth_model_order = pop_models if pop_models is not None else ["uniform", "expgrowthslow", "expgrowthfast", "bottleneck"]

    # Ensure lists
    if isinstance(mutsig_order, str): mutsig_order = [mutsig_order]
    if isinstance(growth_model_order, str): growth_model_order = [growth_model_order]

    ncols = len(growth_model_order)

    # 2 data rows per mutsig + 1 spacer between mutsig groups (except after last)
    title_map = {"expgrowthfast": "Exp-growth fast", "expgrowthslow": "Exp-growth slow",
                 "uniform": "Uniform", "bottleneck": "Bottleneck"}

    single_pop_model = len(growth_model_order) == 1

    if single_pop_model:
        # Layout: cols = mutsigs, rows = [pop, error]
        ncols = len(mutsig_order)
        fig, axes = plt.subplots(2, ncols, figsize=(5 * ncols, 5),
                                 gridspec_kw={"hspace": 0.4, "wspace": 0.4})
        if ncols == 1:
            axes = axes.reshape(2, 1)
    else:
        # Layout: cols = pop models, rows = mutsig groups (pop + error + spacer)
        row_heights = [1, 1, 0.3] * len(mutsig_order)
        row_heights = row_heights[:-1]
        fig = plt.figure(figsize=(5 * ncols, 2.5 * len(row_heights)))
        gs = gridspec.GridSpec(len(row_heights), ncols, height_ratios=row_heights, hspace=0.6, wspace=0.4)

    for i, mutsig in enumerate(mutsig_order):
        for j, growth_model in enumerate(growth_model_order):
            if single_pop_model:
                ax_pop = axes[0, i]
                ax_err = axes[1, i]
            else:
                row_pop = i * 3
                row_err = i * 3 + 1
                ax_pop = fig.add_subplot(gs[row_pop, j])
                ax_err = fig.add_subplot(gs[row_err, j])

            subset_pop = df_population[
                (df_population["mutation_signal"] == mutsig) &
                (df_population["population_model"] == growth_model)]
            subset_err = df_error[
                (df_error["mutation_signal"] == mutsig) &
                (df_error["population_model"] == growth_model)]

            if not subset_pop.empty:
                plot_population_summary_ax(
                    ax=ax_pop,
                    skyline_all_times=subset_pop["skyline_times"].tolist(),
                    skyline_all_medians=subset_pop["skyline_medians"].tolist(),
                    constant_all_estimates=subset_pop["coalescent_median"].tolist(),
                    true_traj=get_true_traj(growth_model, subset_pop.iloc[0]),
                    time_horizon=time_horizon, y_range=pop_y_range, mode=mode)
                ax_pop.set_xlim(time_horizon, 0)

            if not subset_err.empty:
                if plot_type == "scatter":
                    plot_height_errors_over_time_bin(
                        ax=ax_err, df_sub=subset_err, t_max=time_horizon,
                        bins=bins, bin_width=bin_width, error_col=error_col, x_col=x_col,
                        y_range=error_y_range, add_legend=False)
                else:
                    plot_height_errors_by_time_bin(
                        ax=ax_err, df_sub=subset_err, t_max=time_horizon,
                        bins=bins, bin_width=bin_width,
                        error_col=error_col, x_col=x_col, y_range=error_y_range,
                        plot_type=plot_type, add_legend=False,
                        skyline_col="#80557e", constcoal_col="#397398")
                    ax_err.invert_xaxis()

            if single_pop_model:
                ax_pop.set_title(f"{mutsig} mut. signal")
                if i == 0:
                    ax_pop.set_ylabel("Population size")
                    ax_err.set_ylabel(error_col)
            else:
                if j == 0:
                    ax_pop.set_ylabel(f"{mutsig} mut. signal\nPopulation size")
                    ax_err.set_ylabel(error_col)
                if i == 0:
                    ax_pop.set_title(title_map.get(growth_model, growth_model))

    handles = [Line2D([0], [0], color="#80557e", lw=4, label="Skyline"),
               Line2D([0], [0], color="#397398", lw=4, label="Const. Coalescent")]
    fig.legend(handles=handles, title="Model", loc='upper right', bbox_to_anchor=(0.98, 1.02), ncol=2)
    plt.suptitle(f"Population Size & {error_col} ({sampling})", fontsize=16)
    plt.show()


def compute_pairwise_snp_distances(repo_base=None):
    """
    For each SNP alignment file in results/run1/simulated_data/{sampling}/{pop_model}/,
    compute the average pairwise SNP distance between all pairs of isolates.

    Returns a DataFrame with columns: sampling, population_model, tree_index, mean_pairwise_snps
    """
    if repo_base is None:
        repo_base = Path("/Users/mariebecker/Documents/Uni/ETH/RotationStadler/BESP_paper-analyses")

    samplings = ["independenthomochronous", "linearconstant"]
    pop_models = ["expgrowthfast", "expgrowthslow", "uniform", "bottleneck"]

    rows = []
    for sampling in samplings:
        for pop_model in pop_models:
            data_dir = repo_base / "results/run1/simulated_data" / sampling / pop_model
            for fasta_file in sorted(data_dir.glob(f"{pop_model}_*_snps.fasta")):
                # Extract tree index from filename
                match = re.search(r'_(\d+)_snps\.fasta$', fasta_file.name)
                if not match:
                    continue
                tree_index = int(match.group(1))

                # Read sequences as numpy byte array (one row per sequence)
                seqs = [str(rec.seq) for rec in SeqIO.parse(fasta_file, "fasta")]
                if len(seqs) < 2:
                    continue

                mat = np.frombuffer("".join(seqs).encode(), dtype=np.uint8).reshape(len(seqs), -1)
                n = len(seqs)
                total_diffs = sum(int((mat[i] != mat[j]).sum()) for i in range(n) for j in range(i+1, n))
                mean_dist = total_diffs / (n * (n - 1) / 2)

                rows.append({
                    "sampling": sampling,
                    "population_model": pop_model,
                    "tree_index": tree_index,
                    "mean_pairwise_snps": mean_dist
                })

    return pd.DataFrame(rows)


def compute_simulated_tree_lengths(repo_base=None):
    """
    For each simulated tree in results/run1/simulated_data/{sampling}/{pop_model}/{pop_model}.trees,
    compute total branch length and tip branch length.

    Returns a DataFrame with columns:
    sampling, population_model, tree_index, total_branch_length_sim, tip_branch_length_sim
    """
    if repo_base is None:
        repo_base = Path("/Users/mariebecker/Documents/Uni/ETH/RotationStadler/BESP_paper-analyses")

    samplings = ["independenthomochronous", "linearconstant"]
    pop_models = ["expgrowthfast", "expgrowthslow", "uniform", "bottleneck"]

    rows = []
    for sampling in samplings:
        for pop_model in pop_models:
            trees_file = repo_base / "results/run1/simulated_data" / sampling / pop_model / f"{pop_model}.trees"
            if not trees_file.exists():
                continue
            for tree_index, tree in enumerate(Phylo.parse(str(trees_file), "newick")):
                tip_names = {c.name for c in tree.get_terminals()}
                total_bl = sum(c.branch_length for c in tree.find_clades()
                               if c.branch_length is not None)
                tip_bl = sum(c.branch_length for c in tree.find_clades()
                             if c.branch_length is not None and c.name in tip_names)

                # For each internal node: sum of its two direct child branch lengths
                child_pair_bls = []
                for clade in tree.find_clades(order="preorder"):
                    if clade.clades:  # internal node
                        child_sum = sum(c.branch_length for c in clade.clades
                                        if c.branch_length is not None)
                        child_pair_bls.append(child_sum)
                mean_child_pair_bl = np.mean(child_pair_bls) if child_pair_bls else np.nan

                rows.append({
                    "sampling": sampling,
                    "population_model": pop_model,
                    "tree_index": tree_index,
                    "total_branch_length_sim": total_bl,
                    "tip_branch_length_sim": tip_bl,
                    "mean_child_pair_branch_length_sim": mean_child_pair_bl,
                })

    return pd.DataFrame(rows)


def plot_tree_length_by_scenario(df_tree_lengths, metric="total_branch_length_sim",
                                  title=None, y_range=None):
    """
    Violin + strip plot of tree branch length metrics per scenario.

    metric: "total_branch_length_sim" or "tip_branch_length_sim"
    Two colors for the two sampling types (independenthomochronous and linearconstant).
    One subplot per population model (columns) × mutation signal (rows).
    """
    pop_model_order = ["uniform", "expgrowthslow", "expgrowthfast", "bottleneck"]
    title_map = {"expgrowthfast": "Exp-growth fast", "expgrowthslow": "Exp-growth slow",
                 "uniform": "Uniform", "bottleneck": "Bottleneck"}

    sampling_label_map = {"independenthomochronous": "Homochronous", "linearconstant": "Heterochronous"}
    sampling_colors = {"independenthomochronous": "#397398", "linearconstant": "#a6444f"}
    hue_order = ["independenthomochronous", "linearconstant"]

    fig, axes = plt.subplots(1, len(pop_model_order),
                             figsize=(4 * len(pop_model_order), 4),
                             sharey=False, squeeze=False)

    for j, pop_model in enumerate(pop_model_order):
            ax = axes[0, j]
            subset = df_tree_lengths[df_tree_lengths["population_model"] == pop_model]

            if subset.empty:
                ax.axis("off")
                continue

            sns.violinplot(data=subset, y=metric, x="sampling", hue="sampling",
                           order=hue_order, palette=sampling_colors, cut=0, inner="box",
                           density_norm="width", ax=ax, alpha=0.6, legend=False)
            sns.stripplot(data=subset, y=metric, x="sampling", hue="sampling",
                          order=hue_order, palette=sampling_colors, ax=ax,
                          size=3, alpha=0.5, jitter=True, legend=False)
            ax.grid(axis='y', linestyle='--', alpha=0.5, zorder=0)

            if y_range:
                ax.set_ylim(y_range)
            ax.set_xlabel("")
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["Homochronous", "Heterochronous"], fontsize=8, rotation=15, ha="right")
            if j == 0:
                ax.set_ylabel(metric, fontsize=9)
            else:
                ax.set_ylabel("")
            ax.set_title(title_map.get(pop_model, pop_model), fontsize=11)

    handles = [plt.Rectangle((0, 0), 1, 1, color=c, label=sampling_label_map[s])
               for s, c in sampling_colors.items()]
    fig.legend(handles=handles, loc='upper right', bbox_to_anchor=(0.98, 1.02),
               ncol=2, frameon=False)
    plt.suptitle(title or metric, fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.show()


def plot_bl_sim_vs_estimate_grid(
    tree_metrics_combined,
    pop_models,
    mutation_signals,
    tree_index=3,
    n_values=None,
    time_range=None,
    model="skyline",
    hetero_sampling="linearconstant",
    homo_sampling="independenthomochronous",
    sort_by="bl_sim",
    alpha=0.3,
    s=20,
    figsize_per_subplot=(5, 4),
    mode="bl",
    substitution_rate=None,
    genome_length=None,
    sub_range=None,
):
    """
    Plot bl_sim vs bl_estimate for heterochronous and homochronous sampling.

    Columns = population models
    Rows = mutation signals

    Use either:
    - n_values: take the first n rows after sorting by sort_by
    - time_range: tuple (min_height, max_height), filter by height_sim

    Exactly one of n_values or time_range should be set.

    mode="bl": plot raw branch lengths in years.
    mode="expected_substitutions": multiply branch lengths by substitution_rate * genome_length
        to get expected number of substitutions per branch.
    """
    if mode == "expected_substitutions":
        if substitution_rate is None or genome_length is None:
            raise ValueError("substitution_rate and genome_length must be set when mode='expected_substitutions'")
        if not isinstance(genome_length, dict):
            genome_length = {"low": genome_length, "med": genome_length, "high": genome_length}

    if sub_range is not None:
        n_values = None
        time_range = None
    elif n_values is not None and time_range is not None:
        raise ValueError("Set at most one of n_values or time_range.")

    if time_range is not None:
        if not isinstance(time_range, tuple) or len(time_range) != 2:
            raise ValueError("time_range must be a tuple like (0, 100).")

        time_min, time_max = time_range

    n_rows = len(mutation_signals)
    n_cols = len(pop_models)

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(figsize_per_subplot[0] * n_cols, figsize_per_subplot[1] * n_rows),
        sharex=False,
        sharey=False,
        squeeze=False,
    )

    def filter_slice(df, sampling, pop_model, mut_signal):
        subset = tree_metrics_combined[
            (tree_metrics_combined["sampling"] == sampling) &
            (tree_metrics_combined["population_model"] == pop_model) &
            (tree_metrics_combined["mutation_signal"] == mut_signal) &
            (tree_metrics_combined["model"] == model) &
            (tree_metrics_combined["tree_index"] == tree_index)
        ].copy()

        if time_range is not None:
            subset = subset[
                (subset["height_sim"] >= time_min) &
                (subset["height_sim"] <= time_max)
            ].copy()

            subset = subset.sort_values(sort_by).copy()

        else:
            subset = subset.sort_values(sort_by).iloc[:n_values].copy()

        return subset

    for row_idx, mut_signal in enumerate(mutation_signals):
        for col_idx, pop_model in enumerate(pop_models):
            ax = axes[row_idx, col_idx]

            hetero = filter_slice(
                tree_metrics_combined,
                hetero_sampling,
                pop_model,
                mut_signal,
            )

            homo = filter_slice(
                tree_metrics_combined,
                homo_sampling,
                pop_model,
                mut_signal,
            )

            hetero["sampling_type"] = "heterochronous"
            homo["sampling_type"] = "homochronous"

            plot_df = pd.concat([hetero, homo], ignore_index=True)

            if plot_df.empty:
                ax.set_title(f"{pop_model}\n{mut_signal}\n(no data)")
                ax.axis("off")
                continue

            if mode == "expected_substitutions":
                scale = substitution_rate * genome_length[mut_signal]
                plot_df["bl_sim"]      = plot_df["bl_sim"] * scale
                plot_df["bl_estimate"] = plot_df["bl_estimate"] * scale
                xlabel = "Simulated expected substitutions"
                ylabel = "Estimated expected substitutions"
            else:
                xlabel = "Simulated branch length (years)"
                ylabel = "Estimated branch length (years)"

            if sub_range is not None:
                plot_df = plot_df[
                    (plot_df["bl_sim"] >= sub_range[0]) & (plot_df["bl_sim"] <= sub_range[1])
                ]

            color_map = {"heterochronous": "#a6444f", "homochronous": "#397398"}
            for sampling_type, group in plot_df.groupby("sampling_type"):
                ax.scatter(
                    group["bl_sim"],
                    group["bl_estimate"],
                    label=sampling_type,
                    color=color_map[sampling_type],
                    alpha=alpha,
                    s=s,
                )

            min_val = min(plot_df["bl_sim"].min(), plot_df["bl_estimate"].min())
            max_val = max(plot_df["bl_sim"].max(), plot_df["bl_estimate"].max())

            ax.plot([min_val, max_val], [min_val, max_val],
                    linestyle="--", color="black", linewidth=1)

            ax.set_title(f"{pop_model} | {mut_signal}")
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.grid(alpha=0.3)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            labels,
            loc="upper right",
            ncol=2,
            frameon=False,
        )

    if time_range is not None:
        filter_description = f"height_sim in [{time_min}, {time_max}]"
    else:
        filter_description = f"first {n_values} rows sorted by {sort_by}"

    fig.suptitle(
        f"bl_sim vs bl_estimate by sampling type | "
        f"tree_index={tree_index}, {filter_description}",
        y=1.02,
        fontsize=14,
    )

    plt.tight_layout()
    plt.show()

def read_tree_flexible(tree_path, tree_index=0):
    """
    Read a tree from either:
      - a newick .trees file (multiple trees) — picks tree_index (0-based)
      - a nexus .tree or .nex file (single summary tree)
    """
    path = str(tree_path)
    if path.endswith(".trees") and not path.endswith("_summary.tree"):
        trees = list(Phylo.parse(path, "newick"))
        return trees[tree_index]
    else:
        return Phylo.read(path, "nexus")




def coalescent_intensity_all_trees(trees_df):
    """
    For each replicate, compute cumulative coalescent intensity at each internal
    node of the TRUE simulated tree, evaluated under three Ne functions:
        - true trajectory
        - constcoal constant Ne (BEAST estimate)
        - skyline step-function Ne (BEAST estimate)

    All three are integrated along the simulated tree's branch lengths and
    evaluated at the simulated tree's node times. This ensures a fair comparison:
    the only difference is the Ne function, not the tree topology or node times.

    For the skyline, the step function uses BEAST-estimated skyline_times and
    skyline_medians but is integrated up to the true node time (which may be
    shorter or longer than the BEAST tree's root height). Beyond the last skyline
    boundary, the last Ne value is used.
    """
    results = []

    for _, row in trees_df.iterrows():
        pop_model = row["population_model"]
        params = assign_model_params(pop_model)

        # --- true trajectory integral function ---
        if pop_model in ("expgrowthfast", "expgrowthslow"):
            true_mode = "true_exp"
            ci_params_true = {"N0": params["present_pop_size"], "rate": params["growth_rate"]}
        elif pop_model.startswith("bottleneck"):
            true_mode = "true_bottleneck"
            ci_params_true = {
                "Ne_normal":        params["present_pop_size"],
                "Ne_bottleneck":    params["bottleneck_size"],
                "bottleneck_start": params["bottleneck_start"],
                "bottleneck_end":   params["bottleneck_end"],
            }
        else:
            true_mode = "true_uniform"
            ci_params_true = {"Ne": params["present_pop_size"]}

        constcoal_Ne     = row["coalescent_median"]
        skyline_boundaries = [0.0] + list(row["skyline_times"])
        skyline_pop_sizes  = list(row["skyline_medians"])

        meta = {
            "tree_index":       row["tree_index"],
            "population_model": pop_model,
            "mutation_signal":  row["mutation_signal"],
            "sampling":         row["sampling"],
        }

        # --- inline rate helpers (instantaneous 1/Ne at time t) ---
        def _rate_true(t):
            if true_mode == "true_uniform":
                return 1.0 / ci_params_true["Ne"]
            if true_mode == "true_exp":
                return exp(ci_params_true["rate"] * t) / ci_params_true["N0"]
            if true_mode == "true_bottleneck":
                bs, be = ci_params_true["bottleneck_start"], ci_params_true["bottleneck_end"]
                Ne = ci_params_true["Ne_bottleneck"] if bs < t < be else ci_params_true["Ne_normal"]
                return 1.0 / Ne

        def _rate_skyline(t):
            for a, b, Ne in zip(skyline_boundaries[:-1], skyline_boundaries[1:], skyline_pop_sizes):
                if a <= t < b:
                    return 1.0 / Ne
            return 1.0 / skyline_pop_sizes[-1]

        # --- inline integral helpers ---
        def _integral_true(t0, t1):
            if t1 <= t0:
                return 0.0
            if true_mode == "true_uniform":
                return (t1 - t0) / ci_params_true["Ne"]
            if true_mode == "true_exp":
                N0, rate = ci_params_true["N0"], ci_params_true["rate"]
                return (exp(rate * t1) - exp(rate * t0)) / (N0 * rate)
            if true_mode == "true_bottleneck":
                bs, be = ci_params_true["bottleneck_start"], ci_params_true["bottleneck_end"]
                Ne_n, Ne_b = ci_params_true["Ne_normal"], ci_params_true["Ne_bottleneck"]
                total = 0.0
                for a, b, Ne in [(0.0, bs, Ne_n), (bs, be, Ne_b), (be, inf, Ne_n)]:
                    left, right = max(t0, a), min(t1, b)
                    if right > left:
                        total += (right - left) / Ne
                return total

        def _integral_constcoal(t0, t1):
            return (t1 - t0) / constcoal_Ne if t1 > t0 else 0.0

        def _integral_skyline(t0, t1):
            if t1 <= t0:
                return 0.0
            total = 0.0
            for a, b, Ne in zip(skyline_boundaries[:-1], skyline_boundaries[1:], skyline_pop_sizes):
                left, right = max(t0, a), min(t1, b)
                if right > left:
                    total += (right - left) / Ne
            # beyond last boundary: use last Ne value
            if t1 > skyline_boundaries[-1]:
                left = max(t0, skyline_boundaries[-1])
                total += (t1 - left) / skyline_pop_sizes[-1]
            return total

        try:
            tree = read_tree_flexible(row["sim_tree_path"], tree_index=row["tree_index"])
        except Exception as e:
            print(f"Warning: skipping T{row['tree_index']} {pop_model}: {e}")
            continue

        # assign node IDs by preorder
        node_ids = {}
        internal_counter = 0
        for clade in tree.find_clades(order="preorder"):
            if clade.name:
                node_ids[id(clade)] = clade.name
            else:
                node_ids[id(clade)] = f"internal_{internal_counter}"
                internal_counter += 1

        depths      = tree.depths()
        tips        = tree.get_terminals()
        max_tip_depth = max(depths[tip] for tip in tips)
        node_time   = {node: max_tip_depth - depths[node] for node in tree.find_clades()}

        # build event list from simulated tree only
        events = [(node_time[tip], "sample", tip) for tip in tips]
        events += [(node_time[n], "coal", n) for n in tree.get_nonterminals()]
        events.sort(key=lambda x: x[0])

        # single sweep: accumulate H and record rates for all three Ne functions
        H_true = H_constcoal = H_skyline = 0.0
        k = 0
        current_time = 0.0
        ci_true = {}
        ci_constcoal = {}
        ci_skyline = {}
        rates_true = {}
        rates_constcoal = {}
        rates_skyline = {}
        node_times_out = {}

        i = 0
        while i < len(events):
            t = events[i][0]

            if t > current_time and k >= 2:
                coef = comb(k, 2)
                H_true      += coef * _integral_true(current_time, t)
                H_constcoal += coef * _integral_constcoal(current_time, t)
                H_skyline   += coef * _integral_skyline(current_time, t)

            same_time_events = []
            while i < len(events) and events[i][0] == t:
                same_time_events.append(events[i])
                i += 1

            for _, event_type, node in same_time_events:
                if event_type == "coal":
                    node_name = node_ids[id(node)]
                    ci_true[node_name]         = H_true
                    ci_constcoal[node_name]    = H_constcoal
                    ci_skyline[node_name]      = H_skyline
                    rates_true[node_name]      = _rate_true(t)
                    rates_constcoal[node_name] = 1.0 / constcoal_Ne
                    rates_skyline[node_name]   = _rate_skyline(t)
                    node_times_out[node_name]  = t

            for _, event_type, node in same_time_events:
                if event_type == "sample":
                    k += 1
                elif event_type == "coal":
                    k -= max(len(node.clades) - 1, 1)

            current_time = t

        for node in ci_true:
            for model, ci_est, rates_est in [
                ("constcoal", ci_constcoal, rates_constcoal),
                ("skyline",   ci_skyline,   rates_skyline),
            ]:
                results.append({
                    **meta,
                    "model":          model,
                    "node":           node,
                    "node_time":      node_times_out[node],
                    "ci_true":        ci_true[node],
                    "ci_estimated":   ci_est[node],
                    "rate_true":      rates_true[node],
                    "rate_estimated": rates_est[node],
                })

    df = pd.DataFrame(results)
    df["ci_diff"]         = df["ci_estimated"] - df["ci_true"]
    df["ci_rel_error"]    = df["ci_diff"] / df["ci_true"]
    df["ci_abs_diff"]     = df["ci_diff"].abs()
    df["ci_abs_rel_error"] = df["ci_rel_error"].abs()
    df["rate_diff"]         = df["rate_estimated"] - df["rate_true"]
    df["rate_rel_error"]    = df["rate_diff"] / df["rate_true"]
    df["rate_abs_diff"]     = df["rate_diff"].abs()
    df["rate_abs_rel_error"] = df["rate_rel_error"].abs()
    return df


# =============================================================================
# Loaders for compute_errors.py outputs
# =============================================================================

def load_pop_summaries(eval_dir):
    """
    Read all *_pop_summary.pkl files from eval_dir, concatenate, and merge
    with the scenario TSVs (high/med/low.tsv) to add population parameters
    (present_pop_size, growth_rate, bottleneck_size, bottleneck_start,
    bottleneck_end) per replicate.
    """
    eval_dir = Path(eval_dir)

    # Load all pop summaries
    dfs = []
    for pkl in sorted(eval_dir.rglob("*_pop_summary.pkl")):
        dfs.append(pd.read_pickle(pkl))
    if not dfs:
        raise FileNotFoundError(f"No *_pop_summary.pkl files found under {eval_dir}")
    df = pd.concat(dfs, ignore_index=True)
    df = df.drop(columns=["scenario"], errors="ignore")

    # Load all scenario TSVs and concatenate
    pop_params_cols = ["sampling", "population_model", "mutation_signal", "tree_index",
                       "present_pop_size", "growth_rate",
                       "bottleneck_size", "bottleneck_start", "bottleneck_end"]
    tsv_dfs = []
    for tsv in sorted(eval_dir.rglob("*.tsv")):
        if "node_errors" in tsv.name:
            continue
        raw = pd.read_csv(tsv, sep="\t")
        if not all(c in raw.columns for c in pop_params_cols):
            continue
        tsv_dfs.append(raw[pop_params_cols])
    if tsv_dfs:
        params_df = pd.concat(tsv_dfs, ignore_index=True).drop_duplicates()
        df = df.merge(params_df, on=["sampling", "population_model", "mutation_signal", "tree_index"], how="inner")

    return df


def load_node_errors(eval_dir):
    """
    Read all *_node_errors.tsv files from eval_dir and concatenate into one
    DataFrame matching the structure of tree_metrics_combined, with separate
    columns for population_model, mutation_signal, and sampling parsed from
    the scenario string.
    """
    eval_dir = Path(eval_dir)
    dfs = []
    for tsv in sorted(eval_dir.rglob("*_node_errors.tsv")):
        dfs.append(pd.read_csv(tsv, sep="\t"))
    if not dfs:
        raise FileNotFoundError(f"No *_node_errors.tsv files found under {eval_dir}")
    df = pd.concat(dfs, ignore_index=True)

    def parse_scenario(s):
        mutsig = s.split("_")[-1].replace("mutsig", "")
        rest   = s[: s.rfind(f"_{mutsig}mutsig")]
        for sampling in ("independenthomochronous", "linearconstant"):
            if rest.startswith(sampling):
                return sampling, rest[len(sampling) + 1:], mutsig
        return None, rest, mutsig

    parsed = df["scenario"].apply(lambda s: pd.Series(
        parse_scenario(s), index=["sampling", "population_model", "mutation_signal"]
    ))
    df = pd.concat([df.drop(columns=["scenario"]), parsed], axis=1)
    return df


