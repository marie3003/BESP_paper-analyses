from Bio import Nexus, Phylo, SeqIO
from math import comb, exp, inf, ceil, floor, log10
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
    title_map       = {"expgrowthfast": "Exp. Growth (fast)", "expgrowthslow": "Exp. Growth (slow)",
                       "uniform": "Uniform", "bottleneck": "Bottleneck"}

    _cmap_colors = ["#c5dfc0", "#7db87a", "#4e9a8a", "#397398", "#00181E"]
    cmap   = plt.matplotlib.colors.LinearSegmentedColormap.from_list("green_cmap", _cmap_colors)
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
        fig.savefig(save_path, bbox_inches="tight", transparent=True)
    plt.show()


def plot_ci_coverage_heatmap_nodes(tree_metrics_combined, sampling, mode="height",
                                    title=None, figsize=None, save_path=None,
                                    cbar_label="Fraction of nodes"):
    """
    CI coverage contingency heatmap from full tree_metrics_combined (wide format).

    mode="height": node height CI coverage, internal nodes only (tips excluded via NaN)
    mode="bl":     branch length CI coverage, non-root nodes only (root excluded via NaN)
    """
    if mode == "height":
        const_col = "constcoal_height_inside_ci"
        sky_col   = "skyline_height_inside_ci"
        default_title = "Node height CI coverage (internal nodes)"
    elif mode == "bl":
        const_col = "constcoal_bl_inside_ci"
        sky_col   = "skyline_bl_inside_ci"
        default_title = "Branch length CI coverage (non-root nodes)"
    else:
        raise ValueError("mode must be 'height' or 'bl'")

    mutsig_order    = ["low", "med", "high"]
    mutsig_labels   = {"low": "Low mut. signal", "med": "Med. mut. signal", "high": "High mut. signal"}
    pop_model_order = ["uniform", "expgrowthslow", "expgrowthfast", "bottleneck"]
    title_map       = {"expgrowthfast": "Exp. Growth (fast)", "expgrowthslow": "Exp. Growth (slow)",
                       "uniform": "Uniform", "bottleneck": "Bottleneck"}

    _cmap_colors = ["#c5dfc0", "#7db87a", "#4e9a8a", "#397398", "#00181E"]
    cmap_r = plt.matplotlib.colors.LinearSegmentedColormap.from_list("green_cmap", _cmap_colors).reversed()

    df = tree_metrics_combined[tree_metrics_combined["sampling"] == sampling].copy()
    df = df.dropna(subset=[const_col, sky_col])
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
            n = max(len(subset), 1)

            matrix = pd.DataFrame(
                [[subset["only skyline"].sum(), subset["both"].sum()],
                 [subset["neither"].sum(),      subset["only constcoal"].sum()]],
                index=["Skyline ✓", "Skyline ✗"],
                columns=["Constcoal ✗", "Constcoal ✓"],
                dtype=float
            ) / n

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

    suptitle = title if title is not None else f"{default_title} — {sampling}"
    fig.suptitle(suptitle, fontsize=14)
    plt.tight_layout(rect=[0, 0, 0.91, 0.95])
    sm      = plt.cm.ScalarMappable(cmap=cmap_r, norm=plt.Normalize(vmin=0, vmax=100))
    cbar_ax = fig.add_axes([0.91, 0.15, 0.02, 0.65])
    cbar    = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_ticks([0, 25, 50, 75, 100])
    cbar.set_label(cbar_label, fontsize=13)
    cbar.ax.tick_params(labelsize=12)

    if save_path is not None:
        fig.savefig(save_path, bbox_inches="tight", transparent=True)
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
        ax.hlines(const_median, 0, t_max, color=color_const, linestyle=':', label="Constant Estimate", linewidth = 2, zorder = 5)
        ax.fill_between(t_vals, const_lower, const_upper, color=color_const, alpha=0.3, zorder = 4)

        ax.plot(t_vals, skyline_median, color=color_sky, linestyle='--', label="Skyline Estimate", linewidth = 2, zorder = 5)
        ax.fill_between(t_vals, skyline_lower, skyline_upper, color=color_sky, alpha=0.3, zorder = 4)
    elif mode == 'skyline':
        ax.plot(t_vals, skyline_median, color=color_sky, linestyle='--', label="Skyline Estimate", linewidth = 2, zorder = 5)
        ax.fill_between(t_vals, skyline_lower, skyline_upper, color=color_sky, alpha=0.3, zorder = 4)
    elif mode == 'constcoal':
        ax.hlines(const_median, 0, t_max, color=color_const, linestyle=':', label="Constant Estimate", linewidth = 2, zorder = 5)
        ax.fill_between(t_vals, const_lower, const_upper, color=color_const, alpha=0.3, zorder = 4)
    
    if y_range:
        ax.set_ylim(y_range)
    
    ax.set_xlabel("Time before present")
    ax.set_ylabel("Population Size")
    ax.invert_xaxis()

def _add_skyline_changepoint_strip(ax, subset, x_range=None, color_sky="#80557e",
                                    dot_size=25, grey_alpha=0.3, y_frac=0.9):
    """Plot skyline interval change points near the top of ax."""
    all_times_arr = np.array(list(subset["skyline_times"]))[:, :-1]  # drop last (root height)

    ylims = ax.get_ylim()
    y_pos = ylims[0] + (ylims[1] - ylims[0]) * y_frac

    for row in all_times_arr:
        ax.scatter(row, np.full(len(row), y_pos), color="darkgrey", s=dot_size, alpha=grey_alpha, zorder=4, linewidths=0, clip_on=True)

    medians = np.median(all_times_arr, axis=0)
    ax.scatter(medians, np.full(len(medians), y_pos), color=color_sky, s=dot_size * 1.5, zorder=5, linewidths=0, clip_on=True)


def plot_population_summary(path_info_df, sampling, x_range=None, title=None, y_range=None,
                            add_samples=False, mode='both', pop_models=None, mut_signals=None,
                            show_changepoints=True, changepoint_dot_size=25, changepoint_grey_alpha=0.3,
                            x_label=None, y_label=None,
                            trajectory_df=None,
                            show_pop_model_labels=True, show_mutsig_labels=True, show_legend=True,
                            figsize=None,
                            label_fontsize=13, tick_fontsize=12,
                            save_path=None):
    """
    Summary plot per condition (mutation signal × population model) for a given sampling type,
    showing median and 95% CI for skyline and constant-coalescent estimates.
    """
    if pop_models is None:
        pop_models = sorted(path_info_df["population_model"].unique())
    if mut_signals is None:
        mut_signals = ["low", "med", "high"]
    if isinstance(pop_models, str):
        pop_models = [pop_models]
    if isinstance(mut_signals, str):
        mut_signals = [mut_signals]
    ncols, nrows = len(pop_models), len(mut_signals)

    title_map   = {"expgrowthfast": "Exp. Growth (fast)", "expgrowthslow": "Exp. Growth (slow)",
                   "uniform": "Uniform", "bottleneck": "Bottleneck"}
    mutsig_labels = {"low": "Low mut. signal", "med": "Med. mut. signal", "high": "High mut. signal"}

    df = path_info_df[path_info_df["sampling"] == sampling]
    time_horizon = x_range[1] if x_range else 0

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize or (5 * ncols, 3.5 * nrows), squeeze=False)

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

            if trajectory_df is not None:
                traj_sub = trajectory_df[
                    (trajectory_df["sampling"] == sampling) &
                    (trajectory_df["mutation_signal"] == mut_sig) &
                    (trajectory_df["population_model"] == pop_model)
                ]
                if not traj_sub.empty:
                    t_vals = traj_sub["time"].values
                    color_exp  = "#a6444f"
                    color_sky  = "#80557e"
                    color_const = "#397398"
                    t_max = x_range[1] if x_range else t_vals.max()
                    mask = t_vals <= t_max
                    t_plot = t_vals[mask]
                    t_smooth = np.linspace(0, t_max, 1000)
                    N_true = np.array([true_traj(t) for t in t_smooth])
                    ax.plot(t_smooth, N_true, color=color_exp, linewidth=2, label="True Population Size", zorder=6)
                    if mode in ('both', 'skyline'):
                        ax.plot(t_plot, traj_sub["skyline_median"].values[mask],
                                color=color_sky, linestyle='--', label="Skyline Estimate", linewidth=2, zorder=5)
                        ax.fill_between(t_plot,
                                        traj_sub["skyline_hpd_lower"].values[mask],
                                        traj_sub["skyline_hpd_upper"].values[mask],
                                        color=color_sky, alpha=0.3, zorder=4)
                    if mode in ('both', 'constcoal'):
                        ax.plot(t_plot, traj_sub["constcoal_median"].values[mask],
                                color=color_const, linestyle=':', label="Constant Estimate", linewidth=2, zorder=5)
                        ax.fill_between(t_plot,
                                        traj_sub["constcoal_hpd_lower"].values[mask],
                                        traj_sub["constcoal_hpd_upper"].values[mask],
                                        color=color_const, alpha=0.3, zorder=4)
                    if y_range:
                        ax.set_ylim(y_range)
                    ax.invert_xaxis()
            else:
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

            if show_changepoints:
                _add_skyline_changepoint_strip(ax, subset, x_range=x_range,
                                               dot_size=changepoint_dot_size,
                                               grey_alpha=changepoint_grey_alpha)

            ax.tick_params(labelsize=tick_fontsize)
            ax.set_xlabel(x_label or "Time before present", fontsize=label_fontsize)
            ax.set_ylabel(y_label or "Population Size", fontsize=label_fontsize)
            if i == 0 and show_pop_model_labels:
                ax.set_title(title_map.get(pop_model, pop_model), fontsize=15)

    handles, labels = axes[0][0].get_legend_handles_labels()
    if add_samples:
        if mode == 'skyline':
            handles.append(Line2D([0], [0], color='lightgrey', linestyle='--'))
            labels.append('Sample Skyline Trajectories')
        elif mode == 'constcoal':
            handles.append(Line2D([0], [0], color='lightgrey', linestyle=':'))
            labels.append('Sample Constant Trajectories')
    if show_changepoints:
        handles.append(Line2D([0], [0], marker='o', color='w', markerfacecolor='darkgrey',
                               markersize=8, label='Skyline change points (per tree)'))
        labels.append('Skyline change points (per tree)')
        handles.append(Line2D([0], [0], marker='o', color='w', markerfacecolor='#80557e',
                               markersize=10, label='Median change point'))
        labels.append('Median change point')
    if show_legend:
        fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(0.99, 0.99),
                   ncol=len(handles), fontsize=14, frameon=False,
                   handletextpad=0.4, columnspacing=1.2)

    suptitle = title if title is not None else sampling
    plt.suptitle(suptitle, fontsize=17)
    top = 0.96 if suptitle else 1.0
    plt.tight_layout(rect=[0, 0, 1, top])

    # Mutsig row labels left of the y-axis, matching plot_height_vs_popsize_error style
    if show_mutsig_labels:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        inv = fig.transFigure.inverted()
        for i, mut_sig in enumerate(mut_signals):
            ax = axes[i][0]
            ylabel_bb = ax.yaxis.label.get_window_extent(renderer)
            x_fig = inv.transform((ylabel_bb.x0, ylabel_bb.y0 + ylabel_bb.height / 2))[0]
            pos = ax.get_position()
            y_fig = (pos.y0 + pos.y1) / 2
            fig.text(x_fig - 0.01, y_fig, mutsig_labels.get(mut_sig, mut_sig),
                     fontsize=15, ha="center", va="center", rotation=90,
                     transform=fig.transFigure)

    if save_path is not None:
        fig.savefig(save_path, bbox_inches="tight", transparent=True)
    else:
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
              alpha=alpha, label = "Constant Estimate")
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
                         violin=False, legend_inside=False,
                         show_pop_model_labels=True,
                         label_fontsize=14, tick_fontsize=14, legend_fontsize=14,
                         show_diff=False, diff_mode="relative", diff_y_range=None):
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
    diff_color  = "#338778"
    labels      = {"constcoal": "Const. coal.", "skyline": "Skyline"}
    mutsig_order    = ["low", "med", "high"]
    mutsig_labels   = {"low": "Low", "med": "Med.", "high": "High"}
    title_map       = {"expgrowthfast": "Exp. Growth (fast)", "expgrowthslow": "Exp. Growth (slow)",
                       "uniform": "Uniform", "bottleneck": "Bottleneck"}
    pop_model_order = pop_models or ["uniform", "expgrowthslow", "expgrowthfast", "bottleneck"]

    df = root_height_df_wide[root_height_df_wide["sampling"] == sampling].copy()

    # x positions: each mutsig gets a slot of width 1; constcoal left, skyline right
    offsets = {"constcoal": -0.12, "skyline": 0.12}
    diff_offset = 0.36
    x_base  = {ms: i * 1.5 for i, ms in enumerate(mutsig_order)}

    ncols = len(pop_model_order)
    fig, axes = plt.subplots(1, ncols, figsize=figsize or (4.5 * ncols, 4.5))
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

        # third dot: constcoal - skyline difference per tree
        col_cc = f"constcoal_{y_col}"
        col_sk = f"skyline_{y_col}"
        use_twin = show_diff and diff_mode == "relative"
        ax_diff = ax.twinx() if use_twin else None
        if ax_diff is not None:
            ax_diff.tick_params(axis="y", labelsize=tick_fontsize, colors=diff_color)
            ax_diff.spines["right"].set_color(diff_color)
            if diff_y_range:
                ax_diff.set_ylim(diff_y_range)

        if show_diff and col_cc in df.columns and col_sk in df.columns:
            target = ax_diff if use_twin else ax
            for ms in mutsig_order:
                sub = df[(df["population_model"] == pop_model) & (df["mutation_signal"] == ms)].dropna(subset=[col_cc, col_sk])
                diff_vals = ((sub[col_cc] - sub[col_sk]) / sub[col_sk]).values if diff_mode == "relative" else (sub[col_cc] - sub[col_sk]).values
                x_center = x_base[ms] + diff_offset
                jitter = rng.uniform(-0.06, 0.06, size=len(diff_vals))
                target.scatter(x_center + jitter, diff_vals, color="#aaaaaa", s=15,
                               alpha=0.8, zorder=2, linewidths=0)
                if len(diff_vals) >= 2:
                    hpd_lo, hpd_hi = hpd(np.sort(diff_vals))
                    target.plot([x_center, x_center], [hpd_lo, hpd_hi],
                                color=diff_color, linewidth=1.5, zorder=3, solid_capstyle="round")
                    cap = 0.05
                    target.plot([x_center - cap, x_center + cap], [hpd_lo, hpd_lo],
                                color=diff_color, linewidth=1.5, zorder=3)
                    target.plot([x_center - cap, x_center + cap], [hpd_hi, hpd_hi],
                                color=diff_color, linewidth=1.5, zorder=3)
                    if violin:
                        kde = gaussian_kde(diff_vals, bw_method=0.4)
                        ys = np.linspace(hpd_lo, hpd_hi, 200)
                        xs_kde = kde(ys)
                        xs_kde = xs_kde / xs_kde.max() * 0.12
                        target.fill_betweenx(ys, x_center - xs_kde, x_center + xs_kde,
                                             color=diff_color, alpha=0.35, zorder=3, linewidth=0)
                target.scatter(x_center, np.median(diff_vals), color=diff_color, s=75,
                               zorder=4, edgecolors="white", linewidths=0.5)

        tick_positions = [x_base[ms] + (diff_offset / 2 if show_diff else 0) for ms in mutsig_order]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels([mutsig_labels[ms] for ms in mutsig_order], fontsize=tick_fontsize)
        ax.set_xlabel("Mutation signal", fontsize=label_fontsize)
        if show_pop_model_labels:
            ax.set_title(title_map.get(pop_model, pop_model), fontsize=label_fontsize + 1)
        ax.xaxis.grid(False)
        if y_range:
            ax.set_ylim(y_range)
        ax.yaxis.grid(True, color='lightgray', linewidth=0.8, alpha=0.8, zorder=5)
        ax.axhline(0, color='gray', linestyle='--', linewidth=0.8)
        ax.tick_params(labelsize=tick_fontsize, labelleft=True)
        ax.set_ylabel(y_label or y_col, fontsize=label_fontsize + 1, labelpad=2)
        if use_twin:
            ax_diff.set_ylabel("Rel. diff.", fontsize=label_fontsize + 1, labelpad=4, color=diff_color)
    from matplotlib.lines import Line2D as _L2D
    handles = [_L2D([], [], marker="o", linestyle="none", markerfacecolor=c,
                    markeredgecolor="none", markersize=9, label=labels[m])
               for m, c in colors.items()]
    if show_diff and not (show_diff and diff_mode == "relative"):
        handles.append(_L2D([], [], marker="o", linestyle="none", markerfacecolor=diff_color,
                            markeredgecolor="none", markersize=9, label="Diff."))
    _leg_kwargs = dict(handles=handles, frameon=False, fontsize=legend_fontsize,
                       handlelength=0.0, handletextpad=0.6, columnspacing=0.8)
    suptitle = title if title is not None else sampling
    plt.suptitle(suptitle, fontsize=17)
    top = 0.93 if suptitle else 1.0
    plt.tight_layout(rect=[0, 0, 1, top])
    plt.subplots_adjust(wspace=0.35)

    if legend_inside:
        axes[-1].legend(loc="lower right", ncol=2, bbox_to_anchor=(1.0, -0.038), **_leg_kwargs)
    else:
        fig.legend(loc="upper center", bbox_to_anchor=(0.5, 1.01),
                   ncol=3, **_leg_kwargs)

    if save_path is not None:
        fig.savefig(save_path, bbox_inches="tight", transparent=True)
    from IPython.display import display as _display
    _display(fig)
    plt.close(fig)



def plot_height_vs_popsize_error(root_height_df_wide, sampling,
                                  x_col="pop_diff_median", y_col="height_diff_median",
                                  x_range=None, y_range=None, alpha=0.6, s=10,
                                  title=None, x_label=None, y_label=None,
                                  pop_models=None, mutsigs=None,
                                  figsize=None, show_zero_lines=True,
                                  legend_inside=True,
                                  color_by_time=False,
                                  time_col="height_sim",
                                  log_colorbar=False,
                                  cbar_label_constcoal=None,
                                  cbar_label_skyline=None,
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
    title_map     = {"expgrowthfast": "Exp. Growth (fast)", "expgrowthslow": "Exp. Growth (slow)",
                     "uniform": "Uniform", "bottleneck": "Bottleneck"}
    legend_labels = {"constcoal": "Const. coal.", "skyline": "Skyline"}
    colors        = {"constcoal": "#397398", "skyline": "#80557e"}
    markers       = {"constcoal": "o",       "skyline": "^"}
    model_cmaps   = {
        "constcoal": plt.matplotlib.colors.LinearSegmentedColormap.from_list(
                         "constcoal_time", ["#b5d2f2", "#7394c2", "#397398"]),
        "skyline":   plt.matplotlib.colors.LinearSegmentedColormap.from_list(
                         "skyline_time",   ["#a6444f", "#d991b4", "#80557e"]),
    }

    if color_by_time:
        df_all = root_height_df_wide[root_height_df_wide["sampling"] == sampling]
        t_min_pos = max(df_all[time_col].min(), 1)
        t_max     = df_all[time_col].max()
        cnorm = plt.matplotlib.colors.LogNorm(vmin=t_min_pos, vmax=t_max) if log_colorbar \
                else plt.Normalize(vmin=t_min_pos, vmax=t_max)

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

                if color_by_time:
                    t_vals = subset[time_col].clip(lower=t_min_pos)
                    point_colors = [model_cmaps[model](cnorm(t)) for t in t_vals]
                else:
                    point_colors = [color] * len(xm)

                if has_x_ci or has_y_ci:
                    xerr = ([xm - subset[xl_key], subset[xh_key] - xm] if has_x_ci else None)
                    yerr = ([ym - subset[yl_key], subset[yh_key] - ym] if has_y_ci else None)
                    for idx, (xi, yi, c) in enumerate(zip(xm, ym, point_colors)):
                        xe = ([xm.iloc[idx] - subset[xl_key].iloc[idx]],
                              [subset[xh_key].iloc[idx] - xm.iloc[idx]]) if has_x_ci else (None, None)
                        ye = ([ym.iloc[idx] - subset[yl_key].iloc[idx]],
                              [subset[yh_key].iloc[idx] - ym.iloc[idx]]) if has_y_ci else (None, None)
                        eb = ax.errorbar(xi, yi,
                                         xerr=([[xe[0][0]], [xe[1][0]]] if xe[0] is not None else None),
                                         yerr=([[ye[0][0]], [ye[1][0]]] if ye[0] is not None else None),
                                         fmt=markers[model], color=c, ecolor=c,
                                         elinewidth=0.8, capsize=2, markersize=4,
                                         label=legend_labels[model] if idx == 0 else "_nolegend_")
                        eb[0].set_alpha(alpha)
                        for line in eb[1]: line.set_alpha(0.15)
                        for line in eb[2]: line.set_alpha(0.15)
                else:
                    ax.scatter(xm, ym, color=point_colors, marker=markers[model],
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

    if color_by_time:
        fig.subplots_adjust(right=0.86)
        cbar_const = fig.colorbar(plt.cm.ScalarMappable(cmap=model_cmaps["constcoal"], norm=cnorm),
                                   cax=fig.add_axes([0.88, 0.15, 0.015, 0.7]))
        cbar_const.set_label(cbar_label_constcoal or f"True root height — Const. Coalescent",
                             rotation=270, labelpad=18, fontsize=14)
        cbar_const.ax.tick_params(labelsize=13)
        cbar_sky = fig.colorbar(plt.cm.ScalarMappable(cmap=model_cmaps["skyline"], norm=cnorm),
                                cax=fig.add_axes([0.97, 0.15, 0.015, 0.7]))
        cbar_sky.set_label(cbar_label_skyline or f"True root height — Skyline",
                           rotation=270, labelpad=18, fontsize=14)
        cbar_sky.ax.tick_params(labelsize=13)

    if save_path is not None:
        fig.savefig(save_path, bbox_inches="tight", transparent=True)
    plt.show()


def plot_constcoal_vs_skyline_error(root_height_df_wide, sampling,
                                     error_col,
                                     error_label=None,
                                     x_range=None, y_range=None,
                                     alpha=0.5, s=10,
                                     title=None, x_label=None, y_label=None,
                                     pop_models=None, mutsigs=None,
                                     figsize=None,
                                     color_by_time=False,
                                     time_col="height_sim",
                                     log_colorbar=False,
                                     cbar_label=None,
                                     show_pop_model_labels=True,
                                     show_mutsig_labels=True,
                                     label_fontsize=14, tick_fontsize=13,
                                     legend_fontsize=11, cbar_tick_fontsize=10,
                                     save_path=None):
    """
    Scatter plot of constcoal error vs skyline error per replicate.
    Rows = mutation signal, columns = population model.
    Each point is one replicate; error bars show HPD intervals.
    The x=y diagonal is drawn as a black dashed line.

    error_col : str
        Base column name without model prefix, e.g. 'height_diff_median'.
        HPD columns are inferred by replacing '_median' with '_hpd_lower/upper'.
    """
    df = root_height_df_wide[root_height_df_wide["sampling"] == sampling].copy()

    mutsig_order       = mutsigs    if mutsigs    is not None else ["low", "med", "high"]
    growth_model_order = pop_models if pop_models is not None else ["uniform", "expgrowthslow", "expgrowthfast", "bottleneck"]
    if isinstance(mutsig_order, str):       mutsig_order       = [mutsig_order]
    if isinstance(growth_model_order, str): growth_model_order = [growth_model_order]

    mutsig_labels = {"low": "Low mut. signal", "med": "Med. mut. signal", "high": "High mut. signal"}
    title_map     = {"expgrowthfast": "Exp. Growth (fast)", "expgrowthslow": "Exp. Growth (slow)",
                     "uniform": "Uniform", "bottleneck": "Bottleneck"}
    color = "#555555"
    _cmap_colors = ["#c5dfc0", "#7db87a", "#4e9a8a", "#397398", "#00181E"]
    time_cmap = plt.matplotlib.colors.LinearSegmentedColormap.from_list(
                    "green_cmap", _cmap_colors)

    # per-subplot norms computed inside the loop when color_by_time=True

    def _key(model, col):
        full = f"{model}_{col}"
        return full if full in df.columns else None

    med_col  = error_col if "_median" in error_col else f"{error_col}_median"
    lo_col   = med_col.replace("_median", "_hpd_lower")
    hi_col   = med_col.replace("_median", "_hpd_upper")

    nrows, ncols = len(mutsig_order), len(growth_model_order)
    default_figsize = (4 * ncols, 3.5 * nrows)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                             figsize=figsize or default_figsize,
                             squeeze=False)

    for i, mutsig in enumerate(mutsig_order):
        for j, pop_model in enumerate(growth_model_order):
            ax = axes[i, j]
            subset = df[(df["mutation_signal"] == mutsig) & (df["population_model"] == pop_model)]

            if subset.empty:
                ax.axis("off")
                continue

            xm_key = _key("skyline",   med_col)
            ym_key = _key("constcoal", med_col)
            xl_key = _key("skyline",   lo_col)
            xh_key = _key("skyline",   hi_col)
            yl_key = _key("constcoal", lo_col)
            yh_key = _key("constcoal", hi_col)

            if xm_key is None or ym_key is None:
                ax.axis("off")
                continue

            xm = subset[xm_key]
            ym = subset[ym_key]
            has_x_ci = xl_key is not None and xh_key is not None
            has_y_ci = yl_key is not None and yh_key is not None

            if color_by_time:
                t_min_pos = max(subset[time_col].min(), 1)
                t_max     = subset[time_col].max()
                cnorm = plt.matplotlib.colors.LogNorm(vmin=t_min_pos, vmax=t_max) if log_colorbar \
                        else plt.Normalize(vmin=t_min_pos, vmax=t_max)
                t_vals = subset[time_col].clip(lower=t_min_pos)
                point_colors = [time_cmap(cnorm(t)) for t in t_vals]
            else:
                cnorm = None
                point_colors = [color] * len(xm)

            if has_x_ci or has_y_ci:
                for idx, (xi, yi, c) in enumerate(zip(xm, ym, point_colors)):
                    xe = ([[xm.iloc[idx] - subset[xl_key].iloc[idx]],
                            [subset[xh_key].iloc[idx] - xm.iloc[idx]]] if has_x_ci else None)
                    ye = ([[ym.iloc[idx] - subset[yl_key].iloc[idx]],
                            [subset[yh_key].iloc[idx] - ym.iloc[idx]]] if has_y_ci else None)
                    eb = ax.errorbar(xi, yi, xerr=xe, yerr=ye,
                                     fmt="o", color=c, ecolor=c,
                                     elinewidth=0.8, capsize=2, markersize=4, alpha=alpha)
                    for line in eb[1]: line.set_alpha(0.15)
                    for line in eb[2]: line.set_alpha(0.15)
            else:
                ax.scatter(xm, ym, color=point_colors, s=s, alpha=alpha)

            # x = y diagonal
            ax.axline((0, 0), slope=1, color="black", linestyle="--", linewidth=1.0, zorder=1)

            ax.set_aspect("equal")
            ax.grid(True, color="lightgray", linewidth=0.4, alpha=0.5, zorder=0)
            if x_range: ax.set_xlim(x_range)
            if y_range: ax.set_ylim(y_range)
            ax.tick_params(labelsize=tick_fontsize)
            ax.set_xlabel(x_label or f"Skyline — {med_col}", fontsize=label_fontsize)
            ax.set_ylabel(y_label or f"Const. Coalescent — {med_col}", fontsize=label_fontsize, labelpad=-5 if color_by_time else 8)
            if i == 0 and show_pop_model_labels:
                ax.set_title(title_map.get(pop_model, pop_model), fontsize=label_fontsize + 1)

            marker_line = Line2D([0], [0], marker="o", color=color, markersize=5,
                                 linestyle="none", alpha=1.0)
            caplines    = (Line2D([0], [0], color=color, linewidth=1.0),
                           Line2D([0], [0], color=color, linewidth=1.0))
            barlines    = (LineCollection([[[0, -0.3], [0, 0.3]]], colors=[color], linewidths=[0.8]),
                           LineCollection([[[-0.3, 0], [0.3, 0]]], colors=[color], linewidths=[0.8]))
            handle = ErrorbarContainer((marker_line, caplines, barlines),
                                       has_xerr=True, has_yerr=True, label=error_label or med_col)
            ax.legend(handles=[handle], fontsize=legend_fontsize, frameon=False,
                      loc="upper left", borderpad=0.3, handletextpad=0.6,
                      handlelength=0.5, bbox_to_anchor=(0.0, 1.0))

            if color_by_time and cnorm is not None:
                cbar = fig.colorbar(plt.cm.ScalarMappable(cmap=time_cmap, norm=cnorm), ax=ax,
                                    fraction=0.046, pad=0.04, aspect=20)
                if log_colorbar:
                    _cbar_ticks = {
                        "uniform":       ([1e3, 1e4],   ["$10^3$", "$10^4$"]),
                        "bottleneck":    ([1e3, 5e3],   ["$10^3$", "$5\\times10^3$"]),
                        "expgrowthslow": ([3e2, 4e2],   ["$3\\times10^2$", "$4\\times10^2$"]),
                        "expgrowthfast": ([2e2, 2.5e2], ["$2\\times10^2$", "$2.5\\times10^2$"]),
                    }
                    ticks, labels = _cbar_ticks.get(pop_model, ([t_min_pos, t_max], [f"{t_min_pos:.0f}", f"{t_max:.0f}"]))
                    cbar.set_ticks(ticks)
                    cbar.set_ticklabels(labels)
                    cbar.ax.yaxis.set_minor_formatter(plt.matplotlib.ticker.NullFormatter())
                cbar.ax.tick_params(labelsize=cbar_tick_fontsize)
                if i == 0 and j == len(growth_model_order) - 1:
                    cbar.set_label(cbar_label or "True root age", rotation=270,
                                   labelpad=20, fontsize=legend_fontsize)

    suptitle = title if title is not None else sampling
    plt.suptitle(suptitle, fontsize=17)
    top = 0.97 if suptitle else 1.0
    plt.tight_layout(rect=[0, 0, 1, top])

    if show_mutsig_labels:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        inv = fig.transFigure.inverted()
        for i, mutsig in enumerate(mutsig_order):
            ax = axes[i][0]
            ylabel_bb = ax.yaxis.label.get_window_extent(renderer)
            x_fig = inv.transform((ylabel_bb.x0, ylabel_bb.y0 + ylabel_bb.height / 2))[0]
            pos = ax.get_position()
            y_fig = (pos.y0 + pos.y1) / 2
            fig.text(x_fig - (0.025 if color_by_time else 0.01), y_fig, mutsig_labels.get(mutsig, mutsig),
                     fontsize=label_fontsize + 1, ha="center", va="center", rotation=90,
                     transform=fig.transFigure)

    if save_path is not None:
        fig.savefig(save_path, bbox_inches="tight", transparent=True)
    from IPython.display import display as _display
    _display(fig)
    plt.close(fig)


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



def plot_binned_error_timecoloring(time_bin_df, sampling,
                                   x_col="pop_rel_error",
                                   y_col="bl_rel_error",
                                   x_range=None, y_range=None,
                                   show_errorbars=True,
                                   min_nodes_per_bin=10,
                                   log_colorbar=False,
                                   time_col="bin_center",
                                   x_label=None, y_label=None,
                                   cbar_label="Time before present",
                                   cbar_label_constcoal=None, cbar_label_skyline=None,
                                   pop_models=None, mutsigs=None,
                                   figsize=None,
                                   per_tree_df=None,
                                   show_pop_model_labels=True,
                                   show_mutsig_labels=True,
                                   label_fontsize=14, tick_fontsize=13,
                                   legend_fontsize=13, cbar_tick_fontsize=13,
                                   title=None, save_path=None):
    """
    Scatter plot of y_col vs x_col per pre-computed time bin, colored by time.
    Rows = mutation signal, columns = population model. Each point is one bin;
    error bars show the pre-computed HPD. Bins with fewer than min_nodes_per_bin
    nodes are excluded. Two colorbars (constcoal blue, skyline red→purple).
    """
    df = time_bin_df[time_bin_df["sampling"] == sampling].copy()
    if per_tree_df is not None:
        per_tree_df = per_tree_df[per_tree_df["sampling"] == sampling].copy()
        if "bin_center" not in per_tree_df.columns:
            per_tree_df["bin_center"] = (per_tree_df["bin_lower"] + per_tree_df["bin_upper"]) / 2

    mutsig_order       = mutsigs    if mutsigs    is not None else ["low", "med", "high"]
    growth_model_order = pop_models if pop_models is not None else ["uniform", "expgrowthslow", "expgrowthfast", "bottleneck"]
    if isinstance(mutsig_order, str):       mutsig_order       = [mutsig_order]
    if isinstance(growth_model_order, str): growth_model_order = [growth_model_order]

    title_map     = {"expgrowthfast": "Exp. Growth (fast)", "expgrowthslow": "Exp. Growth (slow)",
                     "uniform": "Uniform", "bottleneck": "Bottleneck"}
    mutsig_labels = {"low": "Low mut. signal", "med": "Med. mut. signal", "high": "High mut. signal"}
    model_markers = {"constcoal": "o", "skyline": "^"}
    model_cmaps   = {
        "constcoal": plt.matplotlib.colors.LinearSegmentedColormap.from_list(
                         "constcoal_time", ["#b5d2f2", "#7394c2", "#397398"]),
        "skyline":   plt.matplotlib.colors.LinearSegmentedColormap.from_list(
                         "skyline_time",   ["#f2c4d8", "#d991b4", "#80557e"]),
    }

    x_med = f"{x_col}_median"
    x_lo  = f"{x_col}_hpd_lower"
    x_hi  = f"{x_col}_hpd_upper"
    y_med = f"{y_col}_median"
    y_lo  = f"{y_col}_hpd_lower"
    y_hi  = f"{y_col}_hpd_upper"
    has_err = x_lo in df.columns and y_lo in df.columns

    t_min = df[time_col].min()
    t_max = df[time_col].max()
    t_min_pos = max(t_min, 1)
    cnorm = plt.matplotlib.colors.LogNorm(vmin=t_min_pos, vmax=t_max) if log_colorbar \
            else plt.Normalize(vmin=t_min_pos, vmax=t_max)

    nrows, ncols = len(mutsig_order), len(growth_model_order)
    default_figsize = (4.5 * ncols, 3.5 * nrows)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                             figsize=figsize or default_figsize,
                             squeeze=False)

    for i, mutsig in enumerate(mutsig_order):
        for j, pop_model in enumerate(growth_model_order):
            ax = axes[i, j]
            df_sub = df[(df["mutation_signal"] == mutsig) & (df["population_model"] == pop_model)]

            for model in ("constcoal", "skyline"):
                df_m = df_sub[df_sub["model"] == model].copy()
                df_m = df_m[df_m["n_nodes"] >= min_nodes_per_bin]
                df_m = df_m.dropna(subset=[x_med, y_med])

                if per_tree_df is not None:
                    pt_sub = per_tree_df[
                        (per_tree_df["mutation_signal"] == mutsig) &
                        (per_tree_df["population_model"] == pop_model) &
                        (per_tree_df["model"] == model)
                    ].dropna(subset=[x_med, y_med])
                    bg_marker  = "o" if model == "constcoal" else "^"
                    bg_colors  = [model_cmaps[model](cnorm(max(t, t_min_pos))) for t in pt_sub[time_col]]
                    ax.scatter(pt_sub[x_med], pt_sub[y_med],
                               color=bg_colors, marker=bg_marker, s=8, alpha=0.15, zorder=1, linewidths=0)

                for _, row in df_m.iterrows():
                    t     = max(row[time_col], t_min_pos)
                    color = model_cmaps[model](cnorm(t))
                    xm, ym = row[x_med], row[y_med]

                    if show_errorbars and has_err and per_tree_df is None:
                        ax.errorbar(xm, ym,
                                    xerr=[[xm - row[x_lo]], [row[x_hi] - xm]],
                                    yerr=[[ym - row[y_lo]], [row[y_hi] - ym]],
                                    fmt=model_markers[model], color=color,
                                    ecolor="lightgray", elinewidth=0.8,
                                    capsize=2, markersize=6, alpha=0.85,
                                    zorder=3)
                    else:
                        ax.scatter(xm, ym, color=color,
                                   marker=model_markers[model], s=45,
                                   alpha=1.0, zorder=3,
                                   edgecolors="black", linewidths=0.5)

            ax.axhline(0, color="black", linestyle="--", linewidth=1.5, zorder=1)
            ax.axvline(0, color="black", linestyle="--", linewidth=1.5, zorder=1)
            ax.grid(True, color="lightgray", linewidth=0.4, alpha=0.5, zorder=0)
            if x_range: ax.set_xlim(x_range)
            if y_range: ax.set_ylim(y_range)
            ax.tick_params(labelsize=tick_fontsize)
            ax.set_xlabel(x_label or x_col, fontsize=label_fontsize)
            if i == 0 and show_pop_model_labels:
                ax.set_title(title_map.get(pop_model, pop_model), fontsize=label_fontsize + 1)
            if j == 0:
                ax.set_ylabel(y_label or y_col, fontsize=label_fontsize)
            else:
                ax.set_ylabel("")

    suptitle = title if title is not None else sampling
    plt.suptitle(suptitle, fontsize=17)
    top = 0.97 if suptitle else 1.0
    plt.tight_layout(rect=[0, 0, 0.86, top])

    if show_mutsig_labels:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        inv = fig.transFigure.inverted()
        for i, mutsig in enumerate(mutsig_order):
            ax = axes[i][0]
            ylabel_bb = ax.yaxis.label.get_window_extent(renderer)
            x_fig = inv.transform((ylabel_bb.x0, ylabel_bb.y0 + ylabel_bb.height / 2))[0]
            pos = ax.get_position()
            y_fig = (pos.y0 + pos.y1) / 2
            fig.text(x_fig - 0.01, y_fig, mutsig_labels.get(mutsig, mutsig),
                     fontsize=15, ha="center", va="center", rotation=90,
                     transform=fig.transFigure)

    single = (nrows == 1 and ncols == 1)
    right = 0.76 if single else 0.86
    fig.subplots_adjust(right=right)
    cbar_const = fig.colorbar(plt.cm.ScalarMappable(cmap=model_cmaps["constcoal"], norm=cnorm),
                               cax=fig.add_axes([0.78, 0.15, 0.025, 0.7] if single else [0.88, 0.15, 0.015, 0.7]))
    cbar_sky = fig.colorbar(plt.cm.ScalarMappable(cmap=model_cmaps["skyline"], norm=cnorm),
                             cax=fig.add_axes([0.95, 0.15, 0.025, 0.7] if single else [0.97, 0.15, 0.015, 0.7]))
    if single:
        cbar_const.set_label("Const. coalescent", rotation=270, labelpad=5, fontsize=legend_fontsize)
        cbar_sky.set_label("Skyline", rotation=270, labelpad=5, fontsize=legend_fontsize)
        # shared label above both colorbars
        shared_label = "Time (b.p.)" if cbar_label == "Time before present" else cbar_label
        cbar_sky.ax.figure.text(
            (0.78 + 0.92 + 0.015) / 2 + 0.04, 0.87,
            shared_label, ha="center", va="bottom", fontsize=legend_fontsize,
            transform=fig.transFigure
        )
    else:
        cbar_const.set_label(cbar_label_constcoal or f"{cbar_label} — constcoal", rotation=270, labelpad=18, fontsize=legend_fontsize)
        cbar_sky.set_label(cbar_label_skyline or f"{cbar_label} — skyline", rotation=270, labelpad=18, fontsize=legend_fontsize)
    cbar_const.ax.tick_params(labelsize=cbar_tick_fontsize)
    cbar_sky.ax.tick_params(labelsize=cbar_tick_fontsize)

    if save_path is not None:
        fig.savefig(save_path, bbox_inches="tight", transparent=True)
    from IPython.display import display as _display
    _display(fig)
    plt.close(fig)


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

    title_map = {"expgrowthfast": "Exp. Growth (fast)", "expgrowthslow": "Exp. Growth (slow)",
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

    title_map = {"expgrowthfast": "Exp. Growth (fast)", "expgrowthslow": "Exp. Growth (slow)",
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

    title_map = {"expgrowthfast": "Exp. Growth (fast)", "expgrowthslow": "Exp. Growth (slow)",
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
    title_map = {"expgrowthfast": "Exp. Growth (fast)", "expgrowthslow": "Exp. Growth (slow)",
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



def plot_combined_population_and_error(df_population, df_time_bins, sampling,
                                       pop_y_range=(0, 5000),
                                       error_y_range=(-1, 3),
                                       error_col='bl_rel_error_median',
                                       time_col='bin_lower',
                                       mode='both',
                                       pop_models=None, mutsigs=None,
                                       min_nodes_per_bin=10,
                                       per_tree_df=None,
                                       raw_alpha=0.15,
                                       jitter_width=0.05,
                                       marker_edgecolor="none",
                                       errorbar_alpha=0.4,
                                       x_label="Time before present",
                                       pop_y_label="Population size",
                                       error_y_label=None,
                                       title=None,
                                       figsize=None,
                                       save_path=None,
                                       ylabel_x_positions=None,
                                       error_ylabel_y_offset=-0.08,
                                       sampling_end_time=None,
                                       trajectory_df=None,
                                       show_changepoints=True,
                                       changepoint_dot_size=25,
                                       changepoint_grey_alpha=0.3,
                                       changepoint_y_frac=0.7,
                                       label_fontsize=12, tick_fontsize=12,
                                       legend_fontsize=13,
                                       hspace=0.3, wspace=0.4,
                                       sampling_fontsize=None,
                                       xticks=None,
                                       show_legend=True):
    """
    Combined population size trajectory (top) + pre-binned error (bottom) per mutsig × pop model.

    df_time_bins: output of load_time_bins(), already has bin_lower/bin_upper/bin_center,
                  model column, and pre-computed error medians/HPDs.
    error_col: median column to plot (e.g. 'bl_rel_error_median').
    error_lo/hi_col: HPD columns for error bars; defaults to replacing '_median' with '_hpd_lower/upper'.
    """
    df_population = df_population[df_population["sampling"] == sampling].copy()
    df_bins       = df_time_bins[df_time_bins["sampling"] == sampling].copy()

    error_lo_col = error_col.replace("_median", "_hpd_lower")
    error_hi_col = error_col.replace("_median", "_hpd_upper")

    time_horizon = df_bins[time_col].max()

    if error_y_label is None:
        error_y_label = error_col

    default_mutsig_order       = ["low", "med", "high"]
    default_growth_model_order = ["uniform", "expgrowthslow", "expgrowthfast", "bottleneck"]
    mutsig_order       = mutsigs    if mutsigs    is not None else [m for m in default_mutsig_order       if m in df_bins["mutation_signal"].unique()]
    growth_model_order = pop_models if pop_models is not None else [m for m in default_growth_model_order if m in df_bins["population_model"].unique()]
    if isinstance(mutsig_order, str):       mutsig_order       = [mutsig_order]
    if isinstance(growth_model_order, str): growth_model_order = [growth_model_order]

    mutsig_labels = {"low": "Low mut. signal", "med": "Med. mut. signal", "high": "High mut. signal"}
    title_map     = {"expgrowthfast": "Exp. Growth (fast)", "expgrowthslow": "Exp. Growth (slow)",
                     "uniform": "Uniform", "bottleneck": "Bottleneck"}
    model_colors  = {"constcoal": "#397398", "skyline": "#80557e"}
    model_markers = {"constcoal": "o", "skyline": "^"}
    model_offsets = {"constcoal": 0, "skyline": 0}
    _bin_width = df_bins[time_col].diff().abs().median()

    ncols = len(growth_model_order)
    single_pop_model = ncols == 1

    if single_pop_model:
        ncols_fig = len(mutsig_order)
        fig, axes = plt.subplots(2, ncols_fig,
                                 figsize=figsize or (5 * ncols_fig, 5),
                                 gridspec_kw={"hspace": hspace, "wspace": wspace})
        if ncols_fig == 1:
            axes = axes.reshape(2, 1)
    else:
        row_heights = ([1, 1, 0.3] * len(mutsig_order))[:-1]
        fig = plt.figure(figsize=figsize or (8 * ncols, 2.5 * len(row_heights)))
        gs  = gridspec.GridSpec(len(row_heights), ncols,
                                height_ratios=row_heights, hspace=hspace, wspace=wspace)

    ax_pop_first_col = {}
    ax_err_first_col = {}

    for i, mutsig in enumerate(mutsig_order):
        for j, growth_model in enumerate(growth_model_order):
            if single_pop_model:
                ax_pop = axes[0, i]
                ax_err = axes[1, i]
            else:
                ax_pop = fig.add_subplot(gs[i * 3,     j])
                ax_err = fig.add_subplot(gs[i * 3 + 1, j])

            if j == 0:
                ax_pop_first_col[i] = ax_pop
                ax_err_first_col[i] = ax_err

            # ── Population trajectory ──────────────────────────────────────
            subset_pop = df_population[
                (df_population["mutation_signal"]  == mutsig) &
                (df_population["population_model"] == growth_model)]
            if not subset_pop.empty:
                true_traj = get_true_traj(growth_model, subset_pop.iloc[0])
                color_exp   = "#a6444f"
                color_sky   = "#80557e"
                color_const = "#397398"
                if trajectory_df is not None:
                    traj_sub = trajectory_df[
                        (trajectory_df["sampling"]          == sampling) &
                        (trajectory_df["mutation_signal"]   == mutsig) &
                        (trajectory_df["population_model"]  == growth_model)
                    ]
                    if not traj_sub.empty:
                        t_vals = traj_sub["time"].values
                        mask   = t_vals <= time_horizon
                        t_plot = t_vals[mask]
                        t_smooth = np.linspace(0, time_horizon, 1000)
                        N_true = np.array([true_traj(t) for t in t_smooth])
                        ax_pop.plot(t_smooth, N_true, color=color_exp, linewidth=2,
                                    label="True Population Size", zorder=6)
                        if mode in ('both', 'skyline'):
                            ax_pop.plot(t_plot, traj_sub["skyline_median"].values[mask],
                                        color=color_sky, linestyle='--', label="Skyline Estimate",
                                        linewidth=2, zorder=5)
                            ax_pop.fill_between(t_plot,
                                                traj_sub["skyline_hpd_lower"].values[mask],
                                                traj_sub["skyline_hpd_upper"].values[mask],
                                                color=color_sky, alpha=0.3, zorder=4)
                        if mode in ('both', 'constcoal'):
                            ax_pop.plot(t_plot, traj_sub["constcoal_median"].values[mask],
                                        color=color_const, linestyle=':', label="Constant Estimate",
                                        linewidth=2, zorder=5)
                            ax_pop.fill_between(t_plot,
                                                traj_sub["constcoal_hpd_lower"].values[mask],
                                                traj_sub["constcoal_hpd_upper"].values[mask],
                                                color=color_const, alpha=0.3, zorder=4)
                        if pop_y_range:
                            ax_pop.set_ylim(pop_y_range)
                        ax_pop.invert_xaxis()
                else:
                    plot_population_summary_ax(
                        ax=ax_pop,
                        skyline_all_times=subset_pop["skyline_times"].tolist(),
                        skyline_all_medians=subset_pop["skyline_medians"].tolist(),
                        constant_all_estimates=subset_pop["coalescent_median"].tolist(),
                        true_traj=true_traj,
                        time_horizon=time_horizon, y_range=pop_y_range, mode=mode)
                ax_pop.set_xlim(time_horizon, 0)
            ax_pop.tick_params(labelsize=tick_fontsize)
            ax_pop.set_xlabel(x_label, fontsize=label_fontsize)
            show_ylabel = j == 0
            ax_pop.set_ylabel(pop_y_label if show_ylabel else "", fontsize=label_fontsize)

            if show_changepoints and not subset_pop.empty and "skyline_times" in subset_pop.columns:
                _add_skyline_changepoint_strip(ax_pop, subset_pop,
                                               dot_size=changepoint_dot_size,
                                               grey_alpha=changepoint_grey_alpha,
                                               y_frac=changepoint_y_frac)

            if sampling_end_time is not None:
                ylims = ax_pop.get_ylim()
                y_arrow = ylims[1] * 0.85
                ax_pop.annotate(
                    "", xy=(sampling_end_time, y_arrow), xytext=(0, y_arrow),
                    arrowprops=dict(arrowstyle="-|>", color="grey", lw=1.5),
                    annotation_clip=False,
                )
                ax_pop.text(
                    sampling_end_time / 2, y_arrow + (ylims[1] - ylims[0]) * 0.02, "sampling",
                    color="grey", fontsize=sampling_fontsize if sampling_fontsize is not None else label_fontsize - 2, ha="center", va="bottom",
                )

            # ── Pre-binned error ───────────────────────────────────────────
            subset_bins = df_bins[
                (df_bins["mutation_signal"]  == mutsig) &
                (df_bins["population_model"] == growth_model) &
                (df_bins[time_col] <= time_horizon) &
                (df_bins["n_nodes"]   >= min_nodes_per_bin)
            ].copy()

            has_err = error_lo_col in subset_bins.columns and error_hi_col in subset_bins.columns

            # ── Raw individual node dots (background) ──────────────────────
            if per_tree_df is not None:
                pt_sub_all = per_tree_df[
                    (per_tree_df["sampling"]          == sampling) &
                    (per_tree_df["mutation_signal"]   == mutsig) &
                    (per_tree_df["population_model"]  == growth_model)
                ].copy()
                rng = np.random.default_rng(42)
                for model, color in model_colors.items():
                    pt_sub = pt_sub_all[pt_sub_all["model"] == model].dropna(subset=[error_col, time_col])
                    x_jitter = rng.uniform(-_bin_width * jitter_width, _bin_width * jitter_width, size=len(pt_sub))
                    ax_err.scatter(pt_sub[time_col].values + x_jitter, pt_sub[error_col].values,
                                   color=color, marker=model_markers[model],
                                   s=10, alpha=raw_alpha, zorder=1, linewidths=0)

            for model, color in model_colors.items():
                df_m = subset_bins[subset_bins["model"] == model].dropna(subset=[error_col])
                if df_m.empty:
                    continue
                xvals = df_m[time_col] + model_offsets[model]
                yvals = df_m[error_col]
                if has_err and errorbar_alpha > 0:
                    eb = ax_err.errorbar(
                        xvals, yvals,
                        yerr=[yvals - df_m[error_lo_col], df_m[error_hi_col] - yvals],
                        fmt=model_markers[model], color=color,
                        ecolor=color, elinewidth=1.4, capsize=4,
                        markersize=7, markeredgecolor=marker_edgecolor, markeredgewidth=0.8,
                        alpha=1.0, zorder=3)
                    for line in eb[1]: line.set_alpha(errorbar_alpha)
                    for line in eb[2]: line.set_alpha(errorbar_alpha)
                else:
                    ax_err.scatter(xvals, yvals, color=color,
                                   marker=model_markers[model], s=40,
                                   edgecolors=marker_edgecolor, linewidths=0.8,
                                   alpha=1.0, zorder=3)

            ax_err.axhline(0, color="black", linestyle="--", linewidth=1.5, zorder=1)
            ax_err.grid(True, color="lightgray", linewidth=0.4, alpha=0.5, zorder=0)
            ax_err.invert_xaxis()
            padding = (time_horizon - 0) * 0.05
            ax_err.set_xlim(time_horizon + padding, -padding)
            if error_y_range:
                ax_err.set_ylim(error_y_range)
            ax_err.tick_params(labelsize=tick_fontsize)
            ax_err.set_xlabel(x_label, fontsize=label_fontsize)
            ax_err.set_ylabel(error_y_label if show_ylabel else "", fontsize=label_fontsize)

            # ── Column titles (pop model) ──────────────────────────────────
            if i == 0:
                ax_pop.set_title(title_map.get(growth_model, growth_model), fontsize=label_fontsize + 2)

    handles = [Line2D([0], [0], color=c, lw=2, marker=model_markers[m],
                      markersize=7, label={"constcoal": "Const. coal.", "skyline": "Skyline"}[m])
               for m, c in model_colors.items()]

    if single_pop_model:
        for i, mutsig in enumerate(mutsig_order):
            ax_pop_first_col[i].set_title(mutsig_labels.get(mutsig, mutsig), fontsize=label_fontsize)
        if show_legend:
            fig.legend(handles=handles, loc='upper right', bbox_to_anchor=(0.99, 1.10),
                       ncol=2, frameon=False, fontsize=legend_fontsize)
        # place ylabels as fig.text at a shared x, each centered on its own axis
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        inv = fig.transFigure.inverted()
        for i in ax_pop_first_col:
            if ylabel_x_positions is not None and i < len(ylabel_x_positions):
                x_label_x = ylabel_x_positions[i]
            else:
                col_axes = [ax_pop_first_col[i], ax_err_first_col[i]]
                x_left = min(
                    inv.transform([[ax.get_tightbbox(renderer).x0, 0]])[0][0]
                    for ax in col_axes
                )
                x_label_x = x_left - 0.01
            for ax, label, y_offset in [(ax_pop_first_col[i], pop_y_label, 0.0), (ax_err_first_col[i], error_y_label, error_ylabel_y_offset)]:
                ax.set_ylabel("")
                pos = ax.get_position()
                y_fig = (pos.y0 + pos.y1) / 2 + y_offset
                fig.text(x_label_x, y_fig, label, fontsize=label_fontsize,
                         ha="center", va="center", rotation=90, transform=fig.transFigure)
    else:
        # Mutsig row labels centred between the two axes on the left
        fig.canvas.draw()
        renderer  = fig.canvas.get_renderer()
        inv       = fig.transFigure.inverted()
        # place ylabels as fig.text at a shared x, centered on each axis
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        inv = fig.transFigure.inverted()
        # find the leftmost tick label edge across all left-column axes
        all_left = list(ax_pop_first_col.values()) + list(ax_err_first_col.values())
        x_min = min(
            inv.transform(ax.yaxis.get_ticklabels()[-1].get_window_extent(renderer))[0][0]
            if ax.yaxis.get_ticklabels() else 1.0
            for ax in all_left
        )
        x_label_x = x_min - 0.02
        for i in ax_pop_first_col:
            ax_p = ax_pop_first_col[i]
            ax_e = ax_err_first_col[i]
            ax_p.set_ylabel("")
            ax_e.set_ylabel("")
            pos_p = ax_p.get_position()
            pos_e = ax_e.get_position()
            y_pop = (pos_p.y0 + pos_p.y1) / 2
            y_err = (pos_e.y0 + pos_e.y1) / 2
            if i == 0:
                fig.text(x_label_x, y_pop, pop_y_label, fontsize=label_fontsize,
                         ha="center", va="center", rotation=90, transform=fig.transFigure)
                fig.text(x_label_x, y_err, error_y_label, fontsize=label_fontsize,
                         ha="center", va="center", rotation=90, transform=fig.transFigure)
        for i, mutsig in enumerate(mutsig_order):
            ax_p = ax_pop_first_col[i]
            ax_e = ax_err_first_col[i]
            pos_p = ax_p.get_position()
            pos_e = ax_e.get_position()
            y_fig = (pos_p.y0 + pos_e.y1) / 2
            ylabel_bb = ax_p.yaxis.label.get_window_extent(renderer)
            x_fig = inv.transform((ylabel_bb.x0, 0))[0] - 0.01
            fig.text(x_fig, y_fig, mutsig_labels.get(mutsig, mutsig),
                     fontsize=label_fontsize, ha="center", va="center", rotation=90,
                     transform=fig.transFigure)
        if show_legend:
            fig.legend(handles=handles, loc='upper right', bbox_to_anchor=(0.97, 0.97),
                       ncol=2, frameon=False, fontsize=legend_fontsize)

    if xticks is not None:
        for ax in fig.get_axes():
            ax.set_xticks(xticks)

    suptitle = title if title is not None else f"Population size & {error_col} ({sampling})"
    plt.suptitle(suptitle, fontsize=label_fontsize + 2)
    fig.subplots_adjust(left=0.12, right=0.97, top=0.93, bottom=0.07)

    if save_path is not None:
        fig.savefig(save_path, bbox_inches="tight", transparent=True)
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
    title_map = {"expgrowthfast": "Exp. Growth (fast)", "expgrowthslow": "Exp. Growth (slow)",
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

    bl_est_col = f"{model}_bl_est_median"

    def filter_slice(df, sampling, pop_model, mut_signal):
        subset = tree_metrics_combined[
            (tree_metrics_combined["sampling"] == sampling) &
            (tree_metrics_combined["population_model"] == pop_model) &
            (tree_metrics_combined["mutation_signal"] == mut_signal) &
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
                plot_df["bl_sim"]    = plot_df["bl_sim"] * scale
                plot_df[bl_est_col]  = plot_df[bl_est_col] * scale
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
                    group[bl_est_col],
                    label=sampling_type,
                    color=color_map[sampling_type],
                    alpha=alpha,
                    s=s,
                )

            min_val = min(plot_df["bl_sim"].min(), plot_df[bl_est_col].min())
            max_val = max(plot_df["bl_sim"].max(), plot_df[bl_est_col].max())

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

def load_pop_trajectories(eval_dir):
    """
    Read all *_pop_trajectory.tsv files from eval_dir and concatenate into
    one DataFrame with added sampling, population_model, mutation_signal columns.
    """
    eval_dir = Path(eval_dir)
    dfs = []
    for tsv in sorted(eval_dir.rglob("*_pop_trajectory.tsv")):
        df = pd.read_csv(tsv, sep="\t")
        # parse scenario from filename: {sampling}_{pop_model}_{mutsig}mutsig_pop_trajectory.tsv
        stem = tsv.stem.replace("_pop_trajectory", "")
        mutsig = stem.split("_")[-1].replace("mutsig", "")
        rest = stem[: stem.rfind(f"_{mutsig}mutsig")]
        sampling = None
        for s in ("independenthomochronous", "linearconstant"):
            if rest.startswith(s):
                sampling = s
                pop_model = rest[len(s) + 1:]
                break
        df["sampling"] = sampling
        df["population_model"] = pop_model
        df["mutation_signal"] = mutsig
        dfs.append(df)
    if not dfs:
        raise FileNotFoundError(f"No *_pop_trajectory.tsv files found under {eval_dir}")
    return pd.concat(dfs, ignore_index=True)


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


def load_time_bins(eval_dir, bin_width, cutoff):
    """
    Read all time bin TSV files matching a given bin_width and cutoff from
    eval_dir and concatenate into one DataFrame. sampling, population_model,
    and mutation_signal are parsed from the scenario column.

    Parameters
    ----------
    eval_dir  : str or Path
    bin_width : int   e.g. 10
    cutoff    : int   e.g. 400
    """
    eval_dir = Path(eval_dir)
    pattern  = f"*_time_bins_w{int(bin_width)}_c{int(cutoff)}.tsv"
    dfs = []
    for tsv in sorted(eval_dir.rglob(pattern)):
        dfs.append(pd.read_csv(tsv, sep="\t"))
    if not dfs:
        raise FileNotFoundError(
            f"No files matching '{pattern}' found under {eval_dir}"
        )
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



def load_time_bins_per_tree(eval_dir, bin_width, cutoff):
    """
    Like load_time_bins but reads the per-tree files
    (*_time_bins_w{bin_width}_c{cutoff}_per_tree.tsv).
    """
    eval_dir = Path(eval_dir)
    pattern  = f"*_time_bins_w{int(bin_width)}_c{int(cutoff)}_per_tree.tsv"
    dfs = []
    for tsv in sorted(eval_dir.rglob(pattern)):
        dfs.append(pd.read_csv(tsv, sep="\t"))
    if not dfs:
        raise FileNotFoundError(
            f"No files matching '{pattern}' found under {eval_dir}"
        )
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


def load_time_bins_long_bl(eval_dir, bin_width, cutoff):
    """
    Like load_time_bins but reads the long-branch-length subset files
    (*_time_bins_w{bin_width}_c{cutoff}_long_bl.tsv).
    """
    eval_dir = Path(eval_dir)
    pattern  = f"*_time_bins_w{int(bin_width)}_c{int(cutoff)}_long_bl.tsv"
    dfs = []
    for tsv in sorted(eval_dir.rglob(pattern)):
        dfs.append(pd.read_csv(tsv, sep="\t"))
    if not dfs:
        raise FileNotFoundError(
            f"No files matching '{pattern}' found under {eval_dir}"
        )
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


def compute_pairwise_tip_distances(tree_index, repo_base=None, run="run1",
                                   population_model=None, mutation_signal=None,
                                   sampling=None):
    """
    For one or all successful trees, compute for each pair of tips the distance
    of each tip to their MRCA, from constcoal and skyline summary trees.

    tree_index : int, or 'all' to process every successful tree.
    population_model, mutation_signal, sampling : filter to a single scenario
        (recommended when tree_index='all' to keep runtime manageable).

    Distance of tip X to MRCA = h_mrca - h_tip  (both heights in time before present).
    HPD bounds: conservative — lower = h_mrca_lo - h_tip_hi, upper = h_mrca_hi - h_tip_lo.

    Returns
    -------
    DataFrame with one row per tip pair per scenario, columns:
        sampling, population_model, mutation_signal, tree_index, tip_a, tip_b,
        dist_a_constcoal_median, dist_a_constcoal_hpd_lower, dist_a_constcoal_hpd_upper,
        dist_b_constcoal_median, dist_b_constcoal_hpd_lower, dist_b_constcoal_hpd_upper,
        dist_a_skyline_median,   dist_a_skyline_hpd_lower,   dist_a_skyline_hpd_upper,
        dist_b_skyline_median,   dist_b_skyline_hpd_lower,   dist_b_skyline_hpd_upper,
        dist_a_abs_diff  (constcoal_median - skyline_median),
        dist_a_rel_diff  ((constcoal_median - skyline_median) / skyline_median),
        dist_b_abs_diff,
        dist_b_rel_diff
    """
    if repo_base is None:
        repo_base = Path(__file__).resolve().parent.parent
    repo_base = Path(repo_base)

    eval_dir = repo_base / "results" / run / "evaluation"
    tsv_files = sorted({p for p in eval_dir.rglob("*.tsv")
                        if p.stem in ("low", "med", "high")})

    path_df = pd.concat(
        [pd.read_csv(f, sep="\t") for f in tsv_files],
        ignore_index=True
    )

    if tree_index == "all":
        subset = path_df.copy()
    else:
        subset = path_df[path_df["tree_index"] == tree_index]

    if population_model is not None:
        subset = subset[subset["population_model"] == population_model]
    if mutation_signal is not None:
        subset = subset[subset["mutation_signal"] == mutation_signal]
    if sampling is not None:
        subset = subset[subset["sampling"] == sampling]

    def _parse_height(clade):
        h_med, h_lo, h_hi = 0.0, 0.0, 0.0
        if clade.comment:
            m = re.search(r'height_median=([\d\.eE+-]+)', clade.comment)
            if m: h_med = float(m.group(1))
            m = re.search(r'height_95%_HPD=\{([\d\.eE+-]+),([\d\.eE+-]+)\}', clade.comment)
            if m: h_lo, h_hi = float(m.group(1)), float(m.group(2))
        return h_med, h_lo, h_hi

    # key: (sampling, population_model, mutation_signal, tip_a, tip_b)
    # value: dict accumulating per-model distances
    pair_data = {}

    for _, row in subset.iterrows():
        scenario_key = (row["sampling"], row["population_model"], row["mutation_signal"], row["tree_index"])

        for model, col in (("constcoal", "trees_path_constcoal"),
                           ("skyline",   "trees_path_skyline")):
            raw = row[col]
            if pd.isna(raw):
                continue
            summary_path = repo_base / str(raw).replace(".combined.trees", ".combined_summary.tree")
            if not summary_path.exists():
                continue

            tree = Phylo.read(str(summary_path), "nexus")
            tip_info  = {c.name: _parse_height(c) for c in tree.get_terminals()}
            node_info = {id(c): _parse_height(c)   for c in tree.get_nonterminals()}

            tips = list(tip_info.keys())
            for i_t, tip_a in enumerate(tips):
                for tip_b in tips[i_t + 1:]:
                    mrca = tree.common_ancestor(tip_a, tip_b)
                    mrca_med, mrca_lo, mrca_hi = node_info.get(id(mrca), (0.0, 0.0, 0.0))

                    ha_med, ha_lo, ha_hi = tip_info[tip_a]
                    hb_med, hb_lo, hb_hi = tip_info[tip_b]

                    # distance tip → MRCA = h_mrca - h_tip
                    da_med = mrca_med - ha_med
                    da_lo  = mrca_lo  - ha_hi   # conservative lower
                    da_hi  = mrca_hi  - ha_lo   # conservative upper
                    db_med = mrca_med - hb_med
                    db_lo  = mrca_lo  - hb_hi
                    db_hi  = mrca_hi  - hb_lo

                    key = scenario_key + (tip_a, tip_b)
                    if key not in pair_data:
                        pair_data[key] = {}
                        # sampling times are tree-independent; store once
                        pair_data[key]["tip_a_time"] = ha_med
                        pair_data[key]["tip_b_time"] = hb_med
                    pair_data[key][f"dist_a_{model}_median"]    = da_med
                    pair_data[key][f"dist_a_{model}_hpd_lower"] = da_lo
                    pair_data[key][f"dist_a_{model}_hpd_upper"] = da_hi
                    pair_data[key][f"dist_b_{model}_median"]    = db_med
                    pair_data[key][f"dist_b_{model}_hpd_lower"] = db_lo
                    pair_data[key][f"dist_b_{model}_hpd_upper"] = db_hi

    records = []
    for (sampling, pop_model, mut_sig, t_idx, tip_a, tip_b), vals in pair_data.items():
        rec = {
            "sampling":         sampling,
            "population_model": pop_model,
            "mutation_signal":  mut_sig,
            "tree_index":       t_idx,
            "tip_a":            tip_a,
            "tip_b":            tip_b,
        }
        rec.update(vals)
        # diffs (constcoal - skyline)
        for tip_lbl in ("a", "b"):
            cc = rec.get(f"dist_{tip_lbl}_constcoal_median", float("nan"))
            sk = rec.get(f"dist_{tip_lbl}_skyline_median",   float("nan"))
            rec[f"dist_{tip_lbl}_abs_diff"] = cc - sk
            rec[f"dist_{tip_lbl}_rel_diff"] = (cc - sk) / sk if sk != 0 else float("nan")
        records.append(rec)

    return pd.DataFrame(records)


def plot_pairwise_tip_distances(dist_df, sampling,
                                x_range=None, y_range=None,
                                alpha=0.5, errorbar_alpha=0.3, errorbar_width=0.2,
                                n_subsample=None,
                                pop_models=None, mutsigs=None,
                                figsize=None, title=None,
                                x_label="Skyline — pairwise tip distance",
                                y_label="Const. Coalescent — pairwise tip distance",
                                show_pop_model_labels=True,
                                show_mutsig_labels=True,
                                running_avg=False,
                                window_size=500, window_overlap=250,
                                running_avg_color="#c0392b",
                                diff_plot=False,
                                diff_mode="absolute",
                                diff_axis_range=None,
                                label_fontsize=13, tick_fontsize=13, legend_fontsize=11,
                                save_path=None):
    """
    Scatter plot of skyline (x) vs constcoal (y) pairwise tip-to-MRCA distances,
    with HPD error bars. Rows = mutation signal, columns = pop model.
    Optionally adds running median line and difference panel below each column.
    """
    from matplotlib.gridspec import GridSpec

    mutsig_order       = mutsigs    if mutsigs    is not None else ["low", "med", "high"]
    growth_model_order = pop_models if pop_models is not None else ["uniform", "expgrowthslow", "expgrowthfast", "bottleneck"]
    if isinstance(mutsig_order, str):       mutsig_order       = [mutsig_order]
    if isinstance(growth_model_order, str): growth_model_order = [growth_model_order]

    mutsig_labels = {"low": "Low mut. signal", "med": "Med. mut. signal", "high": "High mut. signal"}
    title_map     = {"expgrowthfast": "Exp. Growth (fast)", "expgrowthslow": "Exp. Growth (slow)",
                     "uniform": "Uniform", "bottleneck": "Bottleneck"}
    color = "#555555"

    base = dist_df[dist_df["sampling"] == sampling].copy()
    tip_a_df = base.rename(columns={
        "dist_a_constcoal_median":    "constcoal_dist_median",
        "dist_a_constcoal_hpd_lower": "constcoal_dist_hpd_lower",
        "dist_a_constcoal_hpd_upper": "constcoal_dist_hpd_upper",
        "dist_a_skyline_median":      "skyline_dist_median",
        "dist_a_skyline_hpd_lower":   "skyline_dist_hpd_lower",
        "dist_a_skyline_hpd_upper":   "skyline_dist_hpd_upper",
    })
    tip_b_df = base.rename(columns={
        "dist_b_constcoal_median":    "constcoal_dist_median",
        "dist_b_constcoal_hpd_lower": "constcoal_dist_hpd_lower",
        "dist_b_constcoal_hpd_upper": "constcoal_dist_hpd_upper",
        "dist_b_skyline_median":      "skyline_dist_median",
        "dist_b_skyline_hpd_lower":   "skyline_dist_hpd_lower",
        "dist_b_skyline_hpd_upper":   "skyline_dist_hpd_upper",
    })
    keep_cols = ["sampling", "population_model", "mutation_signal", "tree_index",
                 "constcoal_dist_median", "constcoal_dist_hpd_lower", "constcoal_dist_hpd_upper",
                 "skyline_dist_median",   "skyline_dist_hpd_lower",   "skyline_dist_hpd_upper"]
    wide = pd.concat([tip_a_df[keep_cols], tip_b_df[keep_cols]], ignore_index=True)

    nrows, ncols = len(mutsig_order), len(growth_model_order)

    # layout: each mutsig row = scatter + optional diff panel
    _scatter_h = (figsize[1] if figsize else 3.5 * nrows)
    _scatter_w = (figsize[0] if figsize else 4 * ncols)
    _diff_h_per_row = 1.5  # inches, added on top of scatter height

    if diff_plot:
        gs_nrows = nrows * 2
        height_ratios = [_scatter_h / nrows, _diff_h_per_row] * nrows
        actual_figsize = (_scatter_w, _scatter_h + _diff_h_per_row * nrows)
    else:
        gs_nrows = nrows
        height_ratios = [1] * nrows
        actual_figsize = (_scatter_w, _scatter_h)

    fig = plt.figure(figsize=actual_figsize)
    gs  = GridSpec(gs_nrows, ncols, figure=fig, height_ratios=height_ratios,
                   hspace=0.22 if diff_plot else 0.3, wspace=0.3)

    axes      = np.empty((nrows, ncols), dtype=object)
    axes_diff = np.empty((nrows, ncols), dtype=object) if diff_plot else None

    for i in range(nrows):
        scatter_row = i * 2 if diff_plot else i
        for j in range(ncols):
            axes[i, j] = fig.add_subplot(gs[scatter_row, j])
            if diff_plot:
                axes_diff[i, j] = fig.add_subplot(gs[scatter_row + 1, j], sharex=axes[i, j])

    # ---- sliding window (x=skyline, y=constcoal) ----
    def _sliding_window(x, y, ws, overlap):
        order = np.argsort(x)
        x_s, y_s = x[order], y[order]
        step = max(1, ws - overlap)
        xs, ys, ds, x_lefts, x_rights = [], [], [], [], []
        for start in range(0, len(x_s) - ws + 1, step):
            sl = slice(start, start + ws)
            xw, yw = x_s[sl], y_s[sl]
            xs.append(np.median(xw));  ys.append(np.median(yw))
            x_lefts.append(xw[0]);     x_rights.append(xw[-1])
            d_vals = (yw - xw) / xw if diff_mode == "relative" else yw - xw
            ds.append(np.median(d_vals))
        return (np.array(xs), np.array(ys), np.array(ds),
                np.array(x_lefts), np.array(x_rights))

    def _step_path(ys, x_lefts, x_rights):
        n = len(ys)
        if n == 0:
            return np.array([]), np.array([])
        splits = (x_rights[:-1] + x_lefts[1:]) / 2
        xp, yp = [x_lefts[0]], [ys[0]]
        for k in range(n - 1):
            xp += [splits[k], splits[k]];  yp += [ys[k], ys[k + 1]]
        xp.append(x_rights[-1]);  yp.append(ys[-1])
        return np.array(xp), np.array(yp)

    for i, mutsig in enumerate(mutsig_order):
        for j, pop_model in enumerate(growth_model_order):
            ax = axes[i, j]
            ax_d = axes_diff[i, j] if diff_plot else None

            cell_full = wide[(wide["mutation_signal"] == mutsig) &
                             (wide["population_model"] == pop_model)].dropna(
                subset=["skyline_dist_median", "constcoal_dist_median"])

            if cell_full.empty:
                ax.axis("off")
                if ax_d is not None: ax_d.axis("off")
                continue

            cell = cell_full.sample(n=n_subsample, random_state=42) if n_subsample is not None and len(cell_full) > n_subsample else cell_full

            xi = cell["skyline_dist_median"].values
            yi = cell["constcoal_dist_median"].values
            xerr = np.array([xi - cell["skyline_dist_hpd_lower"].values,
                             cell["skyline_dist_hpd_upper"].values - xi])
            yerr = np.array([yi - cell["constcoal_dist_hpd_lower"].values,
                             cell["constcoal_dist_hpd_upper"].values - yi])
            eb = ax.errorbar(xi, yi, xerr=xerr, yerr=yerr,
                             fmt="o", color=color, ecolor=color,
                             elinewidth=errorbar_width, capsize=0, markersize=1, alpha=alpha)
            for line in eb[1]: line.set_alpha(errorbar_alpha)
            for line in eb[2]: line.set_alpha(errorbar_alpha)

            if running_avg or diff_plot:
                x_all = cell_full["skyline_dist_median"].values
                y_all = cell_full["constcoal_dist_median"].values
                win_x, win_y, win_d, win_xl, win_xr = _sliding_window(x_all, y_all, window_size, window_overlap)

                if running_avg and len(win_x):
                    xp, yp = _step_path(win_y, win_xl, win_xr)
                    ax.plot(xp, yp, color=running_avg_color, linewidth=0.8, zorder=5)

                if diff_plot and ax_d is not None and len(win_x):
                    bg_x = cell["skyline_dist_median"].values
                    bg_d = (cell["constcoal_dist_median"].values - bg_x) / bg_x if diff_mode == "relative" else cell["constcoal_dist_median"].values - bg_x
                    ax_d.scatter(bg_x, bg_d, s=1, color="#bbbbbb", alpha=0.3, zorder=1, rasterized=True)
                    ax_d.scatter(win_x, win_d, s=4, color=running_avg_color, zorder=3)
                    ax_d.axhline(0, color="black", linewidth=0.8, linestyle="--", zorder=4)
                    ax_d.grid(True, color="lightgray", linewidth=0.4, alpha=0.5, zorder=0)
                    if x_range: ax_d.set_xlim(x_range)
                    if diff_axis_range: ax_d.set_ylim(diff_axis_range)
                    ax_d.tick_params(labelsize=tick_fontsize - 2)
                    ylabel = "Rel. diff." if diff_mode == "relative" else "Diff. (years)"
                    if j == 0: ax_d.set_ylabel(ylabel, fontsize=label_fontsize - 1)

            ax.axline((0, 0), slope=1, color="black", linestyle="--", linewidth=1.0, zorder=1)
            ax.set_aspect("equal")
            ax.grid(True, color="lightgray", linewidth=0.4, alpha=0.5, zorder=0)
            if x_range: ax.set_xlim(x_range)
            if y_range: ax.set_ylim(y_range)
            ax.tick_params(labelsize=tick_fontsize)
            ax.set_xlabel("")
            if diff_plot and ax_d is not None:
                ax_d.set_xlabel("Skyline", fontsize=label_fontsize)
            elif not diff_plot:
                ax.set_xlabel("Skyline", fontsize=label_fontsize)
            if j == 0:
                ax.set_ylabel("Const. Coalescent", fontsize=label_fontsize)
            else:
                ax.set_ylabel("")
            if i == 0 and show_pop_model_labels:
                ax.set_title(title_map.get(pop_model, pop_model), fontsize=label_fontsize + 2)

            if j == 0 and i == 0:
                marker_line = Line2D([0], [0], marker="o", color=color, markersize=4,
                                     linestyle="none", alpha=1.0)
                caplines = (Line2D([0], [0], color=color, linewidth=0.4),
                            Line2D([0], [0], color=color, linewidth=0.4))
                barlines = (LineCollection([[[0, -0.3], [0, 0.3]]], colors=[color], linewidths=[0.4]),
                            LineCollection([[[-0.3, 0], [0.3, 0]]], colors=[color], linewidths=[0.4]))
                handle = ErrorbarContainer((marker_line, caplines, barlines),
                                           has_xerr=True, has_yerr=True,
                                           label="Tip dist. from\npairwise MRCA")
                ax.legend(handles=[handle], fontsize=legend_fontsize, frameon=False,
                          loc="upper left", borderpad=0.1, handletextpad=0.5,
                          handlelength=0.5, bbox_to_anchor=(0.0, 1.0))

    suptitle = title if title is not None else f"Pairwise tip distances — {sampling}"
    fig.suptitle(suptitle, fontsize=17)
    top = 0.97 if suptitle else 1.0
    fig.subplots_adjust(top=top)

    # align diff panel widths to the actual equal-aspect scatter square
    if diff_plot:
        fig.canvas.draw()
        _renderer = fig.canvas.get_renderer()
        _inv = fig.transFigure.inverted()
        for i in range(nrows):
            for j in range(ncols):
                if axes_diff[i, j] is None:
                    continue
                bb = axes[i, j].get_window_extent(_renderer)
                x0_fig, _ = _inv.transform((bb.x0, bb.y0))
                x1_fig, _ = _inv.transform((bb.x1, bb.y1))
                dp = axes_diff[i, j].get_position()
                axes_diff[i, j].set_position([x0_fig, dp.y0, x1_fig - x0_fig, dp.height])

    if show_mutsig_labels:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        inv = fig.transFigure.inverted()
        for i, mutsig in enumerate(mutsig_order):
            ax = axes[i][0]
            ylabel_bb = ax.yaxis.label.get_window_extent(renderer)
            x_fig = inv.transform((ylabel_bb.x0, ylabel_bb.y0 + ylabel_bb.height / 2))[0]
            pos = ax.get_position()
            y_fig = (pos.y0 + pos.y1) / 2
            fig.text(x_fig - 0.01, y_fig, mutsig_labels.get(mutsig, mutsig),
                     fontsize=label_fontsize + 2, ha="center", va="center", rotation=90,
                     transform=fig.transFigure)

    if save_path is not None:
        fig.savefig(save_path, bbox_inches="tight", transparent=True)
    from IPython.display import display as _display
    _display(fig)
    plt.close(fig)


def plot_pairwise_tip_distances_zoom(dist_df, sampling,
                                     x_range=None, y_range=None,
                                     zoom_x_range=None, zoom_y_range=None,
                                     alpha=0.5, errorbar_alpha=0.3, errorbar_width=0.2,
                                     n_subsample=None, n_subsample_zoom=None,
                                     pop_model="bottleneck", mutsig="low",
                                     figsize=(12, 5),
                                     title=None,
                                     x_label="Skyline", y_label="Const. Coalescent",
                                     running_avg=False,
                                     window_size=500, window_overlap=250,
                                     window_size_zoom=None, window_overlap_zoom=None,
                                     running_avg_color="#c0392b",
                                     diff_plot=False,
                                     diff_mode="absolute",
                                     diff_axis_range=None,
                                     label_fontsize=13, tick_fontsize=13, legend_fontsize=11,
                                     save_path=None):
    """
    Two-panel figure: left = full scatter, right = zoomed view of a rectangular window.
    Optionally adds a running median line and a difference panel below the zoom panel.

    diff_mode : 'absolute' for constcoal - skyline, 'relative' for (constcoal - skyline) / skyline.
    """
    from matplotlib.patches import Rectangle
    import matplotlib.colors as _mcolors

    color = "#555555"
    zoom = zoom_x_range is not None

    # ---- data prep ----
    base = dist_df[dist_df["sampling"] == sampling].copy()
    tip_a_df = base.rename(columns={
        "dist_a_constcoal_median":    "constcoal_dist_median",
        "dist_a_constcoal_hpd_lower": "constcoal_dist_hpd_lower",
        "dist_a_constcoal_hpd_upper": "constcoal_dist_hpd_upper",
        "dist_a_skyline_median":      "skyline_dist_median",
        "dist_a_skyline_hpd_lower":   "skyline_dist_hpd_lower",
        "dist_a_skyline_hpd_upper":   "skyline_dist_hpd_upper",
    })
    tip_b_df = base.rename(columns={
        "dist_b_constcoal_median":    "constcoal_dist_median",
        "dist_b_constcoal_hpd_lower": "constcoal_dist_hpd_lower",
        "dist_b_constcoal_hpd_upper": "constcoal_dist_hpd_upper",
        "dist_b_skyline_median":      "skyline_dist_median",
        "dist_b_skyline_hpd_lower":   "skyline_dist_hpd_lower",
        "dist_b_skyline_hpd_upper":   "skyline_dist_hpd_upper",
    })
    keep_cols = ["sampling", "population_model", "mutation_signal", "tree_index",
                 "constcoal_dist_median", "constcoal_dist_hpd_lower", "constcoal_dist_hpd_upper",
                 "skyline_dist_median",   "skyline_dist_hpd_lower",   "skyline_dist_hpd_upper"]
    df_full = pd.concat([tip_a_df[keep_cols], tip_b_df[keep_cols]], ignore_index=True)
    df_full = df_full[(df_full["mutation_signal"] == mutsig) & (df_full["population_model"] == pop_model)]
    df_full = df_full.dropna(subset=["constcoal_dist_median", "skyline_dist_median"])

    df_scatter = df_full.sample(n=n_subsample, random_state=42) if n_subsample is not None and len(df_full) > n_subsample else df_full

    # ---- layout ----
    def _actual_pos(left, bottom, w, h, fw, fh):
        side_in = min(w * fw, h * fh)
        aw, ah = side_in / fw, side_in / fh
        return left + (w - aw) / 2, bottom + (h - ah) / 2, aw, ah

    fw, fh = figsize
    fig = plt.figure(figsize=figsize)

    if zoom:
        ml, mb, mw, mh = 0.10, 0.10, 0.48, 0.82
        ax_main = fig.add_axes([ml, mb, mw, mh])
        al, ab, aw, ah = _actual_pos(ml, mb, mw, mh, fw, fh)
        gap          = 0.06
        zoom_left    = al + aw + gap
        zoom_w       = aw * 0.70
        zoom_h       = ah * 0.70
        zoom_bottom  = ab + ah - zoom_h   # top-aligned with main
        ax_zoom = fig.add_axes([zoom_left, zoom_bottom, zoom_w, zoom_h])
        if diff_plot:
            diff_gap = 0.08
            diff_h   = zoom_bottom - ab - diff_gap
            ax_diff  = fig.add_axes([zoom_left, ab, zoom_w, diff_h])
        else:
            ax_diff = None
    else:
        ax_zoom = None
        ax_diff = None
        ax_main = fig.add_axes([0.10, 0.10, 0.82, 0.82])

    # ---- scatter helper (x=skyline, y=constcoal) ----
    def _scatter(ax_, df_):
        xi  = df_["skyline_dist_median"].values
        yi  = df_["constcoal_dist_median"].values
        xlo = df_["skyline_dist_hpd_lower"].values
        xhi = df_["skyline_dist_hpd_upper"].values
        ylo = df_["constcoal_dist_hpd_lower"].values
        yhi = df_["constcoal_dist_hpd_upper"].values
        from matplotlib.collections import LineCollection as _LC
        x_segs = np.stack([np.column_stack([xlo, yi]), np.column_stack([xhi, yi])], axis=1)
        y_segs = np.stack([np.column_stack([xi, ylo]), np.column_stack([xi, yhi])], axis=1)
        ax_.add_collection(_LC(x_segs, colors="#aaaaaa", linewidths=errorbar_width,
                               alpha=errorbar_alpha, zorder=2, rasterized=True))
        ax_.add_collection(_LC(y_segs, colors="#aaaaaa", linewidths=errorbar_width,
                               alpha=errorbar_alpha, zorder=2, rasterized=True))
        ax_.scatter(xi, yi, s=1, color=color, alpha=alpha, zorder=3, rasterized=True)

    # ---- sliding window stats (x=skyline, y=constcoal) ----
    def _sliding_window(x, y, ws, overlap):
        order  = np.argsort(x)
        x_s, y_s = x[order], y[order]
        step   = max(1, ws - overlap)
        xs, ys, ds, d_los, d_his, x_lefts, x_rights = [], [], [], [], [], [], []
        for start in range(0, len(x_s) - ws + 1, step):
            sl = slice(start, start + ws)
            xw, yw = x_s[sl], y_s[sl]
            xs.append(np.median(xw))
            ys.append(np.median(yw))
            x_lefts.append(xw[0]);  x_rights.append(xw[-1])
            d_vals = (yw - xw) / xw if diff_mode == "relative" else yw - xw
            ds.append(np.median(d_vals))
        return (np.array(xs), np.array(ys), np.array(ds),
                np.array(x_lefts), np.array(x_rights))

    def _step_path(ys, x_lefts, x_rights):
        n = len(ys)
        if n == 0:
            return np.array([]), np.array([])
        splits = (x_rights[:-1] + x_lefts[1:]) / 2
        xp, yp = [x_lefts[0]], [ys[0]]
        for i in range(n - 1):
            xp += [splits[i], splits[i]];  yp += [ys[i], ys[i + 1]]
        xp.append(x_rights[-1]);  yp.append(ys[-1])
        return np.array(xp), np.array(yp)

    # ---- compute window stats on full dataset (x=skyline, y=constcoal) ----
    x_all = df_full["skyline_dist_median"].values
    y_all = df_full["constcoal_dist_median"].values
    if running_avg or diff_plot:
        win_x, win_y, win_d, win_xl, win_xr = _sliding_window(
            x_all, y_all, window_size, window_overlap)

    # ---- main panel ----
    _scatter(ax_main, df_scatter)

    if running_avg:
        xp, yp = _step_path(win_y, win_xl, win_xr)
        ax_main.plot(xp, yp, color=running_avg_color, linewidth=1.5, zorder=5, label="Running median")

    lims = list(x_range) if x_range is not None else None
    if lims is None:
        all_vals = np.concatenate([
            df_scatter["skyline_dist_hpd_lower"].values,   df_scatter["skyline_dist_hpd_upper"].values,
            df_scatter["constcoal_dist_hpd_lower"].values, df_scatter["constcoal_dist_hpd_upper"].values,
        ])
        all_vals = all_vals[np.isfinite(all_vals)]
        lims = [float(all_vals.min()), float(all_vals.max())]

    ax_main.plot(lims, lims, "--", color="black", linewidth=1, zorder=4)
    ax_main.set_xlim(lims)
    ax_main.set_ylim(lims if y_range is None else list(y_range))
    ax_main.set_xlabel(x_label, fontsize=label_fontsize)
    ax_main.set_ylabel(y_label, fontsize=label_fontsize)
    ax_main.tick_params(labelsize=tick_fontsize)
    ax_main.set_aspect("equal")
    if title is not None:
        ax_main.set_title(title, fontsize=label_fontsize + 1)

    marker_line = Line2D([0], [0], marker="o", color=color, markersize=4, linestyle="none", alpha=1.0)
    caplines = (Line2D([0], [0], color=color, linewidth=errorbar_width),
                Line2D([0], [0], color=color, linewidth=errorbar_width))
    barlines = (LineCollection([[[0, -0.3], [0, 0.3]]], colors=[color], linewidths=[errorbar_width]),
                LineCollection([[[-0.3, 0], [0.3, 0]]], colors=[color], linewidths=[errorbar_width]))
    handle = ErrorbarContainer((marker_line, caplines, barlines), has_xerr=True, has_yerr=True,
                                label="Tip dist. from\npairwise MRCA")
    ax_main.legend(handles=[handle], fontsize=legend_fontsize, frameon=False,
                   loc="upper left", borderpad=0.5, handletextpad=0.4)

    # ---- zoom panel ----
    if zoom:
        df_zoom_pool = df_full[
            (df_full["skyline_dist_median"]   >= zoom_x_range[0]) & (df_full["skyline_dist_median"]   <= zoom_x_range[1]) &
            (df_full["constcoal_dist_median"] >= zoom_y_range[0]) & (df_full["constcoal_dist_median"] <= zoom_y_range[1])
        ] if zoom_y_range is not None else df_full[
            (df_full["skyline_dist_median"] >= zoom_x_range[0]) & (df_full["skyline_dist_median"] <= zoom_x_range[1])
        ]
        df_zoom = df_zoom_pool.sample(n=n_subsample_zoom, random_state=43) if n_subsample_zoom is not None and len(df_zoom_pool) > n_subsample_zoom else df_zoom_pool
        _scatter(ax_zoom, df_zoom)

        z_xlim = list(zoom_x_range)
        z_ylim = list(zoom_y_range) if zoom_y_range is not None else z_xlim

        if running_avg or diff_plot:
            ws_z = window_size_zoom   if window_size_zoom   is not None else window_size
            ov_z = window_overlap_zoom if window_overlap_zoom is not None else window_overlap
            zwin_x, zwin_y, zwin_d, zwin_xl, zwin_xr = _sliding_window(
                x_all, y_all, ws_z, ov_z)
            inside = (zwin_xr >= zoom_x_range[0]) & (zwin_xl <= zoom_x_range[1])
            sel = np.where(inside)[0]

            if running_avg and len(sel):
                xp, yp = _step_path(zwin_y[sel], zwin_xl[sel], zwin_xr[sel])
                ax_zoom.plot(xp, yp, color=running_avg_color, linewidth=1.5, zorder=5)

            if diff_plot and ax_diff is not None and len(sel):
                bg_x = df_zoom["skyline_dist_median"].values
                bg_d = (df_zoom["constcoal_dist_median"].values - bg_x) / bg_x if diff_mode == "relative" else df_zoom["constcoal_dist_median"].values - bg_x
                if len(bg_x):
                    ax_diff.scatter(bg_x, bg_d, s=1, color="#bbbbbb", alpha=0.4, zorder=1, rasterized=True)
                ax_diff.scatter(zwin_x[sel], zwin_d[sel], s=3**2, color=running_avg_color, zorder=3)
                ax_diff.axhline(0, color="black", linewidth=0.8, linestyle="--", zorder=4)
                ax_diff.set_xlim(z_xlim)
                if diff_axis_range is not None:
                    ax_diff.set_ylim(diff_axis_range)
                ylabel = "Rel. diff." if diff_mode == "relative" else "Diff. (years)"
                ax_diff.set_ylabel(ylabel, fontsize=label_fontsize - 1, labelpad=-2)
                ax_diff.set_xlabel(x_label, fontsize=label_fontsize - 1)
                ax_diff.tick_params(labelsize=tick_fontsize - 1)

        ax_zoom.plot(z_xlim, z_ylim, "--", color="black", linewidth=1, zorder=4)
        ax_zoom.set_xlim(z_xlim)
        ax_zoom.set_ylim(z_ylim)
        ax_zoom.set_aspect("equal")
        ax_zoom.tick_params(labelsize=tick_fontsize - 1)
        ax_zoom.set_xlabel("")
        ax_zoom.set_ylabel("")

        # rectangle + connecting lines
        zx0, zx1 = zoom_x_range
        zy0, zy1 = (zoom_y_range[0], zoom_y_range[1]) if zoom_y_range else (zx0, zx1)
        zoom_color = "#a6444f"
        _r, _g, _b, _ = _mcolors.to_rgba(zoom_color)
        rect = Rectangle((zx0, zy0), zx1 - zx0, zy1 - zy0,
                          linewidth=1.2, edgecolor=zoom_color,
                          facecolor=(_r, _g, _b, 0.15), linestyle="-", zorder=5)
        ax_main.add_patch(rect)
        for corner_main, corner_zoom in [((zx1, zy1), (zx0, zy1)),
                                          ((zx1, zy0), (zx0, zy0))]:
            con = plt.matplotlib.patches.ConnectionPatch(
                xyA=corner_main, coordsA=ax_main.transData,
                xyB=corner_zoom, coordsB=ax_zoom.transData,
                color=zoom_color, linewidth=0.8, linestyle="--", zorder=5)
            fig.add_artist(con)

    if save_path is not None:
        fig.savefig(save_path, bbox_inches="tight", transparent=True)
    from IPython.display import display as _display
    _display(fig)
    plt.close(fig)


def _ltt_events_from_summary_tree(tree_path):
    """
    Returns sorted internal node heights (time before present) and minimum tip
    height from a BEAST summary tree. Uses height_median annotations.
    """
    tree = Phylo.read(str(tree_path), "nexus")

    def _h(c):
        if c.comment:
            m = re.search(r'height_median=([\d\.eE+-]+)', c.comment)
            if m:
                return float(m.group(1))
        return None

    int_h = [_h(c) for c in tree.get_nonterminals()]
    int_h = [h for h in int_h if h is not None]
    if not int_h:
        return None, None
    tip_h = [_h(c) or 0.0 for c in tree.get_terminals()]
    return sorted(int_h, reverse=True), min(tip_h)


def _ltt_events_from_true_tree(tree_obj):
    """
    Returns sorted internal node heights and minimum tip height from a
    time-calibrated newick tree.
    """
    depths     = tree_obj.depths(unit_branch_lengths=False)
    root_depth = max(depths.values())
    tip_h  = [root_depth - d for c, d in depths.items() if c.is_terminal()]
    int_h  = [root_depth - d for c, d in depths.items() if not c.is_terminal()]
    return sorted(int_h, reverse=True), min(tip_h)


def _events_to_ltt(events, min_tip):
    """Convert sorted branching event heights to step-function arrays."""
    root_h = events[0]
    times  = [root_h]
    counts = [1]
    for h in events:
        times.append(h)
        counts.append(counts[-1] + 1)
    times.append(min_tip)
    counts.append(counts[-1])
    return np.array(times), np.array(counts)


def _interpolate_ltt(times, counts, grid):
    """Evaluate a step-function LTT at grid points (grid decreasing, oldest first)."""
    # times are decreasing (past -> present), counts increasing
    return np.array([counts[np.searchsorted(-times, -t, side="right") - 1]
                     for t in grid], dtype=float)


def compute_ltt(csv_path, sim_dir,
                sampling=None, pop_models=None, mutsigs=None,
                repo_base=None,
                n_grid=1000,
                grid_end=None,
                out_per_tree=None,
                out_median=None):
    """
    Compute LTT curves for constcoal, skyline, and true trees.

    Outputs
    -------
    out_per_tree : TSV with columns
        model, sampling, pop_model, mutsig, tree_index, time, lineages
    out_median : TSV with columns
        model, sampling, pop_model, mutsig, time, median_lineages
        (median over all trees evaluated on a shared grid per scenario)

    Returns (per_tree_df, median_df).
    """
    csv_path = Path(csv_path)
    sim_dir  = Path(sim_dir)
    base     = Path(repo_base) if repo_base else csv_path.parents[2]

    df  = pd.read_csv(csv_path)
    pat = re.compile(
        r'(constcoal|skyline)/([^/]+)/([^/]+)/(\w+mutsig)/\w+\.T(\d+)\.combined\.trees'
    )
    rows = []
    for p in df["trees_path"]:
        m = pat.search(p)
        if m:
            rows.append({"model": m.group(1), "sampling": m.group(2),
                         "pop_model": m.group(3),
                         "mutsig": m.group(4).replace("mutsig", ""),
                         "tree_index": int(m.group(5)), "trees_path": p})
    meta = pd.DataFrame(rows)

    if sampling is not None:
        meta = meta[meta["sampling"] == sampling]
    _sampling = sampling or (meta["sampling"].iloc[0] if not meta.empty else "")

    mutsig_order   = ([mutsigs]    if isinstance(mutsigs, str)    else mutsigs)    or ["low", "med", "high"]
    popmodel_order = ([pop_models] if isinstance(pop_models, str) else pop_models) or ["uniform", "expgrowthslow", "expgrowthfast", "bottleneck"]

    per_tree_rows = []
    median_rows   = []

    for ms in mutsig_order:
        for pm in popmodel_order:
            sim_path = sim_dir / _sampling / pm / f"{pm}.trees"
            if not sim_path.exists():
                print(f"Skipping {pm}/{ms}: sim tree not found at {sim_path}")
                continue
            true_trees = list(Phylo.parse(str(sim_path), "newick"))

            sub = meta[(meta["mutsig"] == ms) & (meta["pop_model"] == pm)]

            # accumulate per-model LTT arrays for median computation
            model_ltt_arrays = {"constcoal": [], "skyline": [], "true": []}
            all_times = {"constcoal": [], "skyline": [], "true": []}

            for t_idx, grp in sub.groupby("tree_index"):
                cc_row = grp[grp["model"] == "constcoal"]
                sk_row = grp[grp["model"] == "skyline"]
                if cc_row.empty or sk_row.empty or t_idx >= len(true_trees):
                    continue

                cc_sum = base / cc_row.iloc[0]["trees_path"].replace(".combined.trees", ".combined_summary.tree")
                sk_sum = base / sk_row.iloc[0]["trees_path"].replace(".combined.trees", ".combined_summary.tree")

                if not cc_sum.exists() or cc_sum.stat().st_size == 0:
                    continue
                if not sk_sum.exists() or sk_sum.stat().st_size == 0:
                    continue

                cc_ev, cc_tip = _ltt_events_from_summary_tree(cc_sum)
                sk_ev, sk_tip = _ltt_events_from_summary_tree(sk_sum)
                tr_ev, tr_tip = _ltt_events_from_true_tree(true_trees[t_idx])

                if cc_ev is None or sk_ev is None:
                    continue

                for label, ev, tip in [("constcoal", cc_ev, cc_tip),
                                        ("skyline",   sk_ev, sk_tip),
                                        ("true",      tr_ev, tr_tip)]:
                    times, counts = _events_to_ltt(ev, tip)
                    for t, c in zip(times, counts):
                        per_tree_rows.append({
                            "model": label, "sampling": _sampling,
                            "pop_model": pm, "mutsig": ms,
                            "tree_index": t_idx, "time": t, "lineages": c,
                        })
                    all_times[label].extend(times.tolist())
                    model_ltt_arrays[label].append((times, counts))

            # compute median on shared grid per model
            for label, ltt_list in model_ltt_arrays.items():
                if not ltt_list:
                    continue
                t_max = grid_end if grid_end is not None else max(t[0][0] for t in ltt_list)
                t_min = min(t[0][-1] for t in ltt_list)
                grid  = np.linspace(t_max, t_min, n_grid)
                matrix = np.vstack([_interpolate_ltt(t, c, grid)
                                    for t, c in ltt_list])
                med = np.median(matrix, axis=0)
                for t, m_val in zip(grid, med):
                    median_rows.append({
                        "model": label, "sampling": _sampling,
                        "pop_model": pm, "mutsig": ms,
                        "time": t, "median_lineages": m_val,
                    })

    per_tree_df = pd.DataFrame(per_tree_rows)
    median_df   = pd.DataFrame(median_rows)

    if out_per_tree:
        per_tree_df.to_csv(out_per_tree, sep="\t", index=False)
        print(f"Saved per-tree LTT to {out_per_tree}")
    if out_median:
        median_df.to_csv(out_median, sep="\t", index=False)
        print(f"Saved median LTT to {out_median}")

    return per_tree_df, median_df


def plot_ltt(per_tree_df, median_df,
             sampling=None, pop_models=None, mutsigs=None,
             alpha=0.15, linewidth=0.5,
             median_linewidth=1.8,
             constcoal_color="#397398",
             skyline_color="#7e4d8b",
             true_color="#a6444f",
             x_range=None, y_range=None,
             figsize=None,
             title=None,
             show_pop_model_labels=True,
             show_mutsig_labels=True,
             sampling_end_time=None,
             label_fontsize=14, tick_fontsize=14, legend_fontsize=14,
             log_y=False,
             save_path=None):
    """
    Plot LTT curves from precomputed compute_ltt() outputs.
    Thin lines = individual trees; thick lines = median across trees.
    """
    colors = {"constcoal": constcoal_color, "skyline": skyline_color, "true": true_color}
    labels = {"constcoal": "Const. coal.", "skyline": "Skyline", "true": "True"}

    if sampling is not None:
        per_tree_df = per_tree_df[per_tree_df["sampling"] == sampling]
        median_df   = median_df[median_df["sampling"] == sampling]

    mutsig_order   = ([mutsigs]    if isinstance(mutsigs, str)    else mutsigs)    or sorted(per_tree_df["mutsig"].unique())
    popmodel_order = ([pop_models] if isinstance(pop_models, str) else pop_models) or sorted(per_tree_df["pop_model"].unique())

    title_map     = {"expgrowthfast": "Exp. Growth (fast)", "expgrowthslow": "Exp. Growth (slow)",
                     "uniform": "Uniform", "bottleneck": "Bottleneck"}
    mutsig_labels = {"low": "Low mut. signal", "med": "Med. mut. signal", "high": "High mut. signal"}

    nrows, ncols = len(mutsig_order), len(popmodel_order)
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=figsize or (4 * ncols, 3.5 * nrows),
                             squeeze=False)

    for i, ms in enumerate(mutsig_order):
        for j, pm in enumerate(popmodel_order):
            ax = axes[i, j]
            pt  = per_tree_df[(per_tree_df["mutsig"] == ms) & (per_tree_df["pop_model"] == pm)]
            med = median_df[(median_df["mutsig"] == ms) & (median_df["pop_model"] == pm)]

            if pt.empty:
                ax.axis("off"); continue

            legend_handles = {}
            for model in ("constcoal", "skyline", "true"):
                color = colors[model]
                pt_m  = pt[pt["model"] == model]
                # thin lines per tree
                for t_idx, grp in pt_m.groupby("tree_index"):
                    grp = grp.sort_values("time", ascending=False)
                    l, = ax.step(grp["time"], grp["lineages"], where="post",
                                 color=color, alpha=alpha, linewidth=linewidth)
                    legend_handles.setdefault(model, l)
                # median line
                med_m = med[med["model"] == model].sort_values("time", ascending=False)
                if not med_m.empty:
                    ax.plot(med_m["time"], med_m["median_lineages"],
                            color=color, linewidth=median_linewidth, zorder=5)

            ax.set_xlabel("Time (before present)", fontsize=label_fontsize)
            ax.set_ylabel("n. lineages", fontsize=label_fontsize)
            ax.invert_xaxis()
            if log_y:
                ax.set_yscale("log")
            ax.grid(True, axis="x", color="lightgray", linewidth=0.4, alpha=0.5)
            ax.tick_params(labelsize=tick_fontsize)
            if x_range: ax.set_xlim(x_range)
            if y_range:
                y_lo = max(y_range[0], 1e-1) if log_y else y_range[0]
                ax.set_ylim(y_lo, y_range[1])

            if sampling_end_time is not None:
                ylims = ax.get_ylim()
                y_arrow = ylims[0] + (ylims[1] - ylims[0]) * 0.07
                ax.annotate(
                    "", xy=(sampling_end_time, y_arrow), xytext=(0, y_arrow),
                    arrowprops=dict(arrowstyle="-|>", color="grey", lw=1.5),
                    annotation_clip=False,
                )
                ax.text(
                    sampling_end_time / 2, y_arrow + (ylims[1] - ylims[0]) * 0.02,
                    "sampling", color="grey", fontsize=label_fontsize - 1, ha="center", va="bottom",
                )

            if legend_handles:
                legend_lines = [
                    plt.Line2D([0], [0], color=colors[k], linewidth=2.5)
                    for k in legend_handles
                ]
                ax.legend(
                    handles=legend_lines,
                    labels=[labels[k] for k in legend_handles],
                    fontsize=legend_fontsize, frameon=False)
            if i == 0 and show_pop_model_labels and len(popmodel_order) > 1:
                ax.set_title(title_map.get(pm, pm), fontsize=label_fontsize + 1)
            if show_mutsig_labels and j == 0 and len(mutsig_order) > 1:
                ax.set_ylabel(f"{mutsig_labels.get(ms, ms)}\nn. lineages", fontsize=label_fontsize)
            else:
                ax.set_ylabel("n. lineages", fontsize=label_fontsize)

    suptitle = title if title is not None else (sampling or "")
    if suptitle and suptitle != "linearconstant":
        plt.suptitle(suptitle, fontsize=17)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
    else:
        plt.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, bbox_inches="tight", transparent=True)
    from IPython.display import display as _display
    _display(fig)
    plt.close(fig)

