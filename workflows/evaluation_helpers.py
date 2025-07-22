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


### CALCULATE NODE HEIGHTS AND BRANCH LENGHTS
def get_branch_info(tree):
    """
    Return a list of:
    (clade_id, branch_length, height, branch_length_CI_lower, branch_length_CI_upper,
     height_CI_lower, height_CI_upper)
    
    - height: distance from root to clade
    - confidence intervals (CIs) extracted from comments, if available
    """
    root_height = max(tree.distance(tree.root, leaf) for leaf in tree.get_terminals())

    node_info = []
    i = 0
    for clade in tree.find_clades(order="preorder"):
        # Assign a clade ID (tip name or internal node ID)
        if clade.name:
            clade_id = clade.name
            internal = False
        else:
            clade_id = f"internal_{i}"
            i += 1
            internal = True

        # Branch length (always from node to parent node)
        bl = clade.branch_length if clade.branch_length is not None else 0.0

        # Node height: distance from root
        if internal:
            height = root_height - tree.distance(tree.root, clade)
        else:
            height = 0

        # Confidence interval parsing
        bl_lower_ci, bl_upper_ci = None, None
        height_lower_ci, height_upper_ci = None, None

        if clade.comment:
            # Branch length CI: e.g. length_range={0.1,0.5}
            bl_match = re.search(r'length_range=\{([\d\.eE+-]+),([\d\.eE+-]+)\}', clade.comment)
            if bl_match:
                bl_lower_ci = float(bl_match.group(1))
                bl_upper_ci = float(bl_match.group(2))
            
            # Height CI: e.g. height_range={0.8,1.2}
            h_match = re.search(r'height_range=\{([\d\.eE+-]+),([\d\.eE+-]+)\}', clade.comment)
            if h_match:
                height_lower_ci = float(h_match.group(1))
                height_upper_ci = float(h_match.group(2))

        node_info.append((
            clade_id,
            bl,
            height,
            bl_lower_ci, bl_upper_ci,
            height_lower_ci, height_upper_ci,
            internal
        ))

    return node_info


def compare_tree_metrics(tree_sim, tree_constcoal):
    # --- Branch lengths and node heights ---
    sim_info = get_branch_info(tree_sim)
    const_info = get_branch_info(tree_constcoal)

    branch_length_df = pd.DataFrame({
        "node": [b1[0] for b1 in sim_info],
        "bl_sim": [b1[1] for b1 in sim_info],
        "bl_constcoal": [b2[1] for b2 in const_info],
        "bl_ci_lower_constcoal": [b2[3] for b2 in const_info],
        "bl_ci_upper_constcoal": [b2[4] for b2 in const_info],
        "bl_inside_ci": [
            b2[3] <= b1[1] <= b2[4] if b2[3] is not None and b2[4] is not None else None
            for b1, b2 in zip(sim_info, const_info)
        ],
        "height_sim": [b1[2] for b1 in sim_info],
        "height_constcoal": [b2[2] for b2 in const_info],
        "height_ci_lower_constcoal": [b2[5] for b2 in const_info],
        "height_ci_upper_constcoal": [b2[6] for b2 in const_info],
        "height_inside_ci": [
            b2[5] <= b1[2] <= b2[6] if b2[5] is not None and b2[5] is not None else None
            for b1, b2 in zip(sim_info, const_info)
        ],
        "internal": [b2[7] for b2 in const_info]
    })

    # Errors for branch lengths
    branch_length_df["bl_diff"] = branch_length_df["bl_constcoal"] - branch_length_df["bl_sim"]
    branch_length_df["bl_relative_error"] = branch_length_df["bl_diff"] / branch_length_df["bl_sim"]
    branch_length_df["bl_abs_relative_error"] = np.abs(branch_length_df["bl_relative_error"])

    # Errors for node heights
    branch_length_df["height_diff"] = branch_length_df["height_constcoal"] - branch_length_df["height_sim"]
    branch_length_df["height_relative_error"] = [(d / s) if is_internal and s != 0 else 0.0 for d, s, is_internal in zip(branch_length_df["height_diff"], branch_length_df["height_sim"], branch_length_df["internal"])]
    #branch_length_df["height_relative_error"] = branch_length_df["height_diff"] / branch_length_df["height_sim"]
    branch_length_df["height_abs_relative_error"] = np.abs(branch_length_df["height_relative_error"])

    return branch_length_df

def tree_metrics_all_trees(path_df):
    # Now process each row
    results = []
    for _, row in path_df.iterrows():
        tree_constcoal = Phylo.read(row["tree_path_constcoal"], "nexus")
        tree_skyline = Phylo.read(row["tree_path_skyline"], "nexus")
        sim_trees = list(Phylo.parse(row["sim_tree_path"], "newick"))
        tree_sim = sim_trees[row["tree_index"]]
        
        eval_df_constcoal = compare_tree_metrics(tree_sim, tree_constcoal)
        eval_df_skyline = compare_tree_metrics(tree_sim, tree_skyline)

        eval_df_constcoal["tree_name"] = Path(row["tree_path_constcoal"]).stem
        eval_df_skyline["tree_name"] = Path(row["tree_path_skyline"]).stem
        eval_df_constcoal["tree_index"] = row["tree_index"]
        eval_df_skyline["tree_index"] = row["tree_index"]
        results.append(eval_df_constcoal)
        results.append(eval_df_skyline)

    # Save combined output
    combined_df = pd.concat(results, ignore_index=True)

    combined_df[["model", "growth_model", "mutsig"]] = combined_df["tree_name"].apply(extract_model_components)

    return combined_df

def root_height_all_trees(path_df):
    # Now process each row
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

    combined_df[["model", "growth_model", "mutsig"]] = combined_df["tree_name"].apply(extract_model_components)

    combined_df['diff_root_height'] = combined_df['estimated_root_height'] - combined_df['sim_root_height']
    combined_df['rel_diff_root_height'] = combined_df['diff_root_height'] / combined_df['sim_root_height']
    combined_df['abs_rel_diff_root_height'] = np.abs(combined_df['rel_diff_root_height'])

    return combined_df


### PLOTTING

def plot_population_trajectories_ax(ax, skyline_times, skyline_means, constant_pop_estimate, 
                                     present_pop_size, growth_rate=None,
                                     color_exp="#a6444f", color_sky="#80557e", color_const="#397398",
                                     alpha=0.4, label_prefix="", first_plot=True, time_horizon=0):

    root_height = skyline_times[-1]
    t_max = root_height if time_horizon == 0 else min(time_horizon, root_height)
    t_vals = np.linspace(0, t_max, 1000)
    N_true = present_pop_size * np.exp(-growth_rate * t_vals)

    # Stepwise skyline
    skyline_start_times = [0.0] + skyline_times[:-1]
    step_times, step_values = [], []
    for start, end, value in zip(skyline_start_times, skyline_times, skyline_means):
        if start > t_max:
            break
        end = min(end, t_max)
        step_times.extend([start, end])
        step_values.extend([value, value])
    
    if first_plot:
        ax.plot(t_vals, N_true, color=color_exp, alpha=1, linewidth= 2, label="True Population Size")
        ax.hlines(constant_pop_estimate, 0, t_max, color=color_const, linestyle=':', 
              alpha=alpha)
        ax.plot(step_times, step_values, drawstyle='steps-post', linestyle='--', color=color_sky, 
            alpha=alpha)
    else:
        ax.hlines(constant_pop_estimate, 0, t_max, color=color_const, linestyle=':', 
              alpha=alpha, label="Constant Pop. Estimate")
        ax.plot(step_times, step_values, drawstyle='steps-post', linestyle='--', color=color_sky, 
            alpha=alpha, label="Skyline Estimate")

    if first_plot:
        ax.set_xlabel("Time before present")
        ax.set_ylabel("Population Size")
        ax.invert_xaxis()



def plot_summary_population_grid(path_info_df, num_groups=10, burnin=0.1, time_horizon = 0):
    """
    Plots all population trajectories for all trees in a subplot grid by
    population model (columns) and mutation signal (rows).
    
    Parameters:
        path_info_df (pd.DataFrame): Table with tree/log paths and metadata.
        present_pop_size (float): N0
        growth_rate (float): Exponential growth rate
        num_groups (int): Number of skyline intervals
        burnin (float or int): Burnin for BEAST logs
    """
    pop_models = ["uniform", "expgrowth_slow", "expgrowth_fast"]
    mut_signals = ["low", "med", "high"]
    ncols, nrows = len(pop_models), len(mut_signals)

    fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 3.5*nrows))

    # Ensure axes is always 2D array
    if nrows == 1:
        axes = np.expand_dims(axes, axis=0)
    if ncols == 1:
        axes = np.expand_dims(axes, axis=1)

    first_plot = True
    previous_tree_index = 0

    for idx, row in path_info_df.iterrows():
        
        tree_path_constcoal = row["tree_path_constcoal"]
        tree_path_skyline = row["tree_path_skyline"]
        log_path_constcoal = row["log_path_constcoal"]
        log_path_skyline = row["log_path_skyline"]
        pop_model = row["population_model"]
        mut_sig = row["mutation_signal"]
        present_pop_size = row["present_pop_size"]
        growth_rate = row["growth_rate"]
        tree_index = row["tree_index"]


        row_idx = mut_signals.index(mut_sig)
        col_idx = pop_models.index(pop_model)
        ax = axes[row_idx][col_idx]

        # Parse data
        tree = Phylo.read(tree_path_skyline, "nexus")
        skyline_times = get_skyline_group_boundaries(tree, num_groups=num_groups)
        skyline_means = get_mean_population_size(log_path_skyline, burnin=burnin, mode="skyline")
        coalescent_mean = get_mean_population_size(log_path_constcoal, burnin=burnin, mode="constcoal")

        if tree_index < previous_tree_index:
            first_plot = True
        previous_tree_index = tree_index

        # Plot into axis
        plot_population_trajectories_ax(
            ax,
            skyline_times=skyline_times,
            skyline_means=skyline_means,
            constant_pop_estimate=coalescent_mean,
            present_pop_size=present_pop_size,
            growth_rate=growth_rate,
            alpha=0.5,  # overlay transparency
            first_plot=first_plot,
            time_horizon=time_horizon,
        )
        first_plot = False

    # Axis labeling
    for i, mut_sig in enumerate(mut_signals):
        axes[i][0].set_ylabel(f"{mut_sig.capitalize()} mut.\nsignal\nPopulation Size")

    for j, pop_model in enumerate(pop_models):
        axes[0][j].set_title(f"{pop_model.replace('_', ' ').capitalize()}", fontsize=12)


    # Add global legend
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.02), ncol=3)
    #fig.suptitle("Population Size Trajectories across Trees", y=1.07)
    plt.show()


def boxplot_branch_length_errors(df_complete, metric = "bl_abs_relative_error", title= "", logscale=True):

    if metric.startswith("bl"):
        df = df_complete.dropna().copy()
    elif metric.startswith("height"):
        df = df_complete[df_complete.internal == True].copy()
    else:
        df = df_complete.copy()

    # Define consistent hue order and colors
    hue_order = ["constcoal", "skyline"]
    palette = ["#a6444f", "#397398"]

    # Order of conditions: 3 growth models × 3 mutation signals
    mutsig_order = ["low", "med", "high"]
    growth_model_order = ["uniform", "expgrowth_slow", "expgrowth_fast"]
    condition_order = [
        f"{g}/{m}" for g in growth_model_order for m in mutsig_order
    ]

    # Combine into condition column
    df["condition"] = df["growth_model"] + "/" + df["mutsig"]

    # Plot
    plt.figure(figsize=(14, 6))

    ax = sns.boxplot(
        data=df,
        x="condition",
        y=metric,
        hue="model",
        palette=palette,
        hue_order=hue_order,
        showfliers=False,
        order=condition_order,
    )

    sns.stripplot(
        data=df,
        x="condition",
        y=metric,
        hue="model",
        palette=palette,
        hue_order=hue_order,
        order=condition_order,
        dodge=True,
        alpha=0.2,
        size=2,
        jitter=True,
        linewidth=0
    )

    # Remove duplicate legends
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles[:len(hue_order)],
        labels[:len(hue_order)],
        title="Model",
        loc='center left',
        bbox_to_anchor=(1.01, 0.5),
        frameon=False
    )

    # Log scale and styling
    if logscale:
        ax.set_yscale("log")
    else:
        ax.set_ylim(-5, 10)
    ax.set_ylabel(f"Error (log scale, {metric})")
    ax.set_xlabel("")
    ax.set_xticks(range(len(condition_order)))
    ax.set_xticklabels(
        [cond.split("/")[1].replace("mutsig", "") for cond in condition_order], rotation=0
    )

    # Add growth model labels below
    for i, growth_model in enumerate(growth_model_order):
        mid = i * 3 + 1  # middle of 3 mutation signal blocks
        label = growth_model.replace("expgrowth_", "exp-growth ").capitalize()
        ax.text(mid, -0.08, label, ha='center', va='top', fontsize=10,
                transform=ax.get_xaxis_transform())

    # Tighten layout
    plt.title(title)
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout(rect=[0, 0.03, 0.95, 1])  # make room for legend + bottom labels
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

def plot_height_errors_by_time_bin(ax, df_sub, y_max, bins=10, error_col="height_abs_relative_error"):
    # Bin simulated heights (time since present)
    df_sub = df_sub.copy()
    df_sub["time_bin"] = pd.cut(df_sub["height_sim"], bins=np.linspace(0, y_max, bins + 1), include_lowest=True)

    # Drop NaNs
    df_sub = df_sub.dropna(subset=["time_bin", error_col])

    # Plot
    sns.boxplot(
        data=df_sub,
        x="time_bin",
        y=error_col,
        hue="model",
        ax=ax,
        showfliers=False,
        palette=["#a6444f", "#397398"]
    )

    '''sns.violinplot(
        data=df_sub,
        x="time_bin",
        y=error_col,
        hue="model",
        ax=ax,
        cut=0,               # Don’t extend violins beyond data range
        density_norm='width',      # Same width for all violins
        inner="quartile",    # Show quartiles
        palette=["#a6444f", "#397398"],
        dodge=True           # Separate violins for each hue
    )'''

    '''sns.stripplot(
        data=df_sub,
        x="time_bin",
        y=error_col,
        hue="model",
        ax = ax,
        palette=["#a6444f", "#397398"],
        dodge=True,
        alpha=0.2,
        size=2,
        jitter=True,
        linewidth=0
    )'''
    
    # Formatting
    #ax.set_ylim(-5, 5)
    ax.set_xlabel("Time before present")
    tick_labels = [f"{int(interval.left)}-{int(interval.right)}" for interval in df_sub["time_bin"].cat.categories]
    ax.set_xticks(range(len(tick_labels)))
    ax.set_xticklabels(tick_labels, rotation=45)
    ax.legend(title="Model", loc="upper right")

def plot_height_error_grid(df, y_max, bins=10, error_col="height_abs_relative_error", title = ""):
    mutsig_order = ["low", "med", "high"]
    growth_model_order = ["uniform", "expgrowth_slow", "expgrowth_fast"]
    
    fig, axes = plt.subplots(nrows=len(mutsig_order), ncols=len(growth_model_order), figsize=(18, 10), sharey=False)

    for i, mutsig in enumerate(mutsig_order):
        for j, growth_model in enumerate(growth_model_order):
            ax = axes[i, j]
            df_sub = df[(df["mutsig"] == mutsig) & (df["growth_model"] == growth_model)]
            plot_height_errors_by_time_bin(ax, df_sub, y_max=y_max, bins=bins, error_col=error_col)
            
            if i == 0:
                ax.set_title(growth_model.replace("expgrowth_", "exp-growth ").capitalize())
            if j == 0:
                ax.set_ylabel(f"{mutsig} mutation signal\n{error_col}")

    plt.suptitle(title, fontsize = 20)
    plt.tight_layout()
    plt.show()

