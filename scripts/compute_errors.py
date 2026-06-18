"""
compute_errors.py

For a given scenario TSV (from process_results), computes per-node error
summaries (median + 95% CI across posterior samples) for each replicate.

Topology is fixed: all posterior trees share the same topology as the
simulated tree, so nodes are matched by preorder position throughout.

Each tree is traversed exactly once. The simulated tree produces a node
dictionary and event list reused for all posterior samples.

Output: one TSV per scenario, one row per (replicate, node).

Usage
-----
    python scripts/compute_errors.py \
        --tsv results/run1/evaluation/independenthomochronous/expgrowthfast/low.tsv \
        --out_dir results/run1/evaluation/independenthomochronous/expgrowthfast/
"""

import argparse
import re
import numpy as np
import pandas as pd
from pathlib import Path
from Bio import Phylo
from io import StringIO
from math import comb


# =============================================================================
# Statistics
# =============================================================================

def summarize(values):
    vals = np.array([v for v in values if not np.isnan(v)], dtype=float)
    if len(vals) == 0:
        return np.nan, np.nan, np.nan
    return float(np.median(vals)), float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


# =============================================================================
# Single tree traversal
# =============================================================================

def traverse_tree(tree):
    """
    Single traversal. Assigns node names (internal_i / tip_i) in place,
    builds node dict and event list.

    Returns:
      node_dict : node_idx -> (height, bl, is_internal, node_id)
      events    : sorted [(time, "sample"/"coal", clade)] for CI sweep
    """
    depths      = tree.depths()
    root_height = max(d for node, d in depths.items() if node.is_terminal())
    node_dict = {}
    events    = []

    for node_idx, clade in enumerate(tree.find_clades(order="preorder")):
        is_int  = not clade.is_terminal()
        node_id = f"internal_{node_idx}" if is_int else f"tip_{node_idx}"
        if is_int:
            clade.name = node_id
        h  = root_height - depths[clade]
        bl = clade.branch_length or 0.0
        node_dict[node_idx] = (h, bl, is_int, node_id)
        n_children = len(clade.clades) if is_int else 1
        events.append((h, "coal" if is_int else "sample", node_idx, n_children))

    events.sort(key=lambda x: x[0])
    return node_dict, events


# =============================================================================
# Population trajectory and integrals
# =============================================================================

def get_true_traj(pop_model, row):
    if pop_model in ("expgrowthfast", "expgrowthslow"):
        N0, r = float(row["present_pop_size"]), float(row["growth_rate"])
        return lambda t: N0 * np.exp(-r * t)
    elif pop_model.startswith("bottleneck"):
        N0, Nb = float(row["present_pop_size"]), float(row["bottleneck_size"])
        bs, be = float(row["bottleneck_start"]),  float(row["bottleneck_end"])
        return lambda t: Nb if bs < t < be else N0
    else:
        return lambda t: float(row["present_pop_size"])


def integral_true(t0, t1, pop_model, row):
    if t1 <= t0:
        return 0.0
    if pop_model in ("expgrowthfast", "expgrowthslow"):
        N0, r = float(row["present_pop_size"]), float(row["growth_rate"])
        return (np.exp(r * t1) - np.exp(r * t0)) / (N0 * r)
    elif pop_model.startswith("bottleneck"):
        N0, Nb = float(row["present_pop_size"]), float(row["bottleneck_size"])
        bs, be = float(row["bottleneck_start"]),  float(row["bottleneck_end"])
        total = 0.0
        for a, b, Ne in [(0.0, bs, N0), (bs, be, Nb), (be, np.inf, N0)]:
            l, r = max(t0, a), min(t1, b)
            if r > l: total += (r - l) / Ne
        return total
    else:
        return (t1 - t0) / float(row["present_pop_size"])


def integral_constcoal(t0, t1, Ne):
    return (t1 - t0) / Ne if t1 > t0 else 0.0


def integral_skyline(t0, t1, bounds, pop_sizes):
    if t1 <= t0:
        return 0.0
    total = 0.0
    for a, b, Ne in zip(bounds[:-1], bounds[1:], pop_sizes):
        l, r = max(t0, a), min(t1, b)
        if r > l: total += (r - l) / Ne
    if t1 > bounds[-1]:
        total += (t1 - max(t0, bounds[-1])) / pop_sizes[-1]
    return total


# =============================================================================
# Coalescent intensity sweep
# =============================================================================

def coalescent_intensity(events, integral_fn):
    """
    Accumulate coalescent intensity H at each internal node.
    Returns dict: node_idx -> H at that coalescent event.
    """
    H, k, t_prev, ci = 0.0, 0, 0.0, {}
    i = 0
    while i < len(events):
        t = events[i][0]
        if t > t_prev and k >= 2:
            H += comb(k, 2) * integral_fn(t_prev, t)
        batch = []
        while i < len(events) and events[i][0] == t:
            batch.append(events[i]); i += 1
        for _, etype, node_idx, _ in batch:
            if etype == "coal": ci[node_idx] = H
        for _, etype, node_idx, n_children in batch:
            if etype == "sample": k += 1
            elif etype == "coal": k -= max(n_children - 1, 1)
        t_prev = t
    return ci


# =============================================================================
# Skyline helpers
# =============================================================================

def skyline_boundaries_from_preorder(node_dict, num_groups):
    times       = sorted(h for h, _, is_int, _ in node_dict.values() if is_int)
    base, r     = divmod(len(times), num_groups)
    sizes       = [base + (1 if i < r else 0) for i in range(num_groups)]
    bounds, idx = [], 0
    for i, s in enumerate(sizes):
        idx += s
        bounds.append(times[idx - 1] if i < num_groups - 1 else max(times))
    return bounds


def skyline_pop_sizes(log_row, num_groups):
    sizes = []
    for i in range(1, num_groups + 1):
        for col in (f"skyline.popSize{i}", f"skyline.popSize.{i}"):
            if col in log_row.index:
                sizes.append(float(log_row[col])); break
    return sizes


def skyline_Ne_at(t, bounds, pop_sizes):
    idx = max(int(np.searchsorted(bounds, t, side="left")) - 1, 0)
    return pop_sizes[min(idx, len(pop_sizes) - 1)]


# =============================================================================
# Per-node error computation
# =============================================================================

METRICS = [
    "pop_diff",    "pop_rel_error",    "pop_abs_rel_error",
    "rate_diff",   "rate_rel_error",   "rate_abs_rel_error",
    "ci_diff",     "ci_rel_error",     "ci_abs_rel_error",
    "height_diff", "height_rel_error", "height_abs_rel_error",
    "bl_diff",     "bl_rel_error",     "bl_abs_rel_error",
]


def compute_errors(h_sim, bl_sim, h_est, bl_est, is_int,
                   Ne_true, Ne_est, ci_true, ci_est):
    def sd(a, b): return a / b if b and b != 0 and not np.isnan(b) else np.nan

    pop_d  = Ne_est - Ne_true
    rt     = 1.0 / Ne_true if Ne_true else np.nan
    re     = 1.0 / Ne_est  if Ne_est  else np.nan
    rate_d = re - rt if not np.isnan(rt) else np.nan
    h_d    = h_est - h_sim if is_int else np.nan
    bl_d   = bl_est - bl_sim
    ci_d   = (ci_est - ci_true) if is_int and not np.isnan(ci_true) and not np.isnan(ci_est) else np.nan

    return {
        "pop_diff":             pop_d,          "pop_rel_error":        sd(pop_d, Ne_true),
        "pop_abs_rel_error":    abs(sd(pop_d, Ne_true)),
        "rate_diff":            rate_d,         "rate_rel_error":       sd(rate_d, rt),
        "rate_abs_rel_error":   abs(sd(rate_d, rt)),
        "ci_diff":              ci_d,           "ci_rel_error":         sd(ci_d, ci_true)      if is_int else np.nan,
        "ci_abs_rel_error":     abs(sd(ci_d, ci_true))  if is_int else np.nan,
        "height_diff":          h_d,            "height_rel_error":     sd(h_d, h_sim)         if is_int else np.nan,
        "height_abs_rel_error": abs(sd(h_d, h_sim))     if is_int else np.nan,
        "bl_diff":              bl_d,           "bl_rel_error":         sd(bl_d, bl_sim),
        "bl_abs_rel_error":     abs(sd(bl_d, bl_sim)),
    }


# =============================================================================
# Time bin schemes
# =============================================================================

def get_bin_schemes(pop_model):
    schemes = {
        "medium5_200": list(range(0, 201, 5)) + [np.inf],
        "medium5_400": list(range(0, 401, 5)) + [np.inf],
        "coarse20": list(range(0, 181, 20)) + [np.inf],
    }
    if pop_model.startswith("bottleneck"):
        schemes["fine3"] = list(range(0, 31, 3)) + [np.inf]
        schemes["fine1"] = list(range(0, 31, 1)) + [np.inf]
    else:
        schemes["medium10"] = list(range(0, 401, 10)) + [np.inf]
    return schemes


def assign_bin(height, edges):
    for i in range(len(edges) - 1):
        if edges[i] <= height < edges[i + 1]:
            return i
    return len(edges) - 2


# =============================================================================
# Main
# =============================================================================

def load_beast_trees(path):
    """Parse BEAST nexus trees, stripping metadata annotations BioPython can't handle."""
    trees = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line.lower().startswith("tree "):
                continue
            # Remove optional [&...] block between tree name and '='
                newick = re.sub(r"\[&[^\]]*\]", "", line)
            eq = newick.find("=")
            if eq == -1:
                continue
            newick = newick[eq + 1:].strip().rstrip(";") + ";"
            trees.append(Phylo.read(StringIO(newick), "newick"))
    return trees


def load_log(path):
    df        = pd.read_csv(path, sep="\t", comment="#")
    state_col = next((c for c in ("State", "state", "Sample", "sample") if c in df.columns), None)
    if state_col is None:
        raise KeyError(f"No state/sample column found in {path}. Columns: {list(df.columns)}")
    return df[pd.to_numeric(df[state_col], errors="coerce").notna()].reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tsv",        required=True)
    parser.add_argument("--out_dir",    required=True)
    parser.add_argument("--num_groups", type=int, default=10)
    args = parser.parse_args()

    out_dir     = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    scenario_df = pd.read_csv(args.tsv, sep="\t")

    if scenario_df.empty:
        print(f"No successful runs in {args.tsv}, skipping.")
        return

    row0        = scenario_df.iloc[0]
    pop_model   = row0["population_model"]
    scenario    = f"{row0['sampling']}_{pop_model}_{row0['mutation_signal']}mutsig"
    bin_schemes = get_bin_schemes(pop_model)
    output_rows = []

    for _, rep_row in scenario_df.iterrows():
        tree_index = int(rep_row["tree_index"])

        # Simulated tree — one traversal, builds everything reused across all samples
        with open(rep_row["sim_tree_path"]) as fh:
            lines = [l.strip() for l in fh if l.strip()]
        sim_dict, events = traverse_tree(Phylo.read(StringIO(lines[tree_index]), "newick"))
        true_traj = get_true_traj(pop_model, rep_row)
        ci_true   = coalescent_intensity(
            events, lambda t0, t1, pm=pop_model, r=rep_row: integral_true(t0, t1, pm, r)
        )

        # Load subsampled files — skip replicate if either model's files are empty (ESS failed)
        model_data = {}
        skip = False
        for model in ("constcoal", "skyline"):
            tp = rep_row[f"trees_path_{model}"].replace(".combined.trees", ".subsampled.trees")
            lp = rep_row[f"log_path_{model}"].replace(".combined.log",     ".subsampled.log")
            if not Path(tp).exists() or Path(tp).stat().st_size == 0:
                skip = True; break
            log   = load_log(lp)
            trees = load_beast_trees(tp)
            model_data[model] = {"log": log, "trees": trees, "n": min(len(trees), len(log))}
        if skip:
            continue

        # Accumulators: model -> node_idx -> metric -> [values]
        acc = {m: {n: {met: [] for met in METRICS} for n in range(len(sim_dict))}
               for m in ("constcoal", "skyline")}

        for model in ("constcoal", "skyline"):
            for s_idx in range(model_data[model]["n"]):
                post_dict, _ = traverse_tree(model_data[model]["trees"][s_idx])
                log_row          = model_data[model]["log"].iloc[s_idx]

                if model == "skyline":
                    sky_b  = skyline_boundaries_from_preorder(post_dict, args.num_groups)
                    sky_Ne = skyline_pop_sizes(log_row, args.num_groups)
                    int_fn = lambda t0, t1, b=sky_b, ne=sky_Ne: integral_skyline(t0, t1, [0.0] + b, ne)
                else:
                    Ne_c   = float(log_row["constant.popSize"])
                    int_fn = lambda t0, t1, Ne=Ne_c: integral_constcoal(t0, t1, Ne)
                    sky_b  = sky_Ne = None

                ci_est = coalescent_intensity(events, int_fn)

                for node_idx, (h_est, bl_est, _, _) in post_dict.items():
                    h_sim, bl_sim, is_int, _ = sim_dict[node_idx]
                    Ne_est = (float(log_row["constant.popSize"]) if model == "constcoal"
                              else skyline_Ne_at(h_sim, sky_b, sky_Ne))
                    obs = compute_errors(
                        h_sim, bl_sim, h_est, bl_est, is_int,
                        true_traj(h_sim), Ne_est,
                        ci_true.get(node_idx, np.nan), ci_est.get(node_idx, np.nan)
                    )
                    for met in METRICS:
                        acc[model][node_idx][met].append(obs[met])

        for node_idx, (h_sim, bl_sim, is_int, node_id) in sim_dict.items():
            row = {
                "scenario": scenario, "tree_index": tree_index,
                "node_idx": node_idx, "node_id":    node_id,
                "internal": is_int,   "height_sim": h_sim, "bl_sim": bl_sim,
            }
            for scheme, edges in bin_schemes.items():
                row[f"bin_{scheme}"] = assign_bin(h_sim, edges)
            for model in ("constcoal", "skyline"):
                for met in METRICS:
                    med, lo, hi = summarize(acc[model][node_idx][met])
                    row[f"{model}_{met}_median"]    = med
                    row[f"{model}_{met}_hpd_lower"] = lo
                    row[f"{model}_{met}_hpd_upper"] = hi
            output_rows.append(row)

        print(f"  Processed replicate T{tree_index} ({len(sim_dict)} nodes)")

    pd.DataFrame(output_rows).to_csv(out_dir / f"{scenario}_node_errors.tsv", sep="\t", index=False)
    print(f"Done: {scenario} — {len(output_rows)} rows")


if __name__ == "__main__":
    main()
