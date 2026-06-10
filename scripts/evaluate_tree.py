"""
Per-tree evaluation script. For one (sampling, popmodel, mutsig, tree_index),
subsamples the combined log and trees to N_SAMPLES aligned posterior samples,
then computes per-sample per-node metrics vs the true simulated tree.

Outputs:
  node_samples.pkl  — one row per sample per internal node per model
  node_summary.pkl  — medians and 95% HPDs per node per model
"""

import argparse, pickle, re, sys
import numpy as np
import pandas as pd
from pathlib import Path
from Bio import Phylo
from io import StringIO

sys.path.insert(0, str(Path(__file__).parent.parent / "workflows"))
from evaluation_helpers import assign_model_params, compute_node_rates


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def read_log(path, n_samples):
    df = pd.read_csv(path, comment="#", sep="\t")
    idx = np.round(np.linspace(0, len(df) - 1, min(n_samples, len(df)))).astype(int)
    return df.iloc[idx].reset_index(drop=True)


def read_trees(path, states):
    """Read only the trees matching the given set of states from a BEAST nexus file."""
    states = set(states)
    translation = {}
    result = {}
    in_translate = False

    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.lower().startswith("translate"):
                in_translate = True
                continue
            if in_translate:
                if line == ";":
                    in_translate = False
                    continue
                k, v = line.rstrip(",;").split(None, 1)
                translation[k] = v
                continue
            m = re.match(r"tree\s+STATE_(\d+)\s*=\s*(.+)", line, re.IGNORECASE)
            if m:
                s = int(m.group(1))
                if s in states:
                    nwk = m.group(2)
                    if translation:
                        pat = re.compile(r"\b(" + "|".join(re.escape(k) for k in translation) + r")\b")
                        nwk = pat.sub(lambda mo: translation[mo.group(0)], nwk)
                    result[s] = nwk

    return result


# ---------------------------------------------------------------------------
# Node info from a single posterior tree
# ---------------------------------------------------------------------------

def node_info(tree):
    depths = tree.depths()
    max_d = max(depths[t] for t in tree.get_terminals())
    rows, ic = [], 0
    for clade in tree.find_clades(order="preorder"):
        if clade.name:
            nid, internal = clade.name, False
        else:
            nid, internal = f"internal_{ic}", True
            ic += 1
        rows.append({
            "node":     nid,
            "height":   max_d - depths[clade],
            "bl":       clade.branch_length or 0.0,
            "internal": internal,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Population size lookup
# ---------------------------------------------------------------------------

def ne_skyline(log_row, t):
    cols = sorted([c for c in log_row.index if c.startswith("skyline.popSize")],
                  key=lambda c: int(re.search(r"\d+$", c).group()))
    root_h = log_row.get("treeModel.rootHeight", 1000)
    bounds = np.linspace(0, root_h, len(cols) + 1)
    for i, col in enumerate(cols):
        if bounds[i] <= t < bounds[i + 1]:
            return log_row[col]
    return log_row[cols[-1]]


def skyline_params(log_row):
    cols = sorted([c for c in log_row.index if c.startswith("skyline.popSize")],
                  key=lambda c: int(re.search(r"\d+$", c).group()))
    root_h = log_row.get("treeModel.rootHeight", 1000)
    bounds = list(np.linspace(0, root_h, len(cols) + 1)[1:])
    return {"time_points": bounds, "population_sizes": [log_row[c] for c in cols]}


# ---------------------------------------------------------------------------
# HPD
# ---------------------------------------------------------------------------

def hpd(vals, frac=0.95):
    v = np.sort(vals[~np.isnan(vals)])
    n = len(v)
    if n == 0:
        return np.nan, np.nan
    w = max(int(np.floor(frac * n)), 1)
    if w >= n:
        return v[0], v[-1]
    widths = v[w:] - v[:n - w]
    i = np.argmin(widths)
    return v[i], v[i + w]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--constcoal_log");   p.add_argument("--constcoal_trees")
    p.add_argument("--skyline_log");     p.add_argument("--skyline_trees")
    p.add_argument("--sim_trees");       p.add_argument("--tree_index", type=int)
    p.add_argument("--pop_model");       p.add_argument("--mutation_signal")
    p.add_argument("--sampling");        p.add_argument("--n_samples", type=int, default=1000)
    p.add_argument("--out_samples");     p.add_argument("--out_summary")
    args = p.parse_args()

    params = assign_model_params(args.pop_model)
    pop_model = args.pop_model

    if pop_model in ("expgrowthfast", "expgrowthslow"):
        true_mode = "true_exp"
        true_params = {"N0": params["present_pop_size"], "rate": params["growth_rate"]}
    elif pop_model.startswith("bottleneck"):
        true_mode = "true_bottleneck"
        true_params = {"Ne_normal": params["present_pop_size"],
                       "Ne_bottleneck": params["bottleneck_size"],
                       "bottleneck_start": params["bottleneck_start"],
                       "bottleneck_end": params["bottleneck_end"]}
    else:
        true_mode = "true_uniform"
        true_params = {"Ne": params["present_pop_size"]}

    # true tree
    sim_tree = list(Phylo.parse(args.sim_trees, "newick"))[args.tree_index]
    true_df = node_info(sim_tree).rename(columns={"height": "height_sim", "bl": "bl_sim"})
    true_node_times = dict(zip(true_df["node"], true_df["height_sim"]))

    # rate_true at true node times
    rate_true_map = compute_node_rates(
        {n: t for n, t in true_node_times.items() if true_df.set_index("node").loc[n, "internal"]},
        true_mode, true_params
    )

    meta = dict(tree_index=args.tree_index, population_model=pop_model,
                mutation_signal=args.mutation_signal, sampling=args.sampling)

    all_rows = []

    for model, log_path, trees_path in [
        ("constcoal", args.constcoal_log, args.constcoal_trees),
        ("skyline",   args.skyline_log,   args.skyline_trees),
    ]:
        print(f"Processing {model}...")
        log_df = read_log(log_path, args.n_samples)
        states = log_df["state"].astype(int).tolist()
        trees_map = read_trees(trees_path, set(states))

        for sample_idx, (_, log_row) in enumerate(log_df.iterrows()):
            state = int(log_row["state"])
            nwk = trees_map.get(state)
            if nwk is None:
                continue
            try:
                tree = Phylo.read(StringIO(nwk), "newick")
            except Exception:
                continue

            ndf = node_info(tree)

            for _, nr in ndf[ndf["internal"]].iterrows():
                nid = nr["node"]
                t_est = nr["height"]

                if model == "skyline":
                    ne_est = ne_skyline(log_row, t_est)
                    sp = skyline_params(log_row)
                    t_true = true_node_times.get(nid, np.nan)
                    rate_est = compute_node_rates({nid: t_true}, "sim_skyline", sp).get(nid, np.nan)
                else:
                    ne_est = log_row.get("constant.popSize", np.nan)
                    rate_est = (1.0 / ne_est) if ne_est else np.nan

                t_true = true_node_times.get(nid, np.nan)
                h_sim  = true_df.set_index("node").loc[nid, "height_sim"] if nid in true_df["node"].values else np.nan
                bl_sim = true_df.set_index("node").loc[nid, "bl_sim"]     if nid in true_df["node"].values else np.nan
                rt     = rate_true_map.get(nid, np.nan)
                ne_true = (1.0 / rt) if rt else np.nan

                all_rows.append({
                    **meta,
                    "model": model, "node": nid, "sample_index": sample_idx,
                    "height": t_est, "bl": nr["bl"],
                    "height_sim": h_sim, "bl_sim": bl_sim,
                    "ne_est": ne_est, "ne_true": ne_true,
                    "rate_est": rate_est, "rate_true": rt,
                    "height_diff":   t_est - h_sim,
                    "height_rel_error": (t_est - h_sim) / h_sim if h_sim else np.nan,
                    "bl_diff":       nr["bl"] - bl_sim,
                    "bl_rel_error":  (nr["bl"] - bl_sim) / bl_sim if bl_sim else np.nan,
                    "ne_diff":       ne_est - ne_true,
                    "ne_rel_error":  (ne_est - ne_true) / ne_true if ne_true else np.nan,
                    "rate_diff":     rate_est - rt,
                    "rate_rel_error":(rate_est - rt) / rt if rt else np.nan,
                })

    samples_df = pd.DataFrame(all_rows)

    # summary: median + 95% HPD per node per model
    metric_cols = ["height", "bl", "height_diff", "height_rel_error",
                   "bl_diff", "bl_rel_error", "ne_est", "ne_true",
                   "rate_est", "rate_true", "ne_diff", "ne_rel_error",
                   "rate_diff", "rate_rel_error"]

    summary_rows = []
    for (node, model), grp in samples_df.groupby(["node", "model"]):
        row = {**meta, "node": node, "model": model}
        for col in metric_cols:
            v = grp[col].values.astype(float)
            row[f"{col}_median"] = np.nanmedian(v)
            lo, hi = hpd(v)
            row[f"{col}_hpd_lower"] = lo
            row[f"{col}_hpd_upper"] = hi
        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)

    print(f"Writing {len(samples_df)} sample rows, {len(summary_df)} summary rows")
    with open(args.out_samples, "wb") as f: pickle.dump(samples_df, f)
    with open(args.out_summary, "wb") as f: pickle.dump(summary_df, f)
    print("Done.")


if __name__ == "__main__":
    main()
