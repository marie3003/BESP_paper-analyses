"""
Aggregate all per-tree node_summary.pkl files into a single tree_metrics_combined.pkl.
Also concatenates all node_samples.pkl files into a full samples dataframe.

Usage:
    python aggregate_tree_metrics.py \
        --summary_dir <dir containing node_summary.pkl files> \
        --out_combined <path for tree_metrics_combined.pkl> \
        --out_samples  <path for tree_samples_all.pkl>  (optional, large)
"""

import argparse
import pickle
import pandas as pd
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary_dir",   required=True)
    parser.add_argument("--out_combined",  required=True)
    parser.add_argument("--out_samples",   default=None)
    args = parser.parse_args()

    summary_dir = Path(args.summary_dir)
    summary_files = sorted(summary_dir.rglob("*.node_summary.pkl"))
    print(f"Found {len(summary_files)} summary files")

    summaries = []
    for f in summary_files:
        with open(f, "rb") as fh:
            summaries.append(pickle.load(fh))
    combined = pd.concat(summaries, ignore_index=True)
    with open(args.out_combined, "wb") as fh:
        pickle.dump(combined, fh)
    print(f"Saved combined summary: {len(combined)} rows → {args.out_combined}")

    if args.out_samples:
        sample_files = sorted(summary_dir.rglob("*.node_samples.pkl"))
        print(f"Found {len(sample_files)} sample files")
        samples = []
        for f in sample_files:
            with open(f, "rb") as fh:
                samples.append(pickle.load(fh))
        all_samples = pd.concat(samples, ignore_index=True)
        with open(args.out_samples, "wb") as fh:
            pickle.dump(all_samples, fh)
        print(f"Saved all samples: {len(all_samples)} rows → {args.out_samples}")


if __name__ == "__main__":
    main()
