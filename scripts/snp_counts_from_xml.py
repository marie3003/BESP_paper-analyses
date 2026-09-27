#!/usr/bin/env python3
"""
Count SNPs (= alignment sequence length) directly from BEAST XML files,
per sampling type / population model / mutation signature setting.

SNP counts are identical between constcoal and skyline XMLs for the same
replicate (same alignment), so only one model is scanned.

Usage:
    python scripts/snp_counts_from_xml.py \
        --beast_dir results/run1/beast_inference \
        --out_dir   results/run1/evaluation \
        --model     constcoal \
        --nreplicates 100
"""
import argparse
import csv
import re
from pathlib import Path

import numpy as np

SAMPLING_TYPES = ["independenthomochronous", "linearconstant"]
POP_MODELS     = ["expgrowthfast", "expgrowthslow", "uniform", "bottleneck"]
MUTSIGS        = ["lowmutsig", "medmutsig", "highmutsig"]

SEQ_RE = re.compile(r'<sequence><taxon idref="[^"]+"/>([ACGTNacgtn\-]+)</sequence>')


def count_snps(xml_path: Path) -> int | None:
    with open(xml_path) as fh:
        for line in fh:
            m = SEQ_RE.search(line)
            if m:
                return len(m.group(1))
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--beast_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--model", default="constcoal")
    ap.add_argument("--nreplicates", type=int, default=100)
    args = ap.parse_args()

    beast_dir = Path(args.beast_dir)
    out_dir   = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_rows = []

    for sampling in SAMPLING_TYPES:
        for popmodel in POP_MODELS:
            for mutsig in MUTSIGS:
                counts = []
                for i in range(args.nreplicates):
                    xml_path = (
                        beast_dir / args.model / sampling / popmodel / mutsig /
                        f"{args.model}_{sampling}_{popmodel}_{mutsig}.T{i}.xml"
                    )
                    if not xml_path.exists():
                        continue
                    n_snps = count_snps(xml_path)
                    if n_snps is not None:
                        counts.append((i, n_snps))

                if not counts:
                    continue

                # per-setting file with all replicate counts
                setting_out = out_dir / f"{sampling}_{popmodel}_{mutsig}_snp_counts.tsv"
                with open(setting_out, "w", newline="") as fh:
                    writer = csv.writer(fh, delimiter="\t")
                    writer.writerow(["replicate", "n_snps"])
                    writer.writerows(counts)

                vals = np.array([c for _, c in counts], dtype=float)
                mean_snps = float(vals.mean())
                # sample variance/SD (ddof=1); NaN for a single replicate, where
                # sample variance is undefined
                var_snps = float(vals.var(ddof=1)) if len(vals) > 1 else float("nan")
                std_snps = float(vals.std(ddof=1))  if len(vals) > 1 else float("nan")
                median_snps = float(np.median(vals))
                q2_5, q25, q75, q97_5 = np.percentile(vals, [2.5, 25, 75, 97.5])
                summary_rows.append(
                    [sampling, popmodel, mutsig, len(counts), mean_snps,
                     var_snps, std_snps, median_snps,
                     float(q2_5), float(q25), float(q75), float(q97_5)]
                )
                print(f"{sampling}/{popmodel}/{mutsig}: "
                      f"n={len(counts)} mean_snps={mean_snps:.2f} std={std_snps:.2f} -> {setting_out}")

    summary_out = out_dir / "snp_counts_summary.tsv"
    with open(summary_out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["sampling", "popmodel", "mutsig", "n_replicates", "mean_n_snps",
                          "var_n_snps", "std_n_snps", "median_n_snps",
                          "q2_5_n_snps", "q25_n_snps", "q75_n_snps", "q97_5_n_snps"])
        writer.writerows(summary_rows)

    print(f"\nSummary written to {summary_out}")


if __name__ == "__main__":
    main()
