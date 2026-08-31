"""Aggregate GSE131907 epithelial-cell scores to independent sample units."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

from common import ensure_parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Per-cell CSV containing Sample, Cell_type, ARID1A and coreBAF")
    parser.add_argument("--sample-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--cell-type", default="Epithelial cells")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data = pd.read_csv(args.input, usecols=["Sample", "Cell_type", "ARID1A", "coreBAF"])
    epithelial = data.loc[data["Cell_type"] == args.cell_type].copy()
    by_sample = epithelial.groupby("Sample", as_index=False).agg(
        n_cells=("ARID1A", "size"), ARID1A_mean=("ARID1A", "mean"), coreBAF_mean=("coreBAF", "mean")
    )
    by_sample["ARID1A_group"] = "high"
    by_sample.loc[by_sample["ARID1A_mean"] <= by_sample["ARID1A_mean"].median(), "ARID1A_group"] = "low"
    rho, p_rho = spearmanr(by_sample["ARID1A_mean"], by_sample["coreBAF_mean"])
    low = by_sample.loc[by_sample.ARID1A_group == "low", "coreBAF_mean"]
    high = by_sample.loc[by_sample.ARID1A_group == "high", "coreBAF_mean"]
    u, p_group = mannwhitneyu(low, high, alternative="two-sided", method="exact")
    summary = pd.DataFrame([{
        "n_cells": len(epithelial), "n_samples": len(by_sample),
        "spearman_rho_ARID1A_vs_coreBAF": rho, "spearman_p": p_rho,
        "median_coreBAF_low_samples": low.median(), "median_coreBAF_high_samples": high.median(),
        "mannwhitney_u": u, "mannwhitney_p": p_group,
    }])
    by_sample.to_csv(ensure_parent(args.sample_output), index=False)
    summary.to_csv(ensure_parent(args.summary_output), index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()

