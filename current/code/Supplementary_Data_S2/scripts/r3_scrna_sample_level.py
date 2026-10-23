"""Sample-level sensitivity analysis for the GSE131907 epithelial cells.

Cells are nested within samples, so cell-level tests overstate precision.  This
analysis aggregates the existing per-cell ARID1A expression and core-BAF score
to one observation per sample before testing their association.
"""

from pathlib import Path

import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr


ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "results_tables"


def main():
    data = pd.read_csv(
        TABLES / "gse131907_harmony.csv",
        usecols=["Sample", "Cell_type", "ARID1A", "coreBAF"],
    )
    epithelial = data.loc[data["Cell_type"] == "Epithelial cells"].copy()
    by_sample = (
        epithelial.groupby("Sample", as_index=False)
        .agg(n_cells=("ARID1A", "size"), ARID1A_mean=("ARID1A", "mean"),
             coreBAF_mean=("coreBAF", "mean"))
    )
    by_sample["ARID1A_group"] = "high"
    by_sample.loc[
        by_sample["ARID1A_mean"] <= by_sample["ARID1A_mean"].median(),
        "ARID1A_group",
    ] = "low"

    rho, p_rho = spearmanr(by_sample["ARID1A_mean"], by_sample["coreBAF_mean"])
    low = by_sample.loc[by_sample.ARID1A_group == "low", "coreBAF_mean"]
    high = by_sample.loc[by_sample.ARID1A_group == "high", "coreBAF_mean"]
    u, p_group = mannwhitneyu(low, high, alternative="two-sided", method="exact")

    summary = pd.DataFrame([
        {
            "n_cells": len(epithelial),
            "n_samples": len(by_sample),
            "spearman_rho_ARID1A_vs_coreBAF": rho,
            "spearman_p": p_rho,
            "median_coreBAF_low_samples": low.median(),
            "median_coreBAF_high_samples": high.median(),
            "mannwhitney_u": u,
            "mannwhitney_p": p_group,
        }
    ])
    by_sample.to_csv(TABLES / "r3_scrna_epithelial_sample_scores.csv", index=False)
    summary.to_csv(TABLES / "r3_scrna_sample_level_summary.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
