"""Recompute the binary and graded structure/co-dependency claim audit."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

from common import ensure_parent


ROOT = Path(__file__).parents[1]
SEED = 20260831


def bootstrap_median_difference(a: np.ndarray, b: np.ndarray, rng: np.random.Generator, n_boot: int) -> np.ndarray:
    estimates = np.empty(n_boot)
    for i in range(n_boot):
        estimates[i] = np.median(rng.choice(a, len(a), replace=True)) - np.median(
            rng.choice(b, len(b), replace=True)
        )
    return np.quantile(estimates, [0.025, 0.975])


def compute(tables: Path, n_boot: int) -> pd.DataFrame:
    merged = pd.read_csv(tables / "structural_vs_codependency_merged.csv")
    lung = pd.read_csv(tables / "genomewide_codependency_lung.csv").set_index("gene")
    pan = pd.read_csv(tables / "genomewide_codependency_pancancer.csv").set_index("gene")
    smarcc1 = pd.DataFrame(
        [{
            "gene": "SMARCC1",
            "n_contacts_4.5A": 186,
            "min_dist_A": 2.58,
            "pancancer_rho": pan.loc["SMARCC1", "rho"],
            "lung_rho": lung.loc["SMARCC1", "rho"],
            "structure": "9A0K",
        }]
    )
    merged["structure"] = "6LTJ"
    analysis = pd.concat([merged, smarcc1], ignore_index=True, sort=False)
    analysis["contacting"] = analysis["n_contacts_4.5A"] > 0
    rng = np.random.default_rng(SEED)
    rows = []
    for cohort, column in (("lung", "lung_rho"), ("pan-cancer", "pancancer_rho")):
        a = analysis.loc[analysis.contacting, column].dropna().to_numpy()
        b = analysis.loc[~analysis.contacting, column].dropna().to_numpy()
        u, pvalue = mannwhitneyu(a, b, alternative="two-sided", method="exact")
        low, high = bootstrap_median_difference(a, b, rng, n_boot)
        rows.append({
            "analysis": "binary_contact_vs_noncontact", "cohort": cohort,
            "n_contacting": len(a), "n_noncontacting": len(b),
            "estimate": np.median(a) - np.median(b), "estimate_type": "median_rho_difference",
            "ci95_low": low, "ci95_high": high, "statistic": u,
            "effect_size": 2 * u / (len(a) * len(b)) - 1,
            "effect_size_type": "Cliffs_delta", "p_value": pvalue,
        })
        same_structure = merged.loc[merged["n_contacts_4.5A"] > 0]
        rho, p_rho = spearmanr(same_structure["n_contacts_4.5A"], same_structure[column])
        rows.append({
            "analysis": "graded_contact_count_within_6LTJ_contacts", "cohort": cohort,
            "n_contacting": len(same_structure), "n_noncontacting": np.nan,
            "estimate": rho, "estimate_type": "Spearman_rho",
            "ci95_low": np.nan, "ci95_high": np.nan, "statistic": np.nan,
            "effect_size": np.nan, "effect_size_type": "not_applicable", "p_value": p_rho,
        })
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tables", type=Path, default=ROOT / "results" / "derived")
    parser.add_argument("--output", default=ROOT / "results" / "recomputed" / "r3_structural_claim_audit.csv")
    parser.add_argument("--bootstrap", type=int, default=20000)
    parser.add_argument("--verify", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = compute(args.tables, args.bootstrap)
    output.to_csv(ensure_parent(args.output), index=False)
    if args.verify:
        expected = pd.read_csv(args.tables / "r3_structural_claim_audit.csv")
        for column in ("estimate", "ci95_low", "ci95_high", "effect_size", "p_value"):
            if not np.allclose(output[column], expected[column], equal_nan=True, rtol=0, atol=1e-12):
                raise AssertionError(f"Archived and recomputed values differ in {column}")
    print(output.to_string(index=False))


if __name__ == "__main__":
    main()

