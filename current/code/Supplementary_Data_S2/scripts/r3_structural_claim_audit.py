"""Audit the structural-contact claim used in the ARID1A/BAF manuscript.

The primary structural result is a binary comparison between ARID1A-contacting
and non-contacting subunits.  This script reports an effect size and bootstrap
interval for that comparison, and separately tests whether contact count shows
a graded relationship with co-dependency.  The graded analysis is restricted
to chains resolved in the same cryo-EM structure (PDB 6LTJ).
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr


ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "results_tables"
SEED = 20260831
N_BOOT = 20_000


def bootstrap_median_difference(contacting, noncontacting, rng):
    contacting = np.asarray(contacting, dtype=float)
    noncontacting = np.asarray(noncontacting, dtype=float)
    estimates = np.empty(N_BOOT)
    for i in range(N_BOOT):
        a = rng.choice(contacting, size=len(contacting), replace=True)
        b = rng.choice(noncontacting, size=len(noncontacting), replace=True)
        estimates[i] = np.median(a) - np.median(b)
    return np.quantile(estimates, [0.025, 0.975])


def main():
    merged = pd.read_csv(TABLES / "structural_vs_codependency_merged.csv")

    # SMARCC1 is not resolved in 6LTJ, but the independent 9A0K integrative
    # model gives a binary contact call. Its Project Score correlations come
    # from the same genome-wide tables as the other subunits.
    lung = pd.read_csv(TABLES / "genomewide_codependency_lung.csv").set_index("gene")
    pan = pd.read_csv(TABLES / "genomewide_codependency_pancancer.csv").set_index("gene")
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
    for cohort, value_col in [("lung", "lung_rho"), ("pan-cancer", "pancancer_rho")]:
        a = analysis.loc[analysis.contacting, value_col].dropna().to_numpy()
        b = analysis.loc[~analysis.contacting, value_col].dropna().to_numpy()
        u, p = mannwhitneyu(a, b, alternative="two-sided", method="exact")
        ci_low, ci_high = bootstrap_median_difference(a, b, rng)
        rows.append({
            "analysis": "binary_contact_vs_noncontact",
            "cohort": cohort,
            "n_contacting": len(a),
            "n_noncontacting": len(b),
            "estimate": np.median(a) - np.median(b),
            "estimate_type": "median_rho_difference",
            "ci95_low": ci_low,
            "ci95_high": ci_high,
            "statistic": u,
            "effect_size": 2 * u / (len(a) * len(b)) - 1,
            "effect_size_type": "Cliffs_delta",
            "p_value": p,
        })

        # Dose-response claim: use only the six positive-contact chains measured
        # within the same structure. Mixing 9A0K and 6LTJ counts is avoided.
        same_structure = merged.loc[merged["n_contacts_4.5A"] > 0]
        rho, p_rho = spearmanr(same_structure["n_contacts_4.5A"], same_structure[value_col])
        rows.append({
            "analysis": "graded_contact_count_within_6LTJ_contacts",
            "cohort": cohort,
            "n_contacting": len(same_structure),
            "n_noncontacting": np.nan,
            "estimate": rho,
            "estimate_type": "Spearman_rho",
            "ci95_low": np.nan,
            "ci95_high": np.nan,
            "statistic": np.nan,
            "effect_size": np.nan,
            "effect_size_type": "not_applicable",
            "p_value": p_rho,
        })

    out = pd.DataFrame(rows)
    out.to_csv(TABLES / "r3_structural_claim_audit.csv", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
