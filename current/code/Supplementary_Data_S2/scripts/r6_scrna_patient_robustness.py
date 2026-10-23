"""Patient-aware robustness analysis for the GSE131907 epithelial result.

The inferential unit is the patient.  Sample means are calculated first and
then given equal weight within patient, so patients contributing more than one
specimen do not receive extra weight.  A direct cell-weighted patient mean is
reported as a sensitivity analysis.  The primary tumour-epithelial definition
combines the source labels ``Malignant cells`` and tS1/tS2/tS3: the former occur
in metastatic brain/lymph-node and tL/B specimens, and the latter in tumour-lung
specimens.  The literal ``Malignant cells`` label and all epithelial cells are
reported as narrower and broader sensitivity definitions, respectively.
"""

from __future__ import annotations

import argparse
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import mannwhitneyu, spearmanr


SEED = 20260902
N_BOOT = 20_000
DEFINITIONS = {
    "strict_malignant": {"Malignant cells"},
    "tumour_associated": {"Malignant cells", "tS1", "tS2", "tS3"},
    "all_epithelial": None,
}


def bootstrap_rho(frame: pd.DataFrame, n_boot: int, seed: int) -> tuple[float, float]:
    values = frame[["ARID1A_mean", "coreBAF_mean"]].to_numpy(float)
    rng = np.random.default_rng(seed)
    estimates = []
    for _ in range(n_boot):
        boot = values[rng.integers(0, len(values), len(values))]
        rho = spearmanr(boot[:, 0], boot[:, 1]).statistic
        if np.isfinite(rho):
            estimates.append(float(rho))
    if len(estimates) < n_boot * 0.95:
        raise RuntimeError("Too many undefined bootstrap correlations")
    return tuple(np.quantile(estimates, [0.025, 0.975]).tolist())


def summarize(frame: pd.DataFrame, definition: str, weighting: str) -> dict:
    rho_test = spearmanr(frame["ARID1A_mean"], frame["coreBAF_mean"])
    median = frame["ARID1A_mean"].median()
    low = frame.loc[frame.ARID1A_mean <= median, "coreBAF_mean"]
    high = frame.loc[frame.ARID1A_mean > median, "coreBAF_mean"]
    mw = mannwhitneyu(low, high, alternative="two-sided", method="exact")
    ci_low, ci_high = bootstrap_rho(
        frame, N_BOOT, SEED + list(DEFINITIONS).index(definition) * 10 + (weighting == "cell_weighted")
    )
    return {
        "definition": definition,
        "patient_weighting": weighting,
        "n_cells": int(frame["n_cells"].sum()),
        "n_samples": int(frame["n_samples"].sum()),
        "n_patients": int(len(frame)),
        "spearman_rho": float(rho_test.statistic),
        "spearman_p": float(rho_test.pvalue),
        "spearman_bootstrap_ci_low": ci_low,
        "spearman_bootstrap_ci_high": ci_high,
        "median_coreBAF_low": float(low.median()),
        "median_coreBAF_high": float(high.median()),
        "mannwhitney_u": float(mw.statistic),
        "mannwhitney_p": float(mw.pvalue),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument(
        "--patient-map", required=True,
        help="Source-publication mapping with Sample and patient_id columns",
    )
    parser.add_argument("--patient-map-provenance", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    usecols = ["Sample", "Cell_type", "Cell_subtype", "ARID1A", "coreBAF"]
    data = pd.read_csv(args.input, usecols=usecols)
    epi = data.loc[data.Cell_type.eq("Epithelial cells")].copy()
    patient_map = pd.read_csv(args.patient_map, usecols=["Sample", "patient_id"])
    if patient_map.Sample.duplicated().any():
        raise ValueError("Patient map contains duplicate sample names")
    epi = epi.merge(patient_map, on="Sample", how="left", validate="many_to_one")
    missing_samples = sorted(epi.loc[epi.patient_id.isna(), "Sample"].unique())
    if missing_samples:
        raise ValueError(f"Samples absent from source-publication patient map: {missing_samples}")
    if epi[["Sample", "ARID1A", "coreBAF", "patient_id"]].isna().any().any():
        raise ValueError("Missing values detected in required quantitative or identifier fields")
    # Missing refined labels remain explicit in the all-epithelial sensitivity
    # analysis and can never enter either malignant-cell definition.
    epi["Cell_subtype"] = epi.Cell_subtype.fillna("Unannotated")

    subtype_counts = (
        epi.groupby("Cell_subtype", observed=True)
        .agg(n_cells=("Sample", "size"), n_samples=("Sample", "nunique"),
             n_patients=("patient_id", "nunique"))
        .reset_index()
        .sort_values("n_cells", ascending=False)
    )
    subtype_counts.to_csv(out / "r6_scrna_subtype_counts.csv", index=False)

    mapping = epi[["patient_id", "Sample"]].drop_duplicates().sort_values(["patient_id", "Sample"])
    mapping.to_csv(out / "r6_scrna_patient_sample_mapping.csv", index=False)

    score_frames = []
    summary_rows = []
    for definition, allowed in DEFINITIONS.items():
        subset = epi if allowed is None else epi.loc[epi.Cell_subtype.isin(allowed)]
        sample = (
            subset.groupby(["patient_id", "Sample"], observed=True)
            .agg(ARID1A_mean=("ARID1A", "mean"), coreBAF_mean=("coreBAF", "mean"),
                 n_cells=("Sample", "size"))
            .reset_index()
        )

        equal_sample = (
            sample.groupby("patient_id", observed=True)
            .agg(ARID1A_mean=("ARID1A_mean", "mean"), coreBAF_mean=("coreBAF_mean", "mean"),
                 n_cells=("n_cells", "sum"), n_samples=("Sample", "nunique"))
            .reset_index()
        )
        cell_weighted = (
            subset.groupby("patient_id", observed=True)
            .agg(ARID1A_mean=("ARID1A", "mean"), coreBAF_mean=("coreBAF", "mean"),
                 n_cells=("Sample", "size"), n_samples=("Sample", "nunique"))
            .reset_index()
        )
        for weighting, frame in [("equal_sample", equal_sample), ("cell_weighted", cell_weighted)]:
            if len(frame) < 10:
                raise ValueError(f"Too few patients for {definition}/{weighting}: {len(frame)}")
            tagged = frame.assign(definition=definition, patient_weighting=weighting)
            score_frames.append(tagged)
            summary_rows.append(summarize(frame, definition, weighting))

    scores = pd.concat(score_frames, ignore_index=True)
    summary = pd.DataFrame(summary_rows)
    scores.to_csv(out / "r6_scrna_patient_scores.csv", index=False)
    summary.to_csv(out / "r6_scrna_patient_summary.csv", index=False)

    patient_map_trace = json.loads(Path(args.patient_map_provenance).read_text(encoding="utf-8"))
    trace = {
        "analysis": "GSE131907 patient-aware epithelial robustness",
        "source_file": str(Path(args.input).resolve()),
        "source_rows": int(len(data)),
        "epithelial_rows": int(len(epi)),
        "epithelial_samples": int(epi.Sample.nunique()),
        "epithelial_patients": int(epi.patient_id.nunique()),
        "patient_mapping": patient_map_trace,
        "primary_definition": "source-labelled tumour epithelial states: Malignant cells + tS1 + tS2 + tS3",
        "sensitivity_definitions": ["literal Malignant cells label", "all epithelial cells"],
        "primary_patient_weighting": "equal weight to specimen means within patient",
        "bootstrap_resamples": N_BOOT,
        "random_seed": SEED,
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
    }
    (out / "r6_scrna_run_configuration.json").write_text(
        json.dumps(trace, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
