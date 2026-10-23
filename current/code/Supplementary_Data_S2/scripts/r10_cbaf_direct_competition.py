#!/usr/bin/env python3
"""Direct cBAF-versus-PBAF/ncBAF specificity tests.

The analysis addresses the inferential gap between "cBAF is significant while
another module is not" and a direct test that their associations differ.  It
uses (i) paired cell-line bootstrap in DepMap 26Q1, (ii) lineage as the unit of
replication for archived within-lineage correlations, and (iii) lineage-
stratified label permutation for ARID1A-state interactions in Project DRIVE.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import platform
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.stats import mannwhitneyu, norm, rankdata, spearmanr, wilcoxon
from statsmodels.stats.multitest import multipletests


SEED = 20260902
N_BOOT = 20_000
N_PERM = 20_000
MODULES = {
    "cBAF_specific": ["ARID1B", "DPF1", "DPF2", "DPF3", "BCL7A", "BCL7B", "BCL7C", "SS18", "SS18L1"],
    "PBAF_specific": ["PBRM1", "ARID2", "BRD7", "PHF10"],
    "ncBAF_specific": ["BRD9", "BICRA", "BICRAL"],
}
COMPARATORS = ["PBAF_specific", "ncBAF_specific"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def rank_normal(values: np.ndarray) -> np.ndarray:
    ranks = rankdata(values, method="average")
    return norm.ppf((ranks - 0.5) / len(ranks))


def spearman_one_sided_greater(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    test = spearmanr(x, y)
    rho = float(test.statistic)
    p = float(test.pvalue / 2 if rho >= 0 else 1 - test.pvalue / 2)
    return rho, p


def depmap_direct_tests(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    rng = np.random.default_rng(SEED)
    for comparator in COMPARATORS:
        cols = ["ARID1A", *MODULES["cBAF_specific"], *MODULES[comparator]]
        use = frame[cols].dropna().copy()
        x = use["ARID1A"].to_numpy(float)
        yc = use[MODULES["cBAF_specific"]].mean(axis=1).to_numpy(float)
        yq = use[MODULES[comparator]].mean(axis=1).to_numpy(float)
        rho_c = float(spearmanr(x, yc).statistic)
        rho_q = float(spearmanr(x, yq).statistic)
        delta = rho_c - rho_q
        contrast = rank_normal(yc) - rank_normal(yq)
        contrast_rho, contrast_p = spearman_one_sided_greater(x, contrast)
        boot = np.empty(N_BOOT)
        n = len(use)
        for i in range(N_BOOT):
            idx = rng.integers(0, n, n)
            boot[i] = (spearmanr(x[idx], yc[idx]).statistic
                       - spearmanr(x[idx], yq[idx]).statistic)
        boot = boot[np.isfinite(boot)]
        ci = np.quantile(boot, [0.025, 0.975])
        p_boot = float((1 + np.sum(boot <= 0)) / (len(boot) + 1))
        rows.append({
            "dataset": "DepMap_26Q1", "cohort": "pan_cancer", "comparator": comparator,
            "n_models_complete_case": n, "rho_ARID1A_cBAF": rho_c,
            "rho_ARID1A_comparator": rho_q, "delta_rho_cBAF_minus_comparator": delta,
            "bootstrap_ci_low": float(ci[0]), "bootstrap_ci_high": float(ci[1]),
            "bootstrap_tail_p_delta_le_zero": p_boot,
            "rho_ARID1A_rank_normalized_module_contrast": contrast_rho,
            "spearman_p_one_sided_contrast_greater": contrast_p,
        })
    result = pd.DataFrame(rows)
    result["holm_p_direct_contrast"] = multipletests(
        result["spearman_p_one_sided_contrast_greater"], method="holm"
    )[1]
    return result


def exact_signflip_mean_greater(values: np.ndarray) -> float:
    values = np.asarray(values, float)
    observed = float(values.mean())
    n = len(values)
    total = 1 << n
    exceed = 0
    chunk = 100_000
    shifts = np.arange(n, dtype=np.uint64)
    for start in range(0, total, chunk):
        masks = np.arange(start, min(start + chunk, total), dtype=np.uint64)[:, None]
        signs = (((masks >> shifts) & 1).astype(np.int8) * 2 - 1)
        exceed += int(np.sum((signs * values).mean(axis=1) >= observed - 1e-15))
    return exceed / total


def lineage_direct_tests(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    long_rows, summary_rows = [], []
    rng = np.random.default_rng(SEED + 100)
    for comparator, column in [("PBAF_specific", "rho_PBAF"), ("ncBAF_specific", "rho_ncBAF")]:
        use = frame[["lineage", "n", "rho_cBAF", column]].dropna().copy()
        use["comparator"] = comparator
        use["delta_rho_cBAF_minus_comparator"] = use["rho_cBAF"] - use[column]
        long_rows.append(use[["lineage", "n", "comparator", "rho_cBAF", column,
                              "delta_rho_cBAF_minus_comparator"]].rename(columns={column: "rho_comparator"}))
        values = use["delta_rho_cBAF_minus_comparator"].to_numpy(float)
        boot = rng.choice(values, size=(N_BOOT, len(values)), replace=True).mean(axis=1)
        ci = np.quantile(boot, [0.025, 0.975])
        summary_rows.append({
            "dataset": "DepMap_24Q4_within_lineage", "comparator": comparator,
            "n_lineages": len(values), "mean_delta_rho": float(values.mean()),
            "median_delta_rho": float(np.median(values)),
            "positive_lineages": int(np.sum(values > 0)),
            "bootstrap_mean_ci_low": float(ci[0]), "bootstrap_mean_ci_high": float(ci[1]),
            "p_exact_signflip_mean_greater": exact_signflip_mean_greater(values),
            "p_wilcoxon_one_sided_greater": float(wilcoxon(values, alternative="greater").pvalue),
        })
    summary = pd.DataFrame(summary_rows)
    summary["holm_p_signflip"] = multipletests(
        summary["p_exact_signflip_mean_greater"], method="holm"
    )[1]
    return pd.concat(long_rows, ignore_index=True), summary


def parse_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.strip().str.lower().eq("true")


def stratified_permutation_p(frame: pd.DataFrame, state: str, contrast: str,
                             observed: float, rng: np.random.Generator) -> float:
    labels = parse_bool(frame[state]).to_numpy()
    values = frame[contrast].to_numpy(float)
    lineages = frame["lineage"].astype(str).to_numpy()
    blocks = [np.flatnonzero(lineages == name) for name in np.unique(lineages)]
    null = np.empty(N_PERM)
    for i in range(N_PERM):
        shuffled = labels.copy()
        for idx in blocks:
            shuffled[idx] = rng.permutation(shuffled[idx])
        null[i] = values[shuffled].mean() - values[~shuffled].mean()
    return float((1 + np.sum(null <= observed)) / (N_PERM + 1))


def drive_direct_tests(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    rng = np.random.default_rng(SEED + 200)
    states = ["all_ARID1A_mutation", "predicted_LOF_mutation"]
    for state, comparator in itertools.product(states, COMPARATORS):
        contrast_name = f"cBAF_minus_{comparator}"
        use = frame[["lineage", state, "cBAF_specific", comparator]].dropna().copy()
        use[contrast_name] = use["cBAF_specific"] - use[comparator]
        labels = parse_bool(use[state])
        altered = use.loc[labels, contrast_name].to_numpy(float)
        other = use.loc[~labels, contrast_name].to_numpy(float)
        observed = float(altered.mean() - other.mean())
        boot = np.empty(N_BOOT)
        for i in range(N_BOOT):
            boot[i] = (rng.choice(altered, len(altered), replace=True).mean()
                       - rng.choice(other, len(other), replace=True).mean())
        ci = np.quantile(boot, [0.025, 0.975])
        mw = mannwhitneyu(altered, other, alternative="less")
        p_perm = stratified_permutation_p(use, state, contrast_name, observed, rng)
        rows.append({
            "dataset": "Project_DRIVE_RNAi", "state": state, "comparator": comparator,
            "n_altered": len(altered), "n_other": len(other),
            "interaction_difference_of_differences": observed,
            "specificity_advantage_positive_supports_cBAF": -observed,
            "bootstrap_ci_low_raw_difference": float(ci[0]),
            "bootstrap_ci_high_raw_difference": float(ci[1]),
            "mannwhitney_p_one_sided_less": float(mw.pvalue),
            "lineage_stratified_permutation_p_one_sided_less": p_perm,
        })
    result = pd.DataFrame(rows)
    result["holm_p_permutation_within_state"] = result.groupby("state", observed=True)[
        "lineage_stratified_permutation_p_one_sided_less"
    ].transform(lambda x: multipletests(x, method="holm")[1])
    return result


def make_figure(depmap: pd.DataFrame, lineage_long: pd.DataFrame,
                lineage_summary: pd.DataFrame, drive: pd.DataFrame, out: Path) -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 7.0, "axes.linewidth": 0.7,
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.facecolor": "white",
    })
    colors = {"PBAF_specific": "#0072B2", "ncBAF_specific": "#D55E00"}
    markers = {"PBAF_specific": "o", "ncBAF_specific": "s"}
    fig, axes = plt.subplots(1, 3, figsize=(190 / 25.4, 70 / 25.4), layout="constrained")

    ax = axes[0]
    for y, row in depmap.reset_index(drop=True).iterrows():
        c = colors[row.comparator]
        ax.errorbar(row.delta_rho_cBAF_minus_comparator, y,
                    xerr=[[row.delta_rho_cBAF_minus_comparator - row.bootstrap_ci_low],
                          [row.bootstrap_ci_high - row.delta_rho_cBAF_minus_comparator]],
                    fmt=markers[row.comparator], color=c, mec="white", mew=0.4, capsize=2)
        ax.text(row.bootstrap_ci_high + 0.015, y,
                f"P={row.holm_p_direct_contrast:.2g}", va="center", fontsize=6, color=c)
    ax.axvline(0, color="#666666", linestyle="--", linewidth=0.8)
    ax.set_yticks(range(len(depmap)), [x.replace("_specific", "") for x in depmap.comparator])
    ax.set_xlabel(r"$\Delta\rho$: cBAF $-$ comparator")
    ax.set_title("DepMap 26Q1\npaired cell lines", loc="left", weight="bold")

    ax = axes[1]
    xmap = {"PBAF_specific": 0, "ncBAF_specific": 1}
    for comparator, group in lineage_long.groupby("comparator", sort=False):
        x = xmap[comparator]
        jitter = np.linspace(-0.10, 0.10, len(group))
        ax.scatter(np.full(len(group), x) + jitter, group.delta_rho_cBAF_minus_comparator,
                   s=13, marker=markers[comparator], facecolor=colors[comparator],
                   edgecolor="white", linewidth=0.3, alpha=0.82)
        row = lineage_summary[lineage_summary.comparator.eq(comparator)].iloc[0]
        ax.errorbar(x, row.mean_delta_rho,
                    yerr=[[row.mean_delta_rho-row.bootstrap_mean_ci_low],
                          [row.bootstrap_mean_ci_high-row.mean_delta_rho]],
                    fmt="D", color="black", mfc="white", ms=5, capsize=3, zorder=5)
        ax.text(x, ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 0.1,
                f"Holm P={row.holm_p_signflip:.2g}", ha="center", va="bottom", fontsize=6)
    ax.axhline(0, color="#666666", linestyle="--", linewidth=0.8)
    ax.set_xticks([0, 1], ["PBAF", "ncBAF"])
    ax.set_ylabel(r"Within-lineage $\Delta\rho$")
    ax.set_title("DepMap 24Q4\nraw lineage estimates", loc="left", weight="bold")

    ax = axes[2]
    positions = {("all_ARID1A_mutation", "PBAF_specific"): 0,
                 ("all_ARID1A_mutation", "ncBAF_specific"): 1,
                 ("predicted_LOF_mutation", "PBAF_specific"): 3,
                 ("predicted_LOF_mutation", "ncBAF_specific"): 4}
    for row in drive.itertuples(index=False):
        x = positions[(row.state, row.comparator)]
        estimate = row.specificity_advantage_positive_supports_cBAF
        low = -row.bootstrap_ci_high_raw_difference
        high = -row.bootstrap_ci_low_raw_difference
        ax.errorbar(x, estimate, yerr=[[estimate-low], [high-estimate]],
                    fmt=markers[row.comparator], color=colors[row.comparator],
                    mec="white", mew=0.4, capsize=2)
    ax.axhline(0, color="#666666", linestyle="--", linewidth=0.8)
    ax.set_xticks([0.5, 3.5], ["Any mutation", "Predicted LoF"])
    ax.set_ylabel("cBAF relative-dependency advantage")
    ax.set_title("Project DRIVE\nstate-by-module interaction", loc="left", weight="bold")
    ax.text(0.02, 0.02, "circle: vs PBAF   square: vs ncBAF",
            transform=ax.transAxes, fontsize=5.8)

    for label, ax in zip("ABC", axes):
        ax.text(-0.18, 1.05, label, transform=ax.transAxes, fontsize=10, weight="bold")
        ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(out / "FigS13_cbaf_specificity_direct_tests.pdf")
    fig.savefig(out / "FigS13_cbaf_specificity_direct_tests.png", dpi=600)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depmap", required=True)
    parser.add_argument("--lineage", required=True)
    parser.add_argument("--drive", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    paths = {k: Path(v).resolve() for k, v in {
        "depmap": args.depmap, "lineage": args.lineage, "drive": args.drive
    }.items()}
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)

    depmap = pd.read_csv(paths["depmap"], index_col=0)
    lineage = pd.read_csv(paths["lineage"])
    drive = pd.read_csv(paths["drive"])
    depmap_result = depmap_direct_tests(depmap)
    lineage_long, lineage_result = lineage_direct_tests(lineage)
    drive_result = drive_direct_tests(drive)

    depmap_result.to_csv(out / "r10_depmap26q1_cbaf_direct_contrasts.csv", index=False)
    lineage_long.to_csv(out / "r10_lineage_cbaf_delta_by_lineage.csv", index=False)
    lineage_result.to_csv(out / "r10_lineage_cbaf_direct_contrasts.csv", index=False)
    drive_result.to_csv(out / "r10_project_drive_cbaf_interactions.csv", index=False)
    make_figure(depmap_result, lineage_long, lineage_result, drive_result, out)

    manifest = {
        "analysis": "direct cBAF specificity contrasts",
        "primary_logic": "test cBAF against PBAF and ncBAF directly; do not infer a difference from separate significance calls",
        "modules": MODULES,
        "source_sha256": {name: sha256(path) for name, path in paths.items()},
        "bootstrap_resamples": N_BOOT, "permutations": N_PERM, "seed": SEED,
        "inference_units": {
            "DepMap_26Q1": "cell line, paired across modules",
            "DepMap_24Q4": "lineage, using archived within-lineage correlations",
            "Project_DRIVE": "cell line with lineage-stratified state-label permutation",
        },
        "interpretation_limit": "specificity of association, not proof of long-term causal dependency",
        "software": {"python": platform.python_version(), "numpy": np.__version__,
                     "pandas": pd.__version__, "scipy": scipy.__version__},
    }
    (out / "r10_cbaf_specificity_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print("\nDepMap direct contrasts\n", depmap_result.to_string(index=False))
    print("\nLineage direct contrasts\n", lineage_result.to_string(index=False))
    print("\nProject DRIVE interactions\n", drive_result.to_string(index=False))


if __name__ == "__main__":
    main()
