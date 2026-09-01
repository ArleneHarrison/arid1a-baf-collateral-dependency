# -*- coding: utf-8 -*-
"""Render Supplementary Figure S8 from server-computed robustness tables."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "derived"
OUT = ROOT / "manuscript" / "figures"

BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"
VERMILLION = "#D55E00"
GREY = "#6B7280"
LIGHT = "#E5E7EB"
DARK = "#111827"


def style_axis(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=0.7, length=3)


def panel_label(ax, label):
    ax.text(-0.12, 1.12, label, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", ha="left")


def main():
    mpl.rcParams.update({
        "font.family": "Arial",
        "font.size": 7,
        "axes.linewidth": 0.7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    })

    paired = pd.read_csv(TABLES / "r5_nested_validation" / "nested_paired_comparison.csv")
    nested_summary = json.loads(
        (TABLES / "r5_nested_validation" / "nested_summary.json").read_text(encoding="utf-8")
    )
    structure = pd.read_csv(
        TABLES / "r5_structure_validation" / "structural_group_sensitivity.csv"
    )
    global_q = pd.read_csv(
        TABLES / "r5_structure_validation" / "structure_global_validation.csv"
    )
    module = pd.read_csv(
        TABLES / "r5_arid1a_null_modules" / "arid1a_low_module_sensitivity.csv"
    )
    control = pd.read_csv(
        TABLES / "r5_arid1a_null_modules" / "smarca4_low_smarca2_positive_control.csv"
    )

    fig, axes = plt.subplots(2, 2, figsize=(190 / 25.4, 175 / 25.4),
                             layout="constrained")
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    # a: paired, symmetric nested-null comparison across curated complexes.
    diff = paired["concordance_difference"].to_numpy()
    bins = np.linspace(-1, 1, 31)
    ax_a.hist(diff, bins=bins, color=BLUE, edgecolor="white", linewidth=0.4)
    ax_a.axvline(0, color=DARK, linestyle="--", linewidth=0.8)
    median = nested_summary["paired_median_difference"]
    ci_low, ci_high = nested_summary["paired_median_difference_bootstrap_95ci"]
    ax_a.axvline(median, color=VERMILLION, linewidth=1.4)
    ymax = ax_a.get_ylim()[1]
    ax_a.plot([ci_low, ci_high], [ymax * 0.90, ymax * 0.90], color=VERMILLION,
              linewidth=2.2, solid_capstyle="butt")
    ax_a.text(0.02, 0.96,
              f"n = {len(diff)} paired complexes\nmedian difference = {median:.2f}\n"
              f"95% bootstrap CI {ci_low:.3f} to {ci_high:.3f}\n"
              r"one-sided Wilcoxon $P=1.03\times10^{-19}$",
              transform=ax_a.transAxes, va="top", ha="left", fontsize=6.4)
    ax_a.set(xlabel="Held-out concordance: curated complex − matched null",
             ylabel="Number of complexes",
             title="Symmetric nested validation")
    style_axis(ax_a)
    panel_label(ax_a, "a")

    # b: deposited map-model quality, plus explicit contact persistence audit.
    g = global_q.set_index("pdbid").loc[["6LTJ", "6LTH"]]
    x = np.arange(2)
    bars = ax_b.bar(x, g["global_Q_same_resolution_percentile"], width=0.58,
                    color=[GREY, GREEN], edgecolor=DARK, linewidth=0.7)
    ax_b.set_ylim(0, 100)
    ax_b.set_xticks(x, ["6LTJ\n3.7 Å", "6LTH\n3.0 Å"])
    ax_b.set_ylabel("Global Q-score percentile\nversus same-resolution entries")
    ax_b.set_title("Independent cryo-EM model-quality audit")
    for bar, row in zip(bars, g.itertuples()):
        y = bar.get_height()
        ax_b.text(bar.get_x() + bar.get_width() / 2, max(y + 3, 4),
                  f"{y:.1f}%\nQ={row.global_Q_score:.3f}",
                  ha="center", va="bottom", fontsize=6.5)
    ax_b.text(0.98, 0.60,
              "6/6 ARID1A contacts retained\n"
              "• in both structures\n"
              "• at 4.0, 4.5 and 5.0 Å\n"
              "• after bottom-decile residue-Q filtering",
              transform=ax_b.transAxes, ha="right", va="top", fontsize=6.3,
              bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=LIGHT, lw=0.7))
    style_axis(ax_b)
    panel_label(ax_b, "b")

    # c: structural-group effect-size sensitivity with bootstrap CIs.
    keep = structure[structure["analysis"].isin([
        "original_binary_partition_7v8",
        "experimental_6LTJ_6LTH_contacts_only_6v8",
        "subunits_physically_present_in_6LTJ_only_6v2",
    ])].copy()
    labels = {
        "original_binary_partition_7v8": "Original partition (7 vs 8)",
        "experimental_6LTJ_6LTH_contacts_only_6v8": "Experimental contacts only (6 vs 8)",
        "subunits_physically_present_in_6LTJ_only_6v2": "Present subunits only (6 vs 2)",
    }
    order_c = [
        "original_binary_partition_7v8",
        "experimental_6LTJ_6LTH_contacts_only_6v8",
        "subunits_physically_present_in_6LTJ_only_6v2",
    ]
    short_labels = ["Original partition\n7 vs 8", "Experimental contacts\n6 vs 8",
                    "Present subunits\n6 vs 2"]
    yc = np.arange(3)[::-1]
    for cohort, color, marker, offset, cohort_label in (
        ("lung", BLUE, "o", 0.11, "Lung"),
        ("pancancer", ORANGE, "s", -0.11, "Pan-cancer"),
    ):
        block = keep.set_index(["analysis", "cohort"])
        for y, analysis in zip(yc, order_c):
            row = block.loc[(analysis, cohort)]
            yy = y + offset
            ax_c.plot([row.bootstrap_95ci_low, row.bootstrap_95ci_high], [yy, yy],
                      color=color, linewidth=1.5)
            ax_c.scatter(row.median_difference, yy, s=24, marker=marker,
                         facecolor=color, edgecolor=DARK, linewidth=0.5, zorder=3,
                         label=cohort_label if y == yc[0] else None)
    ax_c.axvline(0, color=DARK, linestyle="--", linewidth=0.8)
    ax_c.set_yticks(yc, short_labels)
    ax_c.set_xlabel(r"Median $\rho$ difference, 95% bootstrap CI")
    ax_c.set_title("Structure–co-dependency sensitivity")
    ax_c.legend(frameon=False, loc="lower right", fontsize=6.2, handletextpad=0.4)
    style_axis(ax_c)
    panel_label(ax_c, "c")

    # d: primary 10% protein-low module tests plus a pipeline positive control.
    d = module[(module["anchor_definition"] == "ARID1A_protein") &
               (module["bottom_fraction"] == 0.10)].copy()
    order = [
        "ARID1A_contacting_cBAF", "core_shared", "cBAF_specific",
        "PBAF_specific", "ncBAF_specific", "ARID1B_paralogue", "SMARCA2_alt_ATPase",
    ]
    pretty = {
        "ARID1A_contacting_cBAF": "ARID1A-contacting cBAF",
        "core_shared": "Core shared", "cBAF_specific": "cBAF specific",
        "PBAF_specific": "PBAF specific", "ncBAF_specific": "ncBAF specific",
        "ARID1B_paralogue": "ARID1B paralogue", "SMARCA2_alt_ATPase": "SMARCA2 ATPase",
    }
    d = d.set_index("module").loc[order].reset_index()
    pc = control[control["bottom_fraction"] == 0.10].iloc[0]
    rows = []
    for row in d.itertuples():
        rows.append((pretty[row.module], row.mean_difference, row.bootstrap_95ci_low,
                     row.bootstrap_95ci_high, row.BH_q_within_threshold, BLUE, "o"))
    rows.append(("Positive control:\nSMARCA4-low → SMARCA2", pc.mean_difference,
                 pc.bootstrap_95ci_low, pc.bootstrap_95ci_high,
                 pc.BH_q_across_thresholds, VERMILLION, "D"))
    yd = np.arange(len(rows))[::-1]
    for y, (label, est, low, high, qvalue, color, marker) in zip(yd, rows):
        ax_d.plot([low, high], [y, y], color=color, linewidth=1.5)
        ax_d.scatter(est, y, s=24, marker=marker, facecolor=color,
                     edgecolor=DARK, linewidth=0.5, zorder=3)
        ax_d.text(0.58, y, f"q={qvalue:.3g}", ha="right", va="center", fontsize=5.8)
    ax_d.axvline(0, color=DARK, linestyle="--", linewidth=0.8)
    ax_d.set_yticks(yd, [r[0] for r in rows])
    ax_d.set_xlabel("Mean within-lineage dependency-z difference, 95% bootstrap CI\n"
                    "negative = stronger dependency in low-expression group")
    ax_d.set_title("ARID1A-low modules remain null after correction")
    ax_d.set_xlim(-2.3, 0.65)
    style_axis(ax_d)
    panel_label(ax_d, "d")

    pdf = OUT / "FigS8_robustness_validation.pdf"
    png = OUT / "FigS8_robustness_validation.png"
    fig.savefig(pdf, dpi=600, metadata={"Title": "Supplementary Figure S8 robustness validation"})
    fig.savefig(png, dpi=600)
    plt.close(fig)
    manifest = {
        "figure": "Supplementary Figure S8",
        "outputs": [str(pdf.relative_to(ROOT)), str(png.relative_to(ROOT))],
        "source_tables": [
            "results/derived/r5_nested_validation/nested_paired_comparison.csv",
            "results/derived/r5_nested_validation/nested_summary.json",
            "results/derived/r5_structure_validation/structure_global_validation.csv",
            "results/derived/r5_structure_validation/structural_group_sensitivity.csv",
            "results/derived/r5_arid1a_null_modules/arid1a_low_module_sensitivity.csv",
            "results/derived/r5_arid1a_null_modules/smarca4_low_smarca2_positive_control.csv",
        ],
        "uncertainty": "95% percentile bootstrap CI, 20,000 resamples",
        "random_seed": 20260901,
        "missing_data": "No values imputed for plotting; server analysis required complete module scores.",
    }
    (OUT / "FigS8_robustness_validation_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(pdf)
    print(png)


if __name__ == "__main__":
    main()
