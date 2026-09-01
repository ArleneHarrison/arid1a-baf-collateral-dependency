# -*- coding: utf-8 -*-
"""Module-level sensitivity analysis for the ARID1A-low state.

Scientific computation is intended for sxycloud.  Gene effects are z-scored
within lineage before pooling.  Predefined complex modules are then averaged
per cell line and compared between the bottom 5/10/15% of ARID1A abundance and
the remaining lines.  Negative differences mean stronger dependency in the
ARID1A-low group.  The analysis reports effect sizes, bootstrap confidence
intervals and BH-adjusted P values; it does not convert association into a
synthetic-lethal claim.
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats


SEED = 20260901
MODULES = {
    "ARID1A_contacting_cBAF": [
        "SMARCE1", "SMARCB1", "SMARCC1", "SMARCC2", "SMARCD1", "DPF2", "SMARCA4"
    ],
    "core_shared": [
        "SMARCB1", "SMARCC1", "SMARCC2", "SMARCD1", "SMARCD2", "SMARCE1", "ACTL6A", "ACTB"
    ],
    "cBAF_specific": [
        "ARID1B", "DPF1", "DPF2", "DPF3", "BCL7A", "BCL7B", "BCL7C", "SS18", "SS18L1"
    ],
    "PBAF_specific": ["PBRM1", "ARID2", "BRD7", "PHF10"],
    "ncBAF_specific": ["BRD9", "BICRA", "BICRAL"],
    "ARID1B_paralogue": ["ARID1B"],
    "SMARCA2_alt_ATPase": ["SMARCA2"],
}


def bh_adjust(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked)
    out[order] = np.minimum(ranked, 1.0)
    return out


def bootstrap_mean_difference(x, y, rng, n_boot=20000):
    values = np.empty(n_boot)
    for i in range(n_boot):
        xb = x[rng.integers(0, len(x), len(x))]
        yb = y[rng.integers(0, len(y), len(y))]
        values[i] = xb.mean() - yb.mean()
    return np.quantile(values, [0.025, 0.975])


def rank_biserial_from_u(u, n_x, n_y):
    # Positive means x tends to be numerically larger than y.
    return float(2.0 * u / (n_x * n_y) - 1.0)


def load_inputs(bulk, proteomics, gdsc):
    raw = pd.read_csv(
        os.path.join(proteomics, "proteomics_extracted", "Protein_matrix_averaged_20250211.tsv"),
        sep="\t", header=None, low_memory=False,
    )
    symbols = raw.iloc[1, 2:].tolist()
    data = raw.iloc[3:].reset_index(drop=True)
    protein = data.iloc[:, 2:].apply(pd.to_numeric, errors="coerce")
    protein.columns = symbols
    protein.index = data.iloc[:, 1].values
    protein = protein.loc[:, ~protein.columns.duplicated()]
    model_list = pd.read_csv(os.path.join(proteomics, "model_list.csv"), low_memory=False)
    model_to_broad = model_list.dropna(subset=["BROAD_ID"]).set_index("model_id")["BROAD_ID"].to_dict()
    protein.index = [model_to_broad.get(i) for i in protein.index]
    protein = protein[[isinstance(i, str) and i.startswith("ACH-") for i in protein.index]]
    protein = protein[~protein.index.duplicated()]

    dependency = pd.read_csv(os.path.join(bulk, "depmap_CRISPRGeneEffect.csv"), index_col=0)
    dependency.columns = [c.split(" (")[0] for c in dependency.columns]
    dependency = dependency.loc[:, ~dependency.columns.duplicated()]

    models = pd.read_csv(os.path.join(bulk, "depmap_Model.csv"), low_memory=False)
    id_col = "ModelID" if "ModelID" in models.columns else models.columns[0]
    lineage = models.set_index(id_col)["OncotreeLineage"]

    expression_path = os.path.join(gdsc, "CCLE_expression.csv")
    header = pd.read_csv(expression_path, nrows=0).columns.tolist()
    arid1a_col = [c for c in header if c.split(" (")[0] == "ARID1A"][0]
    smarca4_col = [c for c in header if c.split(" (")[0] == "SMARCA4"][0]
    expression = pd.read_csv(
        expression_path, usecols=[header[0], arid1a_col, smarca4_col]
    ).set_index(header[0])
    expression = {
        "ARID1A": expression[arid1a_col],
        "SMARCA4": expression[smarca4_col],
    }
    return protein, expression, dependency, lineage


def within_lineage_zscore(frame, lineage):
    labels = pd.Series({i: lineage.get(i) for i in frame.index})
    out = frame.copy().astype(float)
    for _, indices in labels.dropna().groupby(labels.dropna()).groups.items():
        indices = [i for i in indices if i in out.index]
        if len(indices) < 8:
            out.loc[indices] = np.nan
            continue
        block = out.loc[indices]
        out.loc[indices] = (block - block.mean()) / block.std(ddof=0)
    return out


def analyze(anchor, anchor_label, dependency, lineage):
    genes = sorted({g for members in MODULES.values() for g in members if g in dependency.columns})
    lines = sorted(set(anchor.dropna().index) & set(dependency.index))
    z = within_lineage_zscore(dependency.loc[lines, genes], lineage)
    module_scores = pd.DataFrame(index=z.index)
    member_manifest = []
    for module, members in MODULES.items():
        used = [g for g in members if g in z.columns]
        module_scores[module] = z[used].mean(axis=1, skipna=False)
        member_manifest.append({
            "anchor_definition": anchor_label,
            "module": module,
            "n_members_planned": len(members),
            "n_members_used": len(used),
            "members_used": ";".join(used),
        })

    rows = []
    anchor = anchor.loc[lines]
    for fraction in (0.05, 0.10, 0.15):
        threshold = float(anchor.quantile(fraction))
        low = anchor[anchor <= threshold].index
        rest = anchor[anchor > threshold].index
        block = []
        for module in MODULES:
            x = module_scores.loc[module_scores.index.intersection(low), module].dropna().to_numpy()
            y = module_scores.loc[module_scores.index.intersection(rest), module].dropna().to_numpy()
            u = stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic")
            rng = np.random.default_rng(SEED + int(fraction * 100) + list(MODULES).index(module))
            ci_low, ci_high = bootstrap_mean_difference(x, y, rng)
            block.append({
                "anchor_definition": anchor_label,
                "bottom_fraction": fraction,
                "anchor_threshold": threshold,
                "module": module,
                "n_low": len(x),
                "n_rest": len(y),
                "mean_z_low": float(x.mean()),
                "mean_z_rest": float(y.mean()),
                "mean_difference": float(x.mean() - y.mean()),
                "bootstrap_95ci_low": float(ci_low),
                "bootstrap_95ci_high": float(ci_high),
                "rank_biserial": rank_biserial_from_u(u.statistic, len(x), len(y)),
                "mannwhitney_two_sided_P": float(u.pvalue),
            })
        q = bh_adjust([r["mannwhitney_two_sided_P"] for r in block])
        for row, q_value in zip(block, q):
            row["BH_q_within_threshold"] = float(q_value)
        rows.extend(block)
    return pd.DataFrame(rows), pd.DataFrame(member_manifest)


def analyze_single_pair(anchor, anchor_label, target_gene, dependency, lineage):
    """Run the same lineage-standardised design as an assay-dynamic-range control."""
    lines = sorted(set(anchor.dropna().index) & set(dependency.index))
    z = within_lineage_zscore(dependency.loc[lines, [target_gene]], lineage)[target_gene]
    anchor = anchor.loc[lines]
    rows = []
    for fraction in (0.05, 0.10, 0.15):
        threshold = float(anchor.quantile(fraction))
        low = anchor[anchor <= threshold].index
        rest = anchor[anchor > threshold].index
        x = z.loc[z.index.intersection(low)].dropna().to_numpy()
        y = z.loc[z.index.intersection(rest)].dropna().to_numpy()
        u = stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic")
        rng = np.random.default_rng(SEED + 1000 + int(fraction * 100))
        ci_low, ci_high = bootstrap_mean_difference(x, y, rng)
        rows.append({
            "anchor_definition": anchor_label,
            "target_gene": target_gene,
            "bottom_fraction": fraction,
            "anchor_threshold": threshold,
            "n_low": len(x),
            "n_rest": len(y),
            "mean_z_low": float(x.mean()),
            "mean_z_rest": float(y.mean()),
            "mean_difference": float(x.mean() - y.mean()),
            "bootstrap_95ci_low": float(ci_low),
            "bootstrap_95ci_high": float(ci_high),
            "rank_biserial": rank_biserial_from_u(u.statistic, len(x), len(y)),
            "mannwhitney_two_sided_P": float(u.pvalue),
        })
    result = pd.DataFrame(rows)
    result["BH_q_across_thresholds"] = bh_adjust(result["mannwhitney_two_sided_P"])
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bulk", required=True,
                        help="Directory containing DepMap dependency, expression, mutation and model files")
    parser.add_argument("--proteomics", required=True,
                        help="Directory containing the Sanger proteomics matrix")
    parser.add_argument("--gdsc", required=True,
                        help="Directory containing GDSC cell-line annotation files")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    protein, expression, dependency, lineage = load_inputs(args.bulk, args.proteomics, args.gdsc)
    all_results, manifests = [], []
    for anchor, label in (
        (protein["ARID1A"], "ARID1A_protein"),
        (expression["ARID1A"], "ARID1A_mRNA"),
    ):
        result, manifest = analyze(anchor, label, dependency, lineage)
        all_results.append(result)
        manifests.append(manifest)
    result = pd.concat(all_results, ignore_index=True)
    manifest = pd.concat(manifests, ignore_index=True)
    positive_control = analyze_single_pair(
        expression["SMARCA4"], "SMARCA4_mRNA", "SMARCA2", dependency, lineage
    )
    result.to_csv(os.path.join(args.out, "arid1a_low_module_sensitivity.csv"), index=False)
    manifest.to_csv(os.path.join(args.out, "arid1a_low_module_manifest.csv"), index=False)
    positive_control.to_csv(os.path.join(args.out, "smarca4_low_smarca2_positive_control.csv"), index=False)
    print(result.to_string(index=False))
    print("\nPositive control:\n", positive_control.to_string(index=False))


if __name__ == "__main__":
    main()
