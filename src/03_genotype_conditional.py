"""Compare partner gene effects between ARID1A-mutant and wild-type models."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from common import MODEL_COLUMNS, bh_adjust, ensure_parent, load_gene_effect, load_metadata, select_models


GENE_COLUMNS = ("HugoSymbol", "Gene", "gene", "SYMBOL")
CLASS_COLUMNS = ("VariantClassification", "Variant_Classification", "variant_classification")
SILENT_CLASSES = {"Silent", "Intron", "IGR", "3'UTR", "5'UTR", "RNA", "Flank"}
LOF_CLASSES = {"Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site", "Translation_Start_Site"}


def mutation_models(path: str, gene: str, mode: str) -> set[str]:
    mutations = pd.read_csv(path, low_memory=False)
    gene_column = next((c for c in GENE_COLUMNS if c in mutations.columns), None)
    model_column = next((c for c in MODEL_COLUMNS if c in mutations.columns), None)
    if gene_column is None or model_column is None:
        raise ValueError("Mutation table must contain gene-symbol and model-ID columns")
    selected = mutations[mutations[gene_column].astype(str).str.casefold() == gene.casefold()].copy()
    class_column = next((c for c in CLASS_COLUMNS if c in selected.columns), None)
    if class_column and mode == "lof":
        selected = selected[selected[class_column].isin(LOF_CLASSES)]
    elif class_column and mode == "nonsynonymous":
        selected = selected[~selected[class_column].isin(SILENT_CLASSES)]
    return set(selected[model_column].dropna().astype(str))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-effect", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--mutations", required=True)
    parser.add_argument("--tissue")
    parser.add_argument("--anchor", default="ARID1A")
    parser.add_argument("--mutation-mode", choices=("lof", "nonsynonymous", "any"), default="nonsynonymous")
    parser.add_argument("--modules", default=str(Path(__file__).parents[1] / "config" / "modules.json"))
    parser.add_argument("--min-group-size", type=int, default=6)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    matrix = load_gene_effect(args.gene_effect)
    metadata = load_metadata(args.metadata)
    matrix = select_models(matrix, metadata, args.tissue)
    mutant_ids = mutation_models(args.mutations, args.anchor, args.mutation_mode)
    is_mutant = pd.Series(matrix.index.isin(mutant_ids), index=matrix.index)
    modules = json.loads(Path(args.modules).read_text(encoding="utf-8"))
    genes = sorted({g for values in modules.values() for g in values if g in matrix.columns and g != args.anchor})
    rows = []
    for gene in genes:
        mutant = matrix.loc[is_mutant, gene].dropna()
        wild_type = matrix.loc[~is_mutant, gene].dropna()
        if min(len(mutant), len(wild_type)) < args.min_group_size:
            continue
        statistic, pvalue = mannwhitneyu(mutant, wild_type, alternative="two-sided")
        rows.append(
            {
                "gene": gene,
                "n_mutant": len(mutant),
                "n_wild_type": len(wild_type),
                "mean_mutant": mutant.mean(),
                "mean_wild_type": wild_type.mean(),
                "mean_difference_mutant_minus_wild_type": mutant.mean() - wild_type.mean(),
                "median_mutant": mutant.median(),
                "median_wild_type": wild_type.median(),
                "u_statistic": statistic,
                "p_value": pvalue,
                "common_language_effect": statistic / (len(mutant) * len(wild_type)),
            }
        )
    result = pd.DataFrame(rows)
    if result.empty:
        raise ValueError("No gene passed the requested group-size threshold")
    result["q_value"] = bh_adjust(result["p_value"])
    result.sort_values("p_value").to_csv(ensure_parent(args.output), index=False)
    print(f"Wrote {len(result)} genotype-conditional tests to {args.output}")


if __name__ == "__main__":
    main()

