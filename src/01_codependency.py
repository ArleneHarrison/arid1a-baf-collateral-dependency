"""Genome-wide Spearman co-dependency screen anchored on ARID1A."""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from common import bh_adjust, ensure_parent, load_gene_effect, load_metadata, select_models


def run_screen(matrix: pd.DataFrame, anchor: str) -> pd.DataFrame:
    if anchor not in matrix.columns:
        raise KeyError(f"Anchor {anchor!r} is absent from the gene-effect matrix")
    rows: list[dict[str, float | int | str]] = []
    anchor_values = matrix[anchor]
    for gene in matrix.columns:
        valid = anchor_values.notna() & matrix[gene].notna()
        n = int(valid.sum())
        rho, pvalue = (np.nan, np.nan) if n < 4 else spearmanr(
            anchor_values.loc[valid], matrix.loc[valid, gene]
        )
        rows.append({"anchor": anchor, "gene": gene, "rho": rho, "p_value": pvalue, "n_models": n})
    result = pd.DataFrame(rows)
    result["q_value"] = bh_adjust(result["p_value"])
    result = result.sort_values(["rho", "gene"], ascending=[False, True]).reset_index(drop=True)
    result["rank_descending"] = np.arange(1, len(result) + 1)
    result["percentile_descending"] = 1 - (result["rank_descending"] - 1) / max(len(result) - 1, 1)
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-effect", required=True)
    parser.add_argument("--metadata")
    parser.add_argument("--tissue")
    parser.add_argument("--anchor", default="ARID1A")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    matrix = load_gene_effect(args.gene_effect)
    metadata = load_metadata(args.metadata) if args.metadata else None
    matrix = select_models(matrix, metadata, args.tissue)
    result = run_screen(matrix, args.anchor)
    result.to_csv(ensure_parent(args.output), index=False)
    print(f"Wrote {len(result):,} genes across {len(matrix):,} models to {args.output}")


if __name__ == "__main__":
    main()

