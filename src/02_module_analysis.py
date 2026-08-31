"""Correlate ARID1A dependency with prespecified BAF module scores."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

from common import bh_adjust, ensure_parent, load_gene_effect, load_metadata, select_models


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-effect", required=True)
    parser.add_argument("--metadata")
    parser.add_argument("--tissue")
    parser.add_argument("--anchor", default="ARID1A")
    parser.add_argument("--modules", default=str(Path(__file__).parents[1] / "config" / "modules.json"))
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    matrix = load_gene_effect(args.gene_effect)
    metadata = load_metadata(args.metadata) if args.metadata else None
    matrix = select_models(matrix, metadata, args.tissue)
    if args.anchor not in matrix:
        raise KeyError(f"Anchor {args.anchor!r} is absent from the matrix")
    modules = json.loads(Path(args.modules).read_text(encoding="utf-8"))
    rows = []
    for module, requested in modules.items():
        genes = [g for g in requested if g != args.anchor and g in matrix.columns]
        if not genes:
            continue
        module_score = matrix[genes].mean(axis=1, skipna=True)
        valid = matrix[args.anchor].notna() & module_score.notna()
        rho, pvalue = spearmanr(matrix.loc[valid, args.anchor], module_score.loc[valid])
        rows.append(
            {
                "module": module,
                "anchor": args.anchor,
                "rho": rho,
                "p_value": pvalue,
                "n_models": int(valid.sum()),
                "n_genes_present": len(genes),
                "genes_present": ";".join(genes),
            }
        )
    result = pd.DataFrame(rows)
    result["q_value"] = bh_adjust(result["p_value"])
    result.sort_values("rho", ascending=False).to_csv(ensure_parent(args.output), index=False)
    print(f"Wrote {len(result)} module tests to {args.output}")


if __name__ == "__main__":
    main()

