"""Build a gene co-dependency kNN graph and detect Leiden communities."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.neighbors import NearestNeighbors

from common import ensure_parent, load_gene_effect


SEED = 20260831


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-effect", required=True)
    parser.add_argument("--gene-list", help="Optional one-gene-per-line prespecified selective-gene list")
    parser.add_argument("--anchor", default="ARID1A")
    parser.add_argument("--neighbors", type=int, default=15)
    parser.add_argument("--min-rho", type=float, default=0.15)
    parser.add_argument("--resolution", type=float, default=4.0)
    parser.add_argument("--min-dependent-fraction", type=float, default=0.01)
    parser.add_argument("--max-dependent-fraction", type=float, default=0.80)
    parser.add_argument("--output", required=True)
    parser.add_argument("--anchor-output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    matrix = load_gene_effect(args.gene_effect)
    if args.gene_list:
        requested = [line.strip() for line in Path(args.gene_list).read_text().splitlines() if line.strip()]
        genes = [gene for gene in requested if gene in matrix.columns]
    else:
        dependent_fraction = (matrix < -0.5).mean(axis=0)
        genes = dependent_fraction[
            dependent_fraction.between(args.min_dependent_fraction, args.max_dependent_fraction)
        ].index.tolist()
    if args.anchor not in genes and args.anchor in matrix.columns:
        genes.append(args.anchor)
    if len(genes) <= args.neighbors:
        raise ValueError("Too few selected genes for the requested number of neighbours")

    values = matrix[genes].T.to_numpy(dtype=float)
    row_medians = np.nanmedian(values, axis=1)
    missing = np.where(np.isnan(values))
    values[missing] = row_medians[missing[0]]
    ranked = rankdata(values, axis=1)
    ranked -= ranked.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(ranked, axis=1, keepdims=True)
    ranked = np.divide(ranked, norms, out=np.zeros_like(ranked), where=norms != 0)

    knn = NearestNeighbors(n_neighbors=args.neighbors + 1, metric="cosine", n_jobs=-1).fit(ranked)
    distances, indices = knn.kneighbors(ranked)
    edge_weights: dict[tuple[int, int], float] = {}
    for source, (near, distance) in enumerate(zip(indices[:, 1:], distances[:, 1:])):
        for target, dist in zip(near, distance):
            rho = 1 - float(dist)
            if rho > args.min_rho:
                edge = tuple(sorted((source, int(target))))
                edge_weights[edge] = max(edge_weights.get(edge, -1), rho)

    try:
        import igraph as ig
        import leidenalg
    except ImportError as exc:
        raise SystemExit("Install python-igraph and leidenalg from requirements.txt") from exc
    graph = ig.Graph(n=len(genes), edges=list(edge_weights), directed=False)
    weights = [edge_weights[edge] for edge in edge_weights]
    partition = leidenalg.find_partition(
        graph,
        leidenalg.RBConfigurationVertexPartition,
        weights=weights,
        resolution_parameter=args.resolution,
        seed=SEED,
    )
    membership = pd.DataFrame({"gene": genes, "community": partition.membership})
    membership.to_csv(ensure_parent(args.output), index=False)
    anchor_community = int(membership.loc[membership.gene == args.anchor, "community"].iloc[0])
    membership.loc[membership.community == anchor_community].to_csv(
        ensure_parent(args.anchor_output), index=False
    )
    print(
        f"Built {len(edge_weights):,} edges among {len(genes):,} genes; "
        f"{args.anchor} community {anchor_community} contains "
        f"{(membership.community == anchor_community).sum()} genes"
    )


if __name__ == "__main__":
    main()

