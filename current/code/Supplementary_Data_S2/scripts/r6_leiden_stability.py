"""Parameter- and seed-stability audit of the ARID1A Leiden community."""

from __future__ import annotations

import argparse
import json
import platform
from itertools import product
from pathlib import Path

import igraph as ig
import leidenalg
import numpy as np
import pandas as pd
import scipy
from scipy.stats import hypergeom, rankdata
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors


BASE_SEED = 20260902
REFERENCE_SEED = 20260831
BAF_GENES = {
    "ARID1A", "ARID1B", "DPF1", "DPF2", "DPF3", "BCL7A", "BCL7B", "BCL7C", "SS18",
    "SMARCB1", "SMARCC1", "SMARCC2", "SMARCD1", "SMARCD2", "SMARCE1", "SMARCA4",
    "ACTL6A", "ACTB", "PBRM1", "ARID2", "BRD7", "PHF10", "BRD9", "BICRA", "SMARCA2",
}


def parse_values(text: str, cast):
    return [cast(x.strip()) for x in text.split(",") if x.strip()]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-effect", required=True)
    parser.add_argument(
        "--gene-universe-csv",
        help="Optional CSV with a gene column defining the archived analysis universe",
    )
    parser.add_argument("--out", required=True)
    parser.add_argument("--neighbors", default="10,15,20,30")
    parser.add_argument("--min-rho", default="0.10,0.15,0.20")
    parser.add_argument("--resolution", default="3.0,4.0,5.0")
    parser.add_argument("--seeds", default="20260831,20260902,20260903,20260904,20260905")
    parser.add_argument("--min-dependent-fraction", type=float, default=0.01)
    parser.add_argument("--max-dependent-fraction", type=float, default=0.80)
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    neighbors = parse_values(args.neighbors, int)
    cutoffs = parse_values(args.min_rho, float)
    resolutions = parse_values(args.resolution, float)
    seeds = parse_values(args.seeds, int)
    if not ({15} <= set(neighbors) and {0.15} <= set(cutoffs) and {4.0} <= set(resolutions)
            and REFERENCE_SEED in seeds):
        raise ValueError("The prespecified reference K=15/rho=0.15/resolution=4/seed=20260831 must be included")

    matrix = pd.read_csv(args.gene_effect, index_col=0, low_memory=False)
    matrix.columns = [c.split(" (")[0] for c in matrix.columns]
    matrix = matrix.loc[:, ~matrix.columns.duplicated()]
    if args.gene_universe_csv:
        archived = pd.read_csv(args.gene_universe_csv, usecols=["gene"])
        if archived.gene.duplicated().any():
            raise ValueError("Archived gene universe contains duplicate symbols")
        genes = [gene for gene in archived.gene.astype(str) if gene in matrix.columns]
    else:
        dependent_fraction = (matrix < -0.5).mean(axis=0)
        genes = dependent_fraction[
            dependent_fraction.between(args.min_dependent_fraction, args.max_dependent_fraction)
        ].index.tolist()
    if "ARID1A" not in genes and "ARID1A" in matrix.columns:
        genes.append("ARID1A")
    if "ARID1A" not in genes:
        raise ValueError("ARID1A absent from gene-effect matrix")

    values = matrix[genes].T.to_numpy(dtype=float)
    row_medians = np.nanmedian(values, axis=1)
    missing = np.where(np.isnan(values))
    values[missing] = row_medians[missing[0]]
    ranked = rankdata(values, axis=1)
    ranked -= ranked.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(ranked, axis=1, keepdims=True)
    ranked = np.divide(ranked, norms, out=np.zeros_like(ranked), where=norms != 0)

    max_k = max(neighbors)
    model = NearestNeighbors(n_neighbors=max_k + 1, metric="cosine", n_jobs=-1).fit(ranked)
    distances, indices = model.kneighbors(ranked)
    anchor_idx = genes.index("ARID1A")
    baf_in_universe = BAF_GENES & set(genes)

    partitions: dict[tuple, np.ndarray] = {}
    communities: dict[tuple, set[str]] = {}
    graph_meta: dict[tuple, tuple[int, int]] = {}
    rows = []
    for k, cutoff in product(neighbors, cutoffs):
        edge_weights: dict[tuple[int, int], float] = {}
        for source in range(len(genes)):
            for target, dist in zip(indices[source, 1:k + 1], distances[source, 1:k + 1]):
                rho = 1.0 - float(dist)
                if rho > cutoff:
                    edge = tuple(sorted((source, int(target))))
                    edge_weights[edge] = max(edge_weights.get(edge, -1.0), rho)
        graph = ig.Graph(n=len(genes), edges=list(edge_weights), directed=False)
        weights = list(edge_weights.values())
        graph_meta[(k, cutoff)] = (graph.vcount(), graph.ecount())
        for resolution, seed in product(resolutions, seeds):
            part = leidenalg.find_partition(
                graph,
                leidenalg.RBConfigurationVertexPartition,
                weights=weights,
                resolution_parameter=resolution,
                seed=seed,
            )
            membership = np.asarray(part.membership, dtype=np.int32)
            key = (k, cutoff, resolution, seed)
            partitions[key] = membership
            comm_id = membership[anchor_idx]
            community = {genes[i] for i in np.flatnonzero(membership == comm_id)}
            communities[key] = community

    ref_key = (15, 0.15, 4.0, REFERENCE_SEED)
    ref_partition = partitions[ref_key]
    ref_community = communities[ref_key]
    for key in sorted(partitions):
        k, cutoff, resolution, seed = key
        community = communities[key]
        hits = sorted(community & baf_in_universe)
        union = community | ref_community
        jaccard = len(community & ref_community) / len(union) if union else np.nan
        p_enrich = hypergeom.sf(len(hits) - 1, len(genes), len(baf_in_universe), len(community))
        rows.append({
            "neighbors": k,
            "min_rho": cutoff,
            "resolution": resolution,
            "seed": seed,
            "n_genes": len(genes),
            "n_edges": graph_meta[(k, cutoff)][1],
            "n_communities": int(partitions[key].max() + 1),
            "arid1a_community_size": len(community),
            "n_baf_genes": len(hits),
            "baf_genes": ";".join(hits),
            "baf_hypergeom_p": float(p_enrich),
            "community_jaccard_vs_reference": float(jaccard),
            "partition_ari_vs_reference": float(adjusted_rand_score(ref_partition, partitions[key])),
            "is_reference": int(key == ref_key),
        })

    results = pd.DataFrame(rows)
    results.to_csv(out / "r6_leiden_parameter_stability.csv", index=False)
    membership_rows = []
    for key, community in communities.items():
        k, cutoff, resolution, seed = key
        for gene in sorted(community):
            membership_rows.append({"neighbors": k, "min_rho": cutoff, "resolution": resolution,
                                    "seed": seed, "gene": gene, "is_baf": int(gene in BAF_GENES)})
    membership_df = pd.DataFrame(membership_rows)
    membership_df.to_csv(out / "r6_leiden_arid1a_community_membership.csv", index=False)

    recurrence = []
    for gene in sorted(BAF_GENES):
        flags = np.asarray([gene in c for c in communities.values()], dtype=float)
        recurrence.append({"gene": gene, "n_runs": len(flags), "n_present": int(flags.sum()),
                           "recurrence_fraction": float(flags.mean()),
                           "present_in_reference": int(gene in ref_community)})
    pd.DataFrame(recurrence).to_csv(out / "r6_leiden_baf_recurrence.csv", index=False)

    summary = {
        "analysis": "Leiden stability audit",
        "source_file": str(Path(args.gene_effect).resolve()),
        "gene_universe_source": (
            str(Path(args.gene_universe_csv).resolve()) if args.gene_universe_csv
            else "derived from dependent-fraction thresholds"
        ),
        "matrix_shape": list(matrix.shape),
        "selective_gene_count": len(genes),
        "parameter_runs": len(results),
        "neighbors": neighbors,
        "min_rho": cutoffs,
        "resolutions": resolutions,
        "seeds": seeds,
        "reference": {"key": ref_key, "community_size": len(ref_community),
                      "community_members": sorted(ref_community),
                      "baf_members": sorted(ref_community & BAF_GENES)},
        "median_community_jaccard_vs_reference": float(results.community_jaccard_vs_reference.median()),
        "minimum_community_jaccard_vs_reference": float(results.community_jaccard_vs_reference.min()),
        "median_partition_ari_vs_reference": float(results.partition_ari_vs_reference.median()),
        "random_seed_policy": "explicit seed passed to every Leiden run",
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "igraph": ig.__version__,
        "leidenalg": getattr(leidenalg, "__version__", "unknown"),
    }
    (out / "r6_leiden_run_configuration.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
