#!/usr/bin/env python3
"""Audit long-term ARID1A-conditioned BAF dependency in public screens.

The analysis deliberately separates three questions:
1. Is the assay long enough to measure sustained fitness loss (>=14 days)?
2. Is the contrast isogenic for ARID1A loss/knockdown?
3. Is cBAF stronger than the competing PBAF and ncBAF assemblies?

Higher ``interaction_score`` always means stronger depletion/essentiality in the
ARID1A-deficient condition.  Module enrichment is evaluated against a
size-matched whole-genome random-gene null.  Direct cBAF contrasts use exact
label permutations, so a module is not called specific merely because one
within-module test is significant and another is not.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import platform
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy


SEED = 20260902
N_RANDOM = 100_000

MODULES = {
    "ARID1A-contacting cBAF": [
        "SMARCE1", "SMARCB1", "SMARCC1", "SMARCC2", "SMARCD1", "DPF2", "SMARCA4"
    ],
    "cBAF-specific": [
        "ARID1B", "DPF1", "DPF2", "DPF3", "BCL7A", "BCL7B", "BCL7C", "SS18", "SS18L1"
    ],
    "PBAF-specific": ["PBRM1", "ARID2", "BRD7", "PHF10"],
    "ncBAF-specific": ["BRD9", "BICRA", "BICRAL"],
}

# TKOv3-era symbols used GLTSCR1/GLTSCR1L for BICRA/BICRAL.
ALIASES = {"GLTSCR1": "BICRA", "GLTSCR1L": "BICRAL"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bh(values: pd.Series) -> pd.Series:
    p = values.to_numpy(float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked)
    out[order] = np.minimum(ranked, 1.0)
    return pd.Series(out, index=values.index)


def normalize_symbols(series: pd.Series) -> pd.Series:
    symbols = series.astype(str).str.strip().str.upper()
    return symbols.replace(ALIASES)


def load_u2os(path: Path) -> pd.DataFrame:
    data = pd.read_excel(path, sheet_name="ARID1A KO DMSO vs WT DMSO")
    data = data.dropna(subset=["GENE", "normZ"]).copy()
    data["gene"] = normalize_symbols(data["GENE"])
    # DrugZ: negative normZ means stronger depletion in ARID1A KO.
    data["interaction_score"] = -pd.to_numeric(data["normZ"], errors="coerce")
    data = data.dropna(subset=["interaction_score"])
    return data.groupby("gene", as_index=False)["interaction_score"].mean()


def _rmg1_bf_block(path: Path, sheet: str, ko_offset: int) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name=sheet, header=None)
    wt = raw.iloc[2:, [0, 1]].copy()
    wt.columns = ["gene", "wt_bf"]
    ko = raw.iloc[2:, [ko_offset, ko_offset + 1]].copy()
    ko.columns = ["gene", "ko_bf"]
    wt["gene"] = normalize_symbols(wt["gene"])
    ko["gene"] = normalize_symbols(ko["gene"])
    wt["wt_bf"] = pd.to_numeric(wt["wt_bf"], errors="coerce")
    ko["ko_bf"] = pd.to_numeric(ko["ko_bf"], errors="coerce")
    return wt.dropna().merge(ko.dropna(), on="gene", how="inner")


def load_rmg1(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    rep1 = _rmg1_bf_block(path, "1.5 Replicate 1 PR", 5)
    rep2 = _rmg1_bf_block(path, "1.6 Replicate 2 PR", 6)
    rep1["interaction_score"] = rep1["ko_bf"] - rep1["wt_bf"]
    rep2["interaction_score"] = rep2["ko_bf"] - rep2["wt_bf"]
    merged = rep1[["gene", "interaction_score"]].merge(
        rep2[["gene", "interaction_score"]], on="gene", suffixes=("_rep1", "_rep2")
    )
    merged["interaction_score"] = merged[
        ["interaction_score_rep1", "interaction_score_rep2"]
    ].mean(axis=1)
    return merged[["gene", "interaction_score"]], merged


def load_293a(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    data = pd.read_parquet(path).copy()
    data["gene"] = normalize_symbols(data["GENE"])
    wt_cols = ["WT_1", "WT_2", "WT_3", "WT_4"]
    data["interaction_same_batch"] = data["ARID1A_KD"] - data["WT_3"]
    data["interaction_all_wt_mean"] = data["ARID1A_KD"] - data[wt_cols].mean(axis=1)
    data["interaction_all_wt_max"] = data["ARID1A_KD"] - data[wt_cols].max(axis=1)
    primary = data[["gene", "interaction_same_batch"]].rename(
        columns={"interaction_same_batch": "interaction_score"}
    )
    sensitivity = data[[
        "gene", "interaction_same_batch", "interaction_all_wt_mean", "interaction_all_wt_max"
    ]]
    return primary, sensitivity


def empirical_module_tests(
    data: pd.DataFrame, dataset: str, rng: np.random.Generator
) -> pd.DataFrame:
    scores = data.drop_duplicates("gene").set_index("gene")["interaction_score"].dropna()
    universe = scores.to_numpy(float)
    rows = []
    for module, requested in MODULES.items():
        genes = [gene for gene in requested if gene in scores.index]
        observed = scores.loc[genes].to_numpy(float)
        random_means = np.empty(N_RANDOM, dtype=float)
        for start in range(0, N_RANDOM, 10_000):
            stop = min(start + 10_000, N_RANDOM)
            draws = np.vstack([
                rng.choice(universe, size=len(observed), replace=False)
                for _ in range(stop - start)
            ])
            random_means[start:stop] = draws.mean(axis=1)
        p = (1 + np.count_nonzero(random_means >= observed.mean())) / (N_RANDOM + 1)
        rows.append({
            "dataset": dataset,
            "module": module,
            "requested_genes": len(requested),
            "measured_genes": len(genes),
            "genes_used": ";".join(genes),
            "mean_interaction": observed.mean(),
            "median_interaction": np.median(observed),
            "fraction_positive": np.mean(observed > 0),
            "random_null_mean": random_means.mean(),
            "random_null_sd": random_means.std(ddof=1),
            "standardized_mean": (observed.mean() - random_means.mean()) / random_means.std(ddof=1),
            "empirical_p_one_sided": p,
            "universe_genes": len(universe),
            "random_draws": N_RANDOM,
        })
    out = pd.DataFrame(rows)
    out["bh_q_within_dataset"] = bh(out["empirical_p_one_sided"])
    return out


def exact_mean_difference(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    observed = float(a.mean() - b.mean())
    pooled = np.concatenate([a, b])
    n_a = len(a)
    differences = []
    for indices in itertools.combinations(range(len(pooled)), n_a):
        mask = np.zeros(len(pooled), dtype=bool)
        mask[list(indices)] = True
        differences.append(pooled[mask].mean() - pooled[~mask].mean())
    differences = np.asarray(differences)
    p = np.count_nonzero(differences >= observed - 1e-12) / len(differences)
    return observed, float(p)


def bootstrap_ci(a: np.ndarray, b: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    values = np.empty(50_000)
    for idx in range(len(values)):
        values[idx] = rng.choice(a, len(a), replace=True).mean() - rng.choice(
            b, len(b), replace=True
        ).mean()
    return tuple(np.quantile(values, [0.025, 0.975]))


def direct_tests(data: pd.DataFrame, dataset: str, rng: np.random.Generator) -> pd.DataFrame:
    scores = data.drop_duplicates("gene").set_index("gene")["interaction_score"].dropna()
    rank_scores = scores.rank(pct=True)
    rows = []
    for reference in ["cBAF-specific", "ARID1A-contacting cBAF"]:
        for competitor in ["PBAF-specific", "ncBAF-specific"]:
            genes_a = [g for g in MODULES[reference] if g in scores.index]
            genes_b = [g for g in MODULES[competitor] if g in scores.index]
            a = scores.loc[genes_a].to_numpy(float)
            b = scores.loc[genes_b].to_numpy(float)
            rank_a = rank_scores.loc[genes_a].to_numpy(float)
            rank_b = rank_scores.loc[genes_b].to_numpy(float)
            delta, p = exact_mean_difference(a, b)
            ci_low, ci_high = bootstrap_ci(a, b, rng)
            rank_ci_low, rank_ci_high = bootstrap_ci(rank_a, rank_b, rng)
            rows.append({
                "dataset": dataset,
                "reference_module": reference,
                "competitor_module": competitor,
                "n_reference": len(a),
                "n_competitor": len(b),
                "mean_reference": a.mean(),
                "mean_competitor": b.mean(),
                "delta_mean": delta,
                "bootstrap_ci_low": ci_low,
                "bootstrap_ci_high": ci_high,
                "rank_percentile_delta": rank_a.mean() - rank_b.mean(),
                "rank_bootstrap_ci_low": rank_ci_low,
                "rank_bootstrap_ci_high": rank_ci_high,
                "exact_permutation_p_one_sided": p,
            })
    out = pd.DataFrame(rows)
    out["holm_p_within_dataset"] = np.minimum(
        1.0,
        np.maximum.accumulate(
            np.sort(out["exact_permutation_p_one_sided"].to_numpy())
            * np.arange(len(out), 0, -1)
        )[np.argsort(np.argsort(out["exact_permutation_p_one_sided"].to_numpy()))],
    )
    # Recompute Holm in a transparent order-safe form.
    order = np.argsort(out["exact_permutation_p_one_sided"].to_numpy())
    adjusted_sorted = np.maximum.accumulate(
        out["exact_permutation_p_one_sided"].to_numpy()[order] * np.arange(len(out), 0, -1)
    )
    adjusted = np.empty(len(out))
    adjusted[order] = np.minimum(adjusted_sorted, 1.0)
    out["holm_p_within_dataset"] = adjusted
    return out


def baf_gene_table(datasets: dict[str, pd.DataFrame]) -> pd.DataFrame:
    all_baf = sorted({g for genes in MODULES.values() for g in genes})
    rows = []
    for dataset, data in datasets.items():
        scores = data.drop_duplicates("gene").set_index("gene")["interaction_score"]
        for gene in all_baf:
            if gene not in scores.index:
                continue
            rows.append({
                "dataset": dataset,
                "gene": gene,
                "interaction_score": scores[gene],
                "modules": ";".join(name for name, genes in MODULES.items() if gene in genes),
            })
    return pd.DataFrame(rows)


def complete_ko_consensus(u2os: pd.DataFrame, rmg1: pd.DataFrame) -> pd.DataFrame:
    """Combine independent complete-KO screens by within-screen percentile.

    The rank transform prevents the DrugZ and BAGEL scales from being pooled as
    if they were commensurable.  Each gene receives equal weight in each study.
    """
    left = u2os[["gene", "interaction_score"]].copy()
    right = rmg1[["gene", "interaction_score"]].copy()
    left["rank_u2os"] = left["interaction_score"].rank(pct=True)
    right["rank_rmg1"] = right["interaction_score"].rank(pct=True)
    merged = left[["gene", "rank_u2os"]].merge(
        right[["gene", "rank_rmg1"]], on="gene", how="inner"
    )
    merged["interaction_score"] = merged[["rank_u2os", "rank_rmg1"]].mean(axis=1)
    return merged


def rmg1_replicate_module_summary(replicates: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for module, requested in MODULES.items():
        data = replicates[replicates["gene"].isin(requested)]
        for replicate in ["rep1", "rep2"]:
            column = f"interaction_score_{replicate}"
            rows.append({
                "module": module,
                "replicate": replicate,
                "measured_genes": len(data),
                "mean_interaction": data[column].mean(),
                "median_interaction": data[column].median(),
                "fraction_positive": (data[column] > 0).mean(),
            })
    return pd.DataFrame(rows)


def make_figure(
    module_tests: pd.DataFrame,
    direct: pd.DataFrame,
    genes: pd.DataFrame,
    evidence: pd.DataFrame,
    output_pdf: Path,
) -> None:
    plt.rcParams.update({"font.size": 8, "axes.titlesize": 10, "figure.dpi": 180})
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.2), constrained_layout=True)

    colors = {
        "ARID1A-contacting cBAF": "#005F73",
        "cBAF-specific": "#0A9396",
        "PBAF-specific": "#CA6702",
        "ncBAF-specific": "#9B2226",
    }
    datasets = list(module_tests["dataset"].drop_duplicates())
    x = np.arange(len(datasets))
    offsets = np.linspace(-0.24, 0.24, 4)
    ax = axes[0, 0]
    for offset, module in zip(offsets, MODULES):
        sub = module_tests.set_index(["dataset", "module"]).loc[
            [(d, module) for d in datasets]
        ]
        ax.scatter(x + offset, sub["standardized_mean"], s=46, color=colors[module], label=module)
    ax.axhline(0, color="0.35", lw=0.8)
    ax.set_xticks(x, datasets, rotation=18, ha="right")
    ax.set_ylabel("Mean interaction vs random-gene null (SD)")
    ax.set_title("A  Long-term isogenic module enrichment")
    ax.legend(frameon=False, fontsize=7, ncol=2)

    ax = axes[0, 1]
    direct_cbaf = direct[direct["reference_module"] == "cBAF-specific"].copy()
    labels = []
    for idx, row in direct_cbaf.reset_index(drop=True).iterrows():
        y = idx
        ax.errorbar(
            row["rank_percentile_delta"], y,
            xerr=[[row["rank_percentile_delta"] - row["rank_bootstrap_ci_low"]],
                  [row["rank_bootstrap_ci_high"] - row["rank_percentile_delta"]]],
            fmt="o", color="#005F73" if "PBAF" in row["competitor_module"] else "#9B2226",
            capsize=2,
        )
        labels.append(f"{row['dataset']} vs {row['competitor_module'].replace('-specific','')}")
    ax.axvline(0, color="0.35", lw=0.8)
    ax.set_yticks(range(len(labels)), labels)
    ax.set_xlabel("cBAF-specific mean percentile advantage")
    ax.set_title("B  Direct cBAF competition (95% bootstrap CI)")

    ax = axes[1, 0]
    key_genes = ["ARID1B", "DPF2", "SMARCA4", "SMARCC1", "SMARCB1", "PBRM1", "ARID2", "BRD9", "BICRA"]
    pivot = genes[genes["gene"].isin(key_genes)].pivot(index="gene", columns="dataset", values="interaction_score")
    pivot = pivot.reindex(key_genes)
    ordered_columns = [d for d in datasets if d in pivot.columns]
    pivot = pivot[ordered_columns]
    # Column-wise ranks make scales comparable while preserving direction.
    rank = pivot.rank(axis=0, pct=True)
    image = ax.imshow(rank.to_numpy(), cmap="RdBu_r", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(rank.columns)), rank.columns, rotation=18, ha="right")
    ax.set_yticks(range(len(rank.index)), rank.index)
    ax.set_title("C  BAF-gene interaction percentile within each screen")
    fig.colorbar(image, ax=ax, fraction=0.046, label="Percentile (higher = stronger in ARID1A loss)")

    ax = axes[1, 1]
    cols = ["long_term", "isogenic", "complete_arid1a_loss", "lung_model", "direct_cbaf_positive"]
    matrix = evidence.set_index("dataset")[cols].astype(float)
    ax.imshow(matrix.to_numpy(), cmap="YlGn", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(cols)), ["≥14 d", "Isogenic", "Complete\nARID1A loss", "Lung", "Direct cBAF\npositive"], rotation=20, ha="right")
    ax.set_yticks(range(len(matrix.index)), matrix.index)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix.iat[i, j]
            ax.text(j, i, "Yes" if value == 1 else "No", ha="center", va="center",
                    fontsize=7, color="white" if value == 1 else "black")
    ax.set_title("D  Evidence-grade audit")

    fig.suptitle("Long-term ARID1A-conditioned BAF dependency: public-screen validation", fontsize=13)
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_pdf)
    fig.savefig(output_pdf.with_suffix(".png"), dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--u2os", required=True, type=Path)
    parser.add_argument("--rmg1", required=True, type=Path)
    parser.add_argument("--rpe1", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--figure", required=True, type=Path)
    args = parser.parse_args()

    rng = np.random.default_rng(SEED)
    u2os = load_u2os(args.u2os)
    rmg1, rmg1_replicates = load_rmg1(args.rmg1)
    rpe1, rpe1_sensitivity = load_293a(args.rpe1)
    consensus = complete_ko_consensus(u2os, rmg1)
    datasets = {
        "U2-OS ARID1A-KO, day 14": u2os,
        "RMG-1 ARID1A-KO, day 14": rmg1,
        "293A ARID1A-KD, day 21": rpe1,
        "Complete-KO consensus, day 14": consensus,
    }

    modules = pd.concat(
        [empirical_module_tests(data, name, rng) for name, data in datasets.items()],
        ignore_index=True,
    )
    direct = pd.concat(
        [direct_tests(data, name, rng) for name, data in datasets.items()],
        ignore_index=True,
    )
    genes = baf_gene_table(datasets)

    direct_positive = (
        direct.query("reference_module == 'cBAF-specific'")
        .groupby("dataset")
        .apply(lambda x: bool(((x.delta_mean > 0) & (x.holm_p_within_dataset < 0.05)).all()), include_groups=False)
    )
    evidence = pd.DataFrame([
        {"dataset": "DepMap 26Q1 / Project Score", "long_term": 1, "isogenic": 0,
         "complete_arid1a_loss": 0, "lung_model": 1, "direct_cbaf_positive": 1,
         "interpretation": "multi-week cancer fitness association; not ARID1A-loss causal"},
        {"dataset": "U2-OS ARID1A-KO, day 14", "long_term": 1, "isogenic": 1,
         "complete_arid1a_loss": 1, "lung_model": 0,
         "direct_cbaf_positive": int(direct_positive.get("U2-OS ARID1A-KO, day 14", False)),
         "interpretation": "complete-loss isogenic osteosarcoma screen"},
        {"dataset": "RMG-1 ARID1A-KO, day 14", "long_term": 1, "isogenic": 1,
         "complete_arid1a_loss": 1, "lung_model": 0,
         "direct_cbaf_positive": int(direct_positive.get("RMG-1 ARID1A-KO, day 14", False)),
         "interpretation": "complete-loss isogenic ovarian clear-cell cancer screen"},
        {"dataset": "293A ARID1A-KD, day 21", "long_term": 1, "isogenic": 1,
         "complete_arid1a_loss": 0, "lung_model": 0,
         "direct_cbaf_positive": int(direct_positive.get("293A ARID1A-KD, day 21", False)),
         "interpretation": "residual ARID1A isoform; authors explicitly classify as knockdown"},
    ])

    args.output_dir.mkdir(parents=True, exist_ok=True)
    modules.to_csv(args.output_dir / "r11_long_term_module_tests.csv", index=False)
    direct.to_csv(args.output_dir / "r11_long_term_direct_contrasts.csv", index=False)
    genes.to_csv(args.output_dir / "r11_long_term_baf_gene_scores.csv", index=False)
    evidence.to_csv(args.output_dir / "r11_long_term_evidence_matrix.csv", index=False)
    rmg1_replicates.to_csv(args.output_dir / "r11_rmg1_replicate_interactions.csv", index=False)
    rmg1_replicate_module_summary(rmg1_replicates).to_csv(
        args.output_dir / "r11_rmg1_replicate_module_summary.csv", index=False
    )
    rpe1_sensitivity.to_csv(args.output_dir / "r11_293a_wt_reference_sensitivity.csv", index=False)

    make_figure(modules, direct, genes, evidence, args.figure)
    manifest = {
        "seed": SEED,
        "random_gene_sets_per_module": N_RANDOM,
        "analysis": "long-term (>=14 day) ARID1A-conditioned BAF dependency audit",
        "direction": "higher interaction_score means stronger dependency in ARID1A-deficient cells",
        "primary_tests": {
            "module": "size-matched random-gene empirical null; BH within each dataset",
            "specificity": "exact label permutation for direct cBAF-vs-PBAF/ncBAF mean differences; Holm within dataset",
        },
        "source_sha256": {
            "u2os": sha256(args.u2os), "rmg1": sha256(args.rmg1), "rpe1": sha256(args.rpe1)
        },
        "software": {
            "python": platform.python_version(), "numpy": np.__version__,
            "pandas": pd.__version__, "scipy": scipy.__version__,
        },
        "interpretation_limit": "No public NSCLC isogenic ARID1A-loss genome-wide long-term screen was identified; cross-cancer isogenic evidence cannot establish NSCLC causality.",
    }
    (args.output_dir / "r11_long_term_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print("\nMODULE TESTS\n", modules.to_string(index=False))
    print("\nDIRECT TESTS\n", direct.to_string(index=False))
    print("\nEVIDENCE\n", evidence.to_string(index=False))


if __name__ == "__main__":
    main()
