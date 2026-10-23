"""Targeted DepMap Public 26Q1 replication through the public Breadbox API.

Only the prespecified BAF genes are retrieved.  DepMap 24Q4 remains the
historical analysis; 26Q1 is an additive temporal robustness analysis.  Release
comparisons use identical ModelIDs and identical lineage definitions.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import scipy
from scipy.stats import mannwhitneyu, spearmanr


BASE_URL = "https://depmap.org/portal/breadbox"
DOWNLOAD_META_URL = "https://depmap.org/portal/api/no-captcha/download/files"
USER_AGENT = "ARID1A-BAF-noncommercial-research-replication/1.0"
SEED = 20260902
N_BOOT = 20_000
MODULES = {
    "core_shared": ["SMARCB1", "SMARCC1", "SMARCC2", "SMARCD1", "SMARCD2", "SMARCE1", "ACTL6A", "ACTB"],
    "cBAF_specific": ["ARID1B", "DPF1", "DPF2", "DPF3", "BCL7A", "BCL7B", "BCL7C", "SS18", "SS18L1"],
    "PBAF_specific": ["PBRM1", "ARID2", "BRD7", "PHF10"],
    "ncBAF_specific": ["BRD9", "BICRA", "BICRAL"],
    "SMARCA2_alt_ATPase": ["SMARCA2"],
}


def get_json(session: requests.Session, path: str):
    response = session.get(f"{BASE_URL}{path}", timeout=180)
    response.raise_for_status()
    return response.json()


def find_dataset(datasets: list[dict], given_id: str, release_text: str | None = "26Q1") -> dict:
    matches = [d for d in datasets if d.get("given_id") == given_id]
    if len(matches) != 1:
        raise ValueError(f"Expected one dataset for {given_id}, found {len(matches)}")
    dataset = matches[0]
    if release_text and release_text not in dataset.get("name", ""):
        raise ValueError(f"Dataset {given_id} is not the prespecified {release_text} release: {dataset.get('name')}")
    return dataset


def fetch_matrix(session: requests.Session, dataset: dict, features: list[str]) -> pd.DataFrame:
    response = session.post(
        f"{BASE_URL}/datasets/matrix/{dataset['id']}",
        json={"features": features, "feature_identifier": "label", "samples": None,
              "sample_identifier": "id"},
        timeout=300,
    )
    response.raise_for_status()
    frame = pd.DataFrame(response.json())
    absent = sorted(set(features) - set(frame.columns))
    if absent:
        raise ValueError(f"Features absent from {dataset['name']}: {absent}")
    frame.index.name = "ModelID"
    return frame[features]


def clean_24q4(path: str) -> pd.DataFrame:
    frame = pd.read_csv(path, index_col=0, low_memory=False)
    frame.columns = [c.split(" (")[0] for c in frame.columns]
    return frame.loc[:, ~frame.columns.duplicated()]


def bootstrap_rho(x: np.ndarray, y: np.ndarray, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(x)
    values = []
    for _ in range(N_BOOT):
        idx = rng.integers(0, n, n)
        estimate = spearmanr(x[idx], y[idx]).statistic
        if np.isfinite(estimate):
            values.append(float(estimate))
    return tuple(np.quantile(values, [0.025, 0.975]).tolist())


def module_result(frame: pd.DataFrame, model_ids: list[str], module: str, release: str,
                  cohort: str, seed: int) -> dict:
    genes = MODULES[module]
    data = frame.loc[frame.index.intersection(model_ids), ["ARID1A", *genes]].dropna()
    score = data[genes].mean(axis=1)
    test = spearmanr(data.ARID1A, score)
    ci = bootstrap_rho(data.ARID1A.to_numpy(float), score.to_numpy(float), seed)
    return {
        "release": release, "cohort": cohort, "module": module, "n_models": len(data),
        "n_module_genes": len(genes), "genes": ";".join(genes),
        "spearman_rho": float(test.statistic), "spearman_p": float(test.pvalue),
        "bootstrap_ci_low": ci[0], "bootstrap_ci_high": ci[1],
    }


def paired_release_delta(frame24: pd.DataFrame, frame26: pd.DataFrame, model_ids: list[str],
                         module: str, cohort: str, seed: int) -> dict:
    genes = MODULES[module]
    ids = frame24.index.intersection(frame26.index).intersection(model_ids)
    a = frame24.loc[ids, ["ARID1A", *genes]]
    b = frame26.loc[ids, ["ARID1A", *genes]]
    valid = a.notna().all(axis=1) & b.notna().all(axis=1)
    a, b = a.loc[valid], b.loc[valid]
    a_score, b_score = a[genes].mean(axis=1), b[genes].mean(axis=1)
    rho24 = float(spearmanr(a.ARID1A, a_score).statistic)
    rho26 = float(spearmanr(b.ARID1A, b_score).statistic)
    rng = np.random.default_rng(seed)
    delta = np.empty(N_BOOT)
    for i in range(N_BOOT):
        idx = rng.integers(0, len(a), len(a))
        delta[i] = (spearmanr(b.ARID1A.to_numpy()[idx], b_score.to_numpy()[idx]).statistic
                    - spearmanr(a.ARID1A.to_numpy()[idx], a_score.to_numpy()[idx]).statistic)
    delta = delta[np.isfinite(delta)]
    ci = np.quantile(delta, [0.025, 0.975])
    return {
        "cohort": cohort, "module": module, "n_common_models": len(a),
        "rho_24Q4": rho24, "rho_26Q1": rho26, "delta_26Q1_minus_24Q4": rho26 - rho24,
        "delta_bootstrap_ci_low": float(ci[0]), "delta_bootstrap_ci_high": float(ci[1]),
    }


def bh_adjust(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked)
    out[order] = np.minimum(ranked, 1.0)
    return out


def lineage_zscore(frame: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    joined = frame.join(metadata[["OncotreeLineage"]], how="inner")
    numeric = frame.columns
    for gene in numeric:
        joined[gene] = joined.groupby("OncotreeLineage", observed=True)[gene].transform(
            lambda x: (x - x.mean()) / x.std(ddof=0) if x.notna().sum() >= 5 and x.std(ddof=0) > 0 else np.nan
        )
    return joined


def loss_proxy_results(dep26: pd.DataFrame, meta: pd.DataFrame, damaging: pd.Series,
                       expression: pd.Series, copy_number: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    state = pd.concat([
        damaging.rename("damaging_mutation"), expression.rename("expression_log2tpm1"),
        copy_number.rename("relative_copy_number")
    ], axis=1).join(meta[["OncotreeLineage"]], how="inner")
    state["damaging_mutation"] = state.damaging_mutation.fillna(False).astype(bool)
    state["strict_loss_proxy"] = ((state.relative_copy_number < 0.3) | (state.expression_log2tpm1 < 0.1))
    state["expanded_loss_proxy"] = state.strict_loss_proxy | state.damaging_mutation

    z = lineage_zscore(dep26, meta)
    rows = []
    rng = np.random.default_rng(SEED)
    for proxy in ["strict_loss_proxy", "expanded_loss_proxy"]:
        labels = state[proxy].reindex(z.index)
        for module, genes in MODULES.items():
            score = z[genes].mean(axis=1)
            valid = labels.notna() & score.notna()
            low = score[valid & labels]
            other = score[valid & ~labels]
            if len(low) < 3 or len(other) < 3:
                rows.append({"proxy": proxy, "module": module, "n_loss": len(low), "n_other": len(other),
                             "mean_difference_loss_minus_other": np.nan, "ci_low": np.nan, "ci_high": np.nan,
                             "mannwhitney_p": np.nan})
                continue
            diff = float(low.mean() - other.mean())
            boot = np.empty(N_BOOT)
            la, oa = low.to_numpy(float), other.to_numpy(float)
            for i in range(N_BOOT):
                boot[i] = (rng.choice(la, len(la), replace=True).mean()
                           - rng.choice(oa, len(oa), replace=True).mean())
            test = mannwhitneyu(low, other, alternative="two-sided")
            ci = np.quantile(boot, [0.025, 0.975])
            rows.append({"proxy": proxy, "module": module, "n_loss": len(low), "n_other": len(other),
                         "mean_difference_loss_minus_other": diff, "ci_low": float(ci[0]),
                         "ci_high": float(ci[1]), "mannwhitney_p": float(test.pvalue)})
    results = pd.DataFrame(rows)
    results["mannwhitney_q_within_proxy"] = results.groupby("proxy", observed=True)["mannwhitney_p"].transform(
        lambda x: bh_adjust(x) if x.notna().all() else np.nan
    )
    return state.reset_index(), results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depmap24", required=True)
    parser.add_argument("--model24", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT, "Accept": "application/json"})
    datasets = get_json(session, "/datasets/")
    dep_dataset = find_dataset(datasets, "Chronos_Combined")
    genes = sorted({"ARID1A", *sum(MODULES.values(), [])})
    dep26 = fetch_matrix(session, dep_dataset, genes)
    dep26.to_csv(out / "r6_depmap26q1_baf_subset.csv")

    dep24 = clean_24q4(args.depmap24)
    meta = pd.read_csv(args.model24).set_index("ModelID")
    common = dep24.index.intersection(dep26.index).intersection(meta.index).tolist()
    lung = meta.index[meta.OncotreeLineage.astype(str).str.casefold().eq("lung")].intersection(common).tolist()
    cohorts = {"pan_cancer_common_models": common, "lung_common_models": lung}

    module_rows, delta_rows = [], []
    counter = 0
    for cohort, ids in cohorts.items():
        for module in MODULES:
            counter += 1
            module_rows.append(module_result(dep24, ids, module, "24Q4", cohort, SEED + counter))
            module_rows.append(module_result(dep26, ids, module, "26Q1", cohort, SEED + 100 + counter))
            delta_rows.append(paired_release_delta(dep24, dep26, ids, module, cohort, SEED + 200 + counter))
    pd.DataFrame(module_rows).to_csv(out / "r6_depmap_release_module_correlations.csv", index=False)
    pd.DataFrame(delta_rows).to_csv(out / "r6_depmap_release_paired_deltas.csv", index=False)

    gene_rows = []
    for release, frame in [("24Q4", dep24), ("26Q1", dep26)]:
        for cohort, ids in cohorts.items():
            block = frame.loc[frame.index.intersection(ids), genes]
            for gene in genes:
                if gene == "ARID1A":
                    continue
                valid = block[["ARID1A", gene]].dropna()
                test = spearmanr(valid.ARID1A, valid[gene])
                gene_rows.append({"release": release, "cohort": cohort, "gene": gene,
                                  "n_models": len(valid), "spearman_rho": float(test.statistic),
                                  "spearman_p": float(test.pvalue)})
    pd.DataFrame(gene_rows).to_csv(out / "r6_depmap_release_gene_correlations.csv", index=False)

    mutation_dataset = find_dataset(datasets, "mutations_damaging")
    expression_dataset = find_dataset(datasets, "expression")
    cn_dataset = find_dataset(datasets, "OmicsCNGeneWGS")
    damaging = fetch_matrix(session, mutation_dataset, ["ARID1A"])["ARID1A"]
    expression = fetch_matrix(session, expression_dataset, ["ARID1A"])["ARID1A"]
    copy_number = fetch_matrix(session, cn_dataset, ["ARID1A"])["ARID1A"]
    states, loss_results = loss_proxy_results(dep26, meta, damaging, expression, copy_number)
    states.to_csv(out / "r6_arid1a_loss_proxy_states.csv", index=False)
    loss_results.to_csv(out / "r6_arid1a_loss_proxy_module_results.csv", index=False)

    metadata_response = session.get(DOWNLOAD_META_URL, timeout=180)
    metadata_response.raise_for_status()
    download_meta = pd.read_csv(StringIO(metadata_response.text))
    release_row = download_meta.loc[
        download_meta.release.eq("DepMap Public 26Q1") & download_meta.filename.eq("CRISPRGeneEffect.csv")
    ]
    release_record = release_row.iloc[0].to_dict() if len(release_row) == 1 else {}
    trace = {
        "analysis": "DepMap 26Q1 targeted temporal replication",
        "api_base": BASE_URL,
        "dataset": {k: dep_dataset.get(k) for k in ["id", "name", "given_id", "taiga_id", "units"]},
        "download_catalog_record": release_record,
        "features": genes,
        "n_26q1_models": len(dep26),
        "n_common_24q4_26q1_models": len(common),
        "n_common_lung_models": len(lung),
        "subset_sha256": hashlib.sha256((out / "r6_depmap26q1_baf_subset.csv").read_bytes()).hexdigest(),
        "comparison_policy": "24Q4 and 26Q1 evaluated on identical ModelIDs and 24Q4 lineage labels",
        "loss_proxy_policy": {
            "strict": "relative WGS copy number <0.3 OR ARID1A expression log2(TPM+1)<0.1",
            "expanded": "strict proxy OR DepMap damaging-mutation flag",
            "interpretation": "sensitivity proxy, not proof of biallelic loss or protein absence",
        },
        "bootstrap_resamples": N_BOOT,
        "random_seed": SEED,
        "python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__,
        "scipy": scipy.__version__, "requests": requests.__version__,
    }
    (out / "r6_depmap26q1_run_configuration.json").write_text(
        json.dumps(trace, indent=2, ensure_ascii=False, default=str), encoding="utf-8"
    )
    print(json.dumps(trace, indent=2, ensure_ascii=False, default=str))
    print(pd.DataFrame(delta_rows).to_string(index=False))
    print(loss_results.to_string(index=False))


if __name__ == "__main__":
    main()
