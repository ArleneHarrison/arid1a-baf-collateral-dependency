"""Independent result-level audit; does not refit or select any predictive model."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for part in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(part)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    run = args.run
    manifest = json.loads((run / "run_manifest.json").read_text())
    for name, expected in manifest["outputs_sha256"].items():
        assert sha256(run / name) == expected, name
    spec = json.loads((run / "locked_spec.json").read_text())
    raw = pd.read_csv(spec["gene_effect_path"], index_col=0)
    eligible = pd.read_csv(run / "eligible_models_and_groups.csv")
    ids = eligible.ModelID.tolist()
    raw = raw.loc[ids]
    y = raw["ARID1A (8289)"]
    features = raw.drop(columns=["ARID1A (8289)"]).replace([np.inf, -np.inf], np.nan)
    predictions = pd.read_csv(run / "oof_predictions.csv")
    splits = pd.read_csv(run / "split_manifest.csv")
    selected = pd.read_csv(run / "train_only_feature_manifest.csv")
    metrics = pd.read_csv(run / "fold_metrics.csv")
    summary = pd.read_csv(run / "performance_summary.csv")
    singles = pd.read_csv(run / "single_gene_baseline_selection.csv")
    groups = eligible.set_index("ModelID")
    checks = []
    for scheme, scheme_splits in splits.groupby("scheme"):
        assert len(scheme_splits) == len(ids)
        assert scheme_splits.ModelID.is_unique
        assert set(scheme_splits.ModelID) == set(ids)
        for fold in range(1, 6):
            test_ids = scheme_splits.loc[scheme_splits.test_fold == fold, "ModelID"].tolist()
            train_ids = scheme_splits.loc[scheme_splits.test_fold != fold, "ModelID"].tolist()
            assert not set(groups.loc[train_ids, "patient_group"]) & set(groups.loc[test_ids, "patient_group"])
            if scheme == "patient_lineage_component_5fold":
                assert not set(groups.loc[train_ids, "OncotreeLineage"]) & set(groups.loc[test_ids, "OncotreeLineage"])
            train_features = features.loc[train_ids]
            variance = train_features.var(axis=0, ddof=1)
            ranking = pd.DataFrame({"feature": variance.index, "variance": variance.to_numpy()})
            ranking = ranking.loc[np.isfinite(ranking.variance) & (ranking.variance > 0)]
            ranking = ranking.sort_values(["variance", "feature"], ascending=[False, True]).head(3000)
            actual = selected.loc[(selected.scheme == scheme) & (selected.fold == fold)].sort_values("rank")
            assert ranking.feature.tolist() == actual.feature.tolist()
            np.testing.assert_allclose(ranking.variance, actual.training_variance, rtol=1e-11, atol=1e-13)
            medians = train_features[actual.feature].median()
            np.testing.assert_allclose(medians, actual.training_median, rtol=1e-11, atol=1e-13)
            single = singles.loc[(singles.scheme == scheme) & (singles.fold == fold)].iloc[0]
            imputed = train_features[actual.feature].fillna(medians)
            corrs = imputed.corrwith(y.loc[train_ids])
            assert corrs.abs().idxmax() == single.selected_feature
            np.testing.assert_allclose(corrs[single.selected_feature], single.training_pearson_r, atol=1e-12)
            x_single = imputed[single.selected_feature]
            slope = x_single.cov(y.loc[train_ids]) / x_single.var()
            intercept = y.loc[train_ids].mean() - slope * x_single.mean()
            np.testing.assert_allclose([slope, intercept], [single.slope, single.intercept], atol=1e-12)
            pred_fold = predictions.loc[(predictions.scheme == scheme) & (predictions.fold == fold)]
            for name, part in pred_fold.groupby("model"):
                assert set(part.ModelID) == set(test_ids)
                observed = y.loc[part.ModelID].to_numpy()
                np.testing.assert_allclose(observed, part.observed, atol=1e-14)
                np.testing.assert_allclose(part.training_target_mean, y.loc[train_ids].mean(), atol=1e-14)
                if name == "train_mean":
                    np.testing.assert_allclose(part.predicted, y.loc[train_ids].mean(), atol=1e-14)
                if name == "train_selected_single_gene_ols":
                    x_test = features.loc[part.ModelID, single.selected_feature].fillna(medians[single.selected_feature])
                    np.testing.assert_allclose(slope * x_test + intercept, part.predicted, atol=1e-12)
                residual = observed - part.predicted.to_numpy()
                calculated = np.array([1 - np.square(residual).sum() / np.square(observed - observed.mean()).sum(), np.sqrt(np.square(residual).mean()), np.abs(residual).mean()])
                recorded = metrics.loc[(metrics.scheme == scheme) & (metrics.fold == fold) & (metrics.model == name), ["r2", "rmse", "mae"]].iloc[0].to_numpy(dtype=float)
                np.testing.assert_allclose(calculated, recorded, atol=1e-12)
            checks.append({"scheme": scheme, "fold": fold, "split_group_disjointness": "PASS", "independent_pandas_variance_selection": "PASS", "independent_pandas_training_medians": "PASS", "single_gene_train_selection_and_ols": "PASS", "fold_metrics_manual_formula": "PASS", "observations_match_original_input": "PASS"})
        for name, part in predictions.loc[predictions.scheme == scheme].groupby("model"):
            part = part.set_index("ModelID").loc[ids]
            assert len(part) == len(ids)
            observed = y.to_numpy()
            residual = observed - part.predicted.to_numpy()
            calculated = np.array([1 - np.square(residual).sum() / np.square(observed - observed.mean()).sum(), np.sqrt(np.square(residual).mean()), np.abs(residual).mean()])
            record = summary.loc[(summary.scheme == scheme) & (summary.model == name)].iloc[0]
            np.testing.assert_allclose(calculated, record[["pooled_oof_r2", "pooled_oof_rmse", "pooled_oof_mae"]].to_numpy(dtype=float), atol=1e-12)
            folds = metrics.loc[(metrics.scheme == scheme) & (metrics.model == name)]
            for metric in ("r2", "rmse", "mae"):
                np.testing.assert_allclose(folds[metric].mean(), record["fold_mean_" + metric], atol=1e-12)
                np.testing.assert_allclose(folds[metric].std(ddof=1), record["fold_sd_" + metric], atol=1e-12)
            baseline = predictions.loc[(predictions.scheme == scheme) & (predictions.model == "train_mean")].set_index("ModelID").loc[ids]
            baseline_sse = np.square(observed - baseline.predicted.to_numpy()).sum()
            np.testing.assert_allclose(1 - np.square(residual).sum() / baseline_sse, record.oof_skill_vs_train_mean, atol=1e-12)
    result = {
        "status": "PASS", "new_model_fits": 0, "auditor_source_sha256": sha256(Path(__file__)),
        "run_manifest_sha256": sha256(run / "run_manifest.json"),
        "all_recorded_output_hashes_valid": True,
        "independent_pooled_metrics_and_fold_mean_sd": "PASS", "fold_checks": checks,
        "limitations": "Result-level independent formula checks do not constitute external biological validation; XGBoost/ridge were not refitted and held-out SHAP relies on in-run additivity assertions",
    }
    output = run / "independent_result_audit.json"
    assert not output.exists()
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
