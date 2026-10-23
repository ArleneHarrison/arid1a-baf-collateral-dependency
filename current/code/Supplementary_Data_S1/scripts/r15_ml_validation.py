"""Frozen R15 ARID1A dependency prediction audit. Run scientific work on seeta."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import re
import sys
import time
import warnings
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler
import xgboost as xgb


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, data: object) -> None:
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def gene_symbol(label: str) -> str:
    return re.sub(r"\s+\(\d+\)$", "", str(label))


def train_preprocess(raw_train: np.ndarray, labels: np.ndarray, maximum: int):
    """No validation matrix or labels can enter this fit function."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        variances = np.nanvar(raw_train, axis=0, ddof=1)
    eligible = np.flatnonzero(np.isfinite(variances) & (variances > 0))
    ordered = eligible[np.lexsort((labels[eligible], -variances[eligible]))]
    selected = ordered[:maximum]
    if selected.size == 0:
        raise ValueError("No positive-variance training-fold features")
    imputer = SimpleImputer(strategy="median", keep_empty_features=True)
    transformed = imputer.fit_transform(raw_train[:, selected])
    expected = np.nanmedian(raw_train[:, selected], axis=0)
    np.testing.assert_allclose(imputer.statistics_, expected, rtol=0, atol=0)
    assert np.isfinite(transformed).all()
    return selected, variances[selected], imputer, transformed


def preprocessing_unit_test() -> dict:
    """Synthetic test: huge held-out values cannot alter selection/imputation."""
    train = np.array([[0.0, 1.0, np.nan], [2.0, 1.0, 0.0], [4.0, 1.0, 1.0], [6.0, 1.0, 2.0]])
    labels = np.array(["variable", "constant", "missing"])
    selected, variance, imputer, fitted = train_preprocess(train, labels, 2)
    before = (selected.copy(), imputer.statistics_.copy(), fitted.copy())
    validation_a = np.array([[8.0, 1000.0, np.nan], [10.0, -1000.0, 3.0]])
    validation_b = validation_a.copy()
    validation_b[:, 1] *= 1e10
    imputer.transform(validation_a[:, selected])
    imputer.transform(validation_b[:, selected])
    for observed, expected in zip((selected, imputer.statistics_, fitted), before):
        np.testing.assert_array_equal(observed, expected)
    assert labels[selected].tolist() == ["variable", "missing"]
    return {"train_only_selector": "PASS", "heldout_transform_does_not_mutate_fit": "PASS", "train_median": "PASS"}


def patient_lineage_components(patient: np.ndarray, lineage: np.ndarray):
    parent = list(range(len(patient)))

    def find(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    for grouping in (patient, lineage):
        seen = {}
        for index, group in enumerate(grouping):
            if group in seen:
                union(index, seen[group])
            else:
                seen[group] = index
    roots = np.array([find(i) for i in range(len(patient))])
    return np.array([f"component_{root:04d}" for root in roots])


def measures(observed, predicted):
    return {
        "r2": float(r2_score(observed, predicted)),
        "rmse": float(mean_squared_error(observed, predicted) ** 0.5),
        "mae": float(mean_absolute_error(observed, predicted)),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--expected-spec-sha256", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    start = time.time()
    spec_sha = sha256(args.spec)
    assert spec_sha == args.expected_spec_sha256, "Design changed after lock"
    spec = json.loads(args.spec.read_text(encoding="utf-8"))
    args.out.mkdir(parents=True, exist_ok=False)
    (args.out / "locked_spec.json").write_bytes(args.spec.read_bytes())
    source_path = Path(__file__).resolve()
    manifest = {
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "spec_sha256": spec_sha,
        "script_sha256": sha256(source_path),
        "platform": platform.platform(),
        "python": sys.version,
        "versions": {package: importlib.metadata.version(package) for package in ["numpy", "pandas", "scipy", "scikit-learn", "xgboost"]},
        "inputs": {},
        "compute_location": platform.node(),
    }
    for key in ("gene_effect_path", "metadata_path"):
        path = Path(spec[key])
        manifest["inputs"][key] = {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
    write_json(args.out / "pre_fit_manifest.json", manifest)
    print("FROZEN " + json.dumps(manifest), flush=True)
    tests = preprocessing_unit_test()
    data = pd.read_csv(spec["gene_effect_path"], index_col=0)
    data.index = data.index.astype(str)
    assert data.index.is_unique and data.columns.is_unique
    metadata = pd.read_csv(spec["metadata_path"], dtype=str, keep_default_na=True)
    assert metadata["ModelID"].is_unique
    metadata = metadata.set_index("ModelID")
    target_columns = [column for column in data.columns if gene_symbol(column) == spec["target_symbol"]]
    assert len(target_columns) == 1, target_columns
    target_column = target_columns[0]
    all_y = data[target_column].to_numpy(dtype=float)
    eligible = np.isfinite(all_y)
    initial_n = len(data)
    data = data.loc[eligible]
    y = all_y[eligible]
    ids = data.index.to_numpy(dtype=str)
    feature_frame = data.drop(columns=target_columns)
    labels = feature_frame.columns.to_numpy(dtype=str)
    assert all(gene_symbol(column) != spec["target_symbol"] for column in labels)
    X = feature_frame.to_numpy(dtype=np.float64)
    X[~np.isfinite(X)] = np.nan
    del data, feature_frame
    metadata_present = np.isin(ids, metadata.index)
    meta = metadata.reindex(ids)
    raw_patient = meta["PatientID"].fillna("").str.strip()
    raw_lineage = meta["OncotreeLineage"].fillna("").str.strip()
    patient = np.array([f"patient:{value}" if value else f"missing_patient_model:{model}" for value, model in zip(raw_patient, ids)])
    lineage = np.array([value if value else "UNKNOWN_LINEAGE" for value in raw_lineage])
    component = patient_lineage_components(patient, lineage)
    groups_frame = pd.DataFrame({"ModelID": ids, "PatientID": raw_patient.to_numpy(), "patient_group": patient, "OncotreeLineage": lineage, "lineage_component": component, "metadata_present": metadata_present})
    groups_frame.to_csv(args.out / "eligible_models_and_groups.csv", index=False)
    crosses = groups_frame.groupby("patient_group")["OncotreeLineage"].nunique()
    cross_patient = crosses[crosses > 1].index
    groups_frame.loc[groups_frame.patient_group.isin(cross_patient)].to_csv(args.out / "patients_crossing_lineages.csv", index=False)
    group_summary = {
        "input_rows": initial_n, "eligible_rows": len(y), "nonfinite_target_excluded": int(initial_n - len(y)),
        "target_column": target_column, "n_nontarget_features": len(labels),
        "missing_metadata_rows": int((~metadata_present).sum()), "missing_patient_rows": int((raw_patient == "").sum()),
        "missing_lineage_rows": int((raw_lineage == "").sum()), "unique_patient_groups": len(set(patient)),
        "patient_groups_with_multiple_models": int((pd.Series(patient).value_counts() > 1).sum()),
        "unique_lineages_including_unknown": len(set(lineage)), "patient_lineage_components": len(set(component)),
        "patient_groups_crossing_lineages": len(cross_patient),
        "grouping_limitation": "Missing PatientID uses independent ModelID, not proof of unrelated biological origins; unknown lineage is blocked as one group",
    }
    write_json(args.out / "group_summary.json", group_summary)
    print("GROUPS " + json.dumps(group_summary), flush=True)
    all_metrics, all_predictions, all_splits, selection_rows, audit_rows, summary_rows, single_rows = [], [], [], [], [], [], []
    shap_summaries = []
    schemes = dict(zip(spec["schemes"], (patient, component)))
    model_names = ["xgboost_fixed", "ridge_alpha10", "train_selected_single_gene_ols", "train_mean"]
    completed_xgb_fits = 0
    for scheme, groups in schemes.items():
        assert len(set(groups)) >= 5
        splitter = GroupKFold(n_splits=5, shuffle=True, random_state=spec["seed"])
        oof = {name: np.full(len(y), np.nan) for name in model_names}
        visit_count = np.zeros(len(y), dtype=int)
        shap_abs_sum = np.zeros(len(labels))
        shap_signed_sum = np.zeros(len(labels))
        selected_folds = np.zeros(len(labels), dtype=int)
        shap_success = True
        for fold, (train, test) in enumerate(splitter.split(X, y, groups), start=1):
            fold_start = time.time()
            assert not set(ids[train]).intersection(ids[test])
            assert not set(groups[train]).intersection(groups[test])
            assert not set(patient[train]).intersection(patient[test])
            overlap_lineage = sorted(set(lineage[train]).intersection(lineage[test]))
            if scheme == "patient_lineage_component_5fold":
                assert not overlap_lineage
            visit_count[test] += 1
            selected, variances, imputer, train_x = train_preprocess(X[train], labels, spec["max_features"])
            medians_before = imputer.statistics_.copy()
            test_x = imputer.transform(X[test][:, selected])
            np.testing.assert_array_equal(medians_before, imputer.statistics_)
            assert np.isfinite(test_x).all()
            for index, feature in enumerate(selected):
                selection_rows.append({"scheme": scheme, "fold": fold, "rank": index + 1, "feature": labels[feature], "gene": gene_symbol(labels[feature]), "training_variance": variances[index], "training_median": imputer.statistics_[index]})
            selected_folds[selected] += 1
            scaler = StandardScaler()
            ridge_train = scaler.fit_transform(train_x)
            ridge_test = scaler.transform(test_x)
            np.testing.assert_allclose(scaler.mean_, train_x.mean(axis=0), rtol=1e-12, atol=1e-12)
            ridge_model = Ridge(alpha=spec["ridge"]["alpha"], solver=spec["ridge"]["solver"]).fit(ridge_train, y[train])
            oof["ridge_alpha10"][test] = ridge_model.predict(ridge_test)
            centered_x = train_x - train_x.mean(axis=0)
            centered_y = y[train] - y[train].mean()
            denominator = np.sqrt((centered_x ** 2).sum(axis=0) * (centered_y ** 2).sum())
            correlations = np.divide(centered_x.T @ centered_y, denominator, out=np.zeros_like(denominator), where=denominator > 0)
            single_index = int(np.argmax(np.abs(correlations)))
            single_model = LinearRegression().fit(train_x[:, [single_index]], y[train])
            oof["train_selected_single_gene_ols"][test] = single_model.predict(test_x[:, [single_index]])
            single_rows.append({"scheme": scheme, "fold": fold, "selected_feature": labels[selected[single_index]], "training_pearson_r": float(correlations[single_index]), "slope": float(single_model.coef_[0]), "intercept": float(single_model.intercept_)})
            oof["train_mean"][test] = float(y[train].mean())
            model = xgb.XGBRegressor(**spec["xgboost"])
            model.fit(train_x, y[train])
            completed_xgb_fits += 1
            oof["xgboost_fixed"][test] = model.predict(test_x)
            shap_status, shap_error, shap_additivity_error = "PASS", "", None
            try:
                contributions = model.get_booster().predict(xgb.DMatrix(test_x), pred_contribs=True)
                np.testing.assert_allclose(contributions.sum(axis=1), oof["xgboost_fixed"][test], atol=2e-5, rtol=2e-5)
                shap_additivity_error = float(np.max(np.abs(contributions.sum(axis=1) - oof["xgboost_fixed"][test])))
                assert contributions.shape == (len(test), len(selected) + 1)
                shap_abs_sum[selected] += np.abs(contributions[:, :-1]).sum(axis=0)
                shap_signed_sum[selected] += contributions[:, :-1].sum(axis=0)
            except (xgb.core.XGBoostError, AssertionError, TypeError, ValueError) as exc:
                shap_success = False
                shap_status, shap_error = "UNAVAILABLE", str(exc)
            for name in model_names:
                all_metrics.append({"scheme": scheme, "fold": fold, "model": name, "n_train": len(train), "n_test": len(test), **measures(y[test], oof[name][test])})
                for index in test:
                    all_predictions.append({"scheme": scheme, "fold": fold, "model": name, "ModelID": ids[index], "patient_group": patient[index], "OncotreeLineage": lineage[index], "observed": float(y[index]), "predicted": float(oof[name][index]), "training_target_mean": float(y[train].mean())})
            for index in test:
                all_splits.append({"scheme": scheme, "test_fold": fold, "ModelID": ids[index], "patient_group": patient[index], "OncotreeLineage": lineage[index], "scheme_group": groups[index]})
            audit_rows.append({"scheme": scheme, "fold": fold, "n_train": len(train), "n_test": len(test), "n_selected_features": len(selected), "train_patient_groups": len(set(patient[train])), "test_patient_groups": len(set(patient[test])), "train_lineages": len(set(lineage[train])), "test_lineages": len(set(lineage[test])), "overlapping_patient_groups": 0, "overlapping_scheme_groups": 0, "overlapping_lineages": " | ".join(overlap_lineage), "target_absent_from_features": True, "selection_uses_train_only": True, "imputation_uses_train_only": True, "scaler_uses_train_only": True, "no_test_eval_set_or_early_stopping": True, "heldout_shap_status": shap_status, "heldout_shap_error": shap_error, "heldout_shap_max_additivity_error": shap_additivity_error, "elapsed_seconds": time.time() - fold_start})
            print("FOLD " + json.dumps({"scheme": scheme, "fold": fold, "n_train": len(train), "n_test": len(test), "selected_single": labels[selected[single_index]], "xgb_r2": all_metrics[-4]["r2"], "elapsed_seconds": time.time() - fold_start}), flush=True)
            pd.DataFrame(all_metrics).to_csv(args.out / "fold_metrics.csv", index=False)
            pd.DataFrame(audit_rows).to_csv(args.out / "leakage_audit.csv", index=False)
        assert np.all(visit_count == 1)
        baseline_sse = float(np.square(y - oof["train_mean"]).sum())
        for name in model_names:
            assert np.isfinite(oof[name]).all()
            subset = pd.DataFrame([row for row in all_metrics if row["scheme"] == scheme and row["model"] == name])
            summary_rows.append({"scheme": scheme, "model": name, "n_oof": len(y), **{f"pooled_oof_{key}": value for key, value in measures(y, oof[name]).items()}, **{f"fold_mean_{key}": float(subset[key].mean()) for key in ("r2", "rmse", "mae")}, **{f"fold_sd_{key}": float(subset[key].std(ddof=1)) for key in ("r2", "rmse", "mae")}, "oof_skill_vs_train_mean": float(1 - np.square(y - oof[name]).sum() / baseline_sse)})
        if shap_success:
            for index in np.flatnonzero(selected_folds):
                shap_summaries.append({"scheme": scheme, "feature": labels[index], "gene": gene_symbol(labels[index]), "selected_folds": int(selected_folds[index]), "mean_abs_oof_shap": float(shap_abs_sum[index] / len(y)), "mean_signed_oof_shap": float(shap_signed_sum[index] / len(y)), "n_oof_denominator": len(y)})
        tests[scheme + "_oof_once_per_model"] = "PASS"
    assert completed_xgb_fits == 10
    pd.DataFrame(all_predictions).to_csv(args.out / "oof_predictions.csv", index=False)
    pd.DataFrame(all_splits).to_csv(args.out / "split_manifest.csv", index=False)
    pd.DataFrame(selection_rows).to_csv(args.out / "train_only_feature_manifest.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(args.out / "performance_summary.csv", index=False)
    pd.DataFrame(single_rows).to_csv(args.out / "single_gene_baseline_selection.csv", index=False)
    if shap_summaries:
        pd.DataFrame(shap_summaries).sort_values(["scheme", "mean_abs_oof_shap"], ascending=[True, False]).to_csv(args.out / "heldout_xgboost_shap.csv", index=False)
    write_json(args.out / "assertion_tests.json", {"status": "PASS", "tests": tests, "completed_xgboost_fits": completed_xgb_fits})
    manifest["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    manifest["elapsed_seconds"] = time.time() - start
    manifest["outputs_sha256"] = {path.name: sha256(path) for path in sorted(args.out.iterdir()) if path.is_file()}
    write_json(args.out / "run_manifest.json", manifest)
    print("COMPLETE " + json.dumps(summary_rows), flush=True)


if __name__ == "__main__":
    main()
