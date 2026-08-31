"""Leakage-resistant reconstruction of the exploratory ARID1A SHAP model.

This is a fully specified reconstruction, not the unretained historical driver
that generated ``results/derived/ml_shap_arid1a_feature_importance.csv``.
Feature filtering is unsupervised and prediction quality is estimated only from
out-of-fold predictions.
"""

from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import KFold, cross_val_predict
from xgboost import XGBRegressor

from common import ensure_parent, load_gene_effect, load_metadata, select_models


SEED = 20260831


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-effect", required=True)
    parser.add_argument("--metadata")
    parser.add_argument("--tissue")
    parser.add_argument("--target", default="ARID1A")
    parser.add_argument("--max-features", type=int, default=3000)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--output", required=True, help="SHAP feature-importance CSV")
    parser.add_argument("--metrics", required=True, help="Out-of-fold metrics JSON")
    return parser.parse_args()


def make_model() -> XGBRegressor:
    return XGBRegressor(
        objective="reg:squarederror",
        n_estimators=500,
        max_depth=3,
        learning_rate=0.03,
        subsample=0.8,
        colsample_bytree=0.8,
        reg_lambda=1.0,
        random_state=SEED,
        n_jobs=-1,
    )


def main() -> None:
    args = parse_args()
    matrix = load_gene_effect(args.gene_effect)
    metadata = load_metadata(args.metadata) if args.metadata else None
    matrix = select_models(matrix, metadata, args.tissue)
    if args.target not in matrix.columns:
        raise KeyError(f"Target {args.target!r} is absent from the matrix")

    usable = matrix.loc[matrix[args.target].notna()].copy()
    y = usable.pop(args.target).to_numpy(dtype=float)
    variances = usable.var(axis=0, skipna=True).sort_values(ascending=False)
    features = variances.head(min(args.max_features, len(variances))).index
    X_frame = usable.loc[:, features]
    X = SimpleImputer(strategy="median").fit_transform(X_frame)

    folds = KFold(n_splits=args.folds, shuffle=True, random_state=SEED)
    oof = cross_val_predict(make_model(), X, y, cv=folds, n_jobs=1, method="predict")
    metrics = {
        "provenance": "fully specified reconstruction; historical training driver unavailable",
        "seed": SEED,
        "n_models": len(y),
        "n_features": X.shape[1],
        "n_folds": args.folds,
        "oof_r2": float(r2_score(y, oof)),
        "oof_rmse": float(mean_squared_error(y, oof) ** 0.5),
        "oof_mae": float(mean_absolute_error(y, oof)),
    }

    model = make_model().fit(X, y)
    try:
        import shap
    except ImportError as exc:
        raise SystemExit("Install the pinned requirements to compute SHAP values") from exc
    values = np.asarray(shap.TreeExplainer(model).shap_values(X))
    importance = pd.DataFrame(
        {"feature": features, "mean_abs_shap": np.mean(np.abs(values), axis=0)}
    ).sort_values("mean_abs_shap", ascending=False)
    importance.insert(0, "rank", np.arange(1, len(importance) + 1))
    importance.to_csv(ensure_parent(args.output), index=False)
    metrics_path = ensure_parent(args.metrics)
    metrics_path.write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()

