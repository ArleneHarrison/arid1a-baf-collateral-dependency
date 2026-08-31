"""Shared, deliberately small utilities for the portable workflows."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd


MODEL_COLUMNS = ("ModelID", "DepMap_ID", "DepMapID", "model_id")
LINEAGE_COLUMNS = (
    "OncotreeLineage",
    "OncotreePrimaryDisease",
    "PrimaryDisease",
    "Tissue",
    "lineage",
)


def clean_gene_name(value: object) -> str:
    """Convert DepMap labels such as ``ARID1A (8289)`` to HGNC symbols."""
    return re.sub(r"\s*\([^()]+\)\s*$", "", str(value)).strip()


def load_gene_effect(path: str | Path) -> pd.DataFrame:
    """Load a model-by-gene effect matrix and normalise its gene labels."""
    frame = pd.read_csv(path, index_col=0, low_memory=False)
    frame.index = frame.index.astype(str)
    frame.columns = [clean_gene_name(column) for column in frame.columns]
    frame = frame.apply(pd.to_numeric, errors="coerce")
    if frame.columns.duplicated().any():
        frame = frame.T.groupby(level=0, sort=False).mean().T
    return frame


def load_metadata(path: str | Path) -> pd.DataFrame:
    frame = pd.read_csv(path, low_memory=False)
    model_column = next((c for c in MODEL_COLUMNS if c in frame.columns), None)
    if model_column is None:
        first = frame.columns[0]
        if frame[first].astype(str).str.match(r"^(ACH|SIDM)-").mean() > 0.5:
            model_column = first
        else:
            raise ValueError("Could not identify a model-ID column in metadata")
    return frame.set_index(model_column, drop=False)


def select_models(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame | None = None,
    tissue: str | None = None,
) -> pd.DataFrame:
    if tissue is None:
        return matrix
    if metadata is None:
        raise ValueError("--metadata is required when --tissue is supplied")
    lineage_column = next((c for c in LINEAGE_COLUMNS if c in metadata.columns), None)
    if lineage_column is None:
        raise ValueError("Could not identify a tissue/lineage column in metadata")
    values = metadata[lineage_column].fillna("").astype(str)
    keep = metadata.index[values.str.casefold() == tissue.casefold()].astype(str)
    selected = matrix.loc[matrix.index.intersection(keep)]
    if selected.empty:
        available = sorted(values.dropna().unique())[:20]
        raise ValueError(f"No models matched tissue={tissue!r}; examples: {available}")
    return selected


def bh_adjust(pvalues: np.ndarray | pd.Series) -> np.ndarray:
    """Benjamini–Hochberg adjusted P values, preserving missing entries."""
    values = np.asarray(pvalues, dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    if not valid.any():
        return adjusted
    p = values[valid]
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    restored = np.empty_like(q)
    restored[order] = np.minimum(q, 1.0)
    adjusted[valid] = restored
    return adjusted


def ensure_parent(path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output
