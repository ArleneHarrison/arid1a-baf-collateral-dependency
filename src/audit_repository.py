"""Validate numerical claims and cryptographic integrity of this release."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).parents[1]
MANIFEST = ROOT / "metadata" / "output_manifest.tsv"
REQUIRED = (
    "README.md",
    "CITATION.cff",
    "metadata/data_manifest.tsv",
    "results/derived/r3_structural_claim_audit.csv",
    "results/derived/r3_scrna_sample_level_summary.csv",
    "results/derived/ml_xgboost_recomputed_metrics.json",
    "results/derived/ml_shap_arid1a_feature_importance_RECOMPUTED.csv",
    "manuscript/manuscript.pdf",
    "manuscript/supplementary_material.pdf",
)


TEXT_SUFFIXES = {".cff", ".csv", ".json", ".md", ".py", ".r", ".tex", ".tsv", ".txt", ".yml", ".yaml"}
BUILD_SUFFIXES = {".aux", ".bbl", ".blg", ".fdb_latexmk", ".fls", ".log", ".out"}


def canonical_bytes(path: Path) -> bytes:
    content = path.read_bytes()
    if path.suffix.casefold() in TEXT_SUFFIXES:
        content = content.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return content


def digest(path: Path) -> str:
    return hashlib.sha256(canonical_bytes(path)).hexdigest()


def release_files() -> list[Path]:
    patterns = (
        "README.md", "CITATION.cff", "LICENSE", "requirements.txt", "environment.yml",
        "config/**/*", "metadata/data_manifest.tsv", "metadata/claim_manifest.tsv",
        "results/derived/**/*", "src/**/*", "tests/**/*", "manuscript/**/*",
    )
    files: set[Path] = set()
    for pattern in patterns:
        files.update(
            path for path in ROOT.glob(pattern)
            if path.is_file()
            and "__pycache__" not in path.parts
            and path.suffix.casefold() not in BUILD_SUFFIXES
        )
    return sorted(files)


def write_manifest() -> None:
    rows = [
        {"path": path.relative_to(ROOT).as_posix(), "bytes": len(canonical_bytes(path)), "sha256": digest(path)}
        for path in release_files()
    ]
    pd.DataFrame(rows).to_csv(MANIFEST, sep="\t", index=False)


def verify_manifest() -> None:
    if not MANIFEST.exists():
        raise AssertionError("metadata/output_manifest.tsv is missing; run with --write-manifest")
    recorded = pd.read_csv(MANIFEST, sep="\t")
    for row in recorded.itertuples(index=False):
        path = ROOT / row.path
        if not path.is_file():
            raise AssertionError(f"Manifest file is missing: {row.path}")
        if len(canonical_bytes(path)) != row.bytes or digest(path) != row.sha256:
            raise AssertionError(f"Manifest mismatch: {row.path}")


def verify_claims() -> None:
    structural = pd.read_csv(ROOT / "results" / "derived" / "r3_structural_claim_audit.csv")
    lung = structural.query("analysis == 'binary_contact_vs_noncontact' and cohort == 'lung'").iloc[0]
    pan = structural.query("analysis == 'binary_contact_vs_noncontact' and cohort == 'pan-cancer'").iloc[0]
    graded_lung = structural.query("analysis == 'graded_contact_count_within_6LTJ_contacts' and cohort == 'lung'").iloc[0]
    assert int(lung.n_contacting) == 7 and int(lung.n_noncontacting) == 8
    assert np.isclose(lung.estimate, 0.4355973036073978, atol=1e-12)
    assert np.isclose(lung.p_value, 0.0037296037296037296, atol=1e-15)
    assert np.isclose(pan.estimate, 0.31434597940183, atol=1e-12)
    assert np.isclose(graded_lung.estimate, -0.6571428571428573, atol=1e-12)

    scrna = pd.read_csv(ROOT / "results" / "derived" / "r3_scrna_sample_level_summary.csv").iloc[0]
    assert int(scrna.n_samples) == 48 and int(scrna.n_cells) == 35600
    assert np.isclose(scrna.spearman_rho_ARID1A_vs_coreBAF, 0.41576204950065127, atol=1e-12)
    assert np.isclose(scrna.spearman_p, 0.0032940975407258605, atol=1e-15)
    assert np.isclose(scrna.mannwhitney_p, 0.02569220626226525, atol=1e-15)

    for cohort in ("lung", "pancancer"):
        table = pd.read_csv(ROOT / "results" / "derived" / f"genomewide_codependency_{cohort}.csv")
        assert table.gene.is_unique
        assert table.iloc[0].gene in {"SMARCE1", "SMARCB1"}

    metrics = pd.read_json(
        ROOT / "results" / "derived" / "ml_xgboost_recomputed_metrics.json",
        typ="series",
    )
    shap = pd.read_csv(ROOT / "results" / "derived" / "ml_shap_arid1a_feature_importance_RECOMPUTED.csv")
    assert int(metrics.n_models) == 1178 and int(metrics.n_features) == 3000
    assert np.isclose(metrics.oof_r2, 0.3903418708231198, atol=1e-12)
    assert shap.head(5).feature.tolist() == ["SMARCB1", "SMARCE1", "SMARCC1", "ARID1B", "DPF2"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-manifest", action="store_true")
    args = parser.parse_args()
    missing = [item for item in REQUIRED if not (ROOT / item).is_file()]
    if missing:
        raise AssertionError(f"Required release files are missing: {missing}")
    verify_claims()
    if args.write_manifest:
        write_manifest()
    verify_manifest()
    print("PASS: required files, key numerical claims, and SHA-256 manifest verified")


if __name__ == "__main__":
    main()
