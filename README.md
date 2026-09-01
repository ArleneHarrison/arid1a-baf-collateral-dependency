# ARID1A–BAF collateral dependency — reproducibility repository

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22218065.svg)](https://doi.org/10.5281/zenodo.22218065)

This repository accompanies the manuscript **“Structure-informed co-dependency
analysis localises the collateral dependency of ARID1A-deficient non-small-cell
lung cancer to the subunits ARID1A physically contacts.”** The submission-ready
manuscript, supplement, analysis outputs and auditable source code are versioned
together here.

The central estimand is **CRISPR co-dependency**: correlation between ARID1A
gene-effect scores and those of BAF-complex partners across cancer cell lines.
It is not, by itself, evidence of synthetic lethality or a treatment effect.

## Repository map

- `src/` — portable analysis entry points for co-dependency, module,
  genotype-conditional, XGBoost/SHAP, Leiden, TCGA, structural and single-cell
  analyses.
- `config/` — preregistered complex membership used by the portable workflows.
- `metadata/` — source-data and output manifests.
- `results/derived/` — compact manuscript-level derived tables.
- `results/recomputed/` — regenerated audit outputs (created by the scripts).
- `manuscript/` — exact LaTeX sources and PDFs associated with the release.
- `archive/legacy_scripts/` — historical scripts retained for provenance; these
  are not the recommended workflow entry points.
- `tests/` — numerical and repository-integrity checks.

## Fast verification (no raw downloads)

```bash
python -m pip install -r requirements.txt
python src/audit_repository.py
pytest -q
```

This checks the key structural partition and sample-level single-cell results
against the archived derived tables and validates the release manifest.

## Recompute from public source data

Download the files listed in `metadata/data_manifest.tsv` to `data/raw/` (raw
files are intentionally not redistributed), then run, for example:

```bash
python src/01_codependency.py --gene-effect data/raw/CRISPRGeneEffect.csv \
  --metadata data/raw/Model.csv --tissue Lung \
  --output results/recomputed/codependency_lung.csv

python src/02_module_analysis.py --gene-effect data/raw/CRISPRGeneEffect.csv \
  --metadata data/raw/Model.csv --tissue Lung \
  --output results/recomputed/module_codependency_lung.csv
```

All entry points provide `--help`. DepMap identifiers and gene labels are
normalised defensively, pairwise sample sizes are reported, and multiple-testing
adjustment is performed with Benjamini–Hochberg correction.

## Reproducibility tiers and an important provenance note

- **Exact-result verification:** manuscript-level numerical outputs can be
  checked directly from the compact derived tables included in this release.
- **Core recomputation:** co-dependency, module, genotype, structural and
  sample-level single-cell analyses can be rerun from public source data;
  transparent reconstruction entry points are also provided for Leiden and TCGA.
- **Machine-learning reconstruction:** the original XGBoost training driver
  that produced the archived SHAP table was not retained. A leakage-resistant,
  fully specified reanalysis in `src/04_xgboost_shap.py` nevertheless recovered
  the same five leading BAF features and gave out-of-fold $R^2=0.390$ (archived
  analysis: 0.376). Both outputs and the new metrics are retained.

This distinction preserves the historical output while providing an independently
executable sensitivity analysis rather than overstating its provenance.

## Data access

The release contains no restricted patient-level data. DepMap 24Q4, Project
Score, RCSB PDB and GEO source locations and accessions are listed in
`metadata/data_manifest.tsv`. Users must follow the terms of each upstream
resource. No raw source dataset is relicensed by this repository.

## Citation and licence

Use `CITATION.cff` or the archived release DOI
[`10.5281/zenodo.22218065`](https://doi.org/10.5281/zenodo.22218065). Code is
licensed under MIT (`LICENSE`); manuscript text and third-party data retain
their respective rights. Correspondence: Guohua Fan, `fgh202604@163.com`.
