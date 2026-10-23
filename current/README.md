# Reproducibility materials: v1.1.1 (23 October 2026)

Structure-Guided Proteogenomics Identifies ARID1A-Associated cBAF Protein Co-variation in Lung Cancer

Version DOI: https://doi.org/10.5281/zenodo.23587002
Repository: https://github.com/ArleneHarrison/arid1a-baf-collateral-dependency/tree/v1.1.1

This release tests a fixed structure-selected six-protein cBAF panel. Its central
result is replicated residual protein co-variation, not contact-specific causality,
complex occupancy, or whole-panel synthetic lethality. No new biological experiment
or mass-spectrometry acquisition was performed.

## File and data organization

- Main article: four figures and one table; Supporting Information: twenty figures,
  seven tables and extended methods.
- Data S1: grouped ML, cross-complex controls, HAP1 and supporting workflows.
- Data S2 (three ZIPs): genetics, deposited structure and sustained-fitness sources.
- Data S3 (three ZIPs): exact fixed-panel protein/RNA/CN inputs, v3 scripts,
  all bootstrap arrays and independent numerical audits.
- Data S4 (three ZIPs): frozen structural predictions/mapping, exploratory only.
- Data S5: cross-histology and subunit influence analyses.
- Data S6 (two ZIPs): proteome-wide scan, input NPZs, matched candidates and exact null
  draws; `replay_ready.py` and `independent_audit.py` reproduce/audit these endpoints.
- Data S7: original published Roumeliotis workbook and all iPSC extraction/audit code.
- Data S8: current article sources, figure mapping, original plotted source values,
  TOC source and plotting code, structural display inputs, reference metadata audit.

Each numbered ZIP part is a separate ordinary ZIP, not a split binary archive.
Extract all parts of one dataset into the same empty directory, preserving paths.
Keep S1--S8 in separate directories because their generic README/results names differ.
No controlled-access source files, patient identifiers, credentials, referee reports,
cover letters or account configuration are included in the public release.

## Minimal checks (no new downloads)

```sh
python -m pip install numpy pandas scipy matplotlib pymupdf openpyxl
python current/validate_release.py
cd current
python plot_toc.py --source TOC_source_data.csv --outdir replay_toc
latexmk -pdf -interaction=nonstopmode -halt-on-error manuscript.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error supplementary_information.tex
```

TeX Live 2025, `achemso` 3.13h and BibTeX were used. Figures are deterministic
exports. S3/S6 contain frozen analysis inputs; the data-specific README and CLI
`--help` specify exact replay parameters. Full bootstrap/omics computations should
run on a research server, not be inferred from figure images. Older root-level
scripts/results are legacy provenance, not the v1.1.1 primary estimand.

## Data rights and provenance

MIT applies to original project software, not to third-party source material.
Original text, schematics and derived numeric tables in `current/` are CC BY 4.0.
Public source data, coordinates and reproduced author-provided workbooks retain
their original access/reuse conditions and citations. No source paper's license is
replaced by the repository software license. The Ye et al. paper is CC BY-NC-ND 4.0;
its original Source_Data.xlsx is NOT redistributed here. S8 supplies extracted
factual numeric values, their source DOI and the original workbook SHA256
e772e583ad1ec537111f2ac68a2a667736b361fd96c8906b5f35aa7a0664c086.
The displayed EMD-0974 density is a processed deposited-density projection, not a
raw particle image or new cryo-EM acquisition. Predicted structures are labelled
exploratory and are not experimental validation.

## Change from prior releases

v1.1.1 refocuses the article on a bounded protein co-variation question, includes
independent LUAD/LUSC analyses and their null/sensitivity controls, distinguishes
ARID1B fitness evidence from whole-panel null results, corrects reference metadata,
and supplies ACS-format sources and original TOC artwork. Old tags remain intact.
The archived manuscript is an unpublished manuscript, not a peer-reviewed article.
