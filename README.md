# ARID1A–BAF collateral dependency — analysis scripts

Analysis code accompanying the manuscript *"Structure-informed co-dependency
analysis localises the collateral dependency of ARID1A-deficient non-small
cell lung cancer to the subunits ARID1A physically contacts"* (submitted to
*Briefings in Bioinformatics*).

This repository is private during peer review; it will be made public on
acceptance. Reviewers requiring access should contact the corresponding
author (Guohua Fan, fgh202604@163.com) or the editorial office.

## Contents

46 Python/R scripts covering:
- CRISPR co-dependency screen processing and the five convergent
  hypothesis-agnostic methods (pan-essentiality/lineage controls, replication
  across an independent dataset and 19 tissues)
- Cryo-EM structural contact analysis (BAF complex) and module refinement
  (SWI/SNF, PRC2)
- Self-correlation artefact quantification across anchor–module pairs
- AlphaFold-Multimer ipTM vs. interface-accuracy benchmarking
- Mendelian randomization (MR-Base / Wald ratio)
- Figure generation for the manuscript

No raw patient-level data is included; all source datasets are public and
cited in the manuscript Methods with their accessions.
