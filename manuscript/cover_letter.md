# Cover Letter

**To:** The Editors, *Computational and Structural Biotechnology Journal*

**Re:** Structure-informed CRISPR co-dependency maps the ARID1A-associated BAF module across cancer cell lines

Dear Editors,

We submit this manuscript for consideration as a Research Article in the General section of *Computational and Structural Biotechnology Journal*.

ARID1A is recurrently inactivated in non-small cell lung cancer, but continuous co-dependency and loss-conditioned synthetic lethality are often treated as interchangeable. We distinguish these estimands and use genome-scale CRISPR profiles to ask a more focused question: which parts of the human BAF complex share an essentiality pattern with ARID1A across cancer cell lines?

Five methodologically distinct analyses of the same dependency matrix converge on core-shared and cBAF-associated subunits. The pattern is robust to pan-essentiality and lineage controls, is reproduced across 19 tissues and a second CRISPR processing resource, and aligns with physical contact in the assembled human BAF structure. Contacting subunits show substantially higher co-dependency than non-contacting or mutually exclusive subunits, with a lung median difference of 0.436 (95% bootstrap interval 0.293–0.584; Cliff’s δ=0.86; exact P=0.0037). Importantly, contact count does not show a positive dose–response within resolved contacts, so we use structure as a binary prior rather than overstate it as a quantitative rule.

The broader methodological contribution is a reproducibility-based approach to protein-complex module refinement. Excluding the anchor prevents self-correlation inflation, and the identity of a dilutive member is more reproducible across random half-splits for 227 curated complexes than for size-matched random sets. The manuscript also reports evidence boundaries directly: no genotype-conditional dependency survives multiple-testing correction, sample-level single-cell effects are small, and the study does not claim an ARID1A-loss synthetic-lethal or pharmacological vulnerability.

A public versioned repository accompanies the submission, including a pinned environment, public-source manifest, compact derived tables, executable core workflows and automated numerical checks. A newly specified XGBoost reanalysis independently recovers the same five leading BAF features as the historical SHAP output (out-of-fold $R^2=0.390$), closing the principal machine-learning provenance gap without overstating it as independent biological evidence.

This work fits the journal’s General section because it combines computational functional genomics with experimental macromolecular structures to clarify how components of a multi-protein complex work together. Its principal advance is not another cancer biomarker or public-dataset signature, but a general, testable framework for connecting co-essentiality, structural contact and resampling-based module definition.

This manuscript is original, has not been published, and is not under consideration elsewhere. All authors have read and approved the submitted version and agree to the journal’s publication-ethics policies. The authors declare no competing interests. The work used only public datasets and required no new human or animal experimentation. Generative AI assistance is disclosed in the manuscript and was not used to create or alter primary data.

The manuscript was previously screened by *Briefings in Bioinformatics* and *NAR Genomics and Bioinformatics* and was not sent for external review. It has since been substantially reframed: therapeutic and loss-conditioned claims were removed, statistical sensitivity analyses were added, cell-level pseudoreplication was corrected, and the principal contribution is now protein-complex co-essentiality and module refinement.

Thank you for considering our work.

Yours sincerely,

**Guohua Fan, MD**  
Corresponding author, on behalf of all authors  
Renmin Hospital of Wuhan University  
99 Zhangzhidong Road, Wuchang District, Wuhan, Hubei, China  
fgh202604@163.com  
ORCID: 0009-0005-7470-7863
