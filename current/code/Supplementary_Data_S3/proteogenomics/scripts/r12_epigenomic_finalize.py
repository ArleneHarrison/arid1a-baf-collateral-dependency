#!/usr/bin/env python3
"""Independent numeric audit and provenance for the server-only R12 run."""
import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

root = Path("/root/autodl-tmp/arid1a_baf_r12_epigenomic")
out = root / "outputs"
if platform.system() != "Linux":
    raise SystemExit("Run on the authorized server only")
res = pd.read_csv(out / "r12_rna_primary_interaction.csv")
expr = pd.read_csv(out / "r12_rna_selected_gene_FPKM.csv").set_index("gene")
expr = expr.loc[res.gene]
v = np.log2(expr.iloc[:, 1:].to_numpy(dtype=float) + 1)
direct = v[:, 6:8].mean(axis=1) - v[:, 4:6].mean(axis=1) - v[:, 2:4].mean(axis=1) + v[:, :2].mean(axis=1)
coefficient_error = float(np.max(np.abs(direct - res.interaction_log2_effect)))
q = multipletests(res.p_interaction.to_numpy(), method="fdr_bh")[1]
q_error = float(np.max(np.abs(q - res.q_interaction_BH_all_genes)))
assert coefficient_error < 1e-10
assert q_error < 1e-10
assert len(res.gene) == len(set(res.gene))
modules = pd.read_csv(out / "r12_prespecified_module_gene_response.csv")
module_summary = modules.groupby("module").agg(
    genes_requested=("gene", "size"), genes_tested=("interaction_log2_effect", "count"),
    minimum_gene_q=("q_interaction_BH_all_genes", "min"),
    descriptive_median_interaction=("interaction_log2_effect", "median"))
module_summary["inference"] = "Descriptive summary only; no module-level dependency or specificity claim"
module_summary.to_csv(out / "r12_prespecified_module_response_summary.csv")
diagnostics = pd.read_csv(out / "r12_rna_secondary_diagnostic_contrasts.csv")
diagnostic_summary = diagnostics.assign(significant=diagnostics.q_BH_all_four_secondary_contrasts < 0.05).groupby("contrast").agg(
    genes_tested=("gene", "size"), genes_q05_across_all_secondary=("significant", "sum"))
diagnostic_summary["purpose"] = "Diagnostic known-biology contrast; cannot replace the primary interaction"
diagnostic_summary.to_csv(out / "r12_secondary_contrast_summary.csv")
longest = pd.read_csv(out / "r12_rna_sensitivity_longest_transcript.csv")
representatives = res[["gene", "transcript_id"]].merge(longest[["gene", "transcript_id"]], on="gene", suffixes=("_primary", "_longest"))
changed = int((representatives.transcript_id_primary != representatives.transcript_id_longest).sum())
checks = {"host": platform.node(), "completed_utc": datetime.now(timezone.utc).isoformat(),
          "python_vs_R_interaction_max_absolute_error": coefficient_error,
          "Python_statsmodels_vs_R_BH_max_absolute_error": q_error,
          "unique_gene_rows_confirmed": True,
          "alternate_transcript_rule_changed_representatives": changed,
          "transcript_sensitivity_note": "If zero representatives change, this is a mapping consistency check rather than independent sensitivity evidence.",
          "primary_genes_tested": len(res), "primary_genes_q_lt_005": int((q < 0.05).sum()),
          "module_alias_mapping": {"BICRA": "GLTSCR1", "BICRAL": "GLTSCR1L"},
          "not_count_confirmatory": True,
          "chromatin_formal_interaction_not_run": "Only thresholded per-sample called peaks are deposited; common-interval quantitative counts not available"}
(out / "r12_numeric_audit.json").write_text(json.dumps(checks, ensure_ascii=False, indent=2), encoding="utf-8")
def describe(path):
    return {"file": str(path), "bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
provenance = {"host": platform.node(), "executed_on": "User NVIDIA RTX 5090 Linux server",
              "network_note": "Small public source files relayed via workstation because direct NCBI SSL timed out; all scientific transformations/statistics executed here.",
              "raw_SRA_downloaded": False, "maximum_analysis_threads": 1,
              "started_from_original_study": "GSE101966; DOI 10.7554/eLife.30506",
              "inputs": [describe(p) for p in sorted((root / "inputs").iterdir()) if p.is_file()],
              "scripts": [describe(p) for p in sorted((root / "scripts").iterdir()) if p.is_file()],
              "outputs": [describe(p) for p in sorted(out.iterdir()) if p.is_file() and p.name != "r12_execution_provenance.json"],
              "completed_utc": datetime.now(timezone.utc).isoformat()}
(out / "r12_execution_provenance.json").write_text(json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(checks, ensure_ascii=False, indent=2))
