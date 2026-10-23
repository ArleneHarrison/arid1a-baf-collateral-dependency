#!/usr/bin/env python3
"""Verify the four predeclared CPTAC analyses and produce an unadjusted summary."""
from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd
import argparse

ROOT = Path("/root/autodl-tmp/arid1a_baf_r12_cptac")
OUT = ROOT / "analysis_results"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    global OUT
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', default='analysis_results')
    args = parser.parse_args()
    OUT = ROOT / args.outdir
    report = json.loads((OUT / "run_report.json").read_text())
    assert "finished_utc" in report and len(report["runs"]) == 4
    assert all(r["exit_code"] == 0 for r in report["runs"])
    assert sha(OUT / "frozen/r12_protein_multiomics.py") == report["script_sha256"]
    assert sha(OUT / "frozen/r12_analysis_spec.json") == report["spec_sha256"]
    assert all(sha(OUT / name) == checksum for name, checksum in report["artifact_sha256"].items())
    summaries, checks = [], []
    for cohort, n_base, n_purity in [("LUAD", 109, 108), ("LUSC", 108, 107)]:
        inp = ROOT / "pancan_release/analysis_ready" / cohort
        for model, expected in [("purity", n_purity), ("base", n_base)]:
            target = OUT / f"{cohort}_{model}"
            man = json.loads((target / "r12_coupling_manifest.json").read_text())
            assert man["bootstrap"] == 10000
            assert man["spec_sha256"] == report["spec_sha256"]
            assert all(sha(inp / f"{k}.csv") == v for k, v in man["inputs"].items())
            samples = pd.read_csv(target / "r12_analysis_samples.csv", index_col=0)
            assert len(samples) == expected, (cohort, model, len(samples))
            assert samples.index.is_unique and samples.patient_id.is_unique
            assert samples.sample_type.eq("tumor").all()
            assert samples.cohort.eq(cohort).all()
            if model == "purity":
                assert np.isfinite(samples.purity).all()
            archive = np.load(target / "r12_shared_refit_bootstrap.npz")
            assert archive["values"].shape[0] == 10000
            assert archive["values"].shape[1] == len(archive["genes"])
            assert archive["values"].shape[2] == len(archive["layers"]) == 4
            contrast = pd.read_csv(target / "r12_coupling_contrasts.csv")
            primary = contrast.loc[contrast.primary].copy()
            assert len(primary) == 2
            assert primary.n_bootstrap.eq(10000).all()
            assert primary.evaluable.all()
            assert primary.n_valid_bootstrap.ge(9500).all()
            assert primary.centered_bootstrap_approx_p.between(0, 1).all()
            primary.insert(0, "cohort", cohort)
            primary.insert(1, "model", model)
            primary.insert(2, "n_patients", len(samples))
            primary["final_external_FDR_status"] = "Pending root BH across four purity-model primary tests; within-cohort q is not final"
            summaries.append(primary)
            checks.append({"cohort": cohort, "model": model, "n_patients": len(samples),
                           "bootstrap": 10000, "npz_shape": list(archive["values"].shape),
                           "all_input_hashes_match": True, "sample_integrity": "PASS",
                           "min_valid_primary_bootstraps": int(primary.n_valid_bootstrap.min())})
    summary = pd.concat(summaries, ignore_index=True)
    summary.to_csv(OUT / "r12_cptac_four_runs_primary_summary.csv", index=False)
    if OUT != ROOT / 'analysis_results':
        comparisons = []
        for prior in ['analysis_results', 'analysis_results_v2']:
            if ROOT / prior == OUT:
                continue
            for cohort in ['LUAD', 'LUSC']:
                for model in ['purity', 'base']:
                    name = f'{cohort}_{model}/r12_coupling_contrasts.csv'
                    source = ROOT / prior / name
                    if not source.exists():
                        continue
                    old = pd.read_csv(source, index_col=0)
                    new = pd.read_csv(OUT / name, index_col=0)
                    for test in ['primary_own_adjusted', 'primary_raw_layer_difference']:
                        for field in ['estimate','ci_low','ci_high','centered_bootstrap_approx_p']:
                            comparisons.append({'prior_version':prior,'new_version':args.outdir,'cohort':cohort,
                                                'model':model,'test':test,'field':field,
                                                'prior_value':old.loc[test,field],'new_value':new.loc[test,field],
                                                'new_minus_prior':new.loc[test,field]-old.loc[test,field]})
        pd.DataFrame(comparisons).to_csv(OUT/'r12_cptac_numerical_version_comparison.csv',index=False)
    (OUT / "r12_cptac_result_validation.json").write_text(json.dumps({"checks": checks,
        "all_formal_result_hashes_match": True, "primary_model": "purity",
        "external_multiplicity": "root to pool four tests; no final q computed here"}, indent=2))
    print(summary[["cohort", "model", "n_patients", "test", "estimate", "ci_low", "ci_high",
                   "centered_bootstrap_approx_p", "n_valid_bootstrap"]].to_string(index=False))


if __name__ == "__main__":
    main()
