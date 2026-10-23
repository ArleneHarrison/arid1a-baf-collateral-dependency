#!/usr/bin/env python3
"""Independent statsmodels/scipy verification of fixed CPTAC primary statistics.

Does not import the production statistics script or its residual/FWL utilities.
Reconstructs samples from input CSVs; uses explicit-intercept statsmodels OLS;
independently aggregates saved bootstrap arrays and computes four-test BH.
"""
from datetime import datetime, timezone
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
import scipy
from scipy.stats import spearmanr
import statsmodels
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

ROOT = Path("/root/autodl-tmp/arid1a_baf_r12_cptac")
OUT = ROOT / "analysis_results"
TOL = 1e-10
REPLAY_DRAWS = 50


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def independent_pair(arr, categorical=None):
    """Explicit-intercept/dummy OLS with no centering/FWL helper imported.

    Compute fitted values once per exactly equal design row, then expand.
    This preserves mathematical ties for repeated bootstrap draws without
    rounding/merging distinct designs or changing the fitted model.
    """
    mask = np.isfinite(arr).all(axis=1)
    a = arr[mask]
    dummy = (pd.get_dummies(np.asarray(categorical)[mask], drop_first=True,
                          dtype=float).to_numpy()
             if categorical is not None else np.empty((len(a), 0)))
    extra = list(range(6, a.shape[1]))
    def fit_resid(j, columns):
        x = sm.add_constant(np.column_stack([a[:, columns], dummy]), has_constant="add")
        fit = sm.OLS(a[:, j], x, missing="raise").fit()
        unique, inverse = np.unique(x, axis=0, return_inverse=True)
        resid = a[:, j] - (unique @ fit.params)[inverse]
        if categorical is not None:
            labels = np.asarray(categorical)[mask]
            for label in np.unique(labels):
                ix = np.flatnonzero(labels == label)
                # Group intercept constrains the residual to exact zero when
                # both the outcome AND this model's complete design are constant.
                # Covers singleton lineages and repeats of one sampled cell;
                # never zero a constant-outcome group with varying covariates.
                if np.all(a[ix, j] == a[ix[0], j]) and np.all(x[ix] == x[ix[0]]):
                    resid[ix] = 0.
        return resid, fit.df_resid
    own_a, df_a = fit_resid(0, [2, 4] + extra)
    own_b, df_b = fit_resid(1, [3, 5] + extra)
    joint_a, df_ja = fit_resid(0, [2, 3, 4, 5] + extra)
    joint_b, df_jb = fit_resid(1, [2, 3, 4, 5] + extra)
    if min(df_a, df_b) < 20 or len(a) < 40:
        raise ValueError("Independent model unexpectedly fails frozen feasibility thresholds")
    def rho(x, y):
        return float(spearmanr(x, y).statistic) if min(np.std(x), np.std(y)) > 1e-10 else np.nan
    return {"protein_raw": rho(a[:, 0], a[:, 1]),
            "rna_raw": rho(a[:, 2], a[:, 3]),
            "protein_own_adjusted": rho(own_a, own_b),
            "protein_joint_adjusted": rho(joint_a, joint_b) if min(df_ja, df_jb) >= 20 else np.nan}, len(a)


def boot_contrast(values, contact_ids, background_ids, layer):
    c, b = values[:, contact_ids, layer], values[:, background_ids, layer]
    valid = np.isfinite(c).all(axis=1) & np.isfinite(b).all(axis=1)
    ans = np.full(values.shape[0], np.nan)
    ans[valid] = np.median(c[valid], axis=1) - np.median(b[valid], axis=1)
    return ans


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-directory", default="analysis_results_v2")
    parser.add_argument("--cellline", action="store_true")
    args = parser.parse_args()
    root = ROOT if not args.cellline else Path('/root/autodl-tmp/arid1a_baf_r12_cellline')
    out = root / args.results_directory
    specfile = (out / 'frozen/r12_analysis_spec.json'
                if (out / 'frozen/r12_analysis_spec.json').exists()
                else ROOT / "analysis_results/frozen/r12_analysis_spec.json")
    prefix = "r12_cellline" if args.cellline else "r12_cptac"
    replay_draws = 0 if args.cellline else REPLAY_DRAWS
    spec = json.loads(specfile.read_text())
    contact = spec["contact_experimental"]
    background = spec["primary_baf_background"]
    assert len(contact) == len(background) == 6
    assert not set(contact) & set(background)
    genes = contact + background
    anchor = spec["anchor"]
    audit = {"date_utc": datetime.now(timezone.utc).isoformat(), "tolerance": TOL,
             "spec_sha256": sha(specfile), "script_sha256": sha(Path(__file__)),
             "independent_engine": "statsmodels OLS explicit intercept + own RNA/CNV + WES purity or full-rank lineage dummies; scipy Spearman",
             "does_not_import_production_statistics": True,
             "tie_stability_rule": "OLS predictions computed once per exactly equal design row and expanded; no rounding of distinct rows",
             "exact_group_zero_rule": "If outcome and complete model design both constant within a categorical group, fitted group-intercept constraint implies residual exactly zero; applied explicitly",
             "replay_bootstrap_draws_per_cohort": replay_draws,
             "versions": {"numpy": np.__version__, "pandas": pd.__version__,
                          "statsmodels": statsmodels.__version__, "scipy": scipy.__version__},
             "cohorts": [], "deviations_over_tolerance": []}
    per_gene, combined = [], []
    for cohort in (["Sanger_DepMap24Q4"] if args.cellline else ["LUAD", "LUSC"]):
        inputs = root / "analysis_ready" if args.cellline else root / "pancan_release/analysis_ready" / cohort
        run = out if args.cellline else out / f"{cohort}_purity"
        assays = {m: pd.read_csv(inputs / f"{m}.csv", index_col=0) for m in ["protein", "rna", "cnv"]}
        md = pd.read_csv(inputs / "metadata.csv", index_col=0)
        shared = sorted(set(md.index).intersection(*(set(m.index) for m in assays.values())))
        md = md.loc[shared]
        if args.cellline:
            assert md.lineage.notna().all() and md.lineage.str.strip().ne('').all()
            md = md.loc[~md.patient_id.duplicated(keep='first')]
        else:
            md = md.loc[np.isfinite(pd.to_numeric(md.purity, errors="coerce"))]
        assert md.patient_id.is_unique and md.index.is_unique
        assert md.sample_type.eq("cancer_cell_line" if args.cellline else "tumor").all() and md.cohort.eq(cohort).all()
        saved_samples = pd.read_csv(run / "r12_analysis_samples.csv", index_col=0)
        assert md.index.tolist() == saved_samples.index.tolist(), "Independent sample reconstruction differs"
        n_expected = 694 if args.cellline else (108 if cohort == "LUAD" else 107)
        assert len(md) == n_expected
        recorded = pd.read_csv(run / "r12_coupling_by_gene.csv", index_col=0)
        original_contrasts = pd.read_csv(run / "r12_coupling_contrasts.csv", index_col=0)
        production_manifest = json.loads((run / 'r12_coupling_manifest.json').read_text())
        assert production_manifest['spec_sha256'] == sha(specfile), 'Frozen specification hash differs'
        points, arrays = {}, {}
        for gene in genes:
            arr = np.column_stack([assays["protein"].loc[md.index, anchor],
                                   assays["protein"].loc[md.index, gene],
                                   assays["rna"].loc[md.index, anchor],
                                   assays["rna"].loc[md.index, gene],
                                   assays["cnv"].loc[md.index, anchor],
                                   assays["cnv"].loc[md.index, gene]])
            if not args.cellline:
                arr = np.column_stack([arr, md.purity])
            arrays[gene] = arr
            points[gene], n = independent_pair(arr, md.lineage.to_numpy() if args.cellline else None)
            for layer, estimate in points[gene].items():
                reference = float(recorded.loc[gene, layer])
                delta = abs(estimate - reference) if np.isfinite([estimate,reference]).all() else (0.0 if np.isnan([estimate,reference]).all() else np.inf)
                per_gene.append({"cohort": cohort, "gene": gene, "group": "contact" if gene in contact else "background",
                                 "layer": layer, "n": n, "independent_rho": estimate,
                                 "production_rho": reference, "absolute_difference": delta})
                if delta > TOL:
                    audit["deviations_over_tolerance"].append({"cohort": cohort, "gene": gene, "layer": layer, "difference": delta})
        def point_contrast(layer):
            return np.median([points[g][layer] for g in contact]) - np.median([points[g][layer] for g in background])
        archive_path = run / "r12_shared_refit_bootstrap.npz"
        archive = np.load(archive_path)
        values = archive["values"]
        saved_genes, saved_layers = archive["genes"].tolist(), archive["layers"].tolist()
        assert values.shape[0] == spec["bootstrap_replicates"] == 10000
        ci = [saved_genes.index(g) for g in contact]
        bi = [saved_genes.index(g) for g in background]
        boot_own = boot_contrast(values, ci, bi, saved_layers.index("protein_own_adjusted"))
        boot_raw_p = boot_contrast(values, ci, bi, saved_layers.index("protein_raw"))
        boot_raw_r = boot_contrast(values, ci, bi, saved_layers.index("rna_raw"))
        tests = [("primary_own_adjusted", point_contrast("protein_own_adjusted"), boot_own),
                 ("primary_raw_layer_difference", point_contrast("protein_raw") - point_contrast("rna_raw"), boot_raw_p - boot_raw_r)]
        for test, estimate, boot in tests:
            valid = boot[np.isfinite(boot)]
            assert len(valid) >= int(.95 * 10000)
            bounds = np.quantile(valid, [.025, .975], method="linear")
            pvalue = (np.count_nonzero(np.abs(valid - estimate) >= abs(estimate)) + 1) / (len(valid) + 1)
            check_values = {"estimate": estimate, "ci_low": bounds[0], "ci_high": bounds[1],
                            "centered_bootstrap_approx_p": pvalue}
            deltas = {k: abs(v - float(original_contrasts.loc[test, k])) for k, v in check_values.items()}
            for k, delta in deltas.items():
                if delta > TOL:
                    audit["deviations_over_tolerance"].append({"cohort": cohort, "test": test, "field": k, "difference": delta})
            combined.append({"cohort": cohort if args.cellline else f"CPTAC_{cohort}", "model": "own_RNA_CNV_lineage" if args.cellline else "own_RNA_CNV_WES_purity",
                             "n_patients": len(md), "test": test, **check_values,
                             "n_valid_bootstrap": len(valid), "n_bootstrap": 10000,
                             "contact_genes": ";".join(contact), "background_genes": ";".join(background),
                             "max_difference_from_production": max(deltas.values())})
        # Additional bounded check: independently refit the first 50 shared draws.
        rng = np.random.default_rng(spec["seed"])
        replay_max = 0.0
        for draw in range(replay_draws):
            indices = rng.choice(np.arange(len(md)), size=len(md), replace=True)
            for gene in genes:
                replay, _ = independent_pair(arrays[gene][indices])
                for layer, estimate in replay.items():
                    delta = abs(estimate - values[draw, saved_genes.index(gene), saved_layers.index(layer)])
                    replay_max = max(replay_max, delta)
                    if delta > TOL:
                        audit["deviations_over_tolerance"].append({"cohort": cohort, "draw": draw,
                            "gene": gene, "layer": layer, "difference": delta})
        audit["cohorts"].append({"cohort": cohort, "n_patients": len(md),
             "sample_reconstruction_exact": True, "n_genes_checked": len(genes),
             "point_max_absolute_difference": max(r["absolute_difference"] for r in per_gene if r["cohort"] == cohort),
             "bootstrap_replay_max_absolute_difference": replay_max,
             "source_npz_sha256": sha(archive_path),
             "production_script_sha256": production_manifest['analysis_script_sha256'],
             "input_sha256": {m: sha(inputs / f"{m}.csv") for m in ["protein", "rna", "cnv", "metadata"]}})
    gene_frame = pd.DataFrame(per_gene)
    gene_frame.to_csv(out / f"{prefix}_independent_gene_point_audit.csv", index=False)
    audit["status"] = "FAIL_REQUIRES_INVESTIGATION" if audit["deviations_over_tolerance"] else "PASS"
    if audit["deviations_over_tolerance"]:
        (out / f"{prefix}_independent_statistics_audit.json").write_text(json.dumps(audit, indent=2))
        raise SystemExit("Independent calculation deviates above tolerance; do not reconcile by changing values")
    final = pd.DataFrame(combined)
    assert len(final) == (2 if args.cellline else 4)
    _, qvalues, _, _ = multipletests(final.centered_bootstrap_approx_p.to_numpy(), method="fdr_bh", alpha=.05)
    final["q_BH_two_primary_within_cohort" if args.cellline else "q_external_four_tests"] = qvalues
    final["external_family"] = "Two prespecified discovery primary tests" if args.cellline else "Two prespecified primary tests x LUAD/LUSC purity models (four tests)"
    final["independent_statistics_audit"] = "PASS"
    filename = "r12_cellline_independent_primary_tests.csv" if args.cellline else "r12_cptac_external_primary_tests.csv"
    final.to_csv(out / filename, index=False)
    audit["BH_implementation"] = "statsmodels.stats.multitest.multipletests(method='fdr_bh')"
    audit["BH_family_size"] = len(final)
    audit["results_sha256"] = sha(out / filename)
    audit["gene_audit_sha256"] = sha(out / f"{prefix}_independent_gene_point_audit.csv")
    (out / f"{prefix}_independent_statistics_audit.json").write_text(json.dumps(audit, indent=2))
    print(final[["cohort", "n_patients", "test", "estimate", "ci_low", "ci_high",
                 "centered_bootstrap_approx_p", "q_BH_two_primary_within_cohort" if args.cellline else "q_external_four_tests"]].to_string(index=False))
    print("INDEPENDENT_STATISTICS_AUDIT_PASS", flush=True)


if __name__ == "__main__":
    main()
