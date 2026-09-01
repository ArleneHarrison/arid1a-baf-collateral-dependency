from pathlib import Path
import json

import numpy as np
import pandas as pd


ROOT = Path(__file__).parents[1]


def test_structural_binary_and_graded_claims():
    table = pd.read_csv(ROOT / "results" / "derived" / "r3_structural_claim_audit.csv")
    binary = table.query("analysis == 'binary_contact_vs_noncontact'").set_index("cohort")
    assert binary.loc["lung", "n_contacting"] == 7
    assert binary.loc["lung", "n_noncontacting"] == 8
    assert np.isclose(binary.loc["lung", "estimate"], 0.4355973036073978)
    assert np.isclose(binary.loc["lung", "effect_size"], 0.8571428571428572)
    assert np.isclose(binary.loc["pan-cancer", "estimate"], 0.31434597940183)
    graded = table.query("analysis == 'graded_contact_count_within_6LTJ_contacts'").set_index("cohort")
    assert graded.loc["lung", "p_value"] > 0.05
    assert graded.loc["pan-cancer", "p_value"] > 0.05


def test_sample_level_single_cell_claims():
    summary = pd.read_csv(ROOT / "results" / "derived" / "r3_scrna_sample_level_summary.csv").iloc[0]
    samples = pd.read_csv(ROOT / "results" / "derived" / "r3_scrna_epithelial_sample_scores.csv")
    assert summary.n_samples == len(samples) == 48
    assert samples.Sample.is_unique
    assert samples.ARID1A_group.value_counts().to_dict() == {"high": 24, "low": 24}
    assert np.isclose(summary.spearman_rho_ARID1A_vs_coreBAF, 0.41576204950065127)
    assert summary.spearman_p < 0.01
    assert summary.mannwhitney_p < 0.05


def test_genomewide_tables_are_unique_and_ranked():
    for cohort in ("lung", "pancancer"):
        table = pd.read_csv(ROOT / "results" / "derived" / f"genomewide_codependency_{cohort}.csv")
        assert table.gene.is_unique
        assert table.rank_desc.is_monotonic_increasing
        assert table.rho.is_monotonic_decreasing


def test_recomputed_machine_learning_sensitivity():
    metrics = pd.read_json(ROOT / "results" / "derived" / "ml_xgboost_recomputed_metrics.json", typ="series")
    importance = pd.read_csv(ROOT / "results" / "derived" / "ml_shap_arid1a_feature_importance_RECOMPUTED.csv")
    assert metrics.n_models == 1178
    assert np.isclose(metrics.oof_r2, 0.3903418708231198)
    assert importance.head(5).feature.tolist() == ["SMARCB1", "SMARCE1", "SMARCC1", "ARID1B", "DPF2"]


def test_selection_symmetric_nested_validation():
    summary = json.loads(
        (ROOT / "results" / "derived" / "r5_nested_validation" / "nested_summary.json")
        .read_text(encoding="utf-8")
    )
    assert summary["n_paired_complexes_analyzed"] == 266
    assert summary["n_splits_per_module"] == 30
    assert np.isclose(summary["real_median_concordance"], 0.6)
    assert np.isclose(summary["matched_null_median_concordance"], 1 / 3)
    assert np.isclose(summary["paired_median_difference"], 0.25)
    assert summary["paired_wilcoxon_one_sided_P"] < 1.1e-19


def test_structure_quality_sensitivity():
    summary = json.loads(
        (ROOT / "results" / "derived" / "r5_structure_validation" /
         "structure_validation_summary.json").read_text(encoding="utf-8")
    )
    assert summary["all_six_contacts_persist_4_0_to_5_0A"]
    assert summary["all_six_contacts_persist_after_bottom10pct_Q_filter"]
    global_q = pd.read_csv(
        ROOT / "results" / "derived" / "r5_structure_validation" /
        "structure_global_validation.csv"
    ).set_index("pdbid")
    assert np.isclose(global_q.loc["6LTJ", "global_Q_score"], 0.007)
    assert np.isclose(global_q.loc["6LTH", "global_Q_score"], 0.526)


def test_abundance_defined_nulls_and_positive_control():
    nulls = pd.read_csv(
        ROOT / "results" / "derived" / "r5_arid1a_null_modules" /
        "arid1a_low_module_sensitivity.csv"
    )
    assert len(nulls) == 42
    assert nulls.BH_q_within_threshold.min() > 0.05
    controls = pd.read_csv(
        ROOT / "results" / "derived" / "r5_arid1a_null_modules" /
        "smarca4_low_smarca2_positive_control.csv"
    )
    control = controls.loc[np.isclose(controls.bottom_fraction, 0.10)].iloc[0]
    assert np.isclose(control.mean_difference, -1.01184847648039)
    assert control.BH_q_across_thresholds < 3.1e-9
