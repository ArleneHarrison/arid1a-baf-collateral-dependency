"""Read-only statistical replays; output audit artifacts only in the chosen run.

Reference computation uses explicit lineage dummies and full OLS designs,
not the FWL residualization used by the analysis. Predict exact unique design
rows and expand; do not round residuals or ranks.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import r12_dependency_bridge as bridge


def direct_residual(y, design, group):
    beta = np.linalg.lstsq(design, y, rcond=None)[0]
    unique, inverse = np.unique(design, axis=0, return_inverse=True)
    result = y-(unique @ beta)[inverse]
    # A lineage intercept fits a completely constant within-lineage design and
    # outcome exactly. Enforce that algebraic zero, including singleton groups;
    # do not round other residuals or merge merely nearby observations.
    for label in np.unique(group):
        rows = np.flatnonzero(group == label)
        if np.all(y[rows] == y[rows[0]]) and np.all(design[rows] == design[rows[0]]):
            result[rows] = 0.
    return result


def direct_estimate(values, group):
    _, labels = np.unique(group, return_inverse=True)
    dummy = np.eye(len(np.unique(labels)))[labels]
    standardized = []
    for j in range(6):
        start = 3+3*j
        design = np.column_stack([dummy, values[:, [start+1, start+2]]])
        res = direct_residual(values[:, start], design, group)
        standardized.append((res-res.mean())/res.std(ddof=1))
    module = np.column_stack(standardized).mean(axis=1)
    base = np.column_stack([dummy, values[:, :3]])
    full = np.column_stack([base, module])
    module_res = direct_residual(module, base, group)
    rho, delta_r2 = [], []
    for j in range(2):
        y = values[:, 21+j]
        base_res = direct_residual(y, base, group)
        full_res = direct_residual(y, full, group)
        rho.append(float(spearmanr(module_res, base_res).statistic))
        delta_r2.append(float((base_res@base_res-full_res@full_res)/np.sum((y-y.mean())**2)))
    return np.array(rho), np.array(delta_r2)


def main(root, version):
    spec = json.loads((root/'r12_dependency_bridge_spec.json').read_text())
    current = root/('dependency_bridge_results_'+version)
    values, group, meta, flow, _ = bridge.prepare_inputs(
        root/'analysis_ready', root/'analysis_ready/dependency.csv', spec, 'lineage')
    saved_samples = pd.read_csv(current/'r12_dependency_bridge_samples.csv', index_col=0)
    assert meta.index.tolist() == saved_samples.index.tolist()
    assert meta.analysis_patient_id.nunique() == len(meta)
    observed = pd.read_csv(current/'r12_dependency_bridge_summary.csv').set_index('target').loc[spec['dependency_targets']]
    boot = np.load(current/'r12_dependency_bridge_bootstrap.npz')['rho']
    direct_rho, direct_delta = direct_estimate(values, group)
    rng = np.random.default_rng(spec['seed'])
    strata = [np.flatnonzero(group == k) for k in np.unique(group)]
    replay = []
    for b in range(50):
        index = np.concatenate([rng.choice(s, len(s), replace=True) for s in strata])
        ref, _ = direct_estimate(values[index], group[index])
        for j, target in enumerate(spec['dependency_targets']):
            replay.append(dict(replicate=b+1, target=target, saved_rho=float(boot[b, j]),
                               explicit_OLS_rho=float(ref[j]), absolute_difference=float(abs(boot[b, j]-ref[j]))))
    replay = pd.DataFrame(replay)
    # Scan input-defined, rare bootstrap configurations rather than selecting
    # unusual outcomes: an odd-sized lineage drawn from one source observation.
    # Such a constant block must center to exact zero, not floating-point dust.
    rng = np.random.default_rng(spec['seed'])
    bounds = np.cumsum([0]+[len(s) for s in strata])
    rare = []
    for b in range(len(boot)):
        index = np.concatenate([rng.choice(s, len(s), replace=True) for s in strata])
        constant_blocks = []
        for k, stratum in enumerate(strata):
            block = index[bounds[k]:bounds[k+1]]
            if len(block) > 1 and len(block) % 2 == 1 and np.all(block == block[0]):
                constant_blocks.append((k, block))
        if not constant_blocks:
            continue
        centered, _ = bridge.coupling.within_center(values[index], group[index])
        ref, _ = direct_estimate(values[index], group[index])
        for k, block in constant_blocks:
            max_centered = float(np.max(np.abs(centered[bounds[k]:bounds[k+1]])))
            rare.append(dict(replicate=b+1, source_sample=meta.index[block[0]],
                             lineage_rows=len(block), max_abs_centered_constant_block=max_centered,
                             max_abs_rho_difference_vs_explicit_OLS=float(np.max(np.abs(boot[b]-ref)))))
    rare = pd.DataFrame(rare, columns=['replicate','source_sample','lineage_rows',
                                     'max_abs_centered_constant_block','max_abs_rho_difference_vs_explicit_OLS'])
    rare_rho_error = float(rare.max_abs_rho_difference_vs_explicit_OLS.max()) if len(rare) else 0.
    constant_block_error = float(rare.max_abs_centered_constant_block.max()) if len(rare) else 0.
    point_rho_error = float(np.max(np.abs(observed.rho.to_numpy()-direct_rho)))
    point_r2_error = float(np.max(np.abs(observed.delta_in_sample_r_squared.to_numpy()-direct_delta)))
    replay_error = float(replay.absolute_difference.max())
    new_manifest = json.loads((current/'r12_dependency_bridge_manifest.json').read_text())
    comparisons = {}
    for old_version, directory in [('v1', root/'dependency_bridge_results'),
                                    ('v2', root/'dependency_bridge_results_v2')]:
        if old_version == version:
            continue
        old_manifest = json.loads((directory/'r12_dependency_bridge_manifest.json').read_text())
        old = pd.read_csv(directory/'r12_dependency_bridge_summary.csv').set_index('target').loc[spec['dependency_targets']]
        old_boot = np.load(directory/'r12_dependency_bridge_bootstrap.npz')['rho']
        comparisons[old_version] = dict(
            helper_sha256=old_manifest['helper_script_sha256'],
            same_inputs=old_manifest['input_sha256'] == new_manifest['input_sha256'],
            same_specification=old_manifest['specification_sha256'] == new_manifest['specification_sha256'],
            same_bridge_script=old_manifest['script_sha256'] == new_manifest['script_sha256'],
            bootstrap_max_abs_rho_change=np.max(np.abs(boot-old_boot), axis=0).tolist(),
            bootstrap_changed_above_1e_minus12=np.sum(np.abs(boot-old_boot)>1e-12, axis=0).tolist(),
            point_rho_change=(observed.rho-old.rho).to_list(),
            rho_ci_low_change=(observed.rho_ci_low-old.rho_ci_low).to_list(),
            rho_ci_high_change=(observed.rho_ci_high-old.rho_ci_high).to_list(),
            p_change=(observed.centered_bootstrap_approx_p-old.centered_bootstrap_approx_p).to_list(),
            q_change=(observed.q_BH_two_dependency_targets-old.q_BH_two_dependency_targets).to_list())
    same_inputs = all(c['same_inputs'] for c in comparisons.values())
    same_spec = all(c['same_specification'] for c in comparisons.values())
    same_script = all(c['same_bridge_script'] for c in comparisons.values())
    # Independently recalculate the stored percentile intervals, centered P and
    # the closed-form two-test BH adjustment from every saved bootstrap value.
    ref_ci, ref_p, valid_counts = [], [], []
    for j in range(2):
        finite = boot[np.isfinite(boot[:, j]), j]
        valid_counts.append(len(finite))
        ref_ci.append(np.quantile(finite, [.025, .975]))
        ref_p.append((1+np.count_nonzero(np.abs(finite-observed.rho.iloc[j]) >= abs(observed.rho.iloc[j])))/(1+len(finite)))
    ref_ci, ref_p = np.array(ref_ci), np.array(ref_p)
    order = np.argsort(ref_p)
    ref_q = np.empty(2)
    ref_q[order[1]] = min(1., ref_p[order[1]])
    ref_q[order[0]] = min(1., 2*ref_p[order[0]], ref_q[order[1]])
    ci_error = float(np.max(np.abs(observed[['rho_ci_low','rho_ci_high']].to_numpy()-ref_ci)))
    p_error = float(np.max(np.abs(observed.centered_bootstrap_approx_p.to_numpy()-ref_p)))
    q_error = float(np.max(np.abs(observed.q_BH_two_dependency_targets.to_numpy()-ref_q)))
    counts_correct = bool(np.array_equal(observed.n_valid_bootstrap.to_numpy(), valid_counts)
                          and len(boot) == spec['bootstrap_replicates'])
    helper_matches_run = bridge.coupling.file_hash(bridge.coupling.__file__) == new_manifest['helper_script_sha256']
    report = dict(
        passed=bool(point_rho_error <= 1e-12 and point_r2_error <= 1e-12 and replay_error <= 1e-12
                    and rare_rho_error <= 1e-12 and constant_block_error == 0.
                    and ci_error <= 1e-12 and p_error <= 1e-12 and q_error <= 1e-12
                    and counts_correct and helper_matches_run
                    and same_inputs and same_spec and same_script),
        audited_version=version,
        shared_n=len(meta), unique_analysis_patients=int(meta.analysis_patient_id.nunique()),
        lineages=flow['shared_complete_case_lineages'], same_historical_input_hashes=same_inputs,
        same_historical_specification=same_spec, same_historical_bridge_script=same_script,
        current_helper_matches_production_run=helper_matches_run,
        helper_sha256_current=new_manifest['helper_script_sha256'],
        point_max_absolute_rho_error_vs_explicit_OLS=point_rho_error,
        point_max_absolute_delta_r2_error_vs_explicit_OLS=point_r2_error,
        bootstrap_replays=50, bootstrap_max_absolute_rho_error_vs_explicit_OLS=replay_error,
        rare_odd_constant_block_replays=len(rare),
        rare_odd_constant_block_max_abs_centered_value=constant_block_error,
        rare_odd_constant_block_max_abs_rho_error=rare_rho_error,
        confidence_interval_max_abs_recalculation_error=ci_error,
        centered_p_max_abs_recalculation_error=p_error,
        BH_q_max_abs_recalculation_error=q_error,
        bootstrap_valid_counts=valid_counts, bootstrap_counts_match_summary_and_spec=counts_correct,
        target_order=spec['dependency_targets'],
        current_minus_historical_versions=comparisons,
        audit_script_sha256=bridge.coupling.file_hash(__file__))
    replay.to_csv(current/'r12_dependency_bridge_independent_bootstrap_replay.csv', index=False)
    rare.to_csv(current/'r12_dependency_bridge_rare_constant_block_replay.csv', index=False)
    (current/'r12_dependency_bridge_independent_audit.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2), flush=True)
    if not report['passed']:
        raise AssertionError('Independent explicit-OLS replay audit did not pass the stated tolerance')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', required=True)
    parser.add_argument('--version', choices=['v2','v3'], default='v3')
    args = parser.parse_args()
    main(Path(args.root), args.version)
