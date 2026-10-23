"""Frozen contact-protein module to ARID1A/ARID1B dependency association.

Use on the user's server. Both targets share complete cases and bootstrap rows.
No target-dependent module weighting, feature selection, or predictive claims.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

import r12_protein_multiomics as coupling


def blank_result(target, n, reason):
    return dict(target=target, n=n, rho=np.nan, point_evaluable=False,
                non_evaluable_reason=reason, minimum_module_gene_df=np.nan,
                conditional_model_df=np.nan, full_model_df=np.nan,
                base_in_sample_r_squared=np.nan, full_in_sample_r_squared=np.nan,
                delta_in_sample_r_squared=np.nan)


def estimate_bridge(values, groups, spec):
    """Columns: anchor P/RNA/CNV, six contact P/RNA/CNV triples, two outcomes."""
    genes, targets = spec['contact_experimental'], spec['dependency_targets']
    n = len(values)
    results = [blank_result(t, n, 'insufficient_shared_complete_cases') for t in targets]
    scores = np.full((n, len(genes) + 2), np.nan)
    diagnostics = []
    if n < spec['minimum_shared_n']:
        return results, scores, diagnostics
    if values.shape[1] != 3 + 3*len(genes) + len(targets) or not np.isfinite(values).all():
        raise ValueError('Bridge estimation requires the fixed layout and common finite rows')
    c, ngroups = coupling.within_center(values, groups)
    standardized = []
    for j, gene in enumerate(genes):
        start = 3 + 3*j
        resid, df = coupling.residual(c[:, start], c[:, [start+1, start+2]], ngroups)
        sd = float(np.std(resid, ddof=1))
        evaluable = df >= spec['minimum_residual_df'] and sd > spec['minimum_original_scale_sd']
        diagnostics.append(dict(gene=gene, residual_model_df=df,
                                residual_sample_sd=sd, evaluable=bool(evaluable),
                                weight_after_standardization=1/len(genes)))
        standardized.append((resid-resid.mean())/sd if evaluable else np.full(n, np.nan))
    min_gene_df = min(d['residual_model_df'] for d in diagnostics)
    if not all(d['evaluable'] for d in diagnostics):
        for result in results:
            result.update(non_evaluable_reason='fixed_module_member_insufficient_df_or_zero_residual_variance',
                          minimum_module_gene_df=min_gene_df)
        return results, scores, diagnostics
    standardized = np.column_stack(standardized)
    module = standardized.mean(axis=1)
    centered_module, _ = coupling.within_center(module[:, None], groups)
    anchor_covariates = c[:, :3]
    module_resid, module_df = coupling.residual(centered_module[:, 0], anchor_covariates, ngroups)
    added_module_rank = int(np.std(module_resid) > spec['minimum_original_scale_sd'])
    scores = np.column_stack([standardized, module, module_resid])
    dependency_start = 3 + 3*len(genes)
    for j, result in enumerate(results):
        dep = values[:, dependency_start+j]
        dep_resid, dep_df = coupling.residual(c[:, dependency_start+j], anchor_covariates, ngroups)
        conditional_df = min(module_df, dep_df)
        result.update(minimum_module_gene_df=min_gene_df, conditional_model_df=conditional_df,
                      full_model_df=dep_df-added_module_rank,
                      non_evaluable_reason='conditional_model_insufficient_df_or_degenerate_residuals')
        if conditional_df >= spec['minimum_residual_df']:
            rho = coupling.spearman(module_resid, dep_resid)
            if np.isfinite(rho):
                result.update(rho=rho, point_evaluable=True, non_evaluable_reason='')
        # This increment is descriptive and algebraically matches the nested OLS
        # fit by FWL. It is not cross-validated and has no separate test family.
        if (dep_df-added_module_rank >= spec['minimum_residual_df']
                and np.std(dep) > spec['minimum_original_scale_sd']
                and np.std(module_resid) > spec['minimum_original_scale_sd']):
            sst = float((dep-dep.mean()) @ (dep-dep.mean()))
            base_sse = float(dep_resid @ dep_resid)
            beta = float(module_resid @ dep_resid / (module_resid @ module_resid))
            full_resid = dep_resid - beta*module_resid
            full_sse = float(full_resid @ full_resid)
            delta = (base_sse-full_sse)/sst
            if delta < -1e-10:
                raise ArithmeticError('Nested-model R-squared increment is unexpectedly negative')
            result.update(base_in_sample_r_squared=1-base_sse/sst,
                          full_in_sample_r_squared=1-full_sse/sst,
                          delta_in_sample_r_squared=max(0., delta))
    return results, scores, diagnostics


def prepare_inputs(directory, dependency_path, spec, lineage):
    assays, meta, flow, exclusions = coupling.load_inputs(directory, spec, lineage, [], True)
    dep = pd.read_csv(dependency_path, index_col=0)
    dep.index = dep.index.astype(str)
    if not dep.index.is_unique or not dep.columns.is_unique:
        raise ValueError('Dependency IDs and gene columns must be unique')
    dep = dep.apply(pd.to_numeric, errors='raise')
    missing_columns, labels, columns = [], [], []
    for gene in [spec['anchor']] + spec['contact_experimental']:
        for assay in ['protein', 'rna', 'cnv']:
            label = assay + ':' + gene
            labels.append(label)
            if gene in assays[assay]:
                columns.append(assays[assay][gene].to_numpy(float))
            else:
                missing_columns.append(label)
                columns.append(np.full(len(meta), np.nan))
    for target in spec['dependency_targets']:
        label = 'dependency:' + target
        labels.append(label)
        if target in dep:
            columns.append(dep[target].reindex(meta.index).to_numpy(float))
        else:
            missing_columns.append(label)
            columns.append(np.full(len(meta), np.nan))
    values = np.column_stack(columns)
    finite = np.isfinite(values)
    complete = finite.all(axis=1)
    for j in np.flatnonzero(~complete):
        absent = [labels[k] for k in np.flatnonzero(~finite[j])]
        exclusions.append(dict(sample_id=meta.index[j], reason='nonfinite_required_bridge_measurements:'+';'.join(absent)))
    flow.update(dependency_rows=len(dep), before_shared_complete_case_filter=len(meta),
                shared_complete_case_rows=int(complete.sum()), missing_required_columns=missing_columns,
                complete_case_columns=labels)
    meta = meta.loc[complete].copy()
    if meta.analysis_patient_id.duplicated().any():
        raise ValueError('Shared complete-case rows contain repeated analysis patients')
    groups = pd.factorize(meta[lineage], sort=True)[0]
    flow['shared_complete_case_lineages'] = int(len(np.unique(groups)))
    flow['shared_complete_case_patients'] = int(meta.analysis_patient_id.nunique())
    flow['shared_complete_case_lineage_counts'] = {
        str(k): int(v) for k, v in meta[lineage].value_counts().sort_index().items()}
    return values[complete], groups, meta, flow, exclusions


def run(args):
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    manifest_path = out / 'r12_dependency_bridge_manifest.json'
    if manifest_path.exists():
        raise FileExistsError('Completed output exists; choose a fresh output directory')
    spec_path = Path(args.spec)
    spec = json.loads(spec_path.read_text(encoding='utf-8'))
    source_path = Path(args.source_spec) if args.source_spec else spec_path.with_name(spec['source_spec_file'])
    source = json.loads(source_path.read_text(encoding='utf-8'))
    genes, targets = spec['contact_experimental'], spec['dependency_targets']
    if (len(genes) != 6 or len(set(genes)) != 6 or targets != ['ARID1A', 'ARID1B']
            or set(genes).intersection(targets) or spec['anchor'] in genes):
        raise ValueError('The frozen six-gene module and two dependency targets must remain unchanged')
    if genes != source['contact_experimental'] or spec['anchor'] != source['anchor']:
        raise ValueError('Bridge panel differs from the documented source specification')
    dependency_path = Path(args.dependency) if args.dependency else Path(args.input) / 'dependency.csv'
    values, groups, meta, flow, exclusions = prepare_inputs(Path(args.input), dependency_path, spec, args.lineage_column)
    results, scores, diagnostics = estimate_bridge(values, groups, spec)
    n_boot = args.bootstrap if args.bootstrap is not None else spec['bootstrap_replicates']
    if n_boot < 2:
        raise ValueError('At least two bootstrap replicates are required; final inference uses the frozen count')
    rho_boot = np.full((n_boot, len(targets)), np.nan)
    strata = [np.flatnonzero(groups == k) for k in np.unique(groups)]
    rng = np.random.default_rng(spec['seed'])
    start = time.monotonic()
    executed = 0
    if any(result['point_evaluable'] for result in results):
        for b in range(n_boot):
            index = np.concatenate([rng.choice(s, len(s), replace=True) for s in strata])
            resampled, _, _ = estimate_bridge(values[index], groups[index], spec)
            rho_boot[b] = [r['rho'] for r in resampled]
            executed += 1
            if (b+1) % max(1, min(500, n_boot)) == 0:
                print(json.dumps(dict(bootstrap=b+1, of=n_boot,
                                      elapsed_seconds=round(time.monotonic()-start, 1))), flush=True)
    for j, result in enumerate(results):
        valid = rho_boot[np.isfinite(rho_boot[:, j]), j]
        sufficient = result['point_evaluable'] and len(valid) >= spec['minimum_valid_bootstrap_fraction']*n_boot
        lo, hi = (coupling.finite_ci(rho_boot[:, j], n_boot, spec['minimum_valid_bootstrap_fraction'])
                  if result['point_evaluable'] else (np.nan, np.nan))
        p = ((1+np.sum(np.abs(valid-result['rho']) >= abs(result['rho'])))/(1+len(valid))
             if sufficient else np.nan)
        if result['point_evaluable'] and not sufficient:
            result['non_evaluable_reason'] = 'valid_bootstrap_fraction_below_predeclared_minimum'
        result.update(target_role=spec['target_roles'][result['target']], rho_ci_low=lo, rho_ci_high=hi,
                      centered_bootstrap_approx_p=p, inferential_evaluable=bool(sufficient),
                      n_valid_bootstrap=len(valid), n_requested_bootstrap=n_boot,
                      n_executed_bootstrap=executed,
                      direction='negative_rho_indicates_stronger_dependency_at_higher_module_score')
    summary = pd.DataFrame(results)
    summary['q_BH_two_dependency_targets'] = coupling.bh(summary.centered_bootstrap_approx_p)
    summary.to_csv(out / 'r12_dependency_bridge_summary.csv', index=False)
    pd.DataFrame(diagnostics, columns=['gene','residual_model_df','residual_sample_sd','evaluable',
                                      'weight_after_standardization']).to_csv(out / 'r12_dependency_bridge_gene_diagnostics.csv', index=False)
    score_columns = ['residual_z_' + g for g in genes] + ['fixed_module_score', 'anchor_adjusted_module_score']
    sample_output = meta.join(pd.DataFrame(scores, columns=score_columns, index=meta.index))
    for j, target in enumerate(targets):
        sample_output['dependency_' + target] = values[:, 3+3*len(genes)+j]
    sample_output.to_csv(out / 'r12_dependency_bridge_samples.csv', index_label='sample_id')
    pd.DataFrame(exclusions, columns=['sample_id','reason']).to_csv(out / 'r12_dependency_bridge_exclusions.csv', index=False)
    np.savez_compressed(out / 'r12_dependency_bridge_bootstrap.npz', rho=rho_boot, targets=np.array(targets))
    manifest = dict(
        specification=spec, source_specification=source, input_flow=flow,
        specification_sha256=coupling.file_hash(spec_path), source_specification_sha256=coupling.file_hash(source_path),
        script_sha256=coupling.file_hash(__file__), helper_script_sha256=coupling.file_hash(coupling.__file__),
        input_sha256={**{name: coupling.file_hash(Path(args.input)/(name+'.csv'))
                         for name in ['protein','rna','cnv','metadata']},
                      'dependency': coupling.file_hash(dependency_path)},
        versions={'python':coupling.platform.python_version(), 'numpy':np.__version__,
                  'pandas':pd.__version__, 'scipy':coupling.scipy.__version__},
        lineage_column=args.lineage_column, bootstrap_requested=n_boot, bootstrap_executed=executed,
        elapsed_seconds=time.monotonic()-start,
        status='completed' if summary.inferential_evaluable.any() else 'completed_non_evaluable',
        inference='Two-target fixed-module conditional association family; in-sample R-squared is descriptive only')
    input_manifest = Path(args.input)/'r12_input_manifest.json'
    if input_manifest.exists():
        manifest['input_preparation_manifest_sha256'] = coupling.file_hash(input_manifest)
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(summary.to_string(index=False), flush=True)


def self_test():
    """Synthetic-only checks; execute on the server before any real association."""
    from tempfile import TemporaryDirectory
    spec = dict(anchor='ARID1A', contact_experimental=['g1','g2','g3','g4','g5','g6'],
                dependency_targets=['ARID1A','ARID1B'], minimum_shared_n=60,
                minimum_residual_df=20, minimum_original_scale_sd=1e-10)
    rng = np.random.default_rng(5024)
    group = np.repeat(np.arange(4), 80)
    values = rng.normal(size=(len(group), 23))
    latent = rng.normal(size=len(group))
    for j in range(6):
        k = 3+3*j
        values[:, k] += .5*values[:, k+1]+.2*values[:, k+2]+.7*latent+.4*group
    values[:, 21] = -.9*latent+.3*values[:, 0]+rng.normal(scale=.3, size=len(group))
    values[:, 22] = .6*latent+.2*values[:, 1]+rng.normal(scale=.4, size=len(group))
    result, scores, diag = estimate_bridge(values, group, spec)
    assert result[0]['rho'] < -.5 and result[1]['rho'] > .4
    assert all(d['residual_model_df'] == 314 for d in diag)
    assert all(r['conditional_model_df'] == 313 and r['full_model_df'] == 312 for r in result)
    assert np.allclose(scores[:, :6].mean(axis=0), 0, atol=1e-12)
    assert np.allclose(scores[:, :6].std(axis=0, ddof=1), 1, atol=1e-12)
    # Direct dummy-variable OLS must agree with FWL rho and incremental R-squared.
    dummy = np.eye(4)[group]
    base = np.column_stack([dummy, values[:, :3]])
    full = np.column_stack([base, scores[:, 6]])
    mod_res = scores[:, 6] - base @ np.linalg.lstsq(base, scores[:, 6], rcond=None)[0]
    for j in range(2):
        y = values[:, 21+j]
        base_res = y-base @ np.linalg.lstsq(base, y, rcond=None)[0]
        full_res = y-full @ np.linalg.lstsq(full, y, rcond=None)[0]
        delta = ((base_res@base_res)-(full_res@full_res))/np.sum((y-y.mean())**2)
        assert np.isclose(result[j]['rho'], coupling.spearman(mod_res, base_res), atol=1e-12)
        assert np.isclose(result[j]['delta_in_sample_r_squared'], delta, atol=1e-12)
    changed_targets = values.copy()
    changed_targets[:, 21:] = rng.normal(size=(len(group), 2))
    _, unchanged_scores, _ = estimate_bridge(changed_targets, group, spec)
    assert np.allclose(scores, unchanged_scores, atol=1e-12)
    perm = rng.permutation(len(group))
    perm_result, _, _ = estimate_bridge(values[perm], group[perm], spec)
    assert np.allclose([r['rho'] for r in result], [r['rho'] for r in perm_result])
    constant_target = values.copy()
    constant_target[:, 21] = 1.
    constant_result, _, _ = estimate_bridge(constant_target, group, spec)
    assert np.isnan(constant_result[0]['rho']) and np.isfinite(constant_result[1]['rho'])
    constant_gene = values.copy()
    constant_gene[:, 3] = constant_gene[:, 4]
    no_module, _, _ = estimate_bridge(constant_gene, group, spec)
    assert all(not r['point_evaluable'] for r in no_module)
    short, _, _ = estimate_bridge(values[:20], group[:20], spec)
    assert all(not r['point_evaluable'] for r in short)
    draw = np.concatenate([rng.choice(np.flatnonzero(group == k), 80, replace=True) for k in range(4)])
    resampled, _, _ = estimate_bridge(values[draw], group[draw], spec)
    assert all(np.isfinite(r['rho']) for r in resampled)
    # Input integrity: both outcomes share one complete-case table, and a later
    # same-patient cell line cannot replace the predetermined first cell line.
    ids = ['CELL' + str(i).zfill(3) for i in range(len(group))]
    with TemporaryDirectory(prefix='r12_bridge_selftest_') as temp:
        directory = Path(temp)
        for k, assay in enumerate(['protein', 'rna', 'cnv']):
            positions = [k] + [3+3*j+k for j in range(6)]
            pd.DataFrame(values[:, positions], index=ids, columns=['ARID1A']+spec['contact_experimental']).to_csv(directory/(assay+'.csv'))
        metadata = pd.DataFrame(dict(lineage=['lineage_'+str(x) for x in group],
                                     patient_id=['patient_'+str(i) for i in range(len(group))]), index=ids)
        metadata.loc[ids[1], 'patient_id'] = metadata.loc[ids[0], 'patient_id']
        metadata.loc[ids[3], 'lineage'] = 'Unknown'
        metadata.to_csv(directory/'metadata.csv')
        dependency = pd.DataFrame(values[:, 21:], index=ids, columns=spec['dependency_targets'])
        dependency.loc[ids[0], 'ARID1B'] = np.nan
        dependency.to_csv(directory/'dependency.csv')
        complete_values, complete_groups, complete_meta, flow, exclusions = prepare_inputs(
            directory, directory/'dependency.csv', spec, 'lineage')
        assert len(complete_meta) == 317 and complete_values.shape == (317, 23)
        assert np.isfinite(complete_values).all() and len(complete_groups) == 317
        assert not {ids[0], ids[1], ids[3]}.intersection(complete_meta.index)
        assert flow['shared_complete_case_rows'] == 317
        assert any(row['reason'] == 'same_patient_later_lexical_id' for row in exclusions)
    print('SELF_TEST_PASS: fixed outcome-independent module, standardization, FWL rho/R2, '
          'df, direction, row invariance, one constant target, degenerate member, insufficient n, '
          'duplicate-row refit, shared complete cases, patient selection, unknown lineage')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input')
    parser.add_argument('--dependency')
    parser.add_argument('--spec')
    parser.add_argument('--source-spec')
    parser.add_argument('--outdir')
    parser.add_argument('--lineage-column', default='lineage')
    parser.add_argument('--bootstrap', type=int)
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        self_test()
    else:
        if not all([args.input, args.spec, args.outdir]):
            parser.error('--input, --spec and --outdir are required')
        run(args)
