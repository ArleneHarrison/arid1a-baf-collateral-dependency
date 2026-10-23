"""Fixed-panel protein/RNA/CNV coupling, shared-row refit bootstrap.

Inputs are documented transformed matrices, samples in rows, exact gene symbols
in columns. No automatic log transform, outcome-based feature selection, missing
protein imputation, or silent replacement of non-evaluable panel members.
All computations are intended for the user's servers.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
from pathlib import Path
import time

import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata

LAYERS = ['protein_raw', 'rna_raw', 'protein_own_adjusted', 'protein_joint_adjusted']
UNKNOWN_METADATA_VALUES = {
    '', 'na', 'n/a', 'nan', 'none', 'null', 'unknown', 'not reported',
    'not available', 'not applicable', 'not specified', 'unspecified',
    'unreported', 'missing', '<na>', '-', '--',
}


def file_hash(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(2**20), b''):
            h.update(b)
    return h.hexdigest()


def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    # Fitted zero residuals can contain floating-point dust. Ranking that dust
    # would create a correlation in a layer with no measurable variation.
    if (len(x) < 2 or not np.isfinite(x).all() or not np.isfinite(y).all()
            or np.std(x) <= 1e-10 or np.std(y) <= 1e-10):
        return np.nan
    x = rankdata(x)
    y = rankdata(y)
    x -= x.mean()
    y -= y.mean()
    den = np.sqrt(x @ x * (y @ y))
    return float(x @ y / den) if den > 1e-15 else np.nan


def within_center(matrix, group):
    # One categorical factor, plus intercept. A single level is ordinary centering.
    _, labels = np.unique(group, return_inverse=True)
    sizes = np.bincount(labels)
    centered = matrix.copy()
    for j in range(matrix.shape[1]):
        means = np.bincount(labels, weights=matrix[:, j]) / sizes
        centered[:, j] -= means[labels]
        # Exact constant strata center to mathematical zero. Summation of an
        # odd number of identical floats can otherwise split residual ties.
        minima = np.full(len(sizes), np.inf)
        maxima = np.full(len(sizes), -np.inf)
        np.minimum.at(minima, labels, matrix[:, j])
        np.maximum.at(maxima, labels, matrix[:, j])
        centered[(minima == maxima)[labels], j] = 0.
    return centered, len(sizes)


def residual(y, x, n_groups):
    beta, _, rank, _ = np.linalg.lstsq(x, y, rcond=None)
    # BLAS kernels can differ at the tail of an odd-length array. Predict each
    # exact design row once so repeated bootstrap observations retain exact
    # residual ties. Keep distinct outcomes distinct; do not round observations.
    unique_x, inverse = np.unique(x, axis=0, return_inverse=True)
    prediction = (unique_x @ beta)[inverse]
    return y - prediction, len(y) - n_groups - rank


def pair_estimate(values, groups, min_n, min_df):
    # P_anchor,P_gene,RNA_anchor,RNA_gene,CNV_anchor,CNV_gene,numeric_covariates
    ok = np.isfinite(values).all(axis=1)
    z, g = values[ok], groups[ok]
    n = len(z)
    if n < min_n:
        return np.full(4, np.nan), n, np.full(4, np.nan)
    raw = [spearman(z[:, 0], z[:, 1]), spearman(z[:, 2], z[:, 3])]
    c, ngroups = within_center(z, g)
    extra = list(range(6, z.shape[1]))
    ra, dfa = residual(c[:, 0], c[:, [2, 4] + extra], ngroups)
    rb, dfb = residual(c[:, 1], c[:, [3, 5] + extra], ngroups)
    ja, dfja = residual(c[:, 0], c[:, [2, 3, 4, 5] + extra], ngroups)
    jb, dfjb = residual(c[:, 1], c[:, [2, 3, 4, 5] + extra], ngroups)
    own_df, joint_df = min(dfa, dfb), min(dfja, dfjb)
    own = spearman(ra, rb) if own_df >= min_df else np.nan
    joint = spearman(ja, jb) if joint_df >= min_df else np.nan
    # Raw correlations have no nuisance-regression residual degrees of freedom.
    # A non-evaluable layer must not invalidate a different evaluable layer.
    dfs = np.array([np.nan, np.nan, own_df, joint_df], float)
    return np.array(raw + [own, joint]), n, dfs


def finite_ci(values, total, min_fraction=.95):
    finite = values[np.isfinite(values)]
    if len(finite) == 0 or len(finite) < min_fraction * total:
        return np.nan, np.nan
    return tuple(np.quantile(finite, [.025, .975]))


def bh(pvalues):
    p = np.asarray(pvalues, float)
    q = np.full(len(p), np.nan)
    valid = np.flatnonzero(np.isfinite(p))
    order = valid[np.argsort(p[valid])]
    # Include all prespecified tests, even if one is not evaluable.
    if len(order):
        tmp = p[order] * len(p) / np.arange(1, len(order) + 1)
        q[order] = np.minimum.accumulate(tmp[::-1])[::-1].clip(0, 1)
    return q


def contrast(estimates, genes, contact, background, layer):
    if len(contact) < 1 or len(background) < 1:
        return np.nan
    if any(g not in genes for g in contact + background):
        return np.nan
    a = estimates[[genes.index(g) for g in contact], layer]
    b = estimates[[genes.index(g) for g in background], layer]
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        return np.nan
    return float(np.median(a) - np.median(b))


def layer_contrast(estimates, genes, contact, background, layer):
    if layer == -1:
        return (contrast(estimates, genes, contact, background, 0)
                - contrast(estimates, genes, contact, background, 1))
    return contrast(estimates, genes, contact, background, layer)


def layer_result_fields(estimates, dfs):
    fields = dict(zip(LAYERS, estimates))
    for j, layer in enumerate(LAYERS):
        fields[layer + '_model_residual_df'] = dfs[j]
        fields[layer + '_point_evaluable'] = bool(np.isfinite(estimates[j]))
    return fields


def unknown_metadata(values):
    normalized = (values.fillna('').astype(str).str.strip().str.lower()
                  .str.replace(r'[_\s]+', ' ', regex=True))
    return normalized.isin(UNKNOWN_METADATA_VALUES)


def load_inputs(directory, spec, categorical, numeric, dedup):
    assays = {}
    for name in ['protein', 'rna', 'cnv']:
        f = pd.read_csv(directory / f'{name}.csv', index_col=0)
        if not f.index.is_unique or not f.columns.is_unique:
            raise ValueError(f'Duplicate sample or gene IDs in {name}')
        f.index = f.index.astype(str)
        assays[name] = f.apply(pd.to_numeric, errors='raise')
    meta = pd.read_csv(directory / 'metadata.csv', index_col=0)
    meta.index = meta.index.astype(str)
    if not meta.index.is_unique:
        raise ValueError('Metadata samples must be unique')
    ids = sorted(set(meta.index).intersection(*(set(x.index) for x in assays.values())))
    flow = {'metadata': len(meta), **{k: len(v) for k, v in assays.items()},
            'matched_before_covariate_filter': len(ids)}
    meta = meta.loc[ids].copy()
    required = ([categorical] if categorical else []) + numeric
    for c in required:
        if c not in meta:
            raise ValueError(f'Requested covariate missing: {c}')
    for c in numeric:
        meta[c] = pd.to_numeric(meta[c], errors='coerce')
    exclusions = []
    covariate_issues = {}
    for c in numeric:
        covariate_issues['nonfinite_prespecified_numeric_covariate:' + c] = (
            ~np.isfinite(meta[c].to_numpy(float)))
    if categorical:
        covariate_issues['missing_or_unknown_prespecified_categorical_covariate:' + categorical] = (
            unknown_metadata(meta[categorical]).to_numpy())
        meta[categorical] = meta[categorical].astype(str).str.strip()
    if covariate_issues:
        missing = np.logical_or.reduce(list(covariate_issues.values()))
        for j in np.flatnonzero(missing):
            reasons = [reason for reason, mask in covariate_issues.items() if mask[j]]
            exclusions.append({'sample_id': meta.index[j], 'reason': ';'.join(reasons)})
        flow['covariate_exclusion_counts'] = {
            reason: int(mask.sum()) for reason, mask in covariate_issues.items()}
        meta = meta.loc[~missing]
    if 'patient_id' not in meta:
        raise ValueError('Explicit patient_id mapping is required')
    pid = meta.patient_id.fillna('').astype(str).str.strip()
    unknown = unknown_metadata(pid)
    pid.loc[unknown] = ['unknown_' + x for x in meta.index[unknown]]
    meta['analysis_patient_id'] = pid
    if dedup:
        dup = meta.analysis_patient_id.duplicated(keep='first')
        for i in meta.index[dup]:
            exclusions.append({'sample_id': i, 'reason': 'same_patient_later_lexical_id'})
        meta = meta.loc[~dup]
    flow.update(final_rows=len(meta), unique_known_or_cell_specific_patients=meta.analysis_patient_id.nunique(),
                covariates=required, deduplicate_patient=dedup)
    ids = meta.index.tolist()
    for k in assays:
        assays[k] = assays[k].loc[ids]
    return assays, meta, flow, exclusions


def analyze(args):
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'r12_coupling_manifest.json').exists():
        raise FileExistsError('Completed output exists; use a fresh output directory')
    spec = json.loads(Path(args.spec).read_text(encoding='utf-8'))
    assays, meta, flow, exclusions = load_inputs(Path(args.input), spec, args.categorical,
                                               args.numeric, args.deduplicate_patient)
    if len(meta) < spec['minimum_pair_n']:
        raise ValueError('Insufficient matched samples')
    genes = list(dict.fromkeys(spec['contact_experimental'] + spec['primary_baf_background']
                  + spec['integrative_contact_sensitivity'] + spec['alternative_subunit_descriptive']
                  + spec['physically_present_noncontact'] + spec['unmeasured_background']))
    anchor = spec['anchor']
    for layer, f in assays.items():
        if anchor not in f:
            raise ValueError(f'Anchor not measured in {layer}; no analysis possible')
    groups = pd.factorize(meta[args.categorical], sort=True)[0] if args.categorical else np.zeros(len(meta), int)
    extra = meta[args.numeric].to_numpy(float) if args.numeric else np.empty((len(meta), 0))
    data, rows, point = {}, [], []
    for gene in genes:
        absent = [k for k, f in assays.items() if gene not in f]
        values = np.full((len(meta), 6 + len(args.numeric)), np.nan)
        if not absent:
            values = np.column_stack([assays['protein'][anchor], assays['protein'][gene],
                                      assays['rna'][anchor], assays['rna'][gene],
                                      assays['cnv'][anchor], assays['cnv'][gene], extra])
        data[gene] = values
        estimates, n, dfs = pair_estimate(values, groups, spec['minimum_pair_n'], spec['minimum_residual_df'])
        point.append(estimates)
        detection = (float(np.isfinite(assays['protein'][gene].to_numpy(float)).mean())
                     if gene in assays['protein'] else 0.)
        minimum_df = float(np.nanmin(dfs)) if np.isfinite(dfs).any() else np.nan
        rows.append(dict(gene=gene, n=n, minimum_residual_df=minimum_df,
                         protein_detection_fraction=detection, missing_assays=';'.join(absent),
                         evaluable=bool(np.isfinite(estimates).any()),
                         all_layers_evaluable=bool(np.isfinite(estimates).all()),
                         **layer_result_fields(estimates, dfs)))
    point = np.array(point)
    n_boot = args.bootstrap if args.bootstrap is not None else spec['bootstrap_replicates']
    if n_boot < 2:
        raise ValueError('At least two bootstrap replicates are required; use the frozen count for final inference')
    rng = np.random.default_rng(spec['seed'])
    strata = [np.flatnonzero(groups == k) for k in np.unique(groups)]
    boot = np.full((n_boot, len(genes), 4), np.nan)
    start = time.monotonic()
    active = [i for i, r in enumerate(rows) if r['evaluable']]
    for b in range(n_boot):
        idx = np.concatenate([rng.choice(s, len(s), replace=True) for s in strata])
        gb = groups[idx]
        for i in active:
            boot[b, i], _, _ = pair_estimate(data[genes[i]][idx], gb,
                                            spec['minimum_pair_n'], spec['minimum_residual_df'])
        if (b+1) % max(1, min(500, n_boot)) == 0:
            print(json.dumps({'bootstrap': b+1, 'of': n_boot, 'elapsed_seconds': round(time.monotonic()-start, 1)}), flush=True)
    for i, row in enumerate(rows):
        for j, layer in enumerate(LAYERS):
            lo, hi = (finite_ci(boot[:, i, j], n_boot, spec['minimum_valid_bootstrap_fraction'])
                      if np.isfinite(point[i, j]) else (np.nan, np.nan))
            row[layer + '_ci_low'], row[layer + '_ci_high'] = lo, hi
            row[layer + '_valid_bootstraps'] = int(np.isfinite(boot[:, i, j]).sum())
    pd.DataFrame(rows).to_csv(out / 'r12_coupling_by_gene.csv', index=False)
    contact, background = spec['contact_experimental'], spec['primary_baf_background']
    definitions = [
        ('primary_own_adjusted', contact, background, 2, True),
        ('primary_raw_layer_difference', contact, background, -1, True),
        ('sensitivity_joint_adjusted', contact, background, 3, False),
        ('sensitivity_ACTL6A_only', contact, ['ACTL6A'], 2, False),
        ('strict_physical_noncontact', contact, spec['physically_present_noncontact'], 2, False),
        ('sensitivity_integrative_SMARCC1', contact + spec['integrative_contact_sensitivity'], background, 2, False),
    ]
    high = {r['gene'] for r in rows if r['protein_detection_fraction'] >= .90}
    hc, hb = [g for g in contact if g in high], [g for g in background if g in high]
    definitions.append(('sensitivity_high_detection90', hc, hb, 2, False))
    definition_reasons = {}
    if len(hc) < 3 or len(hb) < 2:
        definition_reasons['sensitivity_high_detection90'] = (
            f'insufficient_high_detection_members:requires_contact>=3_and_background>=2;'
            f'observed_contact={len(hc)};observed_background={len(hb)}')
    for gene in contact + background:
        definitions.append(('sensitivity_leave_out_' + gene, [g for g in contact if g != gene],
                            [g for g in background if g != gene], 2, False))
    summaries = []
    for name, con, bg, layer, primary in definitions:
        def stat(m):
            return layer_contrast(m, genes, con, bg, layer)
        reason = definition_reasons.get(name, '')
        estimate = np.nan if reason else stat(point)
        samples = np.full(n_boot, np.nan) if reason else np.array([stat(x) for x in boot])
        valid = samples[np.isfinite(samples)]
        lo, hi = (finite_ci(samples, n_boot, spec['minimum_valid_bootstrap_fraction'])
                  if np.isfinite(estimate) else (np.nan, np.nan))
        evaluable = np.isfinite(estimate) and len(valid) >= spec['minimum_valid_bootstrap_fraction']*n_boot
        if not reason and not np.isfinite(estimate):
            reason = 'at_least_one_fixed_panel_member_or_required_layer_not_evaluable'
        elif not reason and not evaluable:
            reason = 'valid_bootstrap_fraction_below_predeclared_minimum'
        p = (1+np.sum(np.abs(valid-estimate) >= abs(estimate))) / (1+len(valid)) if evaluable and primary else np.nan
        summaries.append(dict(test=name, primary=primary, estimate=estimate, ci_low=lo, ci_high=hi,
                              centered_bootstrap_approx_p=p, n_valid_bootstrap=len(valid), n_bootstrap=n_boot,
                              evaluable=bool(evaluable), non_evaluable_reason=reason,
                              contact_genes=';'.join(con), background_genes=';'.join(bg)))
    summary = pd.DataFrame(summaries)
    summary['q_BH_two_primary_within_cohort'] = np.nan
    primary_ids = summary.index[summary.primary]
    summary.loc[primary_ids, 'q_BH_two_primary_within_cohort'] = bh(summary.loc[primary_ids, 'centered_bootstrap_approx_p'])
    summary.to_csv(out / 'r12_coupling_contrasts.csv', index=False)
    lolo_status = {'status': 'not_applicable_no_lineage_covariate',
                   'inference': 'Point estimates only; no confidence intervals or p-values'}
    if args.categorical and 'lineage' in args.categorical.lower():
        if len(strata) < 2:
            lolo_status['status'] = 'non_evaluable_fewer_than_two_lineages'
        else:
            lolo_genes, lolo_contrasts = [], []
            panel_ids = [genes.index(g) for g in contact + background]
            for code in np.unique(groups):
                keep = groups != code
                label = str(meta[args.categorical].iloc[np.flatnonzero(~keep)[0]])
                estimates_by_gene, n_by_gene, dfs_by_gene = [], [], []
                for gene in genes:
                    estimates, n, dfs = pair_estimate(data[gene][keep], groups[keep],
                                                     spec['minimum_pair_n'], spec['minimum_residual_df'])
                    estimates_by_gene.append(estimates)
                    n_by_gene.append(n)
                    dfs_by_gene.append(dfs)
                    lolo_genes.append(dict(omitted_lineage=label, omitted_rows=int((~keep).sum()),
                                           remaining_rows=int(keep.sum()), gene=gene, n=n,
                                           **layer_result_fields(estimates, dfs)))
                lolo_point = np.asarray(estimates_by_gene)
                lolo_dfs = np.asarray(dfs_by_gene)
                for test, layer in [('primary_own_adjusted', 2),
                                    ('primary_raw_layer_difference', -1)]:
                    estimate = layer_contrast(lolo_point, genes, contact, background, layer)
                    selected_dfs = lolo_dfs[panel_ids, layer] if layer >= 0 else np.array([])
                    minimum_df = (float(np.min(selected_dfs))
                                  if len(selected_dfs) and np.isfinite(selected_dfs).all() else np.nan)
                    lolo_contrasts.append(dict(
                        omitted_lineage=label, omitted_rows=int((~keep).sum()),
                        remaining_rows=int(keep.sum()), test=test, estimate=estimate,
                        evaluable=bool(np.isfinite(estimate)),
                        non_evaluable_reason='' if np.isfinite(estimate) else 'fixed_panel_not_evaluable',
                        minimum_pair_n=int(np.min(np.asarray(n_by_gene)[panel_ids])),
                        minimum_model_residual_df=minimum_df,
                        contact_genes=';'.join(contact), background_genes=';'.join(background),
                        inference='descriptive_point_estimate_only_no_CI_or_P'))
            pd.DataFrame(lolo_genes).to_csv(out / 'r12_leave_one_lineage_out_by_gene.csv', index=False)
            pd.DataFrame(lolo_contrasts).to_csv(out / 'r12_leave_one_lineage_out_contrasts.csv', index=False)
            lolo_status.update(status='completed_descriptive_point_estimates', n_lineages=len(strata))
    meta.to_csv(out / 'r12_analysis_samples.csv', index_label='sample_id')
    pd.DataFrame(exclusions, columns=['sample_id','reason']).to_csv(out / 'r12_sample_exclusions.csv', index=False)
    np.savez_compressed(out / 'r12_shared_refit_bootstrap.npz', values=boot, genes=np.array(genes), layers=np.array(LAYERS))
    manifest = dict(cohort=args.cohort, flow=flow, specification=spec,
                    spec_sha256=file_hash(args.spec),
                    analysis_script_sha256=file_hash(__file__),
                    inputs={name: file_hash(Path(args.input) / f'{name}.csv') for name in ['protein','rna','cnv','metadata']},
                    versions={'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__},
                    bootstrap=n_boot, elapsed_seconds=time.monotonic()-start,
                    leave_one_lineage_out=lolo_status,
                    layer_df_note='Raw layers have no nuisance-model residual df; own/joint df are assessed separately.',
                    minimum_original_scale_std_for_spearman=1e-10,
                    numerical_revision='v3: unique-design prediction and exact-zero centering of constant strata preserve mathematical ties; no rounding, model changes, or result-selected thresholds',
                    inference_scope='Fixed assayed BAF panel; not an isolated physical-contact effect or protein stability proof',
                    external_multiplicity_note='For clinical replication, collect both primary tests from LUAD and LUSC and apply BH across the four planned tests.')
    (out / 'r12_coupling_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(summary.loc[summary.primary].to_string(index=False), flush=True)


def self_test():
    # Synthetic tests run on server only. FWL must equal explicit dummy OLS.
    rng = np.random.default_rng(4901)
    g = np.repeat(np.arange(4), 60)
    z = rng.normal(size=(len(g), 7))
    z[:, 0] += .8*z[:, 2] + .3*z[:, 4] + .4*g
    z[:, 1] += .7*z[:, 3] - .2*z[:, 5] + .2*g
    z[2, 4] = np.nan
    est, n, dfs = pair_estimate(z, g, 40, 20)
    ok = np.isfinite(z).all(axis=1)
    a = z[ok]
    dummies = np.eye(4)[g[ok]]
    x = np.column_stack([dummies, a[:, [2, 4, 6]]])
    xa = a[:, 0] - x @ np.linalg.lstsq(x, a[:, 0], rcond=None)[0]
    y = np.column_stack([dummies, a[:, [3, 5, 6]]])
    ya = a[:, 1] - y @ np.linalg.lstsq(y, a[:, 1], rcond=None)[0]
    assert np.isclose(est[2], spearman(xa, ya), atol=1e-12)
    assert n == 239 and dfs[2] == 232 and dfs[3] == 230
    assert np.isnan(dfs[:2]).all()
    perm = rng.permutation(len(z))
    swapped, _, _ = pair_estimate(z[perm], g[perm], 40, 20)
    assert np.allclose(est, swapped)
    assert np.allclose(bh([.01,.04]), [.02,.04])
    assert np.allclose(bh([.01, np.nan])[0], .02)
    assert np.isnan(pair_estimate(z[:10], g[:10], 40, 20)[0]).all()
    # Insufficient joint df must not erase an otherwise valid own-adjusted layer.
    limited, _, limited_df = pair_estimate(z, g, 40, 231)
    assert np.isfinite(limited[:3]).all() and np.isnan(limited[3])
    assert limited_df[2] == 232 and limited_df[3] == 230
    # Constant RNA invalidates only its raw correlation, not protein coupling.
    constant_rna = z.copy()
    constant_rna[:, 2] = 1.
    constant_est, _, _ = pair_estimate(constant_rna, g, 40, 20)
    assert np.isnan(constant_est[1]) and np.isfinite(constant_est[[0, 2, 3]]).all()
    # Anchor protein equal to partner RNA is exactly removed by joint adjustment;
    # its floating-point residual must not be ranked into a spurious correlation.
    joint_degenerate = z.copy()
    joint_degenerate[:, 0] = joint_degenerate[:, 3]
    degenerate_est, _, _ = pair_estimate(joint_degenerate, g, 40, 20)
    assert np.isfinite(degenerate_est[2]) and np.isnan(degenerate_est[3])
    assert np.isnan(spearman(np.arange(10)*1e-12, np.arange(10)))
    # Odd-length bootstrap draws exposed a real BLAS tail-kernel tie-splitting
    # issue. Exact duplicate inputs must always have exact equal residuals.
    base_x = rng.normal(size=(67, 5))
    base_y = rng.normal(size=67)
    draw = rng.integers(0, 67, 107)
    repeated_resid, _ = residual(base_y[draw], base_x[draw], 1)
    for source in np.unique(draw):
        assert np.unique(repeated_resid[draw == source]).size == 1
    mixed_y = np.array([1., 2., 3., 4.])
    same_x = np.ones((4, 1))
    distinct_resid, _ = residual(mixed_y, same_x, 1)
    assert np.allclose(np.diff(distinct_resid), 1.)
    constant_blocks = np.array([[.1, .920711640924486]] * 3 + [[.3, .7]] * 5)
    stable_center, ng = within_center(constant_blocks, np.array([0] * 3 + [1] * 5))
    assert ng == 2 and np.array_equal(stable_center, np.zeros_like(stable_center))
    assert unknown_metadata(pd.Series(['Unknown', 'not_reported', 'N/A', None, 'lung'])).tolist() == [True]*4 + [False]
    assert all(np.isnan(x) for x in finite_ci(np.array([]), 0))
    print('SELF_TEST_PASS: FWL=dummy OLS, per-layer df, row invariance, pair masks, BH, '
          'insufficient n, constant RNA, degenerate joint residuals, metadata placeholders')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--input')
    p.add_argument('--outdir')
    p.add_argument('--spec')
    p.add_argument('--cohort', default='unspecified')
    p.add_argument('--categorical', default=None)
    p.add_argument('--numeric', nargs='*', default=[])
    p.add_argument('--deduplicate-patient', action='store_true')
    p.add_argument('--bootstrap', type=int)
    p.add_argument('--self-test', action='store_true')
    args = p.parse_args()
    if args.self_test:
        self_test()
    else:
        if not all([args.input,args.outdir,args.spec]):
            p.error('--input, --outdir and --spec are required')
        analyze(args)
