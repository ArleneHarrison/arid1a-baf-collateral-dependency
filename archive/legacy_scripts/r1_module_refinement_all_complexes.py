# -*- coding: utf-8 -*-
"""How general is module refinement? Complex Portal x DepMap, every complex.

The manuscript demonstrates on two complexes and three anchors that a curated
module can be improved by dropping members that dilute the co-dependency signal,
and lists "generalisation to complexes of other architectures remains untested"
as a limitation. This runs the identical procedure over every human complex in
the EBI Complex Portal that DepMap covers well enough to test.

Procedure, matching the manuscript exactly:
  module score  = mean Chronos gene effect over module members, ANCHOR EXCLUDED
  rho           = Spearman(anchor gene-effect profile, module score) across lines
  leave-one-out = rho recomputed with one further member dropped
A refinement is a drop that raises rho.

Null: the same leave-one-out gain computed for size-matched modules of random
genes, so that "some drop always helps by chance" is quantified rather than
assumed.
"""
import os
import sys
import json

import numpy as np
import pandas as pd
from scipy import stats

BULK = '/data/s01011/ggn_project/02_targets/bulk_resources'
HOME = os.path.expanduser('~')
OUT = os.path.join(HOME, 'arid1a_module_generalisation')
os.makedirs(OUT, exist_ok=True)

MIN_MEMBERS = 6      # need >= 6 so that anchor + module + one dropped is meaningful
MIN_LINES = 200      # cell lines with data for both anchor and module
RNG = np.random.default_rng(0)


def log(*a):
    print(*a, flush=True)


# ---------------------------------------------------------------- dependency
log('reading DepMap Chronos matrix ...')
dep = pd.read_csv(os.path.join(BULK, 'depmap_CRISPRGeneEffect.csv'), index_col=0)
dep.columns = [c.split(' (')[0] for c in dep.columns]
dep = dep.loc[:, ~dep.columns.duplicated()]
log(f'  {dep.shape[0]} cell lines x {dep.shape[1]} genes')

# ---------------------------------------------------------------- complexes
u2s = pd.read_csv(os.path.join(HOME, 'up2sym.tsv'), sep='\t')
u2s = dict(zip(u2s['Entry'], u2s['Gene Names (primary)']))
log(f'  UniProt->symbol entries: {len(u2s)}')

cp = pd.read_csv(os.path.join(HOME, 'cp.tsv'), sep='\t', low_memory=False)
mol_col = [c for c in cp.columns if 'Identifiers' in c][0]
name_col = 'Recommended name'
ac_col = cp.columns[0]

complexes = {}
for _, row in cp.iterrows():
    s = row[mol_col]
    if not isinstance(s, str):
        continue
    accs = {a.split('(')[0].split('-')[0] for a in s.split('|') if a}
    syms = {u2s[a] for a in accs if a in u2s and isinstance(u2s.get(a), str)}
    syms = {g for g in syms if g in dep.columns}
    if len(syms) >= MIN_MEMBERS:
        complexes[row[ac_col]] = (row[name_col], sorted(syms))
log(f'  complexes with >= {MIN_MEMBERS} DepMap-covered subunits: {len(complexes)}')


def rho_of(anchor_vec, members, frame):
    """Spearman of the anchor profile against the mean of `members`."""
    sub = frame[members]
    ok = anchor_vec.notna() & sub.notna().all(axis=1)
    if ok.sum() < MIN_LINES:
        return np.nan, int(ok.sum())
    return stats.spearmanr(anchor_vec[ok], sub[ok].mean(axis=1))[0], int(ok.sum())


rows = []
for i, (ac, (name, members)) in enumerate(complexes.items(), 1):
    if i % 25 == 0:
        log(f'  {i}/{len(complexes)} complexes')
    for anchor in members:
        module = [g for g in members if g != anchor]   # anchor ALWAYS excluded
        if len(module) < MIN_MEMBERS - 1:
            continue
        av = dep[anchor]
        base, n = rho_of(av, module, dep)
        if not np.isfinite(base):
            continue
        best_gain, best_drop = 0.0, None
        for m in module:
            rest = [g for g in module if g != m]
            r, _ = rho_of(av, rest, dep)
            if np.isfinite(r) and (r - base) > best_gain:
                best_gain, best_drop = r - base, m
        rows.append(dict(complex_ac=ac, complex_name=name, anchor=anchor,
                         k_module=len(module), n_lines=n, rho_full=base,
                         best_gain=best_gain, dropped=best_drop))

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, 'module_refinement_all_complexes.csv'), index=False)
log(f'\ntested {len(res)} (complex, anchor) pairs over '
    f'{res.complex_ac.nunique()} complexes')

# ------------------------------------------------------------------- null
log('\nbuilding size-matched random-gene null ...')
null_gain = []
sizes = res.k_module.value_counts()
pool = [g for g in dep.columns if dep[g].notna().sum() > MIN_LINES]
for k, cnt in sizes.items():
    for _ in range(min(cnt, 200)):
        pick = RNG.choice(pool, size=k + 1, replace=False)
        av = dep[pick[0]]
        module = list(pick[1:])
        base, n = rho_of(av, module, dep)
        if not np.isfinite(base):
            continue
        bg = 0.0
        for m in module:
            r, _ = rho_of(av, [g for g in module if g != m], dep)
            if np.isfinite(r) and (r - base) > bg:
                bg = r - base
        null_gain.append(dict(k_module=k, rho_full=base, best_gain=bg))
nul = pd.DataFrame(null_gain)
nul.to_csv(os.path.join(OUT, 'module_refinement_null.csv'), index=False)

# ---------------------------------------------------------------- summary
real = res[res.rho_full > 0.2]          # anchors with a dependency to refine
summary = dict(
    n_pairs=int(len(res)), n_complexes=int(res.complex_ac.nunique()),
    n_pairs_rho_gt_0_2=int(len(real)),
    frac_any_gain=float((res.best_gain > 0).mean()),
    frac_gain_gt_0_05=float((res.best_gain > 0.05).mean()),
    median_gain=float(res.best_gain.median()),
    real_median_gain=float(real.best_gain.median()) if len(real) else None,
    real_frac_gain_gt_0_05=float((real.best_gain > 0.05).mean()) if len(real) else None,
    null_median_gain=float(nul.best_gain.median()),
    null_frac_gain_gt_0_05=float((nul.best_gain > 0.05).mean()),
)
if len(real):
    summary['mannwhitney_real_vs_null_P'] = float(
        stats.mannwhitneyu(real.best_gain, nul.best_gain,
                           alternative='greater')[1])
with open(os.path.join(OUT, 'summary.json'), 'w') as f:
    json.dump(summary, f, indent=2)

log('\n' + '=' * 66)
for k, v in summary.items():
    log(f'  {k:28s} {v}')
log('=' * 66)
log('\ntop refinements among anchors with rho_full > 0.2:')
log(real.nlargest(20, 'best_gain')[
    ['complex_name', 'anchor', 'k_module', 'rho_full', 'best_gain', 'dropped']
].to_string(index=False))
log(f'\nwritten to {OUT}')
