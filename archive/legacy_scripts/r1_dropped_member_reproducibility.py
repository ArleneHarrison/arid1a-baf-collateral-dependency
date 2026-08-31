# -*- coding: utf-8 -*-
"""Is the dilutive subunit a property of the complex, or the argmax of noise?

The first systematic run showed that the SIZE of the leave-one-out gain is not
larger for real complexes than for size-matched random gene sets: taking the best
of k drops is an optimisation, so it gains regardless of biology. That kills
"gain magnitude" as evidence.

What it does not test is whether the SAME member is identified as dilutive when
the data are resampled. Noise-driven argmax should not replicate; a genuinely
dilutive member should. This splits the cell lines in half and asks how often the
two halves nominate the same member, against the same size-matched null.
"""
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

BULK = '/data/s01011/ggn_project/02_targets/bulk_resources'
HOME = os.path.expanduser('~')
OUT = '/data/s01011/arid1a_module_generalisation'
N_SPLIT = 20          # independent half-splits per (complex, anchor)
MIN_HALF = 150        # cell lines required in each half
RNG = np.random.default_rng(7)


def log(*a):
    print(*a, flush=True)


log('reading Chronos matrix ...')
dep = pd.read_csv(os.path.join(BULK, 'depmap_CRISPRGeneEffect.csv'), index_col=0)
dep.columns = [c.split(' (')[0] for c in dep.columns]
dep = dep.loc[:, ~dep.columns.duplicated()]
lines = dep.index.to_numpy()
log(f'  {dep.shape[0]} lines x {dep.shape[1]} genes')

real = pd.read_csv(os.path.join(OUT, 'module_refinement_all_complexes.csv'))

# rebuild membership so each (complex, anchor) row can be re-tested
u2s = pd.read_csv(os.path.join(HOME, 'up2sym.tsv'), sep='\t')
u2s = dict(zip(u2s['Entry'], u2s['Gene Names (primary)']))
cp = pd.read_csv(os.path.join(HOME, 'cp.tsv'), sep='\t', low_memory=False)
mol_col = [c for c in cp.columns if 'Identifiers' in c][0]
members_of = {}
for _, row in cp.iterrows():
    s = row[mol_col]
    if not isinstance(s, str):
        continue
    accs = {a.split('(')[0].split('-')[0] for a in s.split('|') if a}
    syms = sorted({u2s[a] for a in accs
                   if a in u2s and isinstance(u2s.get(a), str)} & set(dep.columns))
    members_of[row[cp.columns[0]]] = syms


def best_drop(anchor_vec, module, idx):
    """Member whose removal most improves the anchor-module correlation."""
    sub = dep.loc[idx, module]
    av = anchor_vec.loc[idx]
    ok = av.notna() & sub.notna().all(axis=1)
    if ok.sum() < MIN_HALF:
        return None
    a, s = av[ok], sub[ok]
    base = stats.spearmanr(a, s.mean(axis=1))[0]
    gains = {}
    for m in module:
        rest = [g for g in module if g != m]
        gains[m] = stats.spearmanr(a, s[rest].mean(axis=1))[0] - base
    return max(gains, key=gains.get)


def concordance(anchor, module):
    """Fraction of half-splits nominating the same member in both halves."""
    hits = tot = 0
    for _ in range(N_SPLIT):
        perm = RNG.permutation(lines)
        h1, h2 = perm[:len(perm) // 2], perm[len(perm) // 2:]
        d1 = best_drop(dep[anchor], module, h1)
        d2 = best_drop(dep[anchor], module, h2)
        if d1 is None or d2 is None:
            continue
        tot += 1
        hits += (d1 == d2)
    return (hits / tot if tot else np.nan), tot


cand = real[(real.rho_full > 0.2) & (real.k_module >= 5)
            & (real.k_module <= 20)].copy()
# one anchor per complex keeps complex variants from dominating the count
cand = cand.sort_values('rho_full', ascending=False).drop_duplicates('complex_ac')
log(f'testing {len(cand)} (complex, anchor) pairs')

rows = []
for i, (_, x) in enumerate(cand.iterrows(), 1):
    if i % 20 == 0:
        log(f'  {i}/{len(cand)}')
    mem = members_of.get(x.complex_ac, [])
    module = [g for g in mem if g != x.anchor]
    if len(module) != x.k_module:
        continue
    c, tot = concordance(x.anchor, module)
    if np.isfinite(c):
        rows.append(dict(complex_name=x.complex_name, anchor=x.anchor,
                         k_module=x.k_module, rho_full=x.rho_full,
                         dropped_full=x.dropped, concordance=c, n_splits=tot,
                         chance=1.0 / x.k_module))

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, 'dropped_member_reproducibility.csv'), index=False)

log('\nnull: size-matched random gene sets ...')
nrows = []
pool = [g for g in dep.columns if dep[g].notna().sum() > 400]
for k, cnt in res.k_module.value_counts().items():
    for _ in range(min(cnt, 60)):
        pick = RNG.choice(pool, size=k + 1, replace=False)
        c, tot = concordance(pick[0], list(pick[1:]))
        if np.isfinite(c):
            nrows.append(dict(k_module=k, concordance=c, chance=1.0 / k))
nul = pd.DataFrame(nrows)
nul.to_csv(os.path.join(OUT, 'dropped_member_reproducibility_null.csv'),
           index=False)

summ = dict(
    n_real=int(len(res)), n_null=int(len(nul)),
    real_median_concordance=float(res.concordance.median()),
    null_median_concordance=float(nul.concordance.median()),
    real_median_chance=float(res.chance.median()),
    real_frac_above_half=float((res.concordance > 0.5).mean()),
    null_frac_above_half=float((nul.concordance > 0.5).mean()),
    mannwhitney_P=float(stats.mannwhitneyu(res.concordance, nul.concordance,
                                           alternative='greater')[1]),
)
with open(os.path.join(OUT, 'reproducibility_summary.json'), 'w') as f:
    json.dump(summ, f, indent=2)
log('\n' + '=' * 62)
for k, v in summ.items():
    log(f'  {k:28s} {v}')
log('=' * 62)
log('\n最可重复的 15 个:')
log(res.nlargest(15, 'concordance')[
    ['complex_name', 'anchor', 'k_module', 'dropped_full',
     'concordance', 'chance']].to_string(index=False))
