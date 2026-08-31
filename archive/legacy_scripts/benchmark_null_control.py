# -*- coding: utf-8 -*-
"""Is the 14/14 discrimination benchmark informative, or near-tautological?

Concern: genes in the same complex co-vary in essentiality as a matter of course
(co-essentiality networks are routinely validated by their recovery of protein
complexes). If so, "anchor correlates with its own module more than with a
mutually exclusive one" may restate that known property rather than test anything.

Three controls:
  A. random genes  -- do non-members also get confidently assigned to a module?
  B. decoy anchors -- genes from OTHER complexes, tested against SWI/SNF modules
  C. effect-size separation -- how far do true anchors sit above the null?
"""
import numpy as np
import pandas as pd
from scipy import stats

MAT = '/data/s01011/ggn_project/02_targets/bulk_resources/depmap_CRISPRGeneEffect.csv'
SEED = 20260801

MODULES = {
    'core':  ['SMARCB1', 'SMARCC1', 'SMARCC2', 'SMARCD1', 'SMARCD2',
              'SMARCE1', 'ACTL6A', 'ACTB'],
    'cBAF':  ['ARID1A', 'ARID1B', 'DPF1', 'DPF2', 'DPF3', 'BCL7A',
              'BCL7B', 'BCL7C', 'SS18', 'SS18L1'],
    'PBAF':  ['PBRM1', 'ARID2', 'BRD7', 'PHF10'],
    'ncBAF': ['BRD9', 'BICRA', 'BICRAL'],
}
EXCLUSIVE = ['cBAF', 'PBAF', 'ncBAF']
TRUE_ANCHORS = {'ARID1A': 'cBAF', 'ARID1B': 'cBAF', 'DPF2': 'cBAF', 'SS18': 'cBAF',
                'PBRM1': 'PBAF', 'ARID2': 'PBAF', 'BRD7': 'PBAF',
                'BRD9': 'ncBAF', 'BICRA': 'ncBAF'}
# decoys: members of unrelated, well-defined complexes
DECOYS = ['EZH2', 'SUZ12', 'EED',            # PRC2
          'SMC1A', 'SMC3', 'RAD21',          # cohesin
          'RPL3', 'RPS3',                    # ribosome
          'PSMA1', 'PSMB6',                  # proteasome
          'POLR2A', 'POLR2B',                # RNA Pol II
          'MED1', 'MED12']                   # Mediator

df = pd.read_csv(MAT, index_col=0, low_memory=False)
df.columns = [c.split(' (')[0] for c in df.columns]
df = df.loc[:, ~df.columns.duplicated()]
ALLMOD = set(sum(MODULES.values(), []))


def profile(anchor):
    out = {}
    for m, mem in MODULES.items():
        mm = [g for g in mem if g in df.columns and g != anchor]
        if not mm:
            out[m] = np.nan
            continue
        s = df[mm].mean(axis=1)
        a = df[anchor]
        ok = a.notna() & s.notna()
        out[m] = stats.spearmanr(a[ok], s[ok])[0] if ok.sum() >= 30 else np.nan
    return out


def margin(prof, home):
    """own exclusive module minus the best competing exclusive module"""
    others = [prof[m] for m in EXCLUSIVE if m != home and not np.isnan(prof[m])]
    return prof[home] - max(others) if others else np.nan


print('=== A. true anchors: own-module margin over best competitor ===')
true_margins = []
for a, home in TRUE_ANCHORS.items():
    p = profile(a)
    mg = margin(p, home)
    true_margins.append(mg)
    print(f'  {a:8s} {home:6s} own={p[home]:+.3f}  margin={mg:+.3f}')
true_margins = np.array(true_margins)
print(f'  median margin = {np.median(true_margins):+.3f}   min = {true_margins.min():+.3f}')

print('\n=== B. random genes: best exclusive-module rho, and would they be "assigned"? ===')
rng = np.random.RandomState(SEED)
pool = [g for g in df.columns if g not in ALLMOD and df[g].notna().sum() > 500]
rand = rng.choice(pool, 300, replace=False)
best_r, assigned = [], 0
for g in rand:
    p = profile(g)
    vals = [p[m] for m in EXCLUSIVE if not np.isnan(p[m])]
    if not vals:
        continue
    b = max(vals)
    best_r.append(b)
    if b > 0.1:
        assigned += 1
best_r = np.array(best_r)
print(f'  n = {len(best_r)} random genes')
print(f'  best exclusive-module rho: median {np.median(best_r):+.3f}, '
      f'95th pct {np.percentile(best_r, 95):+.3f}, max {best_r.max():+.3f}')
print(f'  would pass the rho>0.1 assignment threshold: {assigned}/{len(best_r)} '
      f'({100*assigned/len(best_r):.1f}%)')

print('\n=== C. decoy anchors from other complexes, scored against SWI/SNF modules ===')
print(f'  {"gene":8s}{"core":>9s}{"cBAF":>9s}{"PBAF":>9s}{"ncBAF":>9s}  best exclusive')
dec_best = []
for g in DECOYS:
    if g not in df.columns:
        continue
    p = profile(g)
    vals = {m: p[m] for m in EXCLUSIVE if not np.isnan(p[m])}
    b = max(vals, key=vals.get)
    dec_best.append(vals[b])
    print(f'  {g:8s}' + ''.join(f'{p[m]:>9.3f}' for m in ['core'] + EXCLUSIVE)
          + f'   {b} {vals[b]:+.3f}')
dec_best = np.array(dec_best)

print('\n=== verdict ===')
print(f'  true anchors  : own-module rho median {np.median([profile(a)[h] for a, h in TRUE_ANCHORS.items()]):+.3f}')
print(f'  random genes  : best exclusive rho median {np.median(best_r):+.3f}')
print(f'  decoy anchors : best exclusive rho median {np.median(dec_best):+.3f}')
u, p = stats.mannwhitneyu([profile(a)[h] for a, h in TRUE_ANCHORS.items()], best_r,
                          alternative='greater')
print(f'  true vs random: Mann-Whitney P = {p:.2e}')
print(f'  fraction of random genes exceeding the WEAKEST true anchor '
      f'({min(profile(a)[h] for a, h in TRUE_ANCHORS.items()):.3f}): '
      f'{(best_r > min(profile(a)[h] for a, h in TRUE_ANCHORS.items())).mean()*100:.1f}%')
