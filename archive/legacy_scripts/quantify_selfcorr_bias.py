# -*- coding: utf-8 -*-
"""Quantify the self-correlation inflation in module-score co-dependency analysis.

When an anchor gene is left inside the module whose mean it is correlated against,
1/k of the module score is the anchor itself. This inflates rho, and the inflation
should scale as ~1/k. We measure it empirically across three complexes and a range
of module sizes, and derive the size at which the bias becomes material.
"""
import numpy as np
import pandas as pd
from scipy import stats

MAT = '/data/s01011/ggn_project/02_targets/bulk_resources/depmap_CRISPRGeneEffect.csv'
OUT = '/data/s01011/arid1a-nsclc/results/tables/selfcorrelation_bias.csv'

COMPLEXES = {
    'SWI/SNF': {
        'core':  ['SMARCA4', 'SMARCB1', 'SMARCC1', 'SMARCC2', 'SMARCD1',
                  'SMARCD2', 'SMARCE1', 'ACTL6A', 'ACTB'],
        'cBAF':  ['ARID1A', 'ARID1B', 'DPF1', 'DPF2', 'DPF3', 'BCL7A',
                  'BCL7B', 'BCL7C', 'SS18', 'SS18L1'],
        'PBAF':  ['PBRM1', 'ARID2', 'BRD7', 'PHF10'],
        'ncBAF': ['BRD9', 'BICRA', 'BICRAL'],
    },
    'Cohesin': {
        'ring':    ['SMC1A', 'SMC3', 'RAD21'],
        'STAG':    ['STAG1', 'STAG2'],
        'loader':  ['NIPBL', 'MAU2'],
        'release': ['WAPL', 'PDS5A', 'PDS5B'],
    },
    'PRC2': {
        'core':   ['EZH2', 'EED', 'SUZ12', 'RBBP4', 'RBBP7'],
        'PRC2.1': ['MTF2', 'PHF1', 'PHF19', 'EPOP', 'LCOR'],
        'PRC2.2': ['AEBP2', 'JARID2'],
    },
}

df = pd.read_csv(MAT, index_col=0, low_memory=False)
df.columns = [c.split(' (')[0] for c in df.columns]
df = df.loc[:, ~df.columns.duplicated()]

rows = []
for cname, mods in COMPLEXES.items():
    for mod, members in mods.items():
        mem = [g for g in members if g in df.columns]
        for anchor in mem:                       # anchor is a member of this module
            others = [g for g in mem if g != anchor]
            if len(others) < 1:
                continue
            a = df[anchor]
            s_in = df[mem].mean(axis=1)
            s_ex = df[others].mean(axis=1)
            ok = a.notna() & s_in.notna() & s_ex.notna()
            if ok.sum() < 50:
                continue
            r_in = stats.spearmanr(a[ok], s_in[ok])[0]
            r_ex = stats.spearmanr(a[ok], s_ex[ok])[0]
            rows.append(dict(complex=cname, module=mod, anchor=anchor,
                             k_with_anchor=len(mem), rho_anchor_included=r_in,
                             rho_anchor_excluded=r_ex, inflation=r_in - r_ex,
                             n_lines=int(ok.sum())))

res = pd.DataFrame(rows).sort_values('k_with_anchor')
res.to_csv(OUT, index=False)

print(f'{"complex":9s}{"module":8s}{"anchor":9s}{"k":>3s}{"incl":>8s}{"excl":>8s}{"infl":>8s}')
print('-' * 55)
for _, r in res.iterrows():
    print(f'{r["complex"]:9s}{r["module"]:8s}{r["anchor"]:9s}{r["k_with_anchor"]:>3d}'
          f'{r["rho_anchor_included"]:>8.3f}{r["rho_anchor_excluded"]:>8.3f}'
          f'{r["inflation"]:>8.3f}')

print('\n=== inflation as a function of module size ===')
print(f'{"k":>3s}{"n":>4s}{"mean infl":>11s}{"median":>9s}{"max":>8s}')
for k, g in res.groupby('k_with_anchor'):
    print(f'{k:>3d}{len(g):>4d}{g["inflation"].mean():>11.3f}'
          f'{g["inflation"].median():>9.3f}{g["inflation"].max():>8.3f}')

r, p = stats.spearmanr(res['k_with_anchor'], res['inflation'])
print(f'\ninflation vs module size: Spearman rho = {r:+.3f}, P = {p:.2e}, '
      f'n = {len(res)} anchor-module pairs')
big = res[res['inflation'] > 0.1]
print(f'pairs with inflation > 0.10: {len(big)}/{len(res)} '
      f'(all have k <= {big["k_with_anchor"].max() if len(big) else "n/a"})')
print(f'\nwrote {OUT}')
