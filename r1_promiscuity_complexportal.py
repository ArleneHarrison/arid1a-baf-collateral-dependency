# -*- coding: utf-8 -*-
"""Does subunit promiscuity explain the module refinement?

The manuscript reports this as unresolved: the only promiscuity proxy available
at the time was built from the dependency data itself, which conflates
promiscuity with dependency-signal strength. EBI Complex Portal supplies an
external, independent count of how many annotated complexes each subunit belongs
to, so the test can be made non-circular.

Hypothesis: subunits whose removal IMPROVES the module (positive leave-one-out
delta) are the promiscuous ones, because their dependency profile is driven by
the other complexes they sit in.
"""
import re

import pandas as pd
from scipy import stats

import os
D = os.path.dirname(os.path.abspath(__file__))
T = (r'D:\桌面\_按项目整理\02_生信数据与医学AI项目\07_医学AI_深度学习_多项目孵化'
     r'\医学AI深度学习_导师课题项目群\范教授课题_NSCLC-SEER_分子对接多项目\纯生信'
     r'\ARID1A_BAF_最终项目\results_tables')

# ---- 1. UniProt -> gene symbol, via the Complex Portal expanded column -------
cp = pd.read_csv(D + r'\cp.tsv', sep='\t', low_memory=False)
mol_col = [c for c in cp.columns if 'Identifiers' in c][0]

# count, for each UniProt accession, how many distinct complexes contain it
from collections import Counter
cnt = Counter()
for s in cp[mol_col].dropna():
    accs = {a.split('(')[0].split('-')[0] for a in str(s).split('|') if a}
    for a in accs:
        cnt[a] += 1

# UniProt accessions for the subunits under test (canonical, reviewed)
ACC = {
    'SMARCB1': 'Q12824', 'SMARCC1': 'Q92922', 'SMARCC2': 'Q8TAQ2',
    'SMARCD1': 'Q96GM5', 'SMARCD2': 'Q92925', 'SMARCE1': 'Q969G3',
    'ACTL6A': 'O96019', 'ACTB': 'P60709',
    'EED': 'O75530', 'SUZ12': 'Q15022', 'RBBP4': 'Q09028',
    'RBBP7': 'Q16576', 'EZH2': 'Q15910',
}

lo = pd.read_csv(T + r'\contact_vs_promiscuity.csv', comment='#')
lo = lo.dropna(subset=['loo_delta'])
lo['n_complexes'] = lo.subunit.map(lambda g: cnt.get(ACC.get(g, ''), None))

print('=' * 72)
print('每个亚基在 Complex Portal 中归属的复合体数（外部、独立于依赖数据）')
print('=' * 72)
w = lo.drop_duplicates('subunit')[['subunit', 'n_complexes']].sort_values(
    'n_complexes', ascending=False)
for _, r in w.iterrows():
    print(f'  {r.subunit:10s} {r.n_complexes if pd.notna(r.n_complexes) else "未收录":>8}')

d = lo.dropna(subset=['n_complexes'])
print()
print('=' * 72)
print('检验：复合体归属数 是否预测 留一法增益')
print('=' * 72)
print('（留一法增益 > 0 = 去掉该亚基后模块变好 = 它在稀释信号）')
rho, p = stats.spearmanr(d.n_complexes, d.loo_delta)
print(f'  全部 {len(d)} 个亚基:  Spearman rho = {rho:+.3f}, P = {p:.4f}')

for cx in ('SWI/SNF', 'PRC2'):
    s = d[d.complex == cx]
    if len(s) >= 4:
        r2, p2 = stats.spearmanr(s.n_complexes, s.loo_delta)
        print(f'  {cx:8s} ({len(s)} 个):      Spearman rho = {r2:+.3f}, P = {p2:.4f}')

print()
print('稀释型（delta>0，去掉变好） vs 承载型（delta<0）的复合体归属数:')
dil = d[d.loo_delta > 0]
car = d[d.loo_delta < 0]
print(f'  稀释型 n={len(dil):2d}  中位复合体数 = {dil.n_complexes.median():.1f}  '
      f'({sorted(dil.subunit.unique())})')
print(f'  承载型 n={len(car):2d}  中位复合体数 = {car.n_complexes.median():.1f}  '
      f'({sorted(car.subunit.unique())})')
if len(dil) >= 3 and len(car) >= 3:
    u, pu = stats.mannwhitneyu(dil.n_complexes, car.n_complexes)
    print(f'  Mann-Whitney P = {pu:.4f}')
