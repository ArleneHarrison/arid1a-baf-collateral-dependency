# -*- coding: utf-8 -*-
"""Genetic knockout as the surrogate for an ideal inhibitor.

GDSC cannot test ARID1A-genotype selectivity of the SWI/SNF probes because the
probes are inactive there. CRISPR knockout of the same targets is the closest
available substitute: it asks what a fully effective inhibitor of that subunit
would do. Two questions:
  1. do the drug targets (SMARCA4, SMARCA2) show genotype-conditional dependency?
  2. where do they rank against the other subunits on continuous co-dependency?
"""
import pandas as pd

T = (r'D:\桌面\_按项目整理\02_生信数据与医学AI项目\07_医学AI_深度学习_多项目孵化'
     r'\医学AI深度学习_导师课题项目群\范教授课题_NSCLC-SEER_分子对接多项目\纯生信'
     r'\ARID1A_BAF_最终项目\results_tables')

g = pd.read_csv(T + r'\genotype_based_comparison.csv').dropna(subset=['P'])
s = pd.read_csv(T + r'\structural_vs_codependency_merged.csv')

print('=' * 70)
print('1  基因型条件性依赖：突变型是否比野生型更依赖该亚基')
print('=' * 70)
summ = (g.assign(more=g.direction.eq('mutant_more_dependent'))
          .groupby('partner')
          .agg(n_tissue=('P', 'size'), n_more=('more', 'sum'),
               min_P=('P', 'min'))
          .sort_values(['n_more', 'min_P'], ascending=[False, True]))
print(f'{"亚基":10s}{"组织数":>7s}{"更依赖":>8s}{"最小P":>9s}   注')
for gene, r in summ.iterrows():
    note = ''
    if gene in ('SMARCA4', 'SMARCA2'):
        note = '<- 现有药物靶点'
    if gene == 'ARID1B':
        note = '<- 经典合成致死'
    print(f'{gene:10s}{r.n_tissue:>7.0f}{r.n_more:>8.0f}{r.min_P:>9.3f}   {note}')

sig = g[g.padj_BH_global < 0.05] if 'padj_BH_global' in g else pd.DataFrame()
print(f'\n  全局校正后显著的比较数: {len(sig)}')

print()
print('=' * 70)
print('2  连续型共依赖排序：现有靶点排在第几')
print('=' * 70)
r = s[s.n_contacts_4.rename if False else s.columns[0]].notna()  # keep all rows
t = s[['gene', 'n_contacts_4.5A', 'lung_rho']].copy()
t = t[t['n_contacts_4.5A'] > 0].sort_values('lung_rho', ascending=False)
print(f'{"亚基":10s}{"接触原子对":>10s}{"肺 rho":>9s}   注')
for _, row in t.iterrows():
    note = ''
    if row.gene == 'SMARCA4':
        note = '<- FHD-286 / PFI-3 的靶点'
    print(f'{row.gene:10s}{row["n_contacts_4.5A"]:>10.0f}{row.lung_rho:>9.3f}   {note}')

rank = list(t.gene).index('SMARCA4') + 1
print(f'\n  SMARCA4 在 {len(t)} 个接触亚基中排第 {rank}')
print(f'  排在它前面的: {", ".join(list(t.gene)[:rank-1])}')
