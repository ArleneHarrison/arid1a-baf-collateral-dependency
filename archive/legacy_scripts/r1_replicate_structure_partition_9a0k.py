# -*- coding: utf-8 -*-
"""Does the structure-genetics association replicate in a second structure?

The manuscript's test is Mann-Whitney: co-dependency of the subunits ARID1A
contacts vs the subunits it does not, partitioned by PDB 6LTJ (experimental
cryo-EM, 3.70 A). 9A0K is an independent integrative model of the ENDOGENOUS
cBAF complex (3D-EM + crosslinking-MS + comparative models) that additionally
resolves SMARCC1. Re-partitioning by 9A0K is a genuine external replication of
the partition, on the same dependency data.
"""
import io
import sys

import pandas as pd
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
T = (r'D:\桌面\_按项目整理\02_生信数据与医学AI项目\07_医学AI_深度学习_多项目孵化'
     r'\医学AI深度学习_导师课题项目群\范教授课题_NSCLC-SEER_分子对接多项目\纯生信'
     r'\ARID1A_BAF_最终项目\results_tables')

m = pd.read_csv(T + r'\structural_vs_codependency_merged.csv').set_index('gene')

# SMARCC1 is absent from 6LTJ so it is missing from that table; its dependency
# comes from the genome-wide tables (rank 3 in both tissues).
EXTRA = {'SMARCC1': dict(lung_rho=0.38041910822193736,
                         pancancer_rho=0.3644314674382256)}

# contact counts recomputed from 9A0K with the identical 4.5 A heavy-atom rule
C9 = {'SMARCB1': 748, 'SMARCD1': 401, 'SMARCA4': 374, 'SMARCC2': 233,
      'SMARCC1': 186, 'SMARCE1': 0, 'DPF2': 0, 'ACTL6A': 0, 'ACTB': 0}
C6 = {'SMARCC2': 754, 'SMARCA4': 557, 'SMARCD1': 553, 'SMARCB1': 300,
      'SMARCE1': 200, 'DPF2': 137, 'ACTL6A': 0, 'ACTB': 0}
# PBAF/ncBAF-specific subunits are in neither complex -> non-contacting in both
EXCL = ['PBRM1', 'ARID2', 'BRD7', 'PHF10', 'BRD9', 'BICRA']


def rho(gene, tissue):
    col = 'lung_rho' if tissue == 'lung' else 'pancancer_rho'
    if gene in EXTRA:
        return EXTRA[gene][col]
    return float(m.loc[gene, col])


print('=' * 78)
print('ARID1A 接触谱：实验冷冻电镜 6LTJ  vs  内源复合体整合模型 9A0K')
print('=' * 78)
print(f'{"亚基":10s}{"6LTJ":>8s}{"9A0K":>8s}{"肺 rho":>10s}   一致?')
allg = ['SMARCC2', 'SMARCA4', 'SMARCD1', 'SMARCB1', 'SMARCE1', 'DPF2',
        'SMARCC1', 'ACTL6A', 'ACTB']
agree = 0
for g in allg:
    a, b = C6.get(g), C9.get(g)
    a_s = '未解析' if a is None else str(a)
    if a is None:
        ok = '—(6LTJ无)'
    else:
        same = (a > 0) == (b > 0)
        agree += same
        ok = '一致' if same else '**不一致**'
    print(f'{g:10s}{a_s:>8s}{b:>8d}{rho(g, "lung"):>10.3f}   {ok}')
print(f'\n定性接触/不接触的判定一致率: {agree}/8')

for tissue in ('lung', 'pancancer'):
    print()
    print('=' * 78)
    print(f'Mann-Whitney（与正文同一检验），组织 = {tissue}')
    print('=' * 78)
    for name, cmap, extra_noncontact in (
            ('6LTJ（正文所用）', C6, []),
            ('9A0K（独立结构）', C9, [])):
        con = [g for g in cmap if cmap[g] > 0]
        non = [g for g in cmap if cmap[g] == 0] + EXCL + extra_noncontact
        cv = [rho(g, tissue) for g in con]
        nv = [rho(g, tissue) for g in non]
        u, p = stats.mannwhitneyu(cv, nv, alternative='two-sided')
        print(f'  {name}')
        print(f'    接触 (n={len(con)}) 中位 rho = {pd.Series(cv).median():+.3f}  {sorted(con)}')
        print(f'    非接触 (n={len(non)}) 中位 rho = {pd.Series(nv).median():+.3f}')
        print(f'    U = {u:.1f},  P = {p:.4f}')
