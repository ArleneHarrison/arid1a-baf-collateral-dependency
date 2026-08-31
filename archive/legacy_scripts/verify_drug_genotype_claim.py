# -*- coding: utf-8 -*-
"""Re-verify the manuscript's negative claim about ARID1A-genotype selectivity.

A negative result is only worth reporting if the study could have detected the
effect had it been there. This checks three things in order:
  1. is the claim as written consistent with the table?
  2. how large is the observed difference, in the units of the assay?
  3. what is the smallest difference the design could have detected?
"""
import numpy as np
import pandas as pd
from scipy import stats

T = (r'D:\桌面\_按项目整理\02_生信数据与医学AI项目\07_医学AI_深度学习_多项目孵化'
     r'\医学AI深度学习_导师课题项目群\范教授课题_NSCLC-SEER_分子对接多项目\纯生信'
     r'\ARID1A_BAF_最终项目\results_tables\arid1a_drug_sensitivity_v2_summary.csv')
d = pd.read_csv(T)
d['dAUC'] = d.mean_AUC_mutant - d.mean_AUC_wt

print('=' * 74)
print('1  论文原话："no comparison survives global correction"')
print('=' * 74)
print(f'  比较总数              : {len(d)}  （13 药 × 2 层级）')
print(f'  最小 padj_BH_global   : {d.padj_BH_global.min():.4f}')
print(f'  padj < 0.05 的数量    : {(d.padj_BH_global < 0.05).sum()}')
print(f'  未校正 P < 0.05 的数量: {(d.P < 0.05).sum()}   '
      f'（{d.loc[d.P < 0.05, "drug"].tolist()}）')
print(f'  -> 结论与数据一致: {(d.padj_BH_global < 0.05).sum() == 0}')

print()
print('=' * 74)
print('2  SWI/SNF 直接靶向的两个化合物，效应有多大')
print('=' * 74)
swi = d[d.drug.isin(['PFI-3', 'I-BRD9'])].copy()
print(f'{"层级":11s}{"药物":9s}{"n_mut":>6s}{"n_wt":>6s}'
      f'{"AUC_mut":>9s}{"AUC_wt":>9s}{"差值":>9s}{"P":>8s}')
for _, r in swi.iterrows():
    print(f'{r.label:11s}{r.drug:9s}{r.n_mutant:>6.0f}{r.n_wt:>6.0f}'
          f'{r.mean_AUC_mutant:>9.4f}{r.mean_AUC_wt:>9.4f}'
          f'{r.dAUC:>+9.4f}{r.P:>8.3f}')
print()
print('  注意 AUC 本身：AUC 越接近 1 = 药物对细胞几乎无杀伤。')
print(f'  PFI-3 / I-BRD9 的 AUC 全部在 {swi[["mean_AUC_mutant","mean_AUC_wt"]].values.min():.3f}'
      f'–{swi[["mean_AUC_mutant","mean_AUC_wt"]].values.max():.3f} 之间。')

print()
print('=' * 74)
print('3  检验力：这个设计能查出多大的差异？')
print('=' * 74)
# GDSC AUC 的组内标准差在这类面板里通常 0.06-0.10；取 0.08 做中间估计，
# 并给出 0.06 / 0.10 的边界，避免单点假设
for sd in (0.06, 0.08, 0.10):
    print(f'\n  假定组内 SD = {sd:.2f}')
    for label, n1, n2 in [('lung-only', 8, 180), ('pan-cancer', 84, 890)]:
        # 双样本 t 检验，alpha=0.05 双侧，power=0.80
        za, zb = stats.norm.ppf(0.975), stats.norm.ppf(0.80)
        mdd = (za + zb) * sd * np.sqrt(1 / n1 + 1 / n2)
        print(f'    {label:11s} n={n1:>3d} vs {n2:>3d}   '
              f'80% 检验力下最小可检出差值 = {mdd:.4f} AUC')

print()
print('=' * 74)
print('4  判读')
print('=' * 74)
obs = swi.dAUC.abs().max()
print(f'  SWI/SNF 化合物观察到的最大绝对差值: {obs:.4f} AUC')
print('  肺癌层 ARID1A 突变系仅 n=8 —— 这一层检验力很低，')
print('  但泛癌层 n=84 vs 890，最小可检出差值远小于任何有意义的效应，')
print('  且观察值仍接近 0。两层结论一致。')
