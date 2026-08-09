# -*- coding: utf-8 -*-
"""Audit batch 3, then pool batch 2 + batch 3 (n=6) for the Bliss analysis."""
import csv, os
from collections import defaultdict
import numpy as np
from scipy import stats

B2 = r'D:\桌面\湿实验最终图'
B3 = (r'D:\桌面\_按项目整理\02_生信数据与医学AI项目\07_医学AI_深度学习_多项目孵化'
      r'\医学AI深度学习_导师课题项目群\范教授课题_NSCLC-SEER_分子对接多项目\纯生信'
      r'\ARID1A_BAF_最终项目\湿实验数据包_v2')
GENES = ['ARID1A', 'ARID1B', 'DPF2', 'PBRM1']
BLOCKS = [('A549', 'ARID1B'), ('A549', 'DPF2'), ('A549', 'PBRM1'),
          ('H1299', 'ARID1B'), ('H1299', 'DPF2'), ('H1299', 'PBRM1')]
CBAF = {'ARID1B', 'DPF2'}


def rd(base, fn):
    with open(os.path.join(base, fn), encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))


def bh(ps):
    ps = np.asarray(ps, float); o = np.argsort(ps); q = np.empty_like(ps); prev = 1.0
    for rank, i in enumerate(o[::-1]):
        prev = min(prev, ps[i] * len(ps) / (len(ps) - rank)); q[i] = prev
    return q


print('#' * 76)
print('# A. batch-3 internal consistency')
print('#' * 76)
q3 = rd(B3, '图11B_原始数据_qPCR_Ct值_修正版.csv')
bad = 0
for r in q3:
    ct = [float(r[f'Ct_技术重复{i}']) for i in (1, 2, 3)]
    bad += (abs(np.mean(ct) - float(r['平均Ct'])) > .011
            or abs(float(r['平均Ct']) - float(r['内参Ct_GAPDH']) - float(r['ΔCt'])) > .011
            or abs(2 ** (-float(r['ΔΔCt'])) - float(r['相对表达量(2^-ΔΔCt)'])) > .0015)
print(f'  qPCR arithmetic mismatches: {bad}/{len(q3)}')
c3 = rd(B3, '图11D_原始数据_CCK8_OD值_完整版.csv')
bad = sum(abs(np.mean([float(r[f'OD450_复孔{i}']) for i in (1, 2, 3)]) - float(r['平均OD450'])) > .0011
          for r in c3)
print(f'  CCK-8 arithmetic mismatches: {bad}/{len(c3)}')

qd = defaultdict(list)
for r in q3:
    qd[(r['处理组'], r['检测基因'])].append(float(r['相对表达量(2^-ΔΔCt)']))
sup = {r['检测基因']: r for r in rd(B3, '图11B_统计检验_qPCR敲低效率.csv')}
print('\n  qPCR: my recomputation vs their table')
print(f'  {"gene":8s}{"t mine":>9s}{"t theirs":>10s}{"P Student":>12s}{"P Welch":>11s}{"their P":>11s}')
for g in GENES:
    nc, kd = qd[('siNC', g)], qd[(f'si{g}', g)]
    t, p = stats.ttest_ind(nc, kd, equal_var=True)
    tw, pw = stats.ttest_ind(nc, kd, equal_var=False)
    print(f'  {g:8s}{t:>9.2f}{float(sup[g]["t值"]):>10.2f}{p:>12.2e}{pw:>11.2e}{sup[g]["P值"]:>11s}')
    print(f'  {"":8s}  siNC CV mine {np.std(nc,ddof=1)/np.mean(nc)*100:4.1f}%  theirs {sup[g]["siNC_CV%"]}%'
          f'   KD mine {(1-np.mean(kd)/np.mean(nc))*100:.1f}%  theirs {sup[g]["敲低效率%"]}%')

print('\n  off-target (16 non-target combos):',
      sum(stats.ttest_ind(qd[(t_, g)], qd[('siNC', g)], equal_var=False)[1] < .05
          for t_ in ['siARID1A', 'siARID1B', 'siDPF2', 'siPBRM1'] for g in GENES if f'si{g}' != t_),
      '/16 significant')

print('\n' + '#' * 76)
print('# B. batch-to-batch drift (same experiment, three versions)')
print('#' * 76)


def viab(base, fn, armkey):
    d = defaultdict(dict)
    for r in rd(base, fn):
        d[(r['细胞系'], r['伙伴基因'], r['生物学重复'])][r[armkey]] = float(r['相对活力%'])
    return d


v2 = viab(B2, '图11D_原始数据_CCK8_OD值_完整版.csv', '处理类型')
v3 = viab(B3, '图11D_原始数据_CCK8_OD值_完整版.csv', '处理类型')
v1 = defaultdict(dict)
for r in rd(B2, '图11D_原始数据_CCK8_OD值.csv'):
    v1[(r['细胞系'], r['伙伴基因'], r['生物学重复'])][r['处理组']] = float(r['相对活力%'])

print(f'  {"block":16s}{"double: b1":>12s}{"b2":>8s}{"b3":>8s}   |{"excess over Bliss: b2":>24s}{"b3":>8s}')
for cl, pt in BLOCKS:
    d1 = np.mean([v1[(cl, pt, f'Rep{i}')][f'siARID1A+si{pt}'] for i in (1, 2, 3)])
    r2 = [v2[(cl, pt, f'Rep{i}')] for i in (1, 2, 3)]
    r3 = [v3[(cl, pt, f'Rep{i}')] for i in (1, 2, 3)]
    e2 = np.mean([r['siARID1A'] * r['siPartner'] / 100 - r['double'] for r in r2])
    e3 = np.mean([r['siARID1A'] * r['siPartner'] / 100 - r['double'] for r in r3])
    print(f'  {cl+"/"+pt:16s}{d1:>12.1f}{np.mean([r["double"] for r in r2]):>8.1f}'
          f'{np.mean([r["double"] for r in r3]):>8.1f}   |{e2:>24.1f}{e3:>8.1f}')

print('\n' + '#' * 76)
print('# C. POOLED batch2 + batch3  (n = 6 independent replicates)')
print('#' * 76)
rows, praw, pack = [], [], []
for cl, pt in BLOCKS:
    reps = ([v2[(cl, pt, f'Rep{i}')] for i in (1, 2, 3)] +
            [v3[(cl, pt, f'Rep{i}')] for i in (1, 2, 3)])
    a = np.array([r['siARID1A'] for r in reps])
    b = np.array([r['siPartner'] for r in reps])
    d = np.array([r['double'] for r in reps])
    exc = a * b / 100 - d
    t, p = stats.ttest_1samp(exc, 0)
    pack.append(dict(cl=cl, pt=pt, a=a, b=b, d=d, exc=exc, t=t, p=p))
    praw.append(p)
qs = bh(praw)
print(f'  {"block":16s}{"siA":>12s}{"siP":>12s}{"Bliss":>8s}{"double":>12s}'
      f'{"excess":>9s}{"t":>7s}{"P":>10s}{"q(BH)":>9s}')
for k, qv in zip(pack, qs):
    k['q'] = qv
    print(f'  {k["cl"]+"/"+k["pt"]:16s}{k["a"].mean():>7.1f}±{k["a"].std(ddof=1):<4.1f}'
          f'{k["b"].mean():>7.1f}±{k["b"].std(ddof=1):<4.1f}{(k["a"]*k["b"]/100).mean():>8.1f}'
          f'{k["d"].mean():>7.1f}±{k["d"].std(ddof=1):<4.1f}{k["exc"].mean():>9.1f}'
          f'{k["t"]:>7.2f}{k["p"]:>10.4f}{qv:>9.4f}')

cb = [k['exc'].mean() for k in pack if k['pt'] in CBAF]
pb = [k['exc'].mean() for k in pack if k['pt'] not in CBAF]
print(f'\n  cBAF excess {np.mean(cb):.1f} ± {np.std(cb,ddof=1):.1f}   '
      f'PBAF excess {np.mean(pb):.1f} ± {np.std(pb,ddof=1):.1f}')
allc = np.concatenate([k['exc'] for k in pack if k['pt'] in CBAF])
allp = np.concatenate([k['exc'] for k in pack if k['pt'] not in CBAF])
u, pu = stats.mannwhitneyu(allc, allp, alternative='greater')
print(f'  cBAF vs PBAF at replicate level (n={len(allc)} vs {len(allp)}): '
      f'Mann-Whitney U={u:.0f}, P={pu:.2e}')

print('\n  batch effect check (batch2 vs batch3 excess, per block):')
for cl, pt in BLOCKS:
    e2 = np.array([v2[(cl, pt, f'Rep{i}')]['siARID1A'] * v2[(cl, pt, f'Rep{i}')]['siPartner'] / 100
                   - v2[(cl, pt, f'Rep{i}')]['double'] for i in (1, 2, 3)])
    e3 = np.array([v3[(cl, pt, f'Rep{i}')]['siARID1A'] * v3[(cl, pt, f'Rep{i}')]['siPartner'] / 100
                   - v3[(cl, pt, f'Rep{i}')]['double'] for i in (1, 2, 3)])
    p = stats.ttest_ind(e2, e3, equal_var=False)[1]
    print(f'    {cl+"/"+pt:16s} b2 {e2.mean():5.1f}  b3 {e3.mean():5.1f}   P={p:.3f}'
          + ('   !! batch effect' if p < .05 else ''))

print('\n' + '#' * 76)
print('# D. Panel A WB (knockdown specificity, protein level)')
print('#' * 76)
w = rd(B3, '图11C_PanelA_WB敲低特异性_灰度定量.csv')
gap = {r['泳道']: float(r['净灰度值']) for r in w if r['蛋白'] == 'GAPDH'}
bad = 0
for r in w:
    if r['蛋白'] == 'GAPDH':
        continue
    if abs(float(r['ImageJ_IntDen']) - float(r['背景值']) - float(r['净灰度值'])) > 1.5:
        bad += 1
    if abs(float(r['净灰度值']) / gap[r['泳道']] - float(r['GAPDH归一化'])) > .0006:
        bad += 1
print(f'  arithmetic mismatches: {bad}')
print(f'  {"protein":9s}' + ''.join(f'{l:>12s}' for l in ['siNC', 'siARID1A', 'siARID1B', 'siDPF2', 'siPBRM1']))
for g in GENES:
    vals = {r['泳道']: float(r['相对siNC%']) for r in w if r['蛋白'] == g}
    print(f'  {g:9s}' + ''.join(f'{vals[l]:>12.1f}' for l in
                                ['siNC', 'siARID1A', 'siARID1B', 'siDPF2', 'siPBRM1']))
print('  NOTE: one lane per condition -> n = 1, no statistics possible on this panel.')
