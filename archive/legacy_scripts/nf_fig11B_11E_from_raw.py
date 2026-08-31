# -*- coding: utf-8 -*-
"""Redraw 图11B (qPCR knockdown efficiency + specificity) and 图11E panel B
(WB co-stability densitometry) from the MEASURED raw data.

Every plotted point is a real biological replicate read out of the raw CSV --
unlike assaystyle.dot_bar, which synthesises dots from mean+SEM.

Inputs :  D:\\桌面\\湿实验最终图\\图11B_原始数据_qPCR_Ct值.csv
          D:\\桌面\\湿实验最终图\\图11E_PanelB_WB共稳定性_灰度定量.csv
Outputs:  results_figures/fig11B_qPCR_from_raw.{png,pdf}
          results_figures/fig11E_WB_costability_from_raw.{png,pdf}
          results_tables/T1_qPCR_敲低效率_MEASURED.csv
          results_tables/T3_WB_共稳定性_MEASURED.csv
"""
import csv, os, sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import assaystyle as A

RAW = r'D:\桌面\湿实验最终图'
HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.abspath(os.path.join(HERE, '..'))
TABDIR = os.path.abspath(os.path.join(HERE, '..', '..', 'results_tables'))

A.apply()

GENES = ['ARID1A', 'ARID1B', 'DPF2', 'PBRM1']
TREATS = ['siNC', 'siARID1A', 'siARID1B', 'siDPF2', 'siPBRM1']


def rd(fn):
    with open(fn, encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))


def stars(p):
    return '****' if p < 1e-4 else '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < .05 else 'n.s.'


def bar_with_real_points(ax, x, vals, color, width=0.66, jitter=0.13, seed=0):
    """Bar = mean, error bar = SEM, dots = the actual replicates."""
    vals = np.asarray(vals, float)
    m = vals.mean()
    sem = vals.std(ddof=1) / np.sqrt(len(vals))
    ax.bar(x, m, width=width, color=color, zorder=2, edgecolor='white', linewidth=0.6)
    ax.errorbar(x, m, yerr=sem, ecolor=A.INK, elinewidth=1.1, capsize=3.2,
                capthick=1.1, zorder=4, fmt='none')
    rng = np.random.RandomState(seed)
    off = np.linspace(-jitter, jitter, len(vals)) + rng.uniform(-0.015, 0.015, len(vals))
    ax.scatter(x + off, vals, s=16, color=A.DOT, edgecolor='white',
               linewidth=0.4, zorder=5, alpha=0.9)
    return m, sem


# ----------------------------------------------------------------- figure 11B
q = rd(os.path.join(RAW, '图11B_原始数据_qPCR_Ct值.csv'))
qd = defaultdict(list)
for r in q:
    qd[(r['处理组'], r['检测基因'])].append(float(r['相对表达量(2^-ΔΔCt)']) * 100)

fig, axes = plt.subplots(1, 4, figsize=(13.0, 3.5))
rows = []
for gi, gene in enumerate(GENES):
    ax = axes[gi]
    for ti, tr in enumerate(TREATS):
        on_target = (tr == f'si{gene}')
        col = A.CORE if on_target else (A.NCGY if tr == 'siNC' else A.NEG)
        m, sem = bar_with_real_points(ax, ti, qd[(tr, gene)], col, seed=gi * 10 + ti)
        sd = np.std(qd[(tr, gene)], ddof=1)
        p = (np.nan if tr == 'siNC' else
             stats.ttest_ind(qd[(tr, gene)], qd[('siNC', gene)], equal_var=False)[1])
        rows.append([tr, gene, round(m, 1), round(sd, 1), round(sem, 1), len(qd[(tr, gene)]),
                     '' if tr == 'siNC' else f'{p:.4f}',
                     '' if tr == 'siNC' else stars(p),
                     'on-target' if on_target else ('control' if tr == 'siNC' else 'off-target')])
        if on_target:
            ax.text(ti, m + sem + 7, stars(p), ha='center', va='bottom',
                    fontsize=11, color=A.INK, fontweight='bold')

    ax.set_xticks(range(len(TREATS)))
    ax.set_xticklabels(TREATS, rotation=38, ha='right', fontsize=8.5)
    ax.set_ylim(0, 145)
    ax.set_yticks([0, 25, 50, 75, 100, 125])
    ax.axhline(100, color=A.MUTED, lw=0.7, ls=(0, (4, 3)), zorder=1)
    ax.set_title(f'{gene} mRNA', fontsize=10.5, pad=7)
    if gi == 0:
        ax.set_ylabel('Relative mRNA level\n(% of siNC)', fontsize=9.5)
    A.open_axes(ax)

fig.suptitle('Knockdown efficiency and specificity (RT-qPCR, n = 3 biological replicates)',
             fontsize=11.5, y=1.02)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(FIGDIR, f'fig11B_qPCR_from_raw.{ext}'), bbox_inches='tight')
plt.close(fig)

with open(os.path.join(TABDIR, 'T1_qPCR_敲低效率_MEASURED.csv'), 'w',
          newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(['siRNA处理组', '检测基因', '相对表达量均值(%NC)', 'SD', 'SEM', 'n',
                'P(vs siNC, Welch)', '显著性', '类别'])
    w.writerows(rows)

# ----------------------------------------------------------------- figure 11E
w_ = rd(os.path.join(RAW, '图11E_PanelB_WB共稳定性_灰度定量.csv'))
wd = defaultdict(list)
for r in w_:
    if r['蛋白'] != 'GAPDH':
        wd[(r['蛋白'], r['分组'])].append(float(r['相对siNC%']))

CONTACT = {'SMARCA4': 557, 'SMARCC2': 754, 'SMARCB1': 300, 'DPF2': 137}
ORDER = ['ARID1A', 'SMARCA4', 'SMARCC2', 'SMARCB1', 'SMARCC1', 'DPF2',
         'ARID1B', 'BRD9', 'SMARCA2']

TARGET_CLR = '#7B4EA8'   # the knocked-down gene itself, distinct from both groups

fig, ax = plt.subplots(figsize=(10.6, 4.2))
rows3 = []
for i, pr in enumerate(ORDER):
    nc, kd = wd[(pr, 'siNC')], wd[(pr, 'siARID1A')]
    if pr == 'ARID1A':
        clr = TARGET_CLR
    elif pr in CONTACT or pr == 'SMARCC1':
        clr = A.CORE
    else:
        clr = A.NEG
    bar_with_real_points(ax, i - 0.19, nc, A.NCGY, width=0.35, jitter=0.075, seed=i)
    m, sem = bar_with_real_points(ax, i + 0.19, kd, clr, width=0.35, jitter=0.075, seed=i + 50)
    p = stats.ttest_ind(kd, nc, equal_var=False)[1]
    ax.text(i + 0.19, m + sem + 4, stars(p), ha='center', va='bottom',
            fontsize=9.5, color=A.INK, fontweight='bold')
    rows3.append([pr,
                  CONTACT.get(pr, 'no chain in 6LTJ' if pr == 'SMARCC1'
                              else ('target' if pr == 'ARID1A' else 0)),
                  round(np.mean(nc), 1), round(np.std(nc, ddof=1), 1),
                  round(np.mean(kd), 1), round(np.std(kd, ddof=1), 1),
                  len(kd), f'{p:.4f}', stars(p)])

ax.set_xticks(range(len(ORDER)))
ax.set_xticklabels(ORDER, rotation=30, ha='right', fontsize=9.5)
ax.axhline(100, color=A.MUTED, lw=0.7, ls=(0, (4, 3)), zorder=1)
ax.set_ylim(0, 132)
ax.set_yticks([0, 20, 40, 60, 80, 100, 120])
ax.set_ylabel('Protein level (% of siNC,\nGAPDH-normalised)', fontsize=9.5)
ax.set_title('cBAF subunit co-stability after ARID1A knockdown '
             '(Western blot densitometry, n = 3)', fontsize=11, pad=9)
A.open_axes(ax)

h = [plt.Rectangle((0, 0), 1, 1, color=c)
     for c in (A.NCGY, TARGET_CLR, A.CORE, A.NEG)]
ax.legend(h, ['siNC', 'siARID1A — knockdown target',
              'siARID1A — direct contact in 6LTJ / core subunit',
              'siARID1A — no contact with ARID1A'],
          frameon=False, fontsize=8.6, ncol=4,
          loc='upper center', bbox_to_anchor=(0.5, -0.30), handlelength=1.4,
          columnspacing=1.3)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(FIGDIR, f'fig11E_WB_costability_from_raw.{ext}'),
                bbox_inches='tight')
plt.close(fig)

with open(os.path.join(TABDIR, 'T3_WB_共稳定性_MEASURED.csv'), 'w',
          newline='', encoding='utf-8-sig') as f:
    wr = csv.writer(f)
    wr.writerow(['蛋白', '6LTJ接触原子对数', 'siNC均值%', 'siNC_SD',
                 'siARID1A均值%', 'siARID1A_SD', 'n', 'P(Welch)', '显著性'])
    wr.writerows(rows3)

con = [np.mean(wd[(p, 'siARID1A')]) for p in ['SMARCA4', 'SMARCC2', 'SMARCB1', 'SMARCC1', 'DPF2']]
non = [np.mean(wd[(p, 'siARID1A')]) for p in ['ARID1B', 'BRD9', 'SMARCA2']]
u, pu = stats.mannwhitneyu(con, non, alternative='less')
print(f'contact/core {np.mean(con):.1f}%  vs  non-contact {np.mean(non):.1f}%   '
      f'Mann-Whitney P={pu:.4f}')
print('figures + measured tables written')
