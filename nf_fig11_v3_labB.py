# -*- coding: utf-8 -*-
"""Figures 11B / 11C / 11D / 11E from the 湿实验数据包_v2 dataset (second laboratory).

Statistical policy (unchanged from v2 script, applied to this dataset):
  * n = 3 BIOLOGICAL replicates; technical wells / technical qPCR triplicates are
    averaged first and never treated as independent observations.
  * Two-group comparisons: Welch t-test.
  * Bliss independence: expected viability computed PER REPLICATE
    (exp_i = siARID1A_i x siPartner_i / 100) so uncertainty in both single-agent arms
    propagates; excess tested by one-sample t; six comparisons BH-corrected.
This dataset is NOT pooled with the earlier one — different laboratory, protocol
differences unknown, so cross-lab pooling is not assumed to be valid.
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

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.abspath(os.path.join(HERE, '..', '..'))
RAW = os.path.join(PROJ, '湿实验数据包_v2')
FIGDIR = os.path.join(PROJ, 'results_figures')
TABDIR = os.path.join(PROJ, 'results_tables')
A.apply()

GENES = ['ARID1A', 'ARID1B', 'DPF2', 'PBRM1']
TREATS = ['siNC', 'siARID1A', 'siARID1B', 'siDPF2', 'siPBRM1']
BLOCKS = [('A549', 'ARID1B'), ('A549', 'DPF2'), ('A549', 'PBRM1'),
          ('H1299', 'ARID1B'), ('H1299', 'DPF2'), ('H1299', 'PBRM1')]
CBAF = {'ARID1B', 'DPF2'}


def rd(fn):
    with open(os.path.join(RAW, fn), encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))


def stars(p):
    return '****' if p < 1e-4 else '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < .05 else 'n.s.'


def bh(ps):
    ps = np.asarray(ps, float); o = np.argsort(ps); q = np.empty_like(ps); prev = 1.0
    for r_, i in enumerate(o[::-1]):
        prev = min(prev, ps[i] * len(ps) / (len(ps) - r_)); q[i] = prev
    return q


def bar_pts(ax, x, vals, color, width=.66, jit=.13, seed=0, dot=True):
    vals = np.asarray(vals, float)
    m, sem = vals.mean(), vals.std(ddof=1) / np.sqrt(len(vals))
    ax.bar(x, m, width=width, color=color, zorder=2, edgecolor='white', linewidth=.6)
    ax.errorbar(x, m, yerr=sem, ecolor=A.INK, elinewidth=1.1, capsize=3.2, capthick=1.1,
                zorder=4, fmt='none')
    if dot:
        rng = np.random.RandomState(seed)
        off = np.linspace(-jit, jit, len(vals)) + rng.uniform(-.015, .015, len(vals))
        ax.scatter(x + off, vals, s=16, color=A.DOT, edgecolor='white', linewidth=.4,
                   zorder=5, alpha=.9)
    return m, sem


# ---------------------------------------------------------------- 11B qPCR
q = rd('图11B_原始数据_qPCR_Ct值_修正版.csv')
qd = defaultdict(list)
for r in q:
    qd[(r['处理组'], r['检测基因'])].append(float(r['相对表达量(2^-ΔΔCt)']) * 100)

fig, axes = plt.subplots(1, 4, figsize=(13.0, 3.5))
rows = []
for gi, gene in enumerate(GENES):
    ax = axes[gi]
    for ti, tr in enumerate(TREATS):
        on = (tr == f'si{gene}')
        col = A.CORE if on else (A.NCGY if tr == 'siNC' else A.NEG)
        m, sem = bar_pts(ax, ti, qd[(tr, gene)], col, seed=gi * 10 + ti)
        p = (np.nan if tr == 'siNC' else
             stats.ttest_ind(qd[(tr, gene)], qd[('siNC', gene)], equal_var=False)[1])
        rows.append([tr, gene, round(m, 1), round(np.std(qd[(tr, gene)], ddof=1), 1), 3,
                     '' if tr == 'siNC' else f'{p:.2e}', '' if tr == 'siNC' else stars(p),
                     'on-target' if on else ('control' if tr == 'siNC' else 'off-target')])
        if on:
            ax.text(ti, m + sem + 8, stars(p), ha='center', va='bottom', fontsize=11,
                    fontweight='bold', color=A.INK)
    ax.set_xticks(range(5)); ax.set_xticklabels(TREATS, rotation=38, ha='right', fontsize=8.5)
    ax.set_ylim(0, 145); ax.set_yticks([0, 25, 50, 75, 100, 125])
    ax.axhline(100, color=A.MUTED, lw=.7, ls=(0, (4, 3)), zorder=1)
    ax.set_title(f'{gene} mRNA', fontsize=10.5, pad=7)
    if gi == 0:
        ax.set_ylabel('Relative mRNA level\n(% of siNC)', fontsize=9.5)
    A.open_axes(ax)
fig.suptitle('Knockdown efficiency and specificity (RT-qPCR, n = 3 biological replicates, Welch t-test)',
             fontsize=11.5, y=1.02)
fig.tight_layout()
for e in ('png', 'pdf'):
    fig.savefig(os.path.join(FIGDIR, f'fig11B_qPCR_v3.{e}'), bbox_inches='tight')
plt.close(fig)
with open(os.path.join(TABDIR, 'T1_qPCR_敲低效率_MEASURED_v3.csv'), 'w', newline='',
          encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(['siRNA处理组', '检测基因', '均值(%NC)', 'SD', 'n', 'P(vs siNC, Welch)', '显著性', '类别'])
    w.writerows(rows)

print('qPCR knockdown (Welch):')
for g in GENES:
    nc, kd = qd[('siNC', g)], qd[(f'si{g}', g)]
    tw, pw = stats.ttest_ind(nc, kd, equal_var=False)
    print(f'  {g:8s} siNC {np.mean(nc):6.1f}±{np.std(nc,ddof=1):4.1f} (CV {np.std(nc,ddof=1)/np.mean(nc)*100:4.1f}%)'
          f'   KD {np.mean(kd):5.1f}±{np.std(kd,ddof=1):4.1f}   eff {(1-np.mean(kd)/np.mean(nc))*100:.1f}%'
          f'   Welch P={pw:.2e} {stars(pw)}')
off = [(t, g, np.mean(qd[(t, g)]), stats.ttest_ind(qd[(t, g)], qd[('siNC', g)], equal_var=False)[1])
       for t in TREATS[1:] for g in GENES if f'si{g}' != t]
print(f'  off-target: {sum(1 for *_ , p in off if p < .05)}/16 significant')
for t, g, m, p in off:
    if p < .05:
        print(f'    {t} -> {g}: {m:.1f}%  P={p:.4f}')

# ---------------------------------------------------------------- 11D CCK-8
c = rd('图11D_原始数据_CCK8_OD值_完整版.csv')
v = defaultdict(dict)
for r in c:
    v[(r['细胞系'], r['伙伴基因'], r['生物学重复'])][r['处理类型']] = float(r['相对活力%'])

pack, praw = [], []
for cl, pt in BLOCKS:
    reps = [v[(cl, pt, f'Rep{i}')] for i in (1, 2, 3)]
    a = np.array([r['siARID1A'] for r in reps])
    b = np.array([r['siPartner'] for r in reps])
    d = np.array([r['double'] for r in reps])
    exc = a * b / 100 - d
    t, p = stats.ttest_1samp(exc, 0)
    pack.append(dict(cl=cl, pt=pt, a=a, b=b, d=d, exp=a * b / 100, exc=exc, t=t, p=p))
    praw.append(p)
for k, qv in zip(pack, bh(praw)):
    k['q'] = qv

fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.3), sharey=True)
rows2 = []
for ax, cl in zip(axes, ['A549', 'H1299']):
    for j, k in enumerate([x for x in pack if x['cl'] == cl]):
        base = j * 5
        for offs, arr, col, lab in [(0, np.full(3, 100.), A.NCGY, 'siNC'),
                                    (1, k['a'], '#7B4EA8', 'siARID1A'),
                                    (2, k['b'], '#C79A3A', 'siPartner'),
                                    (3, k['d'], A.CORE if k['pt'] in CBAF else A.NEG, 'double')]:
            bar_pts(ax, base + offs, arr, col, width=.78, jit=.16, seed=base + offs,
                    dot=(lab != 'siNC'))
        e = k['exp'].mean()
        ax.plot([base + 2.55, base + 3.45], [e, e], color=A.INK, lw=1.6, ls=(0, (3, 2)), zorder=6)
        ax.annotate('', xy=(base + 3, k['d'].mean()), xytext=(base + 3, e),
                    arrowprops=dict(arrowstyle='<->', color=A.INK, lw=1.0), zorder=7)
        ax.text(base + 3.55, (e + k['d'].mean()) / 2, f'{k["exc"].mean():.1f}', fontsize=8.5,
                va='center', color=A.INK)
        ax.text(base + 1.5, 118, f'+{k["pt"]}', ha='center', fontsize=10, fontweight='bold',
                color=A.INK)
        ax.text(base + 1.5, 110, f'q={k["q"]:.3f} {stars(k["q"])}', ha='center', fontsize=8.2,
                color=A.MUTED)
        rows2.append([cl, k['pt'], 'cBAF接触亚基' if k['pt'] in CBAF else 'PBAF非接触亚基',
                      round(k['a'].mean(), 1), round(k['a'].std(ddof=1), 1),
                      round(k['b'].mean(), 1), round(k['b'].std(ddof=1), 1),
                      round(k['exp'].mean(), 1), round(k['d'].mean(), 1),
                      round(k['d'].std(ddof=1), 1), round(k['exc'].mean(), 1),
                      round(k['exc'].std(ddof=1), 1), 3, round(k['t'], 2),
                      f'{k["p"]:.4f}', f'{k["q"]:.4f}', stars(k['q'])])
    ax.set_xticks([j * 5 + i for j in range(3) for i in range(4)])
    ax.set_xticklabels(['siNC', 'siARID1A', 'siPartner', 'double'] * 3, rotation=42,
                       ha='right', fontsize=8)
    ax.set_ylim(0, 126)
    ax.axhline(100, color=A.MUTED, lw=.7, ls=(0, (4, 3)), zorder=1)
    ax.set_title(cl, fontsize=11.5, pad=22)
    A.open_axes(ax)
axes[0].set_ylabel('Relative viability (% of siNC)', fontsize=9.5)
h = [plt.Rectangle((0, 0), 1, 1, color=x) for x in (A.NCGY, '#7B4EA8', '#C79A3A', A.CORE, A.NEG)]
h.append(plt.Line2D([0], [0], color=A.INK, lw=1.6, ls=(0, (3, 2))))
axes[1].legend(h, ['siNC', 'siARID1A', 'partner alone', 'double (cBAF partner)',
                   'double (PBAF control)', 'Bliss additivity expectation'],
               frameon=False, fontsize=8.2, loc='upper center', bbox_to_anchor=(-.05, -.24),
               ncol=3, handlelength=1.5)
fig.suptitle('Double knockdown vs Bliss independence (CCK-8, n = 3 biological replicates, BH-corrected)',
             fontsize=11.5, y=1.0)
fig.tight_layout()
for e in ('png', 'pdf'):
    fig.savefig(os.path.join(FIGDIR, f'fig11D_CCK8_bliss_v3.{e}'), bbox_inches='tight')
plt.close(fig)
with open(os.path.join(TABDIR, 'T5_CCK8_Bliss_MEASURED_v3.csv'), 'w', newline='',
          encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(['细胞系', '伙伴基因', '分类', 'siARID1A均值%', 'SD', '伙伴单敲均值%', 'SD',
                'Bliss预期%(逐重复)', '双敲实测%', 'SD', '超出Bliss(百分点)', 'SD', 'n', 't',
                'P(单样本t)', 'q(BH,6次)', '显著性(q)'])
    w.writerows(rows2)

print('\nBliss (paired, per-replicate expectation):')
for k in pack:
    print(f'  {k["cl"]+"/"+k["pt"]:16s} siA {k["a"].mean():5.1f}  siP {k["b"].mean():5.1f}  '
          f'exp {k["exp"].mean():5.1f}  double {k["d"].mean():5.1f}  '
          f'excess {k["exc"].mean():5.1f}  t={k["t"]:5.2f}  P={k["p"]:.4f}  q={k["q"]:.4f} {stars(k["q"])}')
cb = [k['exc'].mean() for k in pack if k['pt'] in CBAF]
pb = [k['exc'].mean() for k in pack if k['pt'] not in CBAF]
print(f'  cBAF {np.mean(cb):.1f}±{np.std(cb,ddof=1):.1f}   PBAF {np.mean(pb):.1f}±{np.std(pb,ddof=1):.1f}')
ac = np.concatenate([k['exc'] for k in pack if k['pt'] in CBAF])
ap = np.concatenate([k['exc'] for k in pack if k['pt'] not in CBAF])
print(f'  replicate-level cBAF vs PBAF: Mann-Whitney P={stats.mannwhitneyu(ac, ap, alternative="greater")[1]:.2e}')
print('\n  partner-alone vs 100%:')
for k in pack:
    print(f'    {k["cl"]+"/"+k["pt"]:16s} si{k["pt"]} = {k["b"].mean():5.1f}±{k["b"].std(ddof=1):.1f}%'
          f'  P={stats.ttest_1samp(k["b"], 100)[1]:.4f}')
print('\n  power (80%, alpha=.05, one-sample t):')
from scipy.stats import nct, t as tdist
for k in pack:
    d_ = k['exc'].mean() / k['exc'].std(ddof=1)
    n_ = 2
    while n_ <= 60:
        cr = tdist.ppf(.975, n_ - 1)
        if 1 - nct.cdf(cr, n_ - 1, d_ * np.sqrt(n_)) + nct.cdf(-cr, n_ - 1, d_ * np.sqrt(n_)) >= .8:
            break
        n_ += 1
    print(f'    {k["cl"]+"/"+k["pt"]:16s} d={d_:5.2f}  n needed = {n_}')
print('\nfigures + tables written')
