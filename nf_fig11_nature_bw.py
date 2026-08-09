# -*- coding: utf-8 -*-
"""Figure 11 (b/c/d/e) — Nature typography, GraphPad-Prism black-and-white bars.

Single data source: 湿实验数据包_v2 (one laboratory only). No cross-laboratory pooling
or comparison anywhere in these figures.

Statistics: n = 3 biological replicates, technical wells averaged first; Welch t-test for
two-group comparisons; Bliss expectation computed per replicate
(exp_i = siARID1A_i x siPartner_i / 100), one-sample t-test on the excess, Benjamini-
Hochberg correction across the six blocks. Panel d also reports the pre-specified primary
contrast (cBAF vs PBAF excess, replicate level, Mann-Whitney).
"""
import csv, os, sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nfprism as S

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.abspath(os.path.join(HERE, '..', '..'))
RAW = os.path.join(PROJ, '湿实验数据包_v2')
FIG = os.path.join(PROJ, 'results_figures')
S.apply()

GENES = ['ARID1A', 'ARID1B', 'DPF2', 'PBRM1']
TREATS = ['siNC', 'siARID1A', 'siARID1B', 'siDPF2', 'siPBRM1']
BLOCKS = [('A549', 'ARID1B'), ('A549', 'DPF2'), ('A549', 'PBRM1'),
          ('H1299', 'ARID1B'), ('H1299', 'DPF2'), ('H1299', 'PBRM1')]
CBAF = {'ARID1B', 'DPF2'}
CONTACT = ['SMARCA4', 'SMARCC2', 'SMARCB1', 'SMARCC1', 'DPF2']
NOCONTACT = ['ARID1B', 'BRD9', 'SMARCA2']
MM = 1 / 25.4


def rd(fn):
    with open(os.path.join(RAW, fn), encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))


def bh(ps):
    ps = np.asarray(ps, float); o = np.argsort(ps); q = np.empty_like(ps); prev = 1.
    for r_, i in enumerate(o[::-1]):
        prev = min(prev, ps[i] * len(ps) / (len(ps) - r_)); q[i] = prev
    return q


def figletter(fig, ax, letter, dx=-0.055, dy=0.030):
    """Panel letter placed in figure coordinates, clear of the y-axis label."""
    bb = ax.get_position()
    fig.text(bb.x0 + dx, bb.y1 + dy, letter, fontsize=9, fontweight='bold',
             ha='left', va='bottom')


def save(fig, name):
    for e in ('png', 'pdf'):
        fig.savefig(os.path.join(FIG, f'{name}.{e}'), dpi=600)
    plt.close(fig)


# ══════════════════════════════════════════════════════════ Fig 11b  qPCR
q = rd('图11B_原始数据_qPCR_Ct值_修正版.csv')
qd = defaultdict(list)
for r in q:
    qd[(r['处理组'], r['检测基因'])].append(float(r['相对表达量(2^-ΔΔCt)']) * 100)

fig = plt.figure(figsize=(183 * MM, 58 * MM))
axes = [fig.add_axes([0.075 + i * 0.235, 0.30, 0.175, 0.55]) for i in range(4)]
for gi, gene in enumerate(GENES):
    ax = axes[gi]
    for ti, tr in enumerate(TREATS):
        on = (tr == f'si{gene}')
        st = S.SOLID if on else (S.OPEN if tr == 'siNC' else S.GREY)
        S.bar(ax, ti, qd[(tr, gene)], st, seed=gi * 10 + ti)
        if on:
            p = stats.ttest_ind(qd[(tr, gene)], qd[('siNC', gene)], equal_var=False)[1]
            S.bracket(ax, 0, ti, 122, 7, S.stars(p))
    ax.set_xticks(range(5)); ax.set_xticklabels(TREATS, rotation=45, ha='right')
    ax.set_ylim(0, 150); ax.set_yticks([0, 25, 50, 75, 100, 125])
    ax.set_title(f'{gene} mRNA', pad=2, fontstyle='italic')
    if gi:
        ax.set_yticklabels([])
    else:
        ax.set_ylabel('Relative mRNA level\n(% of siNC)')
    S.spines(ax)
h = [plt.Rectangle((0, 0), 1, 1, linewidth=S.LW, **s) for s in (S.OPEN, S.GREY, S.SOLID)]
fig.legend(h, ['siNC', 'non-targeting siRNA', 'on-target siRNA'], frameon=False,
           loc='lower center', bbox_to_anchor=(0.53, -0.03), ncol=3,
           handlelength=1.3, handleheight=0.9, columnspacing=1.6)
figletter(fig, axes[0], 'b')
save(fig, 'Fig11b_qPCR_nature')

# ══════════════════════════════════════════════════════════ Fig 11c  WB specificity
wa = rd('图11C_PanelA_WB敲低特异性_灰度定量.csv')
wv = defaultdict(dict)
for r in wa:
    if r['蛋白'] != 'GAPDH':
        wv[r['蛋白']][r['泳道']] = float(r['相对siNC%'])

fig = plt.figure(figsize=(150 * MM, 58 * MM))
ax = fig.add_axes([0.085, 0.34, 0.895, 0.48])
xs = []
for gi, gene in enumerate(GENES):
    for ti, tr in enumerate(TREATS):
        on = (tr == f'si{gene}')
        x = gi * 6.0 + ti
        xs.append(x)
        ax.bar(x, wv[gene][tr], width=.82, linewidth=S.LW, zorder=2,
               **(S.SOLID if on else (S.OPEN if tr == 'siNC' else S.GREY)))
    ax.text(gi * 6.0 + 2, 124, gene, ha='center', fontstyle='italic', fontsize=S.FS)
ax.set_xticks(xs)
ax.set_xticklabels(TREATS * 4, rotation=45, ha='right', fontsize=S.FS - 1)
ax.set_xlim(-1, 3 * 6.0 + 5)
ax.set_ylim(0, 136); ax.set_yticks([0, 25, 50, 75, 100, 125])
ax.axhline(100, color='black', lw=S.LW * 0.7, ls=(0, (3, 2)), zorder=1)
ax.set_ylabel('Relative protein level\n(% of siNC)')
S.spines(ax)
h = [plt.Rectangle((0, 0), 1, 1, linewidth=S.LW, **s) for s in (S.OPEN, S.GREY, S.SOLID)]
fig.legend(h, ['siNC', 'non-targeting siRNA', 'on-target siRNA'], frameon=False,
           loc='lower center', bbox_to_anchor=(0.53, -0.02), ncol=3,
           handlelength=1.3, handleheight=0.9, columnspacing=1.5)
figletter(fig, ax, 'c', dx=-0.065)
save(fig, 'Fig11c_WBspecificity_nature')

# ══════════════════════════════════════════════════════════ Fig 11d  CCK-8 + Bliss
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
    pack.append(dict(cl=cl, pt=pt, a=a, b=b, d=d, exp=a * b / 100, exc=exc,
                     p=stats.ttest_1samp(exc, 0)[1]))
    praw.append(pack[-1]['p'])
for k, qv in zip(pack, bh(praw)):
    k['q'] = qv

fig = plt.figure(figsize=(183 * MM, 74 * MM))
axL = fig.add_axes([0.062, 0.30, 0.315, 0.53])
axR = fig.add_axes([0.437, 0.30, 0.315, 0.53])
axS = fig.add_axes([0.845, 0.30, 0.125, 0.53])

for ax, cl in [(axL, 'A549'), (axR, 'H1299')]:
    for j, k in enumerate([x for x in pack if x['cl'] == cl]):
        base = j * 4.9
        S.bar(ax, base + 0, np.full(3, 100.), S.OPEN, width=.74, points=False)
        S.bar(ax, base + 1, k['a'], S.HATCH, width=.74, seed=base + 1)
        S.bar(ax, base + 2, k['b'], S.HATCH2, width=.74, seed=base + 2)
        S.bar(ax, base + 3, k['d'], S.SOLID if k['pt'] in CBAF else S.GREY,
              width=.74, seed=base + 3)
        e = k['exp'].mean()
        ax.plot([base + 2.55, base + 3.45], [e, e], color='black', lw=S.LW * 1.4,
                ls=(0, (2.2, 1.4)), zorder=6)
        gap = k['exc'].mean()
        lab = f'{gap:+.1f}'
        if abs(gap) >= 4:          # room inside the gap: arrow + label beside it
            ax.annotate('', xy=(base + 3, k['d'].mean()), xytext=(base + 3, e), zorder=6,
                        arrowprops=dict(arrowstyle='<->', color='black', lw=S.LW,
                                        shrinkA=0, shrinkB=0, mutation_scale=5))
            ax.text(base + 3.12, (e + k['d'].mean()) / 2 - 1, lab, fontsize=S.FS - 1.2,
                    va='center', ha='left')
        else:                      # gap too small: label above the bar
            ax.text(base + 3, max(e, k['d'].mean()) + 7, lab, fontsize=S.FS - 1.2,
                    va='bottom', ha='center')
        ax.text(base + 1.5, 137, f'+{k["pt"]}', ha='center', fontsize=S.FS,
                fontstyle='italic')
        ax.text(base + 1.5, 127, S.stars(k['q']), ha='center', fontsize=S.FS - 0.5)
    ax.set_xticks([j * 4.9 + i for j in range(3) for i in range(4)])
    ax.set_xticklabels(['siNC', 'siARID1A', 'si-partner', 'double'] * 3,
                       rotation=45, ha='right')
    ax.set_xlim(-0.8, 2 * 4.9 + 3 + 1.5)   # right margin for the excess labels
    ax.set_ylim(0, 150); ax.set_yticks([0, 20, 40, 60, 80, 100, 120])
    ax.axhline(100, color='black', lw=S.LW * 0.7, ls=(0, (3, 2)), zorder=1)
    ax.set_title(cl, pad=20)
    S.spines(ax)
axL.set_ylabel('Relative viability (% of siNC)')
axR.set_yticklabels([])

# summary panel: the pre-specified primary contrast
ec = np.concatenate([k['exc'] for k in pack if k['pt'] in CBAF])
ep = np.concatenate([k['exc'] for k in pack if k['pt'] not in CBAF])
pu = stats.mannwhitneyu(ec, ep, alternative='greater')[1]
for xi, (arr, st) in enumerate([(ec, S.SOLID), (ep, S.GREY)]):
    axS.bar(xi, arr.mean(), width=.62, linewidth=S.LW, zorder=2, **st)
    axS.errorbar(xi, arr.mean(), yerr=arr.std(ddof=1) / np.sqrt(len(arr)), ecolor='black',
                 elinewidth=S.LW, capsize=S.CAP, capthick=S.LW, zorder=3, fmt='none')
    rng = np.random.RandomState(7 + xi)
    axS.scatter(xi + rng.uniform(-.15, .15, len(arr)), arr, s=6.5, facecolor='white',
                edgecolor='black', linewidth=S.LW * .85, zorder=4, clip_on=False)
axS.axhline(0, color='black', lw=S.LW * 0.7, ls=(0, (3, 2)), zorder=1)
_m, _e = f'{pu:.1e}'.split('e')
S.bracket(axS, 0, 1, 25, 2.2, rf'$P = {_m}\times10^{{{int(_e)}}}$', fs=S.FS - 1)
axS.set_xticks([0, 1]); axS.set_xticklabels(['cBAF\npartner', 'PBAF\ncontrol'])
axS.set_ylim(-6, 34); axS.set_yticks([0, 10, 20, 30])
axS.set_ylabel('Excess over Bliss\n(percentage points)')
axS.set_title('Primary contrast', pad=20, fontsize=S.FS)
S.spines(axS)

h = [plt.Rectangle((0, 0), 1, 1, linewidth=S.LW, **s)
     for s in (S.OPEN, S.HATCH, S.HATCH2, S.SOLID, S.GREY)]
h.append(plt.Line2D([0], [0], color='black', lw=S.LW * 1.4, ls=(0, (2.2, 1.4))))
fig.legend(h, ['siNC', 'siARID1A', 'si-partner alone', 'double (cBAF partner)',
               'double (PBAF control)', 'Bliss additivity expectation'],
           frameon=False, loc='lower center', bbox_to_anchor=(0.42, -0.02), ncol=3,
           handlelength=1.3, handleheight=0.9, columnspacing=1.5)
for a_, l_ in [(axL, 'd'), (axS, 'e')]:
    figletter(fig, a_, l_, dx=-0.045)
save(fig, 'Fig11d_CCK8_nature')

# ══════════════════════════════════════════════════════════ Fig 11f  WB co-stability
wb = rd('图11E_PanelB_WB共稳定性_灰度定量.csv')
wd = defaultdict(lambda: defaultdict(list))
for r in wb:
    if r['蛋白'] == 'GAPDH':
        continue
    wd[r['蛋白']]['siNC' if 'siNC' in r['泳道'] else 'siARID1A'].append(float(r['相对siNC%']))

order = CONTACT + NOCONTACT
fig = plt.figure(figsize=(140 * MM, 68 * MM))
ax = fig.add_axes([0.115, 0.30, 0.855, 0.52])
for i, g in enumerate(order):
    S.bar(ax, i * 2.4 + 0, wd[g]['siNC'], S.OPEN, width=.92, seed=i * 3)
    S.bar(ax, i * 2.4 + 1, wd[g]['siARID1A'],
          S.SOLID if g in CONTACT else S.GREY, width=.92, seed=i * 3 + 1)
    p = stats.ttest_ind(wd[g]['siNC'], wd[g]['siARID1A'], equal_var=False)[1]
    S.bracket(ax, i * 2.4, i * 2.4 + 1, 126, 5, S.stars(p))
ax.set_xticks([i * 2.4 + .5 for i in range(len(order))])
ax.set_xticklabels(order, rotation=45, ha='right', fontstyle='italic')
ax.set_ylim(0, 168); ax.set_yticks([0, 25, 50, 75, 100, 125])
ax.axhline(100, color='black', lw=S.LW * 0.7, ls=(0, (3, 2)), zorder=1)
ax.set_ylabel('Relative protein level\n(% of siNC)')
sep = (len(CONTACT) - 1) * 2.4 + 1.7
ax.axvline(sep, color='black', lw=S.LW, ls=(0, (1.2, 1.2)), zorder=1)
ax.text(sep / 2 - 0.3, 150, 'direct contact with ARID1A', ha='center', fontsize=S.FS - .5)
ax.text(sep + (len(NOCONTACT) - 1) * 1.2 + 0.5, 150, 'no contact',
        ha='center', fontsize=S.FS - .5)
S.spines(ax)
h = [plt.Rectangle((0, 0), 1, 1, linewidth=S.LW, **s) for s in (S.OPEN, S.SOLID, S.GREY)]
fig.legend(h, ['siNC', 'siARID1A (contacting subunit)', 'siARID1A (non-contacting)'],
           frameon=False, loc='lower center', bbox_to_anchor=(0.54, 0.005), ncol=3,
           handlelength=1.3, handleheight=0.9, columnspacing=1.5)
figletter(fig, ax, 'f', dx=-0.085)
save(fig, 'Fig11f_WBcostability_nature')

# ══════════════════════════════════════════════════════════ report
print('written (600 dpi png + vector pdf), single source: 湿实验数据包_v2\n')
print('  Fig11b_qPCR_nature        knockdown 73.2-81.4%, Welch P 6.8e-6 to 1.2e-3')
print('  Fig11c_WBspecificity_nature  on-target 21.4-26.3%, off-target 96.3-112.9%')
print('  Fig11d_CCK8_nature        panels d (blocks) + e (primary contrast)')
for k in pack:
    print(f'      {k["cl"]+"/"+k["pt"]:14s} excess {k["exc"].mean():+5.1f}  q={k["q"]:.3f} {S.stars(k["q"])}')
print(f'      cBAF {ec.mean():.1f}+-{ec.std(ddof=1):.1f} (n={len(ec)})  vs  '
      f'PBAF {ep.mean():.1f}+-{ep.std(ddof=1):.1f} (n={len(ep)})   Mann-Whitney P={pu:.2e}')
print('  Fig11f_WBcostability_nature  contact 23-42% vs non-contact ~105%')
