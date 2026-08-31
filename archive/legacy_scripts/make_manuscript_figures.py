# -*- coding: utf-8 -*-
"""Compose the manuscript figures. All text in English; Nature typography via nfprism.

Figures produced (600 dpi PNG + vector PDF):
  Fig1  ARID1A is a recurrent loss-of-function event in NSCLC
  Fig2  Five unbiased methods converge
  Fig3  Structure explains and refines the module definition
  Fig4  Pan-cancer generalisation and independent replication
  Fig5  Computational structural cross-validation
  Fig6  Alternative explanations excluded
"""
import csv, os, sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nfprism as S

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.abspath(os.path.join(HERE, '..', '..'))
TAB = os.path.join(PROJ, 'results_tables')
OUT = os.path.join(PROJ, 'submission', 'figures')
os.makedirs(OUT, exist_ok=True)
S.apply()
MM = 1 / 25.4
GREY, DARK, MID = '#BFBFBF', 'black', '#7F7F7F'


def rd(fn):
    with open(os.path.join(TAB, fn), encoding='utf-8-sig') as f:
        rows = [r for r in csv.DictReader(f) if not (r.get(list(r)[0]) or '').startswith('#')]
    return rows


def save(fig, name):
    for e in ('png', 'pdf'):
        fig.savefig(os.path.join(OUT, f'{name}.{e}'), dpi=600, bbox_inches='tight')
    plt.close(fig)
    print('  wrote', name)


def panel(fig, ax, letter, dx=-0.055, dy=0.028):
    bb = ax.get_position()
    fig.text(bb.x0 + dx, bb.y1 + dy, letter, fontsize=9, fontweight='bold',
             ha='left', va='bottom')


# ═════════════════════════════════════════════ Figure 1
freq = rd('F1A_arid1a_mutation_frequency.csv')
spec = rd('F1A_arid1a_variant_spectrum.csv')

fig = plt.figure(figsize=(183 * MM, 52 * MM))
axA = fig.add_axes([0.055, 0.24, 0.17, 0.62])
axB = fig.add_axes([0.325, 0.24, 0.24, 0.62])
axC = fig.add_axes([0.680, 0.24, 0.29, 0.62])

# a: mutation frequency
labs = [r['cohort'] for r in freq]
vals = [float(r['freq_pct']) for r in freq]
ns = [f"{r['n_mut_patients']}/{r['n_total']}" for r in freq]
for i, (v, n) in enumerate(zip(vals, ns)):
    axA.bar(i, v, width=.62, facecolor='white' if labs[i] != 'NSCLC' else DARK,
            edgecolor='black', linewidth=S.LW, zorder=2)
    axA.text(i, v + 0.25, f'{v:.2f}%', ha='center', fontsize=S.FS - 1.5)
    axA.text(i, 0.35, n, ha='center', fontsize=S.FS - 2,
             color='white' if labs[i] == 'NSCLC' else 'black')
axA.set_xticks(range(len(labs))); axA.set_xticklabels(labs)
axA.set_ylim(0, 8.4); axA.set_ylabel('ARID1A mutation\nfrequency (%)')
S.spines(axA)

# b: variant spectrum, truncating vs missense
order = ['Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del',
         'Frame_Shift_Ins', 'Splice_Site']
cnt = {r['Variant_Classification']: int(r['N']) for r in spec}
tot = sum(cnt.values())
trunc = tot - cnt.get('Missense_Mutation', 0)
for i, k in enumerate(order):
    v = cnt.get(k, 0)
    axB.bar(i, v, width=.66, facecolor='white' if k == 'Missense_Mutation' else DARK,
            edgecolor='black', linewidth=S.LW, zorder=2)
    axB.text(i, v + 0.6, str(v), ha='center', fontsize=S.FS - 1.5)
axB.set_xticks(range(len(order)))
axB.set_xticklabels([k.replace('_Mutation', '').replace('Frame_Shift', 'Frameshift').replace('_', ' ')
                     for k in order], rotation=38, ha='right', fontsize=S.FS - 1)
axB.set_ylim(0, 40); axB.set_ylabel('Variants (n)')
axB.plot([0.6, 4.4], [37, 37], color='black', lw=S.LW, clip_on=False)
axB.text(2.5, 38, f'truncating {trunc}/{tot} = {100*trunc/tot:.1f}%',
         ha='center', fontsize=S.FS - 1)
S.spines(axB)

# c: DEG summary
rows = rd('fig1b_DEG_rerun_archivedStrat.csv')
lfc, padj = [], []
for r in rows:
    try:
        p = float(r['padj']); l = float(r['log2FoldChange'])
    except (ValueError, TypeError):
        continue
    if p > 0 and np.isfinite(l):
        lfc.append(l); padj.append(-np.log10(p))
lfc, padj = np.array(lfc), np.array(padj)
sig = (padj > -np.log10(0.05)) & (np.abs(lfc) > 1)
axC.scatter(lfc[~sig], padj[~sig], s=1.6, color=GREY, linewidth=0, rasterized=True)
axC.scatter(lfc[sig], padj[sig], s=2.4, color=DARK, linewidth=0, rasterized=True)
axC.axvline(1, color=MID, lw=S.LW * .8, ls=(0, (3, 2)))
axC.axvline(-1, color=MID, lw=S.LW * .8, ls=(0, (3, 2)))
axC.axhline(-np.log10(0.05), color=MID, lw=S.LW * .8, ls=(0, (3, 2)))
axC.set_xlim(-6, 6); axC.set_ylim(0, min(30, padj.max() * 1.05))
axC.set_xlabel('log$_2$ fold change (ARID1A-low vs -high)')
axC.set_ylabel('$-$log$_{10}$ adjusted $P$')
axC.text(.03, .95, f'{int(sig.sum())} DEGs\n421 up / 80 down', transform=axC.transAxes,
         va='top', fontsize=S.FS - 1)
S.spines(axC)
for a, l in ((axA, 'a'), (axB, 'b'), (axC, 'c')):
    panel(fig, a, l)
save(fig, 'Fig1')

# ═════════════════════════════════════════════ Figure 2
gw = rd('genomewide_codependency_lung.csv')[:400]
mod = rd('complex_module_codependency.csv')
pe = rd('panessentiality_specificity_control.csv')
shap = rd('ml_shap_arid1a_feature_importance.csv')[:12]

fig = plt.figure(figsize=(183 * MM, 108 * MM))
ax1 = fig.add_axes([0.065, 0.60, 0.24, 0.33])
ax2 = fig.add_axes([0.405, 0.60, 0.24, 0.33])
ax3 = fig.add_axes([0.735, 0.60, 0.235, 0.33])
ax4 = fig.add_axes([0.065, 0.09, 0.24, 0.33])
ax5 = fig.add_axes([0.405, 0.09, 0.24, 0.33])
ax6 = fig.add_axes([0.735, 0.09, 0.235, 0.33])

SWI = {'SMARCE1', 'SMARCB1', 'SMARCC1', 'DPF2', 'SMARCA4', 'ARID1B', 'SS18'}
rank = np.arange(1, len(gw) + 1)
rho = np.array([float(r['rho']) for r in gw])
names = [r['gene'] for r in gw]
ax1.scatter(rank, rho, s=3, color=GREY, linewidth=0)
# ranks 3-7 sit at nearly the same rho, so labels are fanned out with leader lines
FAN = {'SMARCE1': (60, .010), 'SMARCB1': (66, .008), 'SMARCC1': (86, .012),
       'DPF2': (104, -.005), 'SMARCA4': (128, -.022), 'ARID1B': (150, -.040),
       'SS18': (168, -.058)}
for i, n in enumerate(names):
    if n in SWI:
        ax1.scatter(rank[i], rho[i], s=16, facecolor=DARK, edgecolor='black', zorder=4)
        dx, dy = FAN.get(n, (16, .012))
        ax1.annotate(n, (rank[i], rho[i]), (rank[i] + dx, rho[i] + dy),
                     fontsize=S.FS - 2, color='black', va='center',
                     arrowprops=dict(arrowstyle='-', lw=S.LW * .7, color=MID,
                                     shrinkA=1, shrinkB=3))
ax1.set_xlabel('Genome-wide rank'); ax1.set_ylabel(r'Co-dependency $\rho$ with ARID1A')
ax1.set_xlim(0, 260); ax1.set_title('Genome-scale scan (lung)', fontsize=S.FS, pad=3)
S.spines(ax1)

mods = ['Core_shared', 'cBAF_specific', 'PBAF_specific', 'ncBAF_specific', 'Alt_ATPase']
lab = ['Core\nshared', 'cBAF\nspecific', 'PBAF\nspecific', 'ncBAF\nspecific', 'Alt.\nATPase']
lung = {r['module']: float(r['rho']) for r in mod if r['tissue'] == 'lung'}
panc = {r['module']: float(r['rho']) for r in mod if r['tissue'] == 'pancancer'}
x = np.arange(len(mods))
ax2.bar(x - .19, [lung[m] for m in mods], width=.36, facecolor='white',
        edgecolor='black', linewidth=S.LW, label='Lung', zorder=2)
ax2.bar(x + .19, [panc[m] for m in mods], width=.36, facecolor=DARK,
        edgecolor='black', linewidth=S.LW, label='Pan-cancer', zorder=2)
ax2.axhline(0, color='black', lw=S.LW)
ax2.set_xticks(x); ax2.set_xticklabels(lab, fontsize=S.FS - 1.5)
ax2.set_ylabel(r'Module $\rho$ with ARID1A'); ax2.set_ylim(-.25, .55)
ax2.legend(frameon=False, fontsize=S.FS - 1.5, loc='upper right', handlelength=1.2)
ax2.set_title('Complex-module scoring', fontsize=S.FS, pad=3)
S.spines(ax2)

terms = ['SWI/SNF Complex', 'nBAF Complex', 'EBAFb complex',
         'npBAF Complex', 'BAF complex']
adjp = [1.5e-12, 1.3e-11, 2.2e-10, 5.0e-10, 6.0e-10]
y = np.arange(len(terms))[::-1]
ax3.barh(y, [-np.log10(p) for p in adjp], height=.6, facecolor='white',
         edgecolor='black', linewidth=S.LW, zorder=2)
ax3.set_yticks(y); ax3.set_yticklabels(terms, fontsize=S.FS - 1.5)
ax3.set_xlabel(r'$-$log$_{10}$ adjusted $P$')
ax3.set_title('Prior-free enrichment (top-100)', fontsize=S.FS, pad=3)
S.spines(ax3)

com = [r['gene'] for r in rd('coessentiality_ARID1A_community_genes.csv')]
TAG = SWI | {'ARID1A'}
# the SWI/SNF members are adjacent in the alphabetical list, so they are moved to
# evenly spaced positions around the ring to keep their labels legible
swi_in = [g for g in com if g in TAG]
others = [g for g in com if g not in TAG]
step = len(com) // len(swi_in)
order, oi = [], 0
for i in range(len(com)):
    if i % step == 0 and len(order) // 1 < len(com) and swi_in and i // step < len(swi_in):
        order.append(swi_in[i // step])
    elif oi < len(others):
        order.append(others[oi]); oi += 1
order += others[oi:]
order = order[:len(com)]
ang = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, len(order), endpoint=False)
for i, g in enumerate(order):
    is_swi = g in TAG
    x, y = np.cos(ang[i]), np.sin(ang[i])
    ax4.scatter(x, y, s=26 if is_swi else 7, facecolor=DARK if is_swi else 'white',
                edgecolor='black', linewidth=S.LW * .8, zorder=3)
    if is_swi:
        ax4.annotate(g, (x * 1.22, y * 1.22), fontsize=S.FS - 2.5,
                     ha='center' if abs(x) < .35 else ('left' if x > 0 else 'right'),
                     va='center')
ax4.set_xlim(-1.85, 1.85); ax4.set_ylim(-1.5, 1.5); ax4.set_aspect('equal'); ax4.axis('off')
ax4.set_title(f'Leiden community: {len(com)} genes, {len(swi_in)} SWI/SNF\n'
              r'hypergeometric $P$ = 1.06 $\times$ 10$^{-13}$', fontsize=S.FS - 1, pad=1)

sel = [r for r in pe if r['gene'] in
       ('SMARCB1', 'SMARCE1', 'DPF2', 'SMARCA4', 'ARID1B', 'PBRM1',
        'RPL3', 'RPS3', 'PSMA1', 'EIF3B', 'RRM1')]
gn = [r['gene'] for r in sel]
rr = [float(r['arid1a_codep_rho'] or 0) for r in sel]
cls = [r['essentiality_class'] for r in sel]
ax5.barh(np.arange(len(gn))[::-1], rr, height=.62,
         color=[DARK if c != 'PAN-ESSENTIAL' else GREY for c in cls],
         edgecolor='black', linewidth=S.LW, zorder=2)
ax5.axvline(0, color='black', lw=S.LW)
ax5.set_yticks(np.arange(len(gn))[::-1]); ax5.set_yticklabels(gn, fontsize=S.FS - 1.5)
ax5.set_xlabel(r'$\rho$ with ARID1A'); ax5.set_xlim(-.1, .58)
ax5.set_title('Pan-essentiality control', fontsize=S.FS, pad=3)
h = [Rectangle((0, 0), 1, 1, facecolor=c, edgecolor='black', linewidth=S.LW)
     for c in (DARK, GREY)]
ax5.legend(h, ['BAF subunit', 'pan-essential control'], frameon=False,
           fontsize=S.FS - 2, loc='lower right', handlelength=1.1)
S.spines(ax5)

sg = [r['gene'] for r in shap][::-1]
sv = [float(r['mean_abs_shap']) for r in shap][::-1]
sr = [r['swisnf_role'] for r in shap][::-1]
ax6.barh(np.arange(len(sg)), sv, height=.66,
         color=[DARK if x != 'other' else GREY for x in sr],
         edgecolor='black', linewidth=S.LW, zorder=2)
ax6.set_yticks(np.arange(len(sg))); ax6.set_yticklabels(sg, fontsize=S.FS - 1.5)
ax6.set_xlabel('mean |SHAP| value')
ax6.set_title(r'XGBoost + SHAP (CV $R^2$ = 0.376)', fontsize=S.FS, pad=3)
S.spines(ax6)

for a, l in ((ax1, 'a'), (ax2, 'b'), (ax3, 'c'), (ax4, 'd'), (ax5, 'e'), (ax6, 'f')):
    panel(fig, a, l)
save(fig, 'Fig2')

# ═════════════════════════════════════════════ Figure 3
merged = rd('structural_vs_codependency_merged.csv')
loo = rd('core_module_leave_one_out.csv')
ref = rd('core_contact_module_refinement.csv')
prc2 = rd('prc2_refinement_CORRECTED.csv')

fig = plt.figure(figsize=(183 * MM, 100 * MM))
axA = fig.add_axes([0.070, 0.60, 0.35, 0.33])
axB = fig.add_axes([0.575, 0.60, 0.39, 0.33])
axC = fig.add_axes([0.070, 0.09, 0.35, 0.33])
axD = fig.add_axes([0.575, 0.09, 0.39, 0.33])

ct = np.array([float(r['n_contacts_4.5A']) for r in merged])
rl = np.array([float(r['lung_rho']) for r in merged])
gn = [r['gene'] for r in merged]
# contacting subunits get individual labels; the eight non-contacting ones all sit at
# x = 0 and would overprint, so they are listed once in a single annotation instead
nc_names = []
for i, g in enumerate(gn):
    has = ct[i] > 0
    axA.scatter(ct[i], rl[i], s=34, facecolor=DARK if has else 'white',
                edgecolor='black', linewidth=S.LW, zorder=3)
    if has:
        dy = -.050 if g == 'SMARCD1' else .017
        axA.annotate(g, (ct[i], rl[i]), (ct[i] + 24, rl[i] + dy), fontsize=S.FS - 2)
    else:
        nc_names.append(g)
axA.axhline(0, color=MID, lw=S.LW * .8, ls=(0, (3, 2)))
_ncy = float(np.median([rl[i] for i in range(len(gn)) if ct[i] == 0]))
axA.annotate('no contact: ' + ', '.join(nc_names[:4]) + ',\n' + ', '.join(nc_names[4:]),
             (10, _ncy), (200, -.205), fontsize=S.FS - 2, va='center',
             arrowprops=dict(arrowstyle='-', lw=S.LW * .8, color=MID,
                             shrinkA=2, shrinkB=6))
axA.set_xlabel('Contacting heavy-atom pairs in 6LTJ (< 4.5 A)')
axA.set_ylabel(r'Co-dependency $\rho$ (lung)')
axA.set_xlim(-70, 900)
axA.set_ylim(-.30, .68)
# two-sided, and over the seven contacting subunits (SMARCC1 added from 9A0K, which
# is why it is not among the plotted points -- 6LTJ does not resolve it)
axA.text(.97, .97, 'Mann$-$Whitney (two-sided)\n$P$ = 0.0037 (lung)\n$P$ = 0.0059 (pan-cancer)',
         transform=axA.transAxes, ha='right', va='top', fontsize=S.FS - 1.5)
h = [plt.Line2D([0], [0], marker='o', ls='', markerfacecolor=c, markeredgecolor='black',
                markersize=5) for c in (DARK, 'white')]
axA.legend(h, ['direct contact', 'no contact'], frameon=False, fontsize=S.FS - 1.5,
           loc='upper left', bbox_to_anchor=(0.0, 1.02), handlelength=.8)
S.spines(axA)

ll = [r for r in loo if r['tissue'] == 'lung']
ll = sorted(ll, key=lambda r: float(r['delta_vs_full8']))
gnl = [r['removed_gene'] for r in ll]
dl = [float(r['delta_vs_full8']) for r in ll]
cs = [r['contact_status'] for r in ll]
# the module here is the Table 1 composition, which contains SMARCD2 -- a subunit no
# cBAF structure resolves, so it needs its own category rather than being coloured
# as though its contact status were known
LIGHT = '#EDEDED'


def _col(c):
    if 'no_contact' in c:
        return GREY
    if 'not_resolved' in c:
        return LIGHT
    if 'paralog' in c:
        return 'white'
    return DARK


col = [_col(c) for c in cs]
axB.barh(np.arange(len(gnl)), dl, height=.66, color=col, edgecolor='black',
         linewidth=S.LW, zorder=2)
axB.axvline(0, color='black', lw=S.LW)
axB.set_yticks(np.arange(len(gnl))); axB.set_yticklabels(gnl, fontsize=S.FS - 1.5)
axB.set_xlabel(r'$\Delta\rho$ on removal from core module (lung)')
h = [Rectangle((0, 0), 1, 1, facecolor=c, edgecolor='black', linewidth=S.LW)
     for c in (DARK, 'white', GREY, LIGHT)]
axB.legend(h, ['contacts ARID1A', 'contacts, paralogue-redundant', 'no contact',
               'not resolved'],
           frameon=False, fontsize=S.FS - 2, loc='upper center',
           bbox_to_anchor=(0.5, -0.235), ncol=2, handlelength=1.1, columnspacing=1.0)
axB.set_title('Leave-one-out, SWI/SNF core', fontsize=S.FS, pad=3)
S.spines(axB)

lu = {(r['definition'][:14], r['tissue']): float(r['rho']) for r in ref}
pairs = [('Curated core\n(8 subunits)', 0.475, 0.499),
         ('Contact-restricted\n(6 subunits)', 0.584, 0.533)]
x = np.arange(2)
axC.bar(x - .19, [p[1] for p in pairs], width=.36, facecolor='white',
        edgecolor='black', linewidth=S.LW, label='Lung', zorder=2)
axC.bar(x + .19, [p[2] for p in pairs], width=.36, facecolor=DARK,
        edgecolor='black', linewidth=S.LW, label='Pan-cancer', zorder=2)
for i, p in enumerate(pairs):
    axC.text(i - .19, p[1] + .012, f'{p[1]:.3f}', ha='center', fontsize=S.FS - 2)
    axC.text(i + .19, p[2] + .012, f'{p[2]:.3f}', ha='center', fontsize=S.FS - 2)
axC.set_xticks(x); axC.set_xticklabels([p[0] for p in pairs], fontsize=S.FS - 1)
axC.set_ylabel(r'Core-module $\rho$ with ARID1A'); axC.set_ylim(0, .70)
axC.legend(frameon=False, fontsize=S.FS - 1.5, loc='upper left', handlelength=1.2)
axC.set_title('Structure-refined module (SWI/SNF)', fontsize=S.FS, pad=3)
S.spines(axC)

pr = defaultdict(dict)
for r in prc2:
    pr[r['anchor']][r['definition']] = float(r['rho'])
anchors = ['EZH2', 'SUZ12']
defs = [('curated', 'Curated\ncore'), ('contact_restricted', 'Contact-\nrestricted'),
        ('empirical_optimum', 'Exhaustive\noptimum')]
w = .26
for j, (dk, dl_) in enumerate(defs):
    axD.bar(np.arange(2) + (j - 1) * w, [pr[a][dk] for a in anchors], width=w * .92,
            facecolor=['white', GREY, DARK][j], edgecolor='black', linewidth=S.LW, zorder=2)
    for i, a in enumerate(anchors):
        axD.text(i + (j - 1) * w, pr[a][dk] + .012, f'{pr[a][dk]:.3f}',
                 ha='center', fontsize=S.FS - 2.5)
axD.set_xticks(range(2)); axD.set_xticklabels([f'anchor: {a}' for a in anchors])
axD.set_ylabel(r'Core-module $\rho$'); axD.set_ylim(0, .88)
h = [Rectangle((0, 0), 1, 1, facecolor=c, edgecolor='black', linewidth=S.LW)
     for c in ('white', GREY, DARK)]
axD.legend(h, [d[1].replace('\n', ' ') for d in defs], frameon=False,
           fontsize=S.FS - 2, ncol=3, loc='upper center', handlelength=1.1,
           columnspacing=1.0)
axD.set_title('Refinement generalises to PRC2 (PDB 6C23)', fontsize=S.FS, pad=3)
S.spines(axD)

for a, l in ((axA, 'a'), (axB, 'b'), (axC, 'c'), (axD, 'd')):
    panel(fig, a, l)
save(fig, 'Fig3')

# ═════════════════════════════════════════════ Figure 4
pc = rd('pancancer_partner_summary.csv')
dm = rd('depmap24q4_codependency_replication.csv')

fig = plt.figure(figsize=(183 * MM, 56 * MM))
axA = fig.add_axes([0.075, 0.26, 0.38, 0.60])
axB = fig.add_axes([0.595, 0.26, 0.37, 0.60])

pc = sorted(pc, key=lambda r: -int(r['n_tissues_sig_global_BH']))
gn = [r['partner'] for r in pc]
npos = [int(r['n_tissues_positive']) for r in pc]
nsig = [int(r['n_tissues_sig_global_BH']) for r in pc]
y = np.arange(len(gn))[::-1]
axA.barh(y, npos, height=.62, facecolor='white', edgecolor='black',
         linewidth=S.LW, zorder=2, label='positive direction')
axA.barh(y, nsig, height=.62, facecolor=DARK, edgecolor='black',
         linewidth=S.LW, zorder=3, label='passes global BH')
axA.set_yticks(y); axA.set_yticklabels(gn, fontsize=S.FS - 1)
axA.set_xlabel('Tissues (of 19)'); axA.set_xlim(0, 20)
axA.legend(frameon=False, fontsize=S.FS - 1.5, loc='lower right', handlelength=1.1)
axA.set_title('Within-tissue analysis, 930 lines / 19 tissues', fontsize=S.FS, pad=3)
S.spines(axA)

gn2 = [r['partner'] for r in dm]
ps = [float(r['projectscore_lung_rho'] or 'nan') for r in dm]
d24 = [float(r['depmap24q4_lung_rho']) for r in dm]
ok = [i for i, v in enumerate(ps) if np.isfinite(v)]
axB.scatter([ps[i] for i in ok], [d24[i] for i in ok], s=30, facecolor=DARK,
            edgecolor='black', linewidth=S.LW, zorder=3)
for i in ok:
    axB.annotate(gn2[i], (ps[i], d24[i]), (ps[i] + .018, d24[i] - .028), fontsize=S.FS - 2)
lim = [-.2, .65]
axB.plot(lim, lim, color=MID, lw=S.LW, ls=(0, (3, 2)), zorder=1)
axB.axhline(0, color=MID, lw=S.LW * .7); axB.axvline(0, color=MID, lw=S.LW * .7)
axB.set_xlim(lim); axB.set_ylim(lim)
axB.set_xlabel(r'Project Score $\rho$ (lung)'); axB.set_ylabel(r'DepMap 24Q4 $\rho$ (lung)')
axB.set_title('Independent dataset replication', fontsize=S.FS, pad=3)
S.spines(axB)
for a, l in ((axA, 'a'), (axB, 'b')):
    panel(fig, a, l)
save(fig, 'Fig4')

# ═════════════════════════════════════════════ Figure 5
af = rd('af_interface_vs_cryoem_6partners_RECOMPUTED.csv')
ld = rd('lightdock_200swarm_bestoverall_verified.csv')

fig = plt.figure(figsize=(183 * MM, 56 * MM))
axA = fig.add_axes([0.075, 0.26, 0.36, 0.60])
axB = fig.add_axes([0.585, 0.26, 0.38, 0.60])

pn = [r['partner'] for r in ld]
bj = [float(r['best_jaccard_over_all_poses']) for r in ld]
nn = [int(r['n_near_native_J>=0.2']) for r in ld]
x = np.arange(len(pn))
axA.bar(x, bj, width=.6, facecolor='white', edgecolor='black', linewidth=S.LW, zorder=2)
for i in range(len(pn)):
    axA.text(i, bj[i] + .012, f'{bj[i]:.3f}', ha='center', fontsize=S.FS - 1.5)
    axA.text(i, .02, f'{nn[i]} poses', ha='center', fontsize=S.FS - 2.5)
axA.set_xticks(x); axA.set_xticklabels(pn, rotation=30, ha='right')
axA.set_ylabel('Best interface Jaccard\n(6,000 poses per partner)')
axA.set_ylim(0, .60)
axA.set_title('Blind docking recovers the interface', fontsize=S.FS, pad=3)
S.spines(axA)

ip = [float(r['af_iptm']) for r in af]
jc = [float(r['jaccard']) for r in af]
pn2 = [r['partner'] for r in af]
axB.scatter(ip, jc, s=34, facecolor=DARK, edgecolor='black', linewidth=S.LW, zorder=3)
for i, p in enumerate(pn2):
    axB.annotate(p, (ip[i], jc[i]), (ip[i] + .009, jc[i] + .018), fontsize=S.FS - 2)
axB.set_xlabel('AlphaFold-Multimer ipTM (self-reported confidence)')
axB.set_ylabel('Interface Jaccard vs cryo-EM')
axB.set_xlim(.36, .84); axB.set_ylim(-.05, .72)
axB.text(.03, .95, r'Spearman $\rho$ = $-$0.14, $P$ = 0.79' '\nconfidence does not track accuracy',
         transform=axB.transAxes, va='top', fontsize=S.FS - 1.5)
axB.set_title('ipTM is not a reliability filter', fontsize=S.FS, pad=3)
S.spines(axB)
for a, l in ((axA, 'a'), (axB, 'b')):
    panel(fig, a, l)
save(fig, 'Fig5')

# ═════════════════════════════════════════════ Figure 6
mr = rd('mr_arid1a_lungcancer_results.csv')
prot = rd('proteomics_arid1a_baf_correlation_REBUILT.csv')

fig = plt.figure(figsize=(183 * MM, 54 * MM))
axA = fig.add_axes([0.075, 0.28, 0.36, 0.58])
axB = fig.add_axes([0.590, 0.28, 0.375, 0.58])

lungp = [r for r in prot if r.get('tissue') == 'lung']
lungp = sorted(lungp, key=lambda r: -float(r['rho']))
gp = [r['partner'] for r in lungp]
rp = [float(r['rho']) for r in lungp]
CORE = {'SMARCA4', 'SMARCB1', 'SMARCC1', 'SMARCC2', 'SMARCE1', 'SMARCD1', 'DPF2', 'ARID1B'}
axA.barh(np.arange(len(gp))[::-1], rp, height=.66,
         color=[DARK if g in CORE else GREY for g in gp],
         edgecolor='black', linewidth=S.LW, zorder=2)
axA.axvline(0, color='black', lw=S.LW)
axA.set_yticks(np.arange(len(gp))[::-1]); axA.set_yticklabels(gp, fontsize=S.FS - 1.5)
axA.set_xlabel(r'Protein-level $\rho$ with ARID1A (lung)')
h = [Rectangle((0, 0), 1, 1, facecolor=c, edgecolor='black', linewidth=S.LW)
     for c in (DARK, GREY)]
axA.legend(h, ['core / cBAF', 'PBAF / ncBAF / alt. ATPase'], frameon=False,
           fontsize=S.FS - 2, loc='lower right', handlelength=1.1)
axA.set_title('All BAF subunits correlate positively\nat the protein level',
              fontsize=S.FS - .5, pad=3)
S.spines(axA)

lab = ['McKay 2017\n29,266 / 56,450', 'TRICL 2014\n11,348 / 15,861']
orv = [float(r['or']) for r in mr]
lo = [float(r['or_lci95']) for r in mr]
hi = [float(r['or_uci95']) for r in mr]
pv = [float(r['pval']) for r in mr]
yy = np.arange(len(lab))[::-1]
for i in range(len(lab)):
    axB.plot([lo[i], hi[i]], [yy[i], yy[i]], color='black', lw=S.LW * 1.6,
             solid_capstyle='butt', zorder=3)
    axB.plot([lo[i], lo[i]], [yy[i] - .09, yy[i] + .09], color='black', lw=S.LW * 1.6)
    axB.plot([hi[i], hi[i]], [yy[i] - .09, yy[i] + .09], color='black', lw=S.LW * 1.6)
    axB.scatter(orv[i], yy[i], s=42, marker='s', facecolor=DARK, edgecolor='black',
                linewidth=S.LW, zorder=4)
    axB.text(1.10, yy[i], f'OR {orv[i]:.3f} ({lo[i]:.3f}$-${hi[i]:.3f})\n$P$ = {pv[i]:.3f}',
             va='center', fontsize=S.FS - 1.5)
axB.axvline(1, color=MID, lw=S.LW, ls=(0, (3, 2)))
axB.set_yticks(yy); axB.set_yticklabels(lab, fontsize=S.FS - 1)
axB.set_xlim(.93, 1.26); axB.set_ylim(-.6, len(lab) - .25)
axB.set_xlabel('Odds ratio for lung cancer per unit ARID1A expression')
axB.set_title('Germline MR excludes genetic confounding', fontsize=S.FS - .5, pad=3)
S.spines(axB)
for a, l in ((axA, 'a'), (axB, 'b')):
    panel(fig, a, l)
save(fig, 'Fig6')

print('\nall figures written to', OUT)
