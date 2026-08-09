# -*- coding: utf-8 -*-
"""Redraw the single-cell context figure directly from the recovered source table.

The archived render (fig_scrna_arid1a.png) was made from a filtered subset whose
filter could not be reconstructed -- its per-cell-type counts run about 10% below
gse131907_harmony.csv across every category. Rather than reuse a figure whose
numbers cannot be checked, every panel here is recomputed from the table, so the
counts in the caption are the counts in the file.

Source: /data/s01011/arid1a-nsclc/results/tables/gse131907_harmony.csv (sxycloud),
206,826 Harmony-integrated cells from GSE131907 with a pre-computed core/cBAF
module score (column coreBAF).
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.abspath(os.path.join(HERE, '..', '..'))
TAB = os.path.join(PROJ, 'results_tables')
OUT = os.path.join(PROJ, 'submission', 'figures')
MM = 1 / 25.4

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 7, 'axes.linewidth': 0.75,
    'xtick.major.width': 0.75, 'ytick.major.width': 0.75,
    'xtick.major.size': 2.2, 'ytick.major.size': 2.2, 'pdf.fonttype': 42,
})

WARM = LinearSegmentedColormap.from_list(
    'warm', ['#F2F2F2', '#F6C9A8', '#E8825A', '#C0392B', '#8E1B18'])

TYPE_COLOUR = {
    'T lymphocytes': '#1F4E9C', 'Myeloid cells': '#F08A5D',
    'Epithelial cells': '#B8232E', 'B lymphocytes': '#5BA3D9',
    'NK cells': '#7B5EA7', 'Fibroblasts': '#199A62',
    'MAST cells': '#D9A404', 'Endothelial cells': '#E87BA8',
    'Oligodendrocytes': '#9AA0A6', 'Undetermined': '#C9CCD1',
}

d = pd.read_csv(os.path.join(TAB, 'gse131907_harmony.csv'))
print(f'{len(d):,} cells')

# thin the point cloud for rendering only; every printed number uses all cells
rng = np.random.default_rng(0)
idx = rng.permutation(len(d))
plot = d.iloc[idx[:90000]]

fig = plt.figure(figsize=(183 * MM, 150 * MM))


def scatter(ax, c, cmap=None, vmin=None, vmax=None):
    return ax.scatter(plot.UMAP1, plot.UMAP2, c=c, s=0.35, linewidths=0,
                      cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True)


def strip(ax):
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel('UMAP-1', fontsize=6.5)
    ax.set_ylabel('UMAP-2', fontsize=6.5)
    for s in ax.spines.values():
        s.set_linewidth(0.6)


# ---- A: cell types --------------------------------------------------------
axA = fig.add_axes([0.045, 0.585, 0.40, 0.375])
axA.scatter(plot.UMAP1, plot.UMAP2, s=0.35, linewidths=0, rasterized=True,
            c=[TYPE_COLOUR.get(t, '#C9CCD1') for t in plot.Cell_type])
strip(axA)
axA.set_title('a   Cell types, Harmony-integrated', fontsize=8,
              fontweight='bold', loc='left', pad=4)
for t in ['Epithelial cells', 'T lymphocytes', 'Myeloid cells', 'B lymphocytes']:
    sub = plot[plot.Cell_type == t]
    if len(sub):
        axA.text(sub.UMAP1.median(), sub.UMAP2.median(), t.split()[0],
                 fontsize=6.5, fontweight='bold', ha='center', va='center',
                 bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none',
                           alpha=0.75))

# ---- B: ARID1A expression -------------------------------------------------
axB = fig.add_axes([0.545, 0.585, 0.40, 0.375])
sB = scatter(axB, plot.ARID1A, cmap=WARM)
strip(axB)
axB.set_title('b   ARID1A expression', fontsize=8, fontweight='bold',
              loc='left', pad=4)
cb = fig.colorbar(sB, ax=axB, fraction=0.036, pad=0.02)
cb.ax.tick_params(labelsize=6, width=0.6, length=1.8)
cb.outline.set_linewidth(0.6)

# ---- C: core/cBAF module score -------------------------------------------
axC = fig.add_axes([0.045, 0.105, 0.40, 0.375])
lo_q, hi_q = d.coreBAF.quantile([0.01, 0.99])
sC = scatter(axC, plot.coreBAF, cmap=WARM, vmin=lo_q, vmax=hi_q)
strip(axC)
axC.set_title('c   Core/cBAF module score', fontsize=8, fontweight='bold',
              loc='left', pad=4)
cb = fig.colorbar(sC, ax=axC, fraction=0.036, pad=0.02)
cb.ax.tick_params(labelsize=6, width=0.6, length=1.8)
cb.outline.set_linewidth(0.6)

# ---- D: module score in ARID1A-low vs -high epithelial cells --------------
ep = d[d.Cell_type == 'Epithelial cells']
med = ep.ARID1A.median()
lo = ep[ep.ARID1A <= med].coreBAF.dropna()
hi = ep[ep.ARID1A > med].coreBAF.dropna()
U, p = stats.mannwhitneyu(lo, hi)

axD = fig.add_axes([0.575, 0.105, 0.37, 0.375])
parts = axD.violinplot([lo.values, hi.values], positions=[0, 1], widths=0.75,
                       showextrema=False, showmedians=False)
for body, col in zip(parts['bodies'], ['#B8232E', '#2E74B5']):
    body.set_facecolor(col); body.set_alpha(0.80)
    body.set_edgecolor('black'); body.set_linewidth(0.6)
for x, v in zip([0, 1], [lo, hi]):
    axD.plot(x, v.median(), 'o', ms=4, mfc='white', mec='black', mew=0.75, zorder=5)
axD.set_xticks([0, 1])
axD.set_xticklabels([f'ARID1A-low\n(n = {len(lo):,})', f'ARID1A-high\n(n = {len(hi):,})'],
                    fontsize=6.5)
axD.set_ylabel('core/cBAF module score', fontsize=6.5)
axD.set_ylim(lo_q - 0.05, hi_q + 0.05)
axD.set_title('d   Residual core/cBAF is retained when ARID1A is low',
              fontsize=8, fontweight='bold', loc='left', pad=4)
axD.spines[['top', 'right']].set_visible(False)
axD.text(0.5, 0.03, f'medians {lo.median():+.3f} vs {hi.median():+.3f}',
         transform=axD.transAxes, ha='center', fontsize=6, color='#555555')

# ---- legend ---------------------------------------------------------------
counts = d.Cell_type.value_counts()
handles = [plt.Line2D([], [], marker='o', ls='', ms=3.5,
                      mfc=TYPE_COLOUR.get(t, '#C9CCD1'), mec='none',
                      label=f'{t} ({counts[t]/1000:.0f}k)')
           for t in counts.index if counts[t] >= 1000]
fig.legend(handles=handles, loc='center', bbox_to_anchor=(0.5, 0.535),
           ncol=5, frameon=False, fontsize=6.3, handletextpad=0.4,
           columnspacing=1.2)

for ext in ('pdf', 'png'):
    fig.savefig(os.path.join(OUT, f'Fig_scrna.{ext}'), dpi=450,
                bbox_inches='tight', pad_inches=0.02)
plt.close(fig)

print(f'epithelial {len(ep):,}  low {len(lo):,}  high {len(hi):,}')
print(f'medians {lo.median():.4f} vs {hi.median():.4f}   Mann-Whitney P = {p:.3g}')
print('written Fig_scrna.pdf')
