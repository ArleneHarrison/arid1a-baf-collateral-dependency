# -*- coding: utf-8 -*-
"""Redraw the MD validation figure with the FHD-286 trajectory run to completion.

The archived render (deep_md_validation.png) shows FHD-286 stopping at 9 ns of a
planned 20 ns, which is visible as a truncated trace. The trajectory has now been
run to 20 ns on an RTX 5090 and reanalysed with the identical cpptraj definitions
(same pocket residue list, same fit order, same ligand selection), so the figure
is rebuilt from the new data. PFI-3 is unchanged at 30 ns and reproduces the
published values exactly, which is the control that the pipeline is the same one.

Layout, palette and annotation follow the original figure.
"""
import io
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from scipy.stats import gaussian_kde

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'md20')

ORANGE = '#E08C40'      # PFI-3, sampled from the original legend
CRIMSON = '#C1344E'     # FHD-286
INK = '#161B26'
GREY = '#8A8F98'
DT_PS = 10.0            # DCD written every 10 ps

# MM/GBSA decomposition, kcal/mol. PFI-3 from the archived 30 ns run;
# FHD-286 filled in from the 20 ns rerun.
MMGBSA = {
    'PFI-3': dict(vdw=-35.4301, eel=-18.2489, egb=30.4957, esurf=-4.7499,
                  total=-27.9332, sd=2.3344, ns=30),
    'FHD-286': dict(vdw=-53.3940, eel=-45.1635, egb=72.6827, esurf=-7.5737,
                    total=-33.4484, sd=4.2892, ns=20),
}


def load(name):
    d = np.loadtxt(os.path.join(D, name), comments='#')
    return d[:, 0] * DT_PS / 1000.0, d[:, 1]      # ns, Angstrom


def running(v, w=25):
    """Centred running mean with reflected edges.

    A plain same-mode convolution zero-pads, which drags the first and last w/2
    points toward zero and puts a false cliff at the end of every trace. Reflect
    instead so the ends stay on the data.
    """
    pad = w // 2
    ext = np.concatenate([v[pad:0:-1], v, v[-2:-pad - 2:-1]])
    out = np.convolve(ext, np.ones(w) / w, mode='same')
    return out[pad:pad + len(v)]


def trace(ax, axd, t, v, colour, label):
    ax.plot(t, v, color=colour, alpha=0.22, lw=0.6)
    ax.plot(t, running(v), color=colour, lw=1.9, label=label)
    ax.axhline(v.mean(), color=colour, ls='--', lw=1.1, alpha=0.85)
    ys = np.linspace(v.min(), v.max(), 200)
    dens = gaussian_kde(v)(ys)
    axd.fill_betweenx(ys, 0, dens, color=colour, alpha=0.22)
    axd.plot(dens, ys, color=colour, lw=1.5)
    axd.plot([gaussian_kde(v)(v.mean())[0]], [v.mean()], 'o', color=colour, ms=4)


def build(out_png, out_pdf, bare=False):
    """bare=True drops the title and footnote: the manuscript carries those in the
    LaTeX caption, which is why the archived pipeline cropped them off."""
    fig = plt.figure(figsize=(20.0, 5.2 if bare else 6.2))
    gs = GridSpec(1, 5, figure=fig, width_ratios=[5.0, 0.9, 5.0, 0.9, 6.4],
                  wspace=0.06, left=0.045, right=0.985,
                  top=0.92 if bare else 0.80, bottom=0.13 if bare else 0.27)
    axA, axAd = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    axB, axBd = fig.add_subplot(gs[2]), fig.add_subplot(gs[3])
    axC = fig.add_subplot(gs[4])

    series = [
        ('PFI-3', ORANGE, 'pfi3_rmsd_ligand_pocketfit.dat', 'pfi3_rmsd_pocketbb.dat'),
        ('FHD-286', CRIMSON, 'fhd_rmsd_ligand_pocketfit.dat', 'fhd_rmsd_pocketbb.dat'),
    ]
    for lab, col, flig, fbb in series:
        t, v = load(flig)
        trace(axA, axAd, t, v, col, lab)
        t2, v2 = load(fbb)
        trace(axB, axBd, t2, v2, col, lab)
        print(f'  {lab:8s} ligand {v.mean():.3f} +/- {v.std():.3f} over {t[-1]:.0f} ns'
              f' | pocket bb {v2.mean():.3f} +/- {v2.std():.3f}')

    for ax, ylab in ((axA, 'ligand RMSD (Å)'),
                     (axB, 'binding-site backbone RMSD (Å)')):
        ax.set_xlabel('time (ns)', fontsize=12)
        ax.set_ylabel(ylab, fontsize=12)
        ax.set_ylim(bottom=0)
        ax.grid(axis='y', color='#E3E5E8', lw=0.8)
        ax.set_axisbelow(True)
        ax.legend(frameon=False, fontsize=11, loc='upper left')
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
    for axd in (axAd, axBd):
        axd.set_xticks([])
        axd.set_yticks([])
        for s in axd.spines.values():
            s.set_visible(False)
    axAd.set_ylim(axA.get_ylim())
    axBd.set_ylim(axB.get_ylim())
    axA.text(0.99, 0.02, 'raw + running mean · right: density', transform=axA.transAxes,
             ha='right', va='bottom', fontsize=10, style='italic', color=GREY)

    # ---- panel C, MM/GBSA decomposition
    terms = [('vdW', 'vdw'), ('electrostatic', 'eel'),
             ('polar solv.', 'egb'), ('nonpolar solv.', 'esurf'),
             (r'$\Delta G_{bind}$ (net)', 'total')]
    h, gap = 0.34, 0.92
    for i, (name, key) in enumerate(terms):
        y0 = -i * gap
        for j, (lab, col) in enumerate((('PFI-3', ORANGE), ('FHD-286', CRIMSON))):
            val = MMGBSA[lab][key]
            if val is None:
                continue
            y = y0 - j * h
            axC.barh(y, val, height=h * 0.86, color=col,
                     label=lab if i == 0 else None)
            off = 2.5 if val >= 0 else -2.5
            axC.text(val + off, y, f'{val:+.0f}', va='center',
                     ha='left' if val >= 0 else 'right', fontsize=10.5,
                     color=INK, fontweight='bold' if key == 'total' else 'normal')
    vals = [MMGBSA[l][k] for _, k in terms for l in MMGBSA if MMGBSA[l][k] is not None]
    axC.set_xlim(min(vals) * 1.30, max(vals) * 1.30)     # room for the end labels
    axC.axhspan(-4 * gap - 2 * h, -4 * gap + h * 0.7, color='#F2F3F5', zorder=0)
    axC.set_yticks([-i * gap - h / 2 for i in range(len(terms))])
    axC.set_yticklabels([n for n, _ in terms], fontsize=12)
    axC.axvline(0, color=GREY, lw=1.1)
    axC.set_xlabel('energy (kcal/mol)   ← favourable      unfavourable →',
                   fontsize=11.5)
    axC.legend(frameon=False, fontsize=11, loc='upper right')
    for s in ('top', 'right', 'left'):
        axC.spines[s].set_visible(False)
    axC.tick_params(axis='y', length=0)

    if not bare:
        fig.suptitle('Molecular dynamics: both ligands remain stably bound in the '
                     'residual SWI/SNF pockets',
                     fontsize=17, fontweight='bold', color=INK, x=0.045,
                     ha='left', y=0.955)
    for k, (ax, t) in enumerate(((axA, 'A · Ligand stays in the pocket'),
                                 (axB, 'B · Binding site stays rigid'),
                                 (axC, 'C · MM/GBSA decomposition'))):
        ax.set_title(t, fontsize=14, fontweight='bold', color=INK, loc='left', pad=12)

    foot = ('All-atom explicit-solvent MD (Amber ff14SB / GAFF2 / TIP3P, OpenMM, HMR 4 fs). '
            'PFI-3 · SMARCA2 bromodomain: 30 ns from the 5DKC crystal pose;   '
            'FHD-286 · SMARCA4 ATPase: 20 ns from the AutoDock-Vina docked pose '
            '(110k-atom system).\n'
            'RMSD is fit on binding-site residues to isolate the local pocket from '
            'whole-domain motion. Both ligand RMSDs plateau low (dashed = mean) → the '
            'poses are dynamically stable;\nMM/GBSA (C) is single-trajectory GB (igb=2), '
            '0.15 M salt — relative binding strength, not absolute affinity.')
    if not bare:
        fig.text(0.045, 0.155, foot, fontsize=10.2, color=GREY, style='italic',
                 va='top')

    if out_png:
        fig.savefig(out_png, dpi=300, facecolor='white')
    if out_pdf:
        fig.savefig(out_pdf, facecolor='white')
    plt.close(fig)
    print(f'  written: {out_png or ""} {out_pdf or ""}')


if __name__ == '__main__':
    FIGDIR = (r'D:\桌面\_按项目整理\02_生信数据与医学AI项目\07_医学AI_深度学习_多项目孵化'
              r'\医学AI深度学习_导师课题项目群\范教授课题_NSCLC-SEER_分子对接多项目\纯生信'
              r'\ARID1A_BAF_最终项目\results_figures')
    SUB = (r'D:\桌面\_按项目整理\02_生信数据与医学AI项目\07_医学AI_深度学习_多项目孵化'
           r'\医学AI深度学习_导师课题项目群\范教授课题_NSCLC-SEER_分子对接多项目\纯生信'
           r'\ARID1A_BAF_最终项目\submission\figures')
    # standalone version, with title and footnote, for the figure archive
    build(os.path.join(FIGDIR, 'deep_md_validation.png'), None)
    # manuscript version: title and footnote live in the LaTeX caption
    build(os.path.join(SUB, 'Fig10_md.png'),
          os.path.join(SUB, 'Fig10_md.pdf'), bare=True)
