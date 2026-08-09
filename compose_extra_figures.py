# -*- coding: utf-8 -*-
"""Bring the existing structural / docking / MD renders into the submission set.

Only figures whose analysis is actually described in the manuscript are included.
The sources are already English-labelled and already carry their own internal panel
letters, so each is kept as a self-contained figure rather than being nested inside
a composite (which would produce conflicting double lettering). Baked-in titles and
footnotes are cropped away because the LaTeX caption carries that information.
"""
import os, sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.abspath(os.path.join(HERE, '..', '..'))
SRC = os.path.join(PROJ, 'results_figures')
OUT = os.path.join(PROJ, 'submission', 'figures')
os.makedirs(OUT, exist_ok=True)
MM = 1 / 25.4


def emit(src, dst, width_mm=183, top=0.0, bottom=0.0, left=0.0, right=1.0):
    im = mpimg.imread(os.path.join(SRC, src))
    h, w = im.shape[:2]
    im = im[int(h * top):int(h * (1 - bottom)), int(w * left):int(w * right)]
    W = width_mm * MM
    fig = plt.figure(figsize=(W, W * im.shape[0] / im.shape[1]))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(im, interpolation='lanczos')
    ax.axis('off')
    for e in ('png', 'pdf'):
        fig.savefig(os.path.join(OUT, f'{dst}.{e}'), dpi=450,
                    bbox_inches='tight', pad_inches=0.04)
    plt.close(fig)
    print(f'  {dst:22s} <- {src}')


# cryo-EM render: drop the baked-in title strip and the footnote line
emit('fig4_cryoem_white.png', 'Fig3_cryoem', 150, top=0.055, bottom=0.055)

# molecular docking summary: keeps its internal A-D panels
emit('fig6_docking.png', 'Fig8_docking', 165, top=0.045, bottom=0.055)

# PFI-3 binding pose: drop title, subtitle and the stray corner letter
emit('fig_pocket_pfi3_3dmol.png', 'Fig9_pose', 165, top=0.105, bottom=0.015, left=0.035)

# molecular dynamics validation: keeps its internal A-C panels, drop the method footnote
# bottom crop lands in the blank band at rows 1094-1150 so the x-axis labels
# survive and only the baked-in footnote is removed
emit('deep_md_validation.png', 'Fig10_md', 183, top=0.057, bottom=0.068)

# genotype-stratified reconciliation with paralogue synthetic lethality
# 0.045 cut through the middle of the baked-in title; 0.075 is the blank band
# below it, so the title is removed cleanly (the LaTeX caption carries it)
emit('figD_genotype_helming.png', 'Fig11_genotype', 125, top=0.075)

print('\nwritten to', OUT)
