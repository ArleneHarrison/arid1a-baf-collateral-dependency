# -*- coding: utf-8 -*-
"""Build the supplementary figure set for the submission package.

Every analysis figure produced during the project that is (i) not already a main
figure, (ii) not wet-lab, and (iii) not a Chinese-language duplicate of a main
figure is brought in here. Related sources are stacked into one supplementary
figure, grouped by what they show rather than by when they were made.

Each source figure already carries its own internal panel letters (A, B, C ...).
Stacked units are therefore tagged with roman numerals (I, II ...) so that the
two levels never collide -- the caption refers to "I-A", "II-B" and so on.
Baked-in titles and footnotes are cropped because the LaTeX caption carries that
information, and because several of them are in Chinese.

The AlphaFold gallery is the one exception: its labels are Chinese, so the six
structure renders are re-cut from the source and re-labelled in English from
af_interface_vs_cryoem_6partners_RECOMPUTED.csv.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.abspath(os.path.join(HERE, '..', '..'))
SRC = os.path.join(PROJ, 'results_figures')
TAB = os.path.join(PROJ, 'results_tables')
OUT = os.path.join(PROJ, 'submission', 'figures')
os.makedirs(OUT, exist_ok=True)

MM = 1 / 25.4
TEXT_MM = 162          # \textwidth for a4paper with 2.4 cm margins
plt.rcParams.update({'font.family': 'DejaVu Sans', 'pdf.fonttype': 42})


def load(name, top=0.0, bottom=0.0, left=0.0, right=1.0, base=None):
    """read a source render and crop away its baked-in title / footnote bands"""
    im = mpimg.imread(os.path.join(base or SRC, name))
    h, w = im.shape[:2]
    return im[int(h * top):int(h * (1 - bottom)), int(w * left):int(w * right)]


def built(name, **crop):
    """re-read a figure this script already emitted, for use inside a larger one"""
    return load(name, base=OUT, **crop)


def emit_single(src, dst, width_mm=TEXT_MM, **crop):
    """pass a single source figure through unchanged apart from the crop"""
    im = load(src, **crop)
    W = width_mm * MM
    fig = plt.figure(figsize=(W, W * im.shape[0] / im.shape[1]))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(im, interpolation='lanczos')
    ax.axis('off')
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'{dst}.{ext}'), dpi=400,
                    bbox_inches='tight', pad_inches=0.04)
    plt.close(fig)
    print(f'  {dst:16s} <- {src}')


def stack(units, dst, width_mm=TEXT_MM, gap=0.012, pad_left=0.030):
    """stack cropped panels vertically, tagging each with a roman numeral

    units: list of (roman, image). Heights are allocated in proportion to each
    panel's own aspect ratio so nothing is distorted.
    """
    ars = [im.shape[0] / im.shape[1] for _, im in units]           # h/w per unit
    body = 1.0 - pad_left
    heights = [a * body for a in ars]                              # in width units
    total = sum(heights) + gap * (len(units) - 1)

    W = width_mm * MM
    fig = plt.figure(figsize=(W, W * total))
    y = 1.0
    for (roman, im), hh in zip(units, heights):
        fh = hh / total
        y -= fh
        ax = fig.add_axes([pad_left, y, body, fh])
        ax.imshow(im, interpolation='lanczos')
        ax.axis('off')
        if roman:
            fig.text(0.004, y + fh - 0.004, roman, ha='left', va='top',
                     fontsize=13, fontweight='bold')
        y -= gap / total

    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'{dst}.{ext}'), dpi=400,
                    bbox_inches='tight', pad_inches=0.04)
    plt.close(fig)
    print(f'  {dst:16s} {len(units)} unit(s)')


# ---------------------------------------------------------------------------
# AlphaFold gallery, re-labelled in English
# ---------------------------------------------------------------------------
def af_gallery(dst='FigS3a_af'):
    af = pd.read_csv(os.path.join(TAB, 'af_interface_vs_cryoem_6partners_RECOMPUTED.csv'))
    af = af.set_index('partner')
    order = ['SMARCB1', 'SMARCA4', 'SMARCD1', 'SMARCC2', 'SMARCE1', 'DPF2']

    im = mpimg.imread(os.path.join(SRC, 'deep_af_structures.png'))
    h, w = im.shape[:2]
    # structure renders only -- the Chinese title, sub-labels and legend are dropped
    rows = [(0.085, 0.455), (0.525, 0.895)]
    cols = [(0.02, 0.33), (0.35, 0.66), (0.68, 0.99)]
    tiles = []
    for r0, r1 in rows:
        for c0, c1 in cols:
            tiles.append(im[int(h * r0):int(h * r1), int(w * c0):int(w * c1)])

    W = TEXT_MM * MM
    fig = plt.figure(figsize=(W, W * 0.78))
    for i, (g, tile) in enumerate(zip(order, tiles)):
        r, c = divmod(i, 3)
        ax = fig.add_axes([0.005 + c * 0.332, 0.620 - r * 0.450, 0.322, 0.350])
        ax.imshow(tile, interpolation='lanczos')
        ax.axis('off')
        row = af.loc[g]
        ax.text(0.5, -0.045, f'ARID1A–{g}', transform=ax.transAxes,
                ha='center', va='top', fontsize=8.5, fontweight='bold')
        ax.text(0.5, -0.125, f'ipTM {row.af_iptm:.2f}  ·  interface Jaccard '
                             f'{row.jaccard:.2f}',
                transform=ax.transAxes, ha='center', va='top', fontsize=7,
                color='#555555')
    fig.text(0.5, 0.035, 'grey = ARID1A   ·   colour = core/cBAF partner subunit '
                         '(named under each panel)',
             ha='center', fontsize=7, color='#555555')
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'{dst}.{ext}'), dpi=400,
                    bbox_inches='tight', pad_inches=0.04)
    plt.close(fig)
    print(f'  {dst:16s} 6 AlphaFold models re-labelled in English')
    return mpimg.imread(os.path.join(OUT, f'{dst}.png'))


# ---------------------------------------------------------------------------
# MD pose overlay, re-labelled from the verified RMSD table
# ---------------------------------------------------------------------------
def md_pose(dst='FigS6b_mdpose'):
    """The source render carries hand-set RMSD annotations (1.9 / 2.2 A) that do
    not correspond to any row of md_rmsd_summary_verified.csv. The structures are
    kept; the annotations are re-cut from the verified table, using the
    pocket-fitted metric over the full trajectory -- the same definition as the
    main-text molecular-dynamics figure.
    """
    md = pd.read_csv(os.path.join(TAB, 'md_rmsd_summary_verified.csv'))
    md = md[(md.metric == 'ligand_RMSD_pocketfit') & (md.window == 'full_trajectory')]
    md = md.set_index('ligand')

    # the source's own annotation band ends at 0.183 of the image height; crop
    # past it so the unverified hand-set values cannot survive into the figure
    im = load('deep_md_pose.png', top=0.195, bottom=0.075)
    h, w = im.shape[:2]
    halves = [im[:, :w // 2], im[:, w // 2:]]

    W = TEXT_MM * MM
    fig = plt.figure(figsize=(W, W * 0.44))
    for i, (lig, ns, col) in enumerate([('PFI-3', 30, '#D9822B'),
                                        ('FHD-286', 9, '#C0392B')]):
        r = md.loc[lig]
        ax = fig.add_axes([0.01 + i * 0.50, 0.02, 0.48, 0.87])
        ax.imshow(halves[i], interpolation='lanczos')
        ax.axis('off')
        ax.text(0.0, 1.02, lig, transform=ax.transAxes, fontsize=9,
                fontweight='bold', color=col, va='bottom')
        ax.text(0.215, 1.025, f'ligand RMSD {r["mean"]:.2f} $\\pm$ {r["std"]:.2f} Å '
                              f'over {ns} ns',
                transform=ax.transAxes, fontsize=7.5, va='bottom', color='#333333')
    fig.text(0.5, -0.005, 'initial docked pose (grey) overlaid on the MD-relaxed pose '
                          '(coloured); binding-site residues in green',
             ha='center', fontsize=7, color='#555555')
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'{dst}.{ext}'), dpi=400,
                    bbox_inches='tight', pad_inches=0.04)
    plt.close(fig)
    print(f'  {dst:16s} RMSD re-labelled from md_rmsd_summary_verified.csv')
    return mpimg.imread(os.path.join(OUT, f'{dst}.png'))


# ---------------------------------------------------------------------------
print('building supplementary figures ...')

# S1 -- what the complex looks like and how the five methods converge
stack([('I',  load('fig_hero_baf_schematic.png', top=0.075, bottom=0.010)),
       ('II', load('fig_hero_convergence.png',   top=0.108, bottom=0.055))],
      'FigS1_architecture')

# S2 -- the genome-wide picture and its druggable subset
stack([('I',  load('deep_landscape_waterfall.png', top=0.075, bottom=0.005)),
       ('II', load('deep_druggable_net.png',       top=0.070, bottom=0.055))],
      'FigS2_landscape')

# S3 -- structure-level detail: predicted models and cancer mutations
af = af_gallery()
stack([('I',  af),
       ('II', load('deep_interface_mut.png', top=0.150, bottom=0.010))],
      'FigS3_structure')

# S4 -- pocket occupancy and pose retention under dynamics
stack([('I',  load('fig_pocket_fhd286_3dmol.png', top=0.105, bottom=0.015, left=0.035)),
       ('II', md_pose())],
      'FigS4_pose')

# S5 -- structure-based repurposing screen.
# Source recovered 2026-08-01 from sxycloud:
#   /data/s01011/arid1a-nsclc/results/tables/deep_vs_results.csv
# 6,392 compounds; top hits and the artefact counts (7 bromodomain, 10 ATPase)
# reproduce the figure exactly. INC-280 = capmatinib and R-428 = bemcentinib,
# which is why the figure's trade names do not appear verbatim in the table.
stack([('I',  load('deep_vs_hits.png',    top=0.060, bottom=0.068)),
       ('II', load('deep_vs_gallery.png', top=0.035, bottom=0.030))],
      'Fig_screen')

# scGPT neighbour control. Source recovered from the local spoai backup:
#   F:\spoai服务器医学数据备份\home_z\arid1a_dock\scgpt_perturb_out\
# median percentiles 78.12 / 96.45 / 81.17 / 86.93 match the figure's 78/96/81/87.
emit_single('deep_scgpt_neighbor.png', 'Fig_scgpt', 168, top=0.055, bottom=0.115)

# ---------------------------------------------------------------------------
# Consolidated figures: three groups of same-topic figures merged so the
# manuscript carries fewer, larger figures rather than many adjacent small ones.
# Units keep roman numerals; each unit's own internal lettering is untouched.
# ---------------------------------------------------------------------------

# experimental structure + predicted models + where mutations land
stack([('I',   load('fig4_cryoem_white.png', top=0.055, bottom=0.055)),
       ('II',  built('FigS3a_af.png')),
       ('III', load('deep_interface_mut.png', top=0.150, bottom=0.010))],
      'FigM_structure')

# both ligands in their pockets, and both poses under dynamics
stack([('I',   load('fig_pocket_pfi3_3dmol.png', top=0.105, bottom=0.015, left=0.035)),
       ('II',  load('fig_pocket_fhd286_3dmol.png', top=0.105, bottom=0.015, left=0.035)),
       ('III', built('FigS6b_mdpose.png'))],
      'FigM_poses')

# single-cell context and the scGPT control, which answer the same question
stack([('I',  built('Fig_scrna.png')),
       ('II', load('deep_scgpt_neighbor.png', top=0.055, bottom=0.115))],
      'FigM_singlecell')

print('\nwritten to', OUT)
