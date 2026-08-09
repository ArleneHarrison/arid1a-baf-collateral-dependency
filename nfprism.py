# -*- coding: utf-8 -*-
"""Nature-style typography rendered in a GraphPad-Prism black-and-white bar aesthetic.

Nature conventions taken: Arial/Helvetica sans-serif at regular weight (not everything
bold), 7-8 pt labels, hairline 0.75 pt axes, left+bottom spines only, ticks pointing out,
bold panel letters set flush left above the axes, no gridlines, no chartjunk.

Prism conventions taken: white-filled bars with black outlines for controls, hatched and
solid-black bars for treatments, capped error bars, individual replicate points drawn on
top of every bar, and square significance brackets with asterisks.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm, rcParams

# Prism-style bar fills, black and white only
OPEN = dict(facecolor='white', edgecolor='black', hatch=None)
HATCH = dict(facecolor='white', edgecolor='black', hatch='////')
HATCH2 = dict(facecolor='white', edgecolor='black', hatch='....')
GREY = dict(facecolor='#BFBFBF', edgecolor='black', hatch=None)
SOLID = dict(facecolor='black', edgecolor='black', hatch=None)

LW = 0.75          # hairline, Nature spec
FS = 7             # base font size, Nature spec
CAP = 2.2


def apply():
    for fp in [r'C:\Windows\Fonts\arial.ttf', r'C:\Windows\Fonts\Helvetica.ttf',
               r'C:\Windows\Fonts\segoeui.ttf']:
        if os.path.exists(fp):
            fm.fontManager.addfont(fp)
            rcParams['font.family'] = fm.FontProperties(fname=fp).get_name()
            break
    rcParams.update({
        'font.size': FS,
        'font.weight': 'normal',
        'axes.labelsize': FS, 'axes.titlesize': FS,
        'xtick.labelsize': FS - 0.5, 'ytick.labelsize': FS - 0.5,
        'legend.fontsize': FS - 0.5,
        'axes.labelweight': 'normal', 'axes.titleweight': 'normal',
        'axes.linewidth': LW, 'axes.edgecolor': 'black',
        'xtick.major.width': LW, 'ytick.major.width': LW,
        'xtick.major.size': 2.4, 'ytick.major.size': 2.4,
        'xtick.direction': 'out', 'ytick.direction': 'out',
        'axes.unicode_minus': False,
        'axes.grid': False,
        'figure.facecolor': 'white', 'axes.facecolor': 'white',
        'savefig.facecolor': 'white', 'savefig.bbox': 'tight',
        'pdf.fonttype': 42, 'ps.fonttype': 42,
        'hatch.linewidth': 0.55,
    })


def spines(ax):
    """Nature default: left and bottom only."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_linewidth(LW)
        ax.spines[s].set_color('black')
    ax.tick_params(direction='out', width=LW, length=2.4, colors='black')


def bar(ax, x, vals, style, width=0.62, jitter=0.11, seed=0, points=True,
        ms=6.5, err='sem'):
    """Prism-style bar: mean +- SEM with caps, individual replicates overlaid."""
    vals = np.asarray(vals, float)
    m = vals.mean()
    e = vals.std(ddof=1) / np.sqrt(len(vals)) if err == 'sem' else vals.std(ddof=1)
    ax.bar(x, m, width=width, linewidth=LW, zorder=2, **style)
    ax.errorbar(x, m, yerr=e, ecolor='black', elinewidth=LW, capsize=CAP,
                capthick=LW, zorder=3, fmt='none')
    if points and len(vals) > 1:
        rng = np.random.RandomState(int(round(seed)) % (2 ** 31))
        off = np.linspace(-jitter, jitter, len(vals)) + rng.uniform(-.012, .012, len(vals))
        ax.scatter(x + off, vals, s=ms, facecolor='white', edgecolor='black',
                   linewidth=LW * 0.85, zorder=4, clip_on=False)
    return m, e


def bracket(ax, x1, x2, y, h, text='*', fs=None, dy=0.0):
    """Square significance bracket with centred asterisk(s)."""
    fs = fs or FS + 1
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=LW, color='black',
            solid_capstyle='butt', clip_on=False, zorder=5)
    va, pad = ('bottom', h * 0.15) if text != 'n.s.' else ('bottom', h * 0.35)
    ax.text((x1 + x2) / 2, y + h + pad + dy, text, ha='center', va=va,
            fontsize=fs if text != 'n.s.' else FS - 0.5, color='black', zorder=5)


def panel(ax, letter, x=-0.16, y=1.06, size=9):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=size, fontweight='bold',
            ha='left', va='bottom', color='black')


def stars(p):
    return ('****' if p < 1e-4 else '***' if p < 1e-3 else
            '**' if p < 1e-2 else '*' if p < .05 else 'n.s.')
