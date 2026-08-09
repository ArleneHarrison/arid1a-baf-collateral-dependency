"""经典GraphPad Prism风格柱状图——全封闭坐标轴+黑白对比柱+括线星号+全粗体字
模仿参考图：白色空心柱(对照) vs 黑色实心柱(处理)，误差棒上下端有小帽，
组间显著性用"匚"形括线+居中星号，坐标轴四边全部保留，字体全部加粗。
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm, rcParams
from matplotlib.patches import Rectangle

def apply():
    for fp in [r"C:\Windows\Fonts\msyhbd.ttc", r"C:\Windows\Fonts\msyh.ttc", r"C:\Windows\Fonts\simhei.ttf"]:
        if os.path.exists(fp):
            fm.fontManager.addfont(fp)
            rcParams["font.family"] = fm.FontProperties(fname=fp).get_name()
            break
    rcParams.update({
        "axes.unicode_minus": False,
        "font.size": 12,
        "font.weight": "bold",
        "axes.labelweight": "bold",
        "axes.titleweight": "bold",
        "axes.linewidth": 1.6,
        "axes.edgecolor": "black",
        "xtick.major.width": 1.6, "ytick.major.width": 1.6,
        "xtick.major.size": 5, "ytick.major.size": 5,
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "pdf.fonttype": 42,
    })

def boxed(ax):
    for s in ["top", "right", "left", "bottom"]:
        ax.spines[s].set_visible(True); ax.spines[s].set_linewidth(1.6); ax.spines[s].set_color("black")
    ax.tick_params(direction="out", width=1.6, length=5)

def bar2(ax, x, mean, sem, filled, width=0.62, cap=6):
    fc = "black" if filled else "white"
    ax.bar(x, mean, width=width, facecolor=fc, edgecolor="black", linewidth=1.6, zorder=3)
    ax.errorbar(x, mean, yerr=sem, ecolor="black", elinewidth=1.6, capsize=cap/2,
                capthick=1.6, zorder=4, fmt="none")

def bracket(ax, x1, x2, y, h, text="*", fontsize=15, lw=1.6):
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=lw, color="black", solid_capstyle="butt", zorder=5)
    ax.text((x1+x2)/2, y+h*1.05, text, ha="center", va="bottom", fontsize=fontsize,
            fontweight="bold", color="black", zorder=5)

def panel(ax, letter, x=-0.30, y=1.12, size=22):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=size, fontweight="bold",
             ha="left", va="bottom", color="black")
