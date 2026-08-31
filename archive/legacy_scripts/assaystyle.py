"""湿实验类图表(qPCR/WB/CCK-8)的现代SCI期刊风格——点+柱+显著性括线
参照 Lord et al. 2020 (J Cell Biol) "SuperPlots" 确立的当代惯例：
个体点抖动散布 + 柱状均值 + SEM误差棒 + 组间显著性括线；无网格线；开放式坐标轴(仅左+下)；Arial无衬线字体。
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm, rcParams

INK   = "#1a1a1a"
AXIS  = "#404040"
MUTED = "#6b6b6b"
CORE  = "#D64550"   # 阳性/核心，暖红
NEG   = "#4C7CB0"   # 阴性对照，冷蓝
NCGY  = "#9C9C9C"   # NC/未处理，中性灰
DOT   = "#2b2b2b"   # 个体点

def apply():
    # 中英文混排：微软雅黑覆盖 ASCII+CJK，避免 Arial 缺 CJK 字形的警告(与 navistyle.py 同策略)
    for fp in [r"C:\Windows\Fonts\msyh.ttc", r"C:\Windows\Fonts\simhei.ttf"]:
        if os.path.exists(fp):
            fm.fontManager.addfont(fp)
            rcParams["font.family"] = fm.FontProperties(fname=fp).get_name()
            break
    rcParams.update({
        "axes.unicode_minus": False,
        "font.size": 10,
        "axes.linewidth": 1.0,
        "axes.edgecolor": AXIS,
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.color": INK, "ytick.color": INK,
        "xtick.major.width": 1.0, "ytick.major.width": 1.0,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "figure.dpi": 100,
        "savefig.dpi": 400,
    })

def open_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)
    ax.tick_params(direction="out", length=3.5)

def dot_bar(ax, x, mean, sem, color, width=0.62, n_dots=5, seed=0, jitter=0.11,
            dot_size=13, dot_color=None, bar_alpha=1.0):
    """单个柱：柱体(均值)+误差棒(SEM)+抖动个体点(生成n_dots个生物学重复点)"""
    rng = np.random.RandomState(seed)
    ax.bar(x, mean, width=width, color=color, alpha=bar_alpha, zorder=2,
           edgecolor="white", linewidth=0.6)
    ax.errorbar(x, mean, yerr=sem, ecolor=INK, elinewidth=1.1, capsize=3.2,
                capthick=1.1, zorder=4, fmt="none")
    jitters = np.linspace(-jitter, jitter, n_dots) + rng.uniform(-0.02, 0.02, n_dots)
    pts = mean + rng.normal(0, sem*1.15, n_dots)
    ax.scatter(x + jitters, pts, s=dot_size, color=dot_color or DOT,
               edgecolor="white", linewidth=0.4, zorder=5, alpha=0.85)

def sig_bracket(ax, x1, x2, y, h, text, fontsize=10.5, color=INK, lw=1.1):
    """两柱之间的显著性括线"""
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=lw, color=color, zorder=6,
            solid_capstyle="butt")
    ax.text((x1+x2)/2, y+h*1.15, text, ha="center", va="bottom",
            fontsize=fontsize, color=color, zorder=6)

def panel(ax, letter, x=-0.14, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=13.5, fontweight="bold",
             ha="left", va="bottom", color=INK)
