"""
ARID1A-BAF 论文图件统一样式（Nature/SCI 顶刊规范）

配色经 dataviz 验证器校验通过（light, surface #ffffff, categorical 4 slots）：
  最差相邻 正常视觉 ΔE 20.8（底线 15）｜色盲 deutan 15.3 / tritan 9.6（目标 ≥8）
语义沿用全文既有约定（红=核心），但把原先的橙色替换为黄色——
原配色 红#C0392B 与 橙#E8A33D 的正常视觉 ΔE 仅 7.1，远低于 15 底线，属真实可读性缺陷。
"""
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm, rcParams

# ---------- 语义配色（经校验） ----------
CORE   = "#e34948"   # 核心共享亚基 / 主结论
CBAF   = "#eda100"   # cBAF 特异亚基
PBAF   = "#2a78d6"   # PBAF 特异亚基
NCBAF  = "#1baf7a"   # ncBAF 特异亚基
# 中性（非分类语义，用于 null / 对照 / 被解释掉的项）
NULL   = "#898781"   # muted
NULL_L = "#d6d5d0"   # 更浅的中性，用于"未校正/背景"对比

# ---------- 墨色与图面 ----------
INK      = "#0b0b0b"   # 主要文字
INK2     = "#52514e"   # 次要文字
MUTED    = "#898781"   # 轴标签 / 脚注
GRID     = "#e1e0d9"   # 发丝网格
AXIS     = "#c3c2b7"   # 轴线 / 基线
SURFACE  = "#ffffff"

# 顺序色（单色相 blue，用于连续量）
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]


def apply(chinese=True):
    """套用全局 rcParams。中文用微软雅黑；正式英文稿应改为 Arial/Helvetica。"""
    if chinese:
        for fp in [r"C:\Windows\Fonts\msyh.ttc", r"C:\Windows\Fonts\simhei.ttf"]:
            if os.path.exists(fp):
                fm.fontManager.addfont(fp)
                rcParams["font.family"] = fm.FontProperties(fname=fp).get_name()
                break
    else:
        rcParams["font.family"] = "Arial"
    rcParams.update({
        "axes.unicode_minus": False,
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        # 极简边框：只留左下
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": AXIS,
        "axes.linewidth": 0.8,
        # 轴与刻度：recessive
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelcolor": INK2,
        "ytick.labelcolor": INK2,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.direction": "out",
        "ytick.direction": "out",
        # 发丝网格，默认关闭，需要时按轴打开
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "axes.grid": False,
        # 文字
        "axes.labelcolor": INK2,
        "axes.titlecolor": INK,
        "text.color": INK,
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 10.5,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8,
        "legend.frameon": False,
        # 线与marker
        "lines.linewidth": 1.8,
        "lines.markersize": 5,
        "figure.dpi": 200,
        "savefig.dpi": 400,          # 顶刊要求 ≥300 dpi
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,          # 字体可编辑（投稿要求）
        "ps.fonttype": 42,
    })


def title(fig, main, sub=None, y=0.975, ysub=0.932, size=13.5, subsize=8.6):
    """统一标题层级：主标题加粗主墨色，副标题次要墨色。"""
    fig.text(0.5, y, main, ha="center", va="top", fontsize=size,
             fontweight="bold", color=INK)
    if sub:
        fig.text(0.5, ysub, sub, ha="center", va="top", fontsize=subsize,
                 color=INK2, linespacing=1.45)


def panel(ax, letter, text, pad=10):
    """Panel 标签：字母加粗 + 说明，左对齐（Nature 惯例）。"""
    ax.set_title(f"{letter}  {text}", fontsize=10, fontweight="bold",
                 color=INK, loc="left", pad=pad)


def footnote(fig, text, y=0.045, size=7.8):
    """脚注：不用灰底色块，改用细上分隔线 + 次要墨色（更接近期刊排版）。"""
    fig.text(0.5, y, text, ha="center", va="center", fontsize=size,
             color=INK2, linespacing=1.5)


def bar_gap(bars, lw=1.4):
    """相邻/堆叠柱之间留出图面色描边，形成 2px 视觉间隙。"""
    for b in bars:
        b.set_edgecolor(SURFACE)
        b.set_linewidth(lw)


def despine_zero(ax, axis="y"):
    """零基线用轴色发丝线，替代粗黑线。"""
    if axis == "y":
        ax.axhline(0, color=AXIS, lw=0.9, zorder=1)
    else:
        ax.axvline(0, color=AXIS, lw=0.9, zorder=1)


def hairline_grid(ax, axis="y"):
    ax.grid(axis=axis, color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
