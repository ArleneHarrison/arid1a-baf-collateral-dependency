"""图1C-2 重画：GSE131907恶性上皮细胞中当前湿实验靶点(ARID1A/ARID1B/DPF2/PBRM1)的检出率
真实数据(本会话服务器重新提取)。已自查：原图的"ARID1A-low组内伙伴基因检出率"对比经验证被测序深度
严重混杂(控制nGene后rho从0.15-0.18降到0.02-0.03)，故本图不再使用该对比，只呈现干净的整体检出率。
"""
import csv, numpy as np, matplotlib.pyplot as plt
import navistyle as ns
ns.apply()

rows = list(csv.DictReader(open("../results_tables/fig1c2_recompute_detection.csv", encoding="utf-8-sig")))
genes = [r["gene"] for r in rows]
pct = [float(r["pct_detected"]) for r in rows]
n_det = [int(r["n_detected"]) for r in rows]
n_total = int(rows[0]["n_total"])

fig = plt.figure(figsize=(7.6, 5.2))
ns.title(fig, "湿实验靶点在患者肿瘤恶性上皮细胞中的表达",
    f"GSE131907 肺腺癌单细胞（恶性上皮细胞 n={n_total:,}）｜检出率 = 该基因 UMI>0 的细胞占比（真实数据）")
ax = fig.add_axes([0.15, 0.20, 0.80, 0.56])
cols = [ns.CORE if g == "ARID1A" else ns.PBAF for g in genes]
bars = ax.bar(range(len(genes)), pct, color=cols, width=0.58, zorder=3)
ns.bar_gap(bars)
for i, (p, n) in enumerate(zip(pct, n_det)):
    ax.text(i, p + 1.2, f"{p:.1f}%\n(n={n:,})", ha="center", va="bottom", fontsize=8.6, color=ns.INK2)
ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, fontsize=11, fontweight="bold")
ax.set_ylabel("检出率（该基因 UMI>0 的细胞占比）")
ax.set_ylim(0, max(pct)*1.35)
ns.hairline_grid(ax)
ns.footnote(fig,
    f"ARID1B(23.2%)/DPF2(17.6%)/PBRM1(23.9%)——本计划书最终选定的湿实验靶点——在患者恶性上皮细胞中的检出率"
    f"与 ARID1A 本身(42.4%)属同一量级，说明这些靶点确实在肿瘤细胞中表达，为湿实验的临床相关性提供组织层面的合理性。\n"
    "诚实标注：这是描述性的\"靶点存在性/表达语境\"证据，不是依赖性论断。原计划按 ARID1A 检出与否分层比较伙伴基因检出率，"
    "但自查发现该分层与测序深度严重混杂(控制深度后关联基本消失)，故本图不采用该对比，仅呈现整体检出率。",
    y=0.045, size=7.6)
fig.savefig("fig1c2_v2_wetlab_targets.png")
print("saved fig1c2_v2_wetlab_targets.png")
