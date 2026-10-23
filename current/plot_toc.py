"""Regenerate the TOC graphic directly from archived numeric estimates.

Usage: python plot_toc.py --source TOC_source_data.csv --outdir .
No model-generated pixels, synthetic observations or statistical refitting.
"""
import argparse
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,default=Path('TOC_source_data.csv'));p.add_argument('--outdir',type=Path,default=Path('.'));a=p.parse_args()
    d=pd.read_csv(a.source).set_index('cohort');a.outdir.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig=plt.figure(figsize=(3.25,1.75),facecolor='white');ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
    ax.text(.5,.95,'ARID1A-associated cBAF protein co-variation',ha='center',va='top',fontsize=7.3,weight='bold')
    ax.add_patch(FancyBboxPatch((.025,.23),.31,.56,boxstyle='round,pad=.01,rounding_size=.025',fc='#EEF5F9',ec='#0072B2',lw=.7))
    ax.text(.18,.735,'Six structure-selected\ncBAF proteins',ha='center',va='center',fontsize=5.9,weight='bold')
    for i,g in enumerate(['SMARCC2','SMARCA4','SMARCD1','SMARCB1','SMARCE1','DPF2']):ax.text(.18,.63-i*.061,g,ha='center',fontsize=6.2,color='#0072B2')
    ax.add_patch(FancyArrowPatch((.35,.51),(.435,.51),arrowstyle='-|>',mutation_scale=12,lw=.9,color='.3'))
    ax.text(.71,.79,'Independent tumour\nproteogenomes',ha='center',va='top',fontsize=6.2,weight='bold')
    ax.text(.71,.655,'RNA, copy number and purity adjusted',ha='center',fontsize=5.5)
    plot=fig.add_axes([.62,.30,.32,.29])
    for name,y,col in [('LUAD',1,'#0072B2'),('LUSC',0,'#B34D00')]:
        r=d.loc[name];plot.errorbar(r.estimate,y,xerr=[[r.estimate-r.ci_low],[r.ci_high-r.estimate]],fmt='o',ms=3,color=col,capsize=2,lw=.9)
    plot.set(xlim=(-.02,.43),ylim=(-.6,1.6),yticks=[1,0],yticklabels=[f'LUAD, n={int(d.loc["LUAD","n"])}',f'LUSC, n={int(d.loc["LUSC","n"])}'],xticks=[0,.2,.4])
    plot.tick_params(labelsize=5.5,length=2,pad=2);plot.axvline(0,lw=.5,ls='--',color='.6');plot.spines[['top','right','left']].set_visible(False)
    plot.set_xlabel('Selected minus BAF-background\ncorrelation (95% interval)',fontsize=5.5,labelpad=2)
    ax.text(.5,.045,'Protein co-variation does not establish synthetic lethality',ha='center',fontsize=5.8,color='.25')
    for suffix in ['pdf','svg','tif','png']:
        kw={'dpi':600,'facecolor':'white'}
        if suffix=='tif':kw['pil_kwargs']={'compression':'tiff_lzw'}
        fig.savefig(a.outdir/f'TOC_graphic.{suffix}',**kw)
    plt.close(fig)

if __name__=='__main__':main()
