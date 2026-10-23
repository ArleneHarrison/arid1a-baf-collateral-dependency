"""Descriptive re-analysis of published iPSC KO protein/RNA data, not a new experiment."""
from pathlib import Path
import argparse,json,hashlib,platform
import numpy as np
import pandas as pd
CONTACT=['SMARCC2','SMARCA4','SMARCD1','SMARCB1','SMARCE1','DPF2']
BACKGROUND=['ACTL6A','PBRM1','ARID2','BRD7','PHF10','BRD9']
ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(exist_ok=True)
d=pd.read_excel(a.input,header=2);assert d.columns[21]=='WT cloneA.1' and d.columns[23]=='ARID1A KO cloneA.1'
rows=[]
for gene in ['ARID1A']+CONTACT+BACKGROUND:
    hits=d[d['Gene Name']==gene];assert len(hits)==1,(gene,len(hits));r=hits.iloc[0]
    wt=pd.to_numeric(r.iloc[21:23]).to_numpy(dtype=float)
    for query,k,prefix in [('ARID1A',23,'  ARID1A'),('ARID2',25,' ARID2'),('PBRM1',27,' PBRM1')]:
        ko=pd.to_numeric(r.iloc[k:k+2]).to_numpy(dtype=float);ok=np.isfinite(wt).all() and np.isfinite(ko).all() and np.min(np.r_[wt,ko])>0
        fc=float(np.log2(ko.mean()/wt.mean())) if ok else np.nan
        ratios=np.log2(ko[:,None]/wt[None,:]) if ok else np.full((2,2),np.nan)
        rows.append(dict(gene=gene,panel='anchor' if gene=='ARID1A' else 'contact' if gene in CONTACT else 'background',query=query,protein_log2FC=fc,WT_A=wt[0],WT_B=wt[1],KO_A=ko[0],KO_B=ko[1],all_four_ratio_min=np.nanmin(ratios),all_four_ratio_max=np.nanmax(ratios),RNA_log2FC=r['log2FoldChange'+prefix],RNA_source_padj=r['padj'+prefix],source_omnibus_ANOVA_p=r['ANOVA p value'],source_omnibus_FDR_hit=r['permutaiton based FDR<0.1'],n_WT_protein_columns=2,n_KO_protein_columns=2,protein_evaluable=ok))
result=pd.DataFrame(rows);result.to_csv(a.out/'fixed_panel_KO_results.csv',index=False)
summary={'host':platform.node(),'source_sha256':hashlib.sha256(a.input.read_bytes()).hexdigest(),'article':'10.1016/j.celrep.2017.08.010','endpoint':'Published author-normalized protein log2 mean(KO)/mean(WT); source DESeq2 RNA log2FC and adjusted P retained','source_sample_units':'Two labelled clones per protein condition; RNA source methods report six biological replicates per KO and four WT. These are not twelve independent proteins or new lung cancer models.','no_new_gene_or_module_P_values':True,'source_protein_ANOVA_not_ARID1A_pairwise':True,'no_added_imputation':True,'contrasts':[]}
for query,block in result.groupby('query',sort=False):
    c=block[block.panel=='contact'];b=block[block.panel=='background']
    summary['contrasts'].append({'query':query,'contact_protein_median':float(c.protein_log2FC.median()),'background_protein_median':float(b.protein_log2FC.median()),'contact_minus_background':float(c.protein_log2FC.median()-b.protein_log2FC.median()),'contact_RNA_median':float(c.RNA_log2FC.median()),'contact_protein_down_count':int((c.protein_log2FC<0).sum()),'contact_source_RNA_padj_below_005_count':int((c.RNA_source_padj<.05).sum()),'contact_coverage':int(c.protein_evaluable.sum()),'background_coverage':int(b.protein_evaluable.sum())})
    b_ex=b[b.gene!=query]
    summary['contrasts'][-1].update({'background_without_query_n':len(b_ex),'background_without_query_protein_median':float(b_ex.protein_log2FC.median()),'contrast_without_query':float(c.protein_log2FC.median()-b_ex.protein_log2FC.median())})
summary['source_preprocessing_limitation']='Published methods include peptide-channel minimum-intensity handling and loading normalization. No new imputation was performed; upstream imputation cannot be undone from this table.'
(a.out/'ipsc_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
