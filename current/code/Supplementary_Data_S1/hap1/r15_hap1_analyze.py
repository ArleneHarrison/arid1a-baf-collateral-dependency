"""Frozen-panel HAP1 ARID1A GI pilot; execute on seeta, not local workstation."""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import re
import time
import zipfile
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import openpyxl
import scipy
from scipy.stats import fisher_exact, spearmanr
from statsmodels.stats.multitest import multipletests

NS='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''): h.update(b)
    return h.hexdigest()

def save(path,data):
    Path(path).write_text(json.dumps(data,indent=2,ensure_ascii=False,default=str),encoding='utf-8')

def colletters(n):
    s=''
    while n:
        n,r=divmod(n-1,26); s=chr(65+r)+s
    return s

def selected_xml_columns(file,sheet,wanted):
    """Select named columns while streaming the original XLSX worksheet XML.

    Openpyxl is used independently later to audit extracted fixed-gene values.
    Unselected wide matrix values are never accumulated in memory.
    """
    book=openpyxl.load_workbook(file,read_only=True,data_only=True)
    ws=book[sheet]; header=next(ws.iter_rows(min_row=1,max_row=1,values_only=True))
    assert all(x in header for x in wanted)
    positions={colletters(header.index(x)+1):x for x in wanted}
    xmlpath=ws._worksheet_path.lstrip('/'); book.close()
    rows=[]
    with zipfile.ZipFile(file) as z:
        strings=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            with z.open('xl/sharedStrings.xml') as f:
                for _,el in ET.iterparse(f,events=('end',)):
                    if el.tag==NS+'si':
                        strings.append(''.join(x.text or '' for x in el.iter(NS+'t'))); el.clear()
        with z.open(xmlpath) as f:
            for event,el in ET.iterparse(f,events=('end',)):
                if el.tag!=NS+'row': continue
                if int(el.attrib.get('r','0'))<=1:
                    el.clear(); continue
                row={k:None for k in wanted}
                for cell in el:
                    addr=cell.attrib.get('r',''); col=re.sub('[0-9]','',addr)
                    if col not in positions: continue
                    key=positions[col]; v=cell.find(NS+'v'); kind=cell.attrib.get('t')
                    if kind=='s' and v is not None: row[key]=strings[int(v.text)]
                    elif kind=='inlineStr': row[key]=''.join(t.text or '' for t in cell.iter(NS+'t'))
                    elif v is not None:
                        try: row[key]=float(v.text)
                        except (ValueError,TypeError): row[key]=v.text
                rows.append(row); el.clear()
    frame=pd.DataFrame(rows,columns=wanted)
    return frame

def exact_contrast(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float); joined=np.r_[a,b]
    assert np.isfinite(joined).all()
    observed=a.mean()-b.mean(); total=float(joined.sum()); n=len(a)
    values=[]
    for idx in itertools.combinations(range(len(joined)),n):
        x=float(joined[list(idx)].sum())
        values.append(x/n-(total-x)/len(b))
    values=np.array(values)
    exceed=int(np.count_nonzero(abs(values)>=abs(observed)-1e-12))
    return {'mean_difference_negative_qGI':observed,'mean_group_A_negative_qGI':a.mean(),
            'mean_group_B_negative_qGI':b.mean(),'n_permutations':len(values),
            'n_as_or_more_extreme':exceed,'p_two_sided_exact':exceed/len(values),
            'median_difference_negative_qGI_descriptive':float(np.median(a)-np.median(b))}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--out',required=True); args=ap.parse_args()
    out=Path(args.out); raw=out/'raw'; spec=json.loads((out/'R15_HAP1_FROZEN_SPEC.json').read_text())
    sets=spec['fixed_sets']; aliases=spec['aliases']; fixed=sorted(set(itertools.chain.from_iterable(sets.values())))
    assert sha(raw/'File_S4.xlsx')=='9df16c51264e9bfc8bba9620ccb0eb5cdd308a5d053c1fcebe0152ca07ea57b6'
    anchor_column='ARID1A_019_min'
    extracted=[]
    for sheet,newname in [('qGI_scores','author_qGI'),('qGI_FDR','author_FDR'),('queries_lfc','primary_query_LFC')]:
        print('EXTRACTING',sheet,flush=True)
        f=selected_xml_columns(raw/'File_S4.xlsx',sheet,['Gene.symbol',anchor_column])
        f=f.rename(columns={'Gene.symbol':'source_gene',anchor_column:newname})
        f['gene']=f['source_gene'].replace(aliases)
        assert not f['gene'].duplicated().any()
        extracted.append(f)
    data=extracted[0]
    for f in extracted[1:]: data=data.merge(f.drop(columns='source_gene'),on='gene',how='outer',validate='one_to_one')
    data['negative_qGI_dependency_score']=-pd.to_numeric(data['author_qGI'],errors='coerce')
    data['author_FDR']=pd.to_numeric(data['author_FDR'],errors='coerce')
    data['finite_scored']=np.isfinite(data['negative_qGI_dependency_score']) & np.isfinite(data['author_FDR'])
    data['standard_negative_hit']=(data['author_qGI'] < -.3)&(data['author_FDR']<.1)&data['finite_scored']
    data['standard_positive_hit']=(data['author_qGI'] > .3)&(data['author_FDR']<.1)&data['finite_scored']
    data['strict_negative_hit']=(data['author_qGI'] < -.6)&(data['author_FDR']<.01)&data['finite_scored']
    data['fixed_set_membership']=data['gene'].map(lambda g:';'.join(k for k,v in sets.items() if g in v))
    data.to_csv(out/'S4_ARID1A_all_gene_scores.csv',index=False)
    by_gene=data.set_index('gene')
    fixedframe=by_gene.reindex(fixed).reset_index()
    fixedframe.to_csv(out/'S4_ARID1A_fixed_panel_scores.csv',index=False)
    print('PRIMARY_GENE_EXTRACTION_COMPLETE',len(data),flush=True)

    validation=pd.read_excel(raw/'File_S6.xlsx',sheet_name='Validation_queries _LFC_dLFC')
    validation=validation.rename(columns={'Gene.symbol':'source_gene'})
    validation['gene']=validation['source_gene'].replace(aliases)
    assert not validation['gene'].duplicated().any()
    validation['negative_dLFC_dependency_score']=-validation['ARID1A.T19_dLFC']
    validation['dLFC_recalculated_from_author_LFC']=validation['ARID1A.T19_LFC']-validation['WT.T19_LFC']
    validation['dLFC_arithmetic_abs_error']=abs(validation['ARID1A.T19_dLFC']-validation['dLFC_recalculated_from_author_LFC'])
    validation=validation.merge(data[['gene','author_qGI','author_FDR','standard_negative_hit','standard_positive_hit']],on='gene',how='left',validate='one_to_one')
    validation.to_csv(out/'S6_ARID1A_all_available_gene_scores.csv',index=False)
    vfixed=validation.set_index('gene').reindex(fixed).reset_index()
    vfixed.to_csv(out/'S6_ARID1A_fixed_panel_available_cases.csv',index=False)
    vmap=validation.set_index('gene')['negative_dLFC_dependency_score']
    coverage=[]
    for name,members in sets.items():
        primary_present=[g for g in members if g in by_gene.index and by_gene.loc[g,'finite_scored']]
        valid_present=[g for g in members if g in vmap.index and np.isfinite(vmap.loc[g])]
        coverage.append({'set':name,'n_fixed':len(members),'fixed_members':';'.join(members),
            'n_primary_finite':len(primary_present),'primary_missing':';'.join(g for g in members if g not in primary_present),
            'n_S6_available':len(valid_present),'S6_missing':';'.join(g for g in members if g not in valid_present),
            'S6_available_genes':';'.join(valid_present),
            'S6_available_mean_negative_dLFC_descriptive':float(vmap.loc[valid_present].mean()) if valid_present else np.nan,
            'S6_full_fixed_set_evaluable':len(valid_present)==len(members)})
    pd.DataFrame(coverage).to_csv(out/'fixed_panel_coverage.csv',index=False)

    contrasts=[]
    for a,b in spec['primary_contrasts']:
        row={'group_A':a,'group_B':b,'n_A_fixed':len(sets[a]),'n_B_fixed':len(sets[b]),
             'exchangeability_unit':'gene, not independent HAP1 clones/patients'}
        values_a=by_gene.reindex(sets[a]); values_b=by_gene.reindex(sets[b])
        evaluable=values_a['finite_scored'].fillna(False).all() and values_b['finite_scored'].fillna(False).all()
        row['evaluable_full_fixed_panel']=bool(evaluable)
        if evaluable: row.update(exact_contrast(values_a['negative_qGI_dependency_score'],values_b['negative_qGI_dependency_score']))
        else: row['p_two_sided_exact']=np.nan
        contrasts.append(row)
    c=pd.DataFrame(contrasts); p=c['p_two_sided_exact'].fillna(1).to_numpy()
    c['Holm_p_three_prespecified']=multipletests(p,method='holm')[1]
    c.loc[~c['evaluable_full_fixed_panel'],'Holm_p_three_prespecified']=np.nan
    c['supports_greater_group_A_dependency']=(c['mean_difference_negative_qGI']>0)&(c['Holm_p_three_prespecified']<.05)
    c.to_csv(out/'primary_direct_contrasts.csv',index=False)

    measured=data[data.finite_scored & (data.gene!='ARID1A')].set_index('gene')
    enrichment=[]
    for name,genes in sets.items():
        row={'set':name,'n_fixed':len(genes),'background_actual_scored_genes':len(measured),
             'background_negative_hits':int(measured.standard_negative_hit.sum()),'test':'one-sided competitive Fisher exact'}
        if all(g in measured.index for g in genes):
            n_hit=int(measured.loc[genes,'standard_negative_hit'].sum()); rem=measured.drop(index=genes)
            n_other_hit=int(rem.standard_negative_hit.sum())
            odds,p=fisher_exact([[n_hit,len(genes)-n_hit],[n_other_hit,len(rem)-n_other_hit]],alternative='greater')
            row.update({'n_set_negative_hits':n_hit,'odds_ratio':odds,'p':p})
        else: row.update({'n_set_negative_hits':np.nan,'odds_ratio':np.nan,'p':np.nan})
        enrichment.append(row)
    e=pd.DataFrame(enrichment); e['BH_q_five_prespecified']=multipletests(e.p.fillna(1),method='fdr_bh')[1]
    e.loc[e.p.isna(),'BH_q_five_prespecified']=np.nan
    e.to_csv(out/'secondary_negative_hit_enrichment.csv',index=False)

    s1=pd.read_excel(raw/'File_S1.xlsx',sheet_name='Query and screen summary')
    s1_fixed=s1[s1['Query Gene'].isin(fixed+['ARID1A','SMARCC1','ACTB'])]
    s1_fixed.to_csv(out/'S1_BAF_query_screen_metadata.csv',index=False)
    arid_meta=s1[s1['Query Gene']=='ARID1A']
    assert len(arid_meta)==1 and arid_meta.iloc[0]['GIN_ID']=='GIN019'

    print('AUDITING_SELECTED_GUIDE_COUNTS',flush=True)
    guides=pd.read_excel(raw/'File_S6.xlsx',sheet_name='Tiling_library_reads')
    guides['gene']=guides['Gene.symbol'].replace(aliases)
    guide_fixed=guides[guides.gene.isin(fixed)].copy()
    t0cols=[x for x in guides.columns if isinstance(x,str) and x.endswith('.T0')]
    guide_fixed['passes_all_query_T0_min30']=guide_fixed[t0cols].notna().all(axis=1)&(guide_fixed[t0cols]>=30).all(axis=1)
    guide_fixed['passes_ARID1A_WT_T0_min30']=guide_fixed[['WT.T0','ARID1A.T0']].notna().all(axis=1)&(guide_fixed[['WT.T0','ARID1A.T0']]>=30).all(axis=1)
    guide_fixed.to_csv(out/'S6_fixed_panel_guide_counts_raw.csv',index=False)
    guide_summary=guide_fixed.groupby('gene').agg(n_guides=('Guide.seq','size'),n_unique_guides=('Guide.seq','nunique'),n_all_query_T0_pass=('passes_all_query_T0_min30','sum'),n_pair_T0_pass=('passes_ARID1A_WT_T0_min30','sum')).reindex(fixed)
    guide_summary.to_csv(out/'S6_fixed_panel_guide_coverage.csv')

    overlap=validation[validation.gene.isin(fixed)].dropna(subset=['author_qGI','ARID1A.T19_dLFC'])
    match_nonzero=overlap[(overlap.author_qGI!=0)&(overlap['ARID1A.T19_dLFC']!=0)]
    sanity={
      'primary_gene_rows':len(data),'primary_finite_gene_rows':int(data.finite_scored.sum()),
      'primary_unique_gene_symbols':int(data.gene.nunique()),'primary_background_excluding_anchor':len(measured),
      'primary_author_negative_hits':int(measured.standard_negative_hit.sum()),
      'primary_author_positive_hits':int(measured.standard_positive_hit.sum()),
      'validation_gene_rows':len(validation),'validation_unique_genes':int(validation.gene.nunique()),
      'validation_raw_guides':len(guides),'validation_raw_unique_genes':int(guides.gene.nunique()),
      'raw_only_validation_genes':sorted(set(guides.gene)-set(validation.gene)),
      'validation_max_dLFC_arithmetic_abs_error':float(validation.dLFC_arithmetic_abs_error.max()),
      'validation_fixed_overlap_n':len(overlap),'validation_fixed_same_sign_n':int(((match_nonzero.author_qGI>0)==(match_nonzero['ARID1A.T19_dLFC']>0)).sum()),
      'validation_fixed_nonzero_n':len(match_nonzero),
      'validation_fixed_spearman_descriptive':float(spearmanr(overlap.author_qGI,overlap['ARID1A.T19_dLFC']).statistic) if len(overlap)>2 else None,
      'S1_ARID1A_rows':arid_meta.to_dict('records'),
      'dLFC_sign_rule':'negative author qGI and negative validation mutant-minus-WT dLFC both mean greater depletion under ARID1A loss',
      'qGI_is_not_raw_dLFC':'Author qGI includes nonlinear fitness and network-wide corrections; no exact equality to raw query-minus-WT LFC is assumed.',
      'author_FDR_provenance':'Retained directly from complete qGI_FDR matrix, not recalculated on BAF subset. Author moderated t tests use 2-4 independent guides and multiple WT-control contrasts; those are not multiple ARID1A KO clones.',
    }
    save(out/'data_and_direction_audit.json',sanity)
    go=bool(c.supports_greater_group_A_dependency.any())
    report={
      'pilot_decision':'LIMITED_GO_SUPPLEMENTAL_ONLY' if go else 'NO_GO_FOR_NEW_DEPENDENCY_SPECIFICITY_CLAIM',
      'primary_positive_contrasts':c.loc[c.supports_greater_group_A_dependency,['group_A','group_B']].to_dict('records'),
      'full_independent_module_validation_available':bool(all(x['S6_full_fixed_set_evaluable'] for x in coverage)),
      'message':'One HAP1 clone is external real genetic-perturbation evidence, not NSCLC mechanism or independent biological replication. Hit-selected S6 is available-case guide validation, not unbiased genome-wide replication.',
      'novelty_boundary':'The 2026 source study already analysed GI vs co-essentiality and complex enrichment. This fixed ARID1A structural-panel reanalysis is not a new genome-wide discovery or new algorithm.',
      'frozen_spec_sha256':sha(out/'R15_HAP1_FROZEN_SPEC.json'),
      'analysis_script_sha256':sha(__file__), 'run_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
      'runtime':{'host':platform.node(),'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__},
    }
    save(out/'R15_HAP1_PILOT_DECISION.json',report)
    print(c.to_string(index=False),flush=True)
    print(pd.DataFrame(coverage).to_string(index=False),flush=True)
    print(json.dumps(sanity,indent=2,default=str),flush=True)
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__': main()
