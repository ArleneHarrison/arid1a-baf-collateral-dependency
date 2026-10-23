"""Independent cell-line NPZ reaggregation and input-defined rare-tie replay.

No production statistical helpers are imported. Scan all 10,000 shared draws
for exactly constant outcome+design within lineage, then replay the first 50
draws plus the first occurrence of every distinct input pattern with explicit
statsmodels dummy OLS. No rank rounding or association-dependent selection.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import numpy as np
import pandas as pd
from r12_cptac_independent_statistics_audit import independent_pair

ROOT=Path('/root/autodl-tmp/arid1a_baf_r12_cellline')
OUT=ROOT/'analysis_results_v3'
SPEC=Path('/root/autodl-tmp/arid1a_baf_r12_cptac/analysis_results_v3/frozen/r12_analysis_spec.json')
TOL=1e-10

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def same(x,y):
    return (np.isnan(x) and np.isnan(y)) or (np.isfinite(x) and np.isfinite(y) and abs(x-y)<=TOL)

def main():
    spec=json.loads(SPEC.read_text())
    assays={m:pd.read_csv(ROOT/'analysis_ready'/f'{m}.csv',index_col=0) for m in ['protein','rna','cnv']}
    md=pd.read_csv(ROOT/'analysis_ready/metadata.csv',index_col=0).sort_index()
    assert md.patient_id.is_unique and len(md)==694
    saved_samples=pd.read_csv(OUT/'r12_analysis_samples.csv',index_col=0)
    assert saved_samples.index.tolist()==md.index.tolist()
    archive=np.load(OUT/'r12_shared_refit_bootstrap.npz')
    boot=archive['values']
    genes,layers=archive['genes'].tolist(),archive['layers'].tolist()
    assert boot.shape==(10000,17,4)
    gene_summary=pd.read_csv(OUT/'r12_coupling_by_gene.csv',index_col=0)
    contrast_summary=pd.read_csv(OUT/'r12_coupling_contrasts.csv',index_col=0)
    arrays={}
    points=np.full((len(genes),4),np.nan)
    failures=[]
    comparisons=[]
    for gi,gene in enumerate(genes):
        if any(gene not in v for v in assays.values()):
            continue
        a=np.column_stack([assays['protein'].loc[md.index,'ARID1A'],assays['protein'].loc[md.index,gene],
                           assays['rna'].loc[md.index,'ARID1A'],assays['rna'].loc[md.index,gene],
                           assays['cnv'].loc[md.index,'ARID1A'],assays['cnv'].loc[md.index,gene]])
        arrays[gene]=a
        if np.isfinite(a).all(axis=1).sum()>=40:
            estimates,_=independent_pair(a,md.lineage.to_numpy())
            points[gi]=[estimates[layer] for layer in layers]
    for gi,gene in enumerate(genes):
        for li,layer in enumerate(layers):
            valid=boot[:,gi,li][np.isfinite(boot[:,gi,li])]
            bounds=np.quantile(valid,[.025,.975]) if len(valid)>=9500 and np.isfinite(points[gi,li]) else [np.nan,np.nan]
            fields={layer:points[gi,li],layer+'_ci_low':bounds[0],layer+'_ci_high':bounds[1],layer+'_valid_bootstraps':len(valid)}
            for field,value in fields.items():
                expected=float(gene_summary.loc[gene,field])
                ok=same(value,expected)
                comparisons.append({'table':'by_gene','key':gene,'field':field,'independent':value,'saved':expected,'pass':ok})
                if not ok:
                    failures.append(comparisons[-1])
    def contrast(values,contact,background,layer):
        if not contact or not background or any(g not in genes for g in contact+background):
            return np.full(values.shape[0],np.nan)
        c=values[:,[genes.index(g) for g in contact],layer]
        b=values[:,[genes.index(g) for g in background],layer]
        out=np.full(values.shape[0],np.nan)
        ok=np.isfinite(c).all(axis=1)&np.isfinite(b).all(axis=1)
        out[ok]=np.median(c[ok],axis=1)-np.median(b[ok],axis=1)
        return out
    for test,row in contrast_summary.iterrows():
        contact=str(row.contact_genes).split(';') if pd.notna(row.contact_genes) else []
        background=str(row.background_genes).split(';') if pd.notna(row.background_genes) else []
        blocked=test=='sensitivity_high_detection90' and (len(contact)<3 or len(background)<2)
        layer=3 if test=='sensitivity_joint_adjusted' else 2
        if blocked:
            values=np.full(10000,np.nan)
            point=np.nan
        elif test=='primary_raw_layer_difference':
            values=contrast(boot,contact,background,0)-contrast(boot,contact,background,1)
            point=float((contrast(points[None,:,:],contact,background,0)-contrast(points[None,:,:],contact,background,1))[0])
        else:
            values=contrast(boot,contact,background,layer)
            point=float(contrast(points[None,:,:],contact,background,layer)[0])
        valid=values[np.isfinite(values)]
        bounds=np.quantile(valid,[.025,.975]) if len(valid)>=9500 and np.isfinite(point) else [np.nan,np.nan]
        p=(1+np.count_nonzero(abs(valid-point)>=abs(point)))/(1+len(valid)) if row.primary and np.isfinite(point) and len(valid)>=9500 else np.nan
        for field,value in {'estimate':point,'ci_low':bounds[0],'ci_high':bounds[1],
                            'centered_bootstrap_approx_p':p,'n_valid_bootstrap':len(valid)}.items():
            expected=float(row[field])
            ok=same(value,expected)
            comparisons.append({'table':'contrasts','key':test,'field':field,'independent':value,'saved':expected,'pass':ok})
            if not ok:
                failures.append(comparisons[-1])
    pd.DataFrame(comparisons).to_csv(OUT/'r12_cellline_all_npz_summary_audit.csv',index=False)
    print('ALL_NPZ_SUMMARIES_CHECKED',len(comparisons),'failures',len(failures),flush=True)

    # Replay selections are determined solely by exact input configurations.
    primary=spec['contact_experimental']+spec['primary_baf_background']
    labels=md.lineage.to_numpy()
    lineage_names=sorted(set(labels))
    strata=[np.flatnonzero(labels==label) for label in lineage_names]
    complete={g:np.isfinite(arrays[g]).all(axis=1) for g in primary}
    model_columns={'own_anchor':[0,2,4],'own_partner':[1,3,5],
                   'joint_anchor':[0,2,3,4,5],'joint_partner':[1,2,3,4,5]}
    patterns={}
    events=[]
    selected={(draw,gene) for draw in range(50) for gene in primary}
    rng=np.random.default_rng(spec['seed'])
    draws=[]
    for draw in range(10000):
        sampled=[rng.choice(s,len(s),replace=True) for s in strata]
        draws.append(np.concatenate(sampled))
        for gene in primary:
            for label,source in zip(lineage_names,sampled):
                source=source[complete[gene][source]]
                if len(source)<2:
                    continue
                constant=np.all(arrays[gene][source]==arrays[gene][source[0]],axis=0)
                if constant.sum()<3:
                    continue
                source_key=tuple(np.unique(source).tolist())
                for model,columns in model_columns.items():
                    if not constant[columns].all():
                        continue
                    key=(gene,label,model,source_key,len(source))
                    if key not in patterns:
                        patterns[key]={'pattern_id':len(patterns),'gene':gene,'lineage':label,'model':model,
                                       'source_ids':';'.join(md.index[list(source_key)]),'n_repeated_rows':len(source),
                                       'first_draw':draw,'n_events':0}
                        selected.add((draw,gene))
                    patterns[key]['n_events']+=1
                    events.append((draw,patterns[key]['pattern_id']))
        if (draw+1)%2000==0:
            print('SCANNED',draw+1,'patterns',len(patterns),flush=True)
    pd.DataFrame(patterns.values()).to_csv(OUT/'r12_cellline_constant_input_patterns.csv',index=False)
    np.savez_compressed(OUT/'r12_cellline_constant_input_events.npz',draw_pattern_pairs=np.asarray(events,dtype=np.int32))
    replay=[]
    for draw,gene in sorted(selected):
        ix=draws[draw]
        estimates,_=independent_pair(arrays[gene][ix],labels[ix])
        for layer,value in estimates.items():
            saved=float(boot[draw,genes.index(gene),layers.index(layer)])
            ok=same(value,saved)
            record={'draw':draw,'gene':gene,'layer':layer,'independent_rho':value,'saved_rho':saved,
                    'absolute_difference':abs(value-saved),'pass':ok}
            replay.append(record)
            if not ok:
                failures.append(record)
    pd.DataFrame(replay).to_csv(OUT/'r12_cellline_input_pattern_OLS_replay.csv',index=False)
    report={'date_utc':datetime.now(timezone.utc).isoformat(),'status':'PASS' if not failures else 'FAIL_REQUIRES_INVESTIGATION',
            'production_helpers_imported':False,'model':'Explicit statsmodels OLS intercept+lineage dummies+own or common RNA/CNV',
            'n_all_npz_summary_fields':len(comparisons),'n_genes_in_summary':len(genes),'n_contrasts':len(contrast_summary),
            'all_draws_input_scan':10000,'scan_rule':'At least two pair-complete rows within lineage, with EXACTLY equal y and all model-specific numeric X; constants tested independently for both own and both joint models',
            'replay_selection':'First 50 shared draws for all 12 primary genes plus first occurrence of every gene/lineage/model/sourceIDs/repeat-count pattern; defined from inputs, not association values',
            'n_constant_patterns':len(patterns),'n_constant_model_group_events':len(events),
            'n_selected_gene_draw_pairs':len(selected),'n_OLS_rho_replayed':len(replay),
            'maximum_replay_absolute_difference':max(r['absolute_difference'] for r in replay),
            'failures':failures,'spec_sha256':sha(SPEC),'original_npz_sha256':sha(OUT/'r12_shared_refit_bootstrap.npz'),
            'script_sha256':sha(Path(__file__)),
            'output_sha256':{name:sha(OUT/name) for name in ['r12_cellline_all_npz_summary_audit.csv','r12_cellline_constant_input_patterns.csv',
                           'r12_cellline_constant_input_events.npz','r12_cellline_input_pattern_OLS_replay.csv']}}
    (OUT/'r12_cellline_extended_numerical_audit.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({k:v for k,v in report.items() if k not in ['failures','output_sha256']}),flush=True)
    if failures:
        raise SystemExit('Independent audit discrepancies require investigation; never reconcile numbers')

if __name__=='__main__':
    main()
