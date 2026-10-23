"""Exploratory follow-up of archived, shared-row bootstrap distributions.
Run on the user's server. Does not refit or alter the original primary analyses.
"""
from pathlib import Path
import argparse, hashlib, json, platform
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

CONTACT=['SMARCC2','SMARCA4','SMARCD1','SMARCB1','SMARCE1','DPF2']
BACKGROUND=['ACTL6A','PBRM1','ARID2','BRD7','PHF10','BRD9']
PANEL=CONTACT+BACKGROUND
SEED=20260908

def holm(p):
    p=np.asarray(p); order=np.argsort(p); q=np.empty(len(p))
    q[order]=np.minimum(1,np.maximum.accumulate(p[order]*(len(p)-np.arange(len(p)))))
    return q

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    out=a.output;out.mkdir(parents=True,exist_ok=True)
    results={};hashes={};loo=[]
    for cohort in ['LUAD','LUSC']:
        p=a.input/cohort
        by=pd.read_csv(p/'r12_coupling_by_gene.csv').set_index('gene')
        z=np.load(p/'r12_shared_refit_bootstrap.npz')
        genes=z['genes'].tolist(); layers=z['layers'].tolist()
        ix=[genes.index(g) for g in PANEL]; li=layers.index('protein_own_adjusted')
        draws=z['values'][:,ix,li]; point=by.loc[PANEL,'protein_own_adjusted'].to_numpy()
        assert draws.shape==(10000,12) and np.isfinite(draws).all()
        assert by.loc[PANEL,'n'].eq({'LUAD':108,'LUSC':107}[cohort]).all()
        results[cohort]={'point':point,'draws':draws}
        for f in p.iterdir():hashes[cohort+'/'+f.name]=hashlib.sha256(f.read_bytes()).hexdigest()
        for removed in ['none']+PANEL:
            ci=[i for i,g in enumerate(PANEL) if g in CONTACT and g!=removed]
            bi=[i for i,g in enumerate(PANEL) if g in BACKGROUND and g!=removed]
            value=float(np.median(point[ci])-np.median(point[bi]))
            boot=np.median(draws[:,ci],axis=1)-np.median(draws[:,bi],axis=1)
            lo,hi=np.quantile(boot,[.025,.975])
            loo.append(dict(cohort=cohort,removed=removed,estimate=value,ci_low=lo,ci_high=hi,n_bootstrap=len(boot),interpretation='descriptive nominal interval; no multiplicity-adjusted inference'))
        # Audit archived point estimates and percentile intervals against the source tables.
        expected=pd.read_csv(p/'r12_coupling_contrasts.csv').set_index('test').loc['primary_own_adjusted']
        full=next(x for x in loo if x['cohort']==cohort and x['removed']=='none')
        for key in ['estimate','ci_low','ci_high']:
            assert abs(full[key]-expected[key])<1e-12, (cohort,key,full[key],expected[key])
    # Intersection-union positive coupling: both cohorts must individually support rho > 0.
    rows=[]
    for j,g in enumerate(PANEL):
        row={'gene':g,'panel':'contact' if g in CONTACT else 'background'}
        for c in results:
            obs=results[c]['point'][j]; d=results[c]['draws'][:,j]
            p=(1+np.count_nonzero(d-obs>=obs))/(len(d)+1)
            row[c+'_rho']=obs;row[c+'_one_sided_centered_bootstrap_p']=p
        row['intersection_union_p']=max(row[c+'_one_sided_centered_bootstrap_p'] for c in results)
        rows.append(row)
    genes=pd.DataFrame(rows); genes['replicability_p_Holm12']=holm(genes['intersection_union_p'])
    genes['positive_in_both']= (genes.LUAD_rho>0)&(genes.LUSC_rho>0)
    genes['replicated_positive_exploratory']=genes.positive_in_both&(genes.replicability_p_Holm12<.05)
    # Independent cohorts: shuffle one whole bootstrap-row vector, retaining within-cohort dependence.
    a0=results['LUAD'];b0=results['LUSC'];rng=np.random.default_rng(SEED)
    permutation=rng.permutation(10000)
    da=np.median(a0['draws'][:,:6],axis=1)-np.median(a0['draws'][:,6:],axis=1)
    db=np.median(b0['draws'][permutation,:6],axis=1)-np.median(b0['draws'][permutation,6:],axis=1)
    oa=np.median(a0['point'][:6])-np.median(a0['point'][6:]);ob=np.median(b0['point'][:6])-np.median(b0['point'][6:])
    diff=float(oa-ob); bootdiff=da-db
    lo,hi=np.quantile(bootdiff,[.025,.975]); pv=(1+np.count_nonzero(abs(bootdiff-diff)>=abs(diff)))/(len(bootdiff)+1)
    hetero={'contrast':'LUAD minus LUSC contact-background median residual correlation contrast','estimate':diff,'ci_low':float(lo),'ci_high':float(hi),'two_sided_centered_bootstrap_p':float(pv),'n_bootstrap':10000,'seed':SEED,'interpretation':'exploratory difference test; nonsignificance does not establish equivalence'}
    pd.DataFrame(loo).to_csv(out/'leave_one_subunit_out.csv',index=False)
    genes.to_csv(out/'gene_replicability.csv',index=False)
    pd.DataFrame([hetero]).to_csv(out/'histology_difference.csv',index=False)
    np.savez_compressed(out/'derived_bootstrap.npz',histology_difference=bootdiff,LUSC_permutation=permutation)
    summary={'status':'exploratory analysis added after prior primary results were known','new_patients':0,'original_primary_results_unchanged':True,'archived_primary_point_and_CI_replay':'passed within 1e-12','host':platform.node(),'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'seed':SEED,'input_sha256':hashes,'histology_difference':hetero,'replicated_genes':genes.loc[genes.replicated_positive_exploratory,'gene'].tolist(),'leave_one_out':{c:{'point_min':float(min(r['estimate'] for r in loo if r['cohort']==c and r['removed']!='none')),'point_max':float(max(r['estimate'] for r in loo if r['cohort']==c and r['removed']!='none')),'minimum_nominal_lower_CI':float(min(r['ci_low'] for r in loo if r['cohort']==c and r['removed']!='none'))} for c in results},'gene_profile_spearman_descriptive':float(spearmanr(a0['point'],b0['point']).statistic)}
    (out/'analysis_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
