"""Server-run exploratory complete-case proteome analysis. All tests are retained."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import argparse,json,hashlib,platform,sys
import numpy as np
import pandas as pd
from scipy.stats import rankdata

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bh(p):
    p=np.asarray(p);i=np.argsort(p);out=np.empty(len(p));out[i]=np.minimum(1,np.minimum.accumulate((p[i]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return out
def res(y,x):
    x=x-x.mean(0);y=y-y.mean()
    beta=np.linalg.lstsq(x,y,rcond=None)[0]
    ux,inv=np.unique(x,axis=0,return_inverse=True)
    return y-(ux@beta)[inv]
def corr(x,y):
    if min(x.std(),y.std())<=1e-10:return np.nan
    a=rankdata(x);b=rankdata(y);a-=a.mean();b-=b.mean();return float(a@b/np.sqrt((a@a)*(b@b)))
def panel_est(P,R,C,U,cols,extra=None,ix=None):
    if ix is None:ix=np.arange(len(U))
    ext=U[:,None] if extra is None else np.column_stack([U,extra])
    E=np.column_stack([res(P[ix,j],np.column_stack([R[ix,j],C[ix,j],ext[ix]])) for j in cols])
    vals=np.array([corr(E[:,0],E[:,j]) for j in range(1,len(cols))])
    return vals,float(np.median(vals[:6])-np.median(vals[6:]))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--spec',type=Path,required=True);args=ap.parse_args()
    spec=json.loads(args.spec.read_text());out=args.out;out.mkdir(parents=True,exist_ok=False)
    (out/'analysis_plan.json').write_text(json.dumps(spec,indent=2))
    contact=spec['contact'];background=spec['background'];panel=['ARID1A']+contact+background
    nullrows=[];sens=[];influence=[];summaries={};hashes={}
    for k,cohort in enumerate(spec['cohorts']):
        root=args.root/'pancan_release';meta_path=root/'analysis_ready'/cohort/'metadata.csv'
        meta=pd.read_csv(meta_path).dropna(subset=['purity']).set_index('sample_id')
        assert len(meta)=={'LUAD':108,'LUSC':107}[cohort] and meta.patient_id.is_unique
        mats={}
        for layer in ['protein','rna','cnv']:
            p=root/'audit'/f'{cohort}_{layer}_tumor_gene_symbol.tsv.gz'
            mats[layer]=pd.read_csv(p,sep='\t',index_col=0).loc[:,meta.index].T
            assert mats[layer].columns.is_unique
            hashes[str(p)]=digest(p)
        hashes[str(meta_path)]=digest(meta_path)
        genes=sorted(set.intersection(*(set(d.columns) for d in mats.values())))
        complete=np.logical_and.reduce([np.isfinite(d[genes].to_numpy()).all(0) for d in mats.values()])
        genes=np.asarray(genes)[complete].tolist()
        P,R,C=[mats[t][genes].to_numpy() for t in ['protein','rna','cnv']];U=meta.purity.to_numpy()
        valid=(P.std(0)>1e-10)&(R.std(0)>1e-10)
        genes=np.array(genes)[valid].tolist();P,R,C=P[:,valid],R[:,valid],C[:,valid]
        assert set(panel)<=set(genes)
        E=np.column_stack([res(P[:,j],np.column_stack([R[:,j],C[:,j],U])) for j in range(len(genes))])
        valid=E.std(0)>1e-10;genes=np.array(genes)[valid].tolist();P,R,C,E=P[:,valid],R[:,valid],C[:,valid],E[:,valid]
        ai=genes.index('ARID1A');rho=np.array([corr(E[:,ai],E[:,j]) for j in range(len(genes))])
        cols=[genes.index(g) for g in panel]
        original=pd.read_csv(args.root/'analysis_results_v3'/f'{cohort}_purity'/'r12_coupling_by_gene.csv').set_index('gene')
        error=max(abs(rho[genes.index(g)]-original.loc[g,'protein_own_adjusted']) for g in panel[1:]);assert error<1e-12,error
        exclude=set(spec['exclude_from_null_and_PC']);cand=np.array([j for j,g in enumerate(genes) if g not in exclude])
        features=np.column_stack([P.mean(0),P.std(0,ddof=1),R.mean(0),R.std(0,ddof=1)])
        f=np.column_stack([rankdata(features[:,j])/len(genes) for j in range(4)])
        pools=[];matches=[]
        for g in contact:
            j=genes.index(g);ds=((f[cand]-f[j])**2).sum(1);pool=cand[np.argsort(ds,kind='stable')[:100]];pools.append(pool)
            for j2 in pool:matches.append(dict(contact=g,candidate=genes[j2],distance=float(np.linalg.norm(f[j2]-f[j])),**{f'percentile_delta_{a}':float(f[j2,a]-f[j,a]) for a in range(4)}))
        d=out/cohort;d.mkdir()
        scan=pd.DataFrame({'gene':genes,'rho':rho,'protein_mean':features[:,0],'protein_sd':features[:,1],'rna_mean':features[:,2],'rna_sd':features[:,3]})
        scan=scan[scan.gene!='ARID1A'].sort_values(['rho','gene'],ascending=[False,True]);scan['rank']=rankdata(-scan.rho,method='min').astype(int);scan['rank_percent']=scan['rank']/len(scan)*100
        scan['panel']=np.select([scan.gene.isin(contact),scan.gene.isin(background),scan.gene.isin(exclude)],['contact','BAF background','other BAF'],default='other')
        scan.to_csv(d/'proteome_scan.csv',index=False)
        pd.DataFrame(matches).to_csv(d/'matched_candidates.csv',index=False)
        obs=float(np.median(rho[[genes.index(g) for g in contact]]));rng=np.random.default_rng(spec['seed']+k)
        nulldraws={}
        for typ in ['size_only','matched']:
            idx=[]
            for b in range(20000):
                if typ=='size_only':ids=rng.choice(cand,6,replace=False)
                else:
                    ids=np.array([rng.choice(pool) for pool in pools])
                    while len(set(ids))<6:ids=np.array([rng.choice(pool) for pool in pools])
                idx.append(ids)
            idx=np.array(idx);v=np.median(rho[idx],axis=1);nulldraws[typ]=v
            nullrows.append(dict(cohort=cohort,null=typ,observed_contact_median=obs,null_median=float(np.median(v)),null_low=float(np.quantile(v,.025)),null_high=float(np.quantile(v,.975)),competitive_p=(1+sum(v>=obs))/(len(v)+1),universe=len(scan),non_BAF_candidates=len(cand)))
            np.savez_compressed(d/f'{typ}_null_draws.npz',indices=idx,gene_order=np.array(genes),statistics=v)
        # Global PCs are a deliberately conservative biological/technical sensitivity, not proof of batch removal.
        Z=E[:,cand]/E[:,cand].std(0,ddof=1)
        u,s,_=np.linalg.svd(Z,full_matrices=False);pcs=u[:,:5]*s[:5]
        pcvar=(s**2)/(s**2).sum()
        pcgene=[]
        for nPC in [0,1,3,5]:
            vals,value=panel_est(P,R,C,U,cols,extra=pcs[:,:nPC] if nPC else None)
            sens.append(dict(cohort=cohort,n_PC=nPC,estimate=value,interpretation='descriptive point estimate; PCs contain both biological and technical variation'))
            for g,v in zip(panel[1:],vals):pcgene.append(dict(gene=g,n_PC=nPC,rho=v))
        for i,patient in enumerate(meta.patient_id):
            vals,value=panel_est(P,R,C,U,cols,ix=np.delete(np.arange(len(U)),i))
            influence.append(dict(cohort=cohort,omitted_patient=patient,estimate=value))
        pd.DataFrame(pcgene).to_csv(d/'PC_gene_sensitivity.csv',index=False)
        pd.DataFrame(pcs,index=meta.index,columns=[f'PC{i+1}' for i in range(5)]).to_csv(d/'global_PC_scores.csv')
        # Compact, exact complete-case numeric input snapshot for reproducibility.
        np.savez_compressed(d/'analysis_input.npz',P=P,R=R,C=C,purity=U,genes=np.array(genes),samples=meta.index.to_numpy(dtype=str),patients=meta.patient_id.to_numpy(dtype=str))
        qc=pd.DataFrame({'sample_id':meta.index,'n_protein_observed':mats['protein'].notna().sum(1),'protein_median':mats['protein'].median(1),'protein_IQR':mats['protein'].quantile(.75,axis=1)-mats['protein'].quantile(.25,axis=1),'purity':U});qc.to_csv(d/'processed_matrix_QC.csv',index=False)
        summaries[cohort]={'patients':len(meta),'ranked_proteins':len(scan),'complete_non_BAF':len(cand),'primary_max_absolute_replay_error':error,'contact_ranks':scan[scan.panel=='contact'][['gene','rho','rank','rank_percent']].to_dict('records'),'PC_variance_fraction':pcvar[:5].tolist()}
        print(cohort,json.dumps(summaries[cohort]),flush=True)
    tests=pd.DataFrame(nullrows);tests['BH_q_four_tests']=bh(tests.competitive_p);tests.to_csv(out/'competitive_tests.csv',index=False)
    pd.DataFrame(sens).to_csv(out/'PC_sensitivity.csv',index=False);pd.DataFrame(influence).to_csv(out/'patient_influence.csv',index=False)
    summary={'host':platform.node(),'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'spec_sha256':digest(args.spec),'script_sha256':digest(Path(__file__)),'input_sha256':hashes,'cohorts':summaries,'tests':tests.to_dict('records'),'PC_sensitivity':sens,'new_patients':0,'status':'exploratory; all specified endpoints retained'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps({'tests':summary['tests'],'PC_sensitivity':sens}),flush=True)

if __name__=='__main__':main()
