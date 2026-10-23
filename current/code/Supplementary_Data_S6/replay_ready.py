"""Replay full follow-up analysis from the compact processed-matrix snapshot in S6.
Extract both S6 ZIPs into one directory; python replay_ready.py --results results.
This verifies downstream statistics, not raw mass-spectrometry preprocessing.
"""
from proteome_wide_analysis import *

ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True);a=ap.parse_args()
spec=json.loads((a.results/'analysis_plan.json').read_text());report={};testrows=[]
tests=pd.read_csv(a.results/'competitive_tests.csv')
pcold=pd.read_csv(a.results/'PC_sensitivity.csv');infold=pd.read_csv(a.results/'patient_influence.csv')
for k,cohort in enumerate(spec['cohorts']):
    d=a.results/cohort;z=np.load(d/'analysis_input.npz');P,R,C,U=z['P'],z['R'],z['C'],z['purity'];genes=z['genes'].tolist()
    E=np.column_stack([res(P[:,j],np.column_stack([R[:,j],C[:,j],U])) for j in range(len(genes))])
    ai=genes.index('ARID1A');rho=np.array([corr(E[:,ai],E[:,j]) for j in range(len(genes))])
    scan=pd.read_csv(d/'proteome_scan.csv').set_index('gene')
    rhoerr=max(abs(rho[j]-scan.loc[g,'rho']) for j,g in enumerate(genes) if g!='ARID1A');assert rhoerr<1e-12
    cand=np.array([j for j,g in enumerate(genes) if g not in spec['exclude_from_null_and_PC']])
    features=np.column_stack([P.mean(0),P.std(0,ddof=1),R.mean(0),R.std(0,ddof=1)])
    f=np.column_stack([rankdata(features[:,j])/len(genes) for j in range(4)])
    pools=[];oldmatches=pd.read_csv(d/'matched_candidates.csv')
    for g in spec['contact']:
        distances=((f[cand]-f[genes.index(g)])**2).sum(1);pool=cand[np.argsort(distances,kind='stable')[:100]];pools.append(pool)
        assert [genes[j] for j in pool]==oldmatches[oldmatches.contact==g].candidate.tolist()
    rng=np.random.default_rng(spec['seed']+k);obs=float(np.median(rho[[genes.index(g) for g in spec['contact']]]))
    for typ in ['size_only','matched']:
        old=np.load(d/f'{typ}_null_draws.npz');idx=[]
        for b in range(20000):
            if typ=='size_only':ids=rng.choice(cand,6,replace=False)
            else:
                ids=np.array([rng.choice(p) for p in pools])
                while len(set(ids))<6:ids=np.array([rng.choice(p) for p in pools])
            idx.append(ids)
        idx=np.array(idx);assert np.array_equal(idx,old['indices'])
        values=np.median(rho[idx],axis=1);assert np.allclose(values,old['statistics'],atol=1e-12)
        p=(1+sum(values>=obs))/20001;expected=tests[(tests.cohort==cohort)&(tests['null']==typ)].iloc[0]
        assert abs(p-expected.competitive_p)<1e-12;testrows.append(p)
    cols=[genes.index(g) for g in ['ARID1A']+spec['contact']+spec['background']]
    Z=E[:,cand]/E[:,cand].std(0,ddof=1);u,s,_=np.linalg.svd(Z,full_matrices=False);pcs=u[:,:5]*s[:5]
    pcerr=[];inferr=[]
    for n in [0,1,3,5]:
        _,value=panel_est(P,R,C,U,cols,extra=pcs[:,:n] if n else None)
        old=pcold[(pcold.cohort==cohort)&(pcold.n_PC==n)].estimate.iloc[0];pcerr.append(abs(value-old))
    for i,patient in enumerate(z['patients']):
        _,value=panel_est(P,R,C,U,cols,ix=np.delete(np.arange(len(U)),i))
        old=infold[(infold.cohort==cohort)&(infold.omitted_patient==patient)].estimate.iloc[0];inferr.append(abs(value-old))
    assert max(pcerr+inferr)<1e-12
    report[cohort]={'correlation_max_error':rhoerr,'PC_contrast_max_error':max(pcerr),'patient_omission_max_error':max(inferr),'all_40000_seeded_panels_exact':True,'all_600_matched_candidates_exact':True}
assert np.allclose(bh(testrows),tests.BH_q_four_tests,atol=1e-12)
(a.results/'full_replay_audit.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
