"""Independent intercept-model audit; replay archived null draws without rerandomization."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy.linalg import lstsq

ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True);a=ap.parse_args();report={}
for cohort in ['LUAD','LUSC']:
    root=a.results/cohort;z=np.load(root/'analysis_input.npz');P,R,C,U=z['P'],z['R'],z['C'],z['purity'];genes=z['genes'].tolist()
    e=[]
    for j,g in enumerate(genes):
        x=np.column_stack([np.ones(len(U)),R[:,j],C[:,j],U]);beta=lstsq(x,P[:,j],lapack_driver='gelsy')[0];e.append(P[:,j]-x@beta)
    ai=genes.index('ARID1A');rho=np.array([spearmanr(e[ai],v).statistic for v in e]);scan=pd.read_csv(root/'proteome_scan.csv').set_index('gene')
    err=max(abs(scan.loc[g,'rho']-rho[j]) for j,g in enumerate(genes) if g!='ARID1A');assert err<1e-12,err
    for typ in ['size_only','matched']:
        draws=np.load(root/f'{typ}_null_draws.npz');assert len(draws['statistics'])==20000
        idx=draws['indices'];assert all(len(set(x))==6 for x in idx)
        assert np.allclose(np.median(rho[idx],axis=1),draws['statistics'],atol=1e-12)
    report[cohort]={'ranked_proteins':len(scan),'independent_QR_intercept_max_rho_difference':err,'null_draw_replay':'40000 sets passed; six distinct genes per set'}
(a.results/'independent_audit.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
