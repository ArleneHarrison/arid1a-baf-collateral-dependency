# -*- coding: utf-8 -*-
"""Positive control for the pooled null-state test.

The pooled test returns nothing significant for any BAF subunit. Before that can
be reported as a real negative, the same pipeline -- same lines, same
within-tissue standardisation, same statistics -- must be shown to recover a
signal that is known to be there. Two controls:

  A. The continuous axis. Correlating the ARID1A dependency profile against each
     subunit's, within tissue, must reproduce the manuscript's result: core/cBAF
     strongly positive, PBAF/ncBAF null.
  B. A planted effect. Genes known to be selectively essential in a definable
     subset (here, lines with the corresponding paralogue lost) must come out.

If A reproduces and the null-state test still returns nothing, the negative is
about the biology of the conditional axis, not about the pipeline.
"""
import io
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
BULK = '/data/s01011/ggn_project/02_targets/bulk_resources'
P = '/data/s01011/arid1a-nsclc/data/raw/sanger_proteomics'
G = '/data/s01011/00_Raw_Data/GDSC/GDSC_DrugSensitivity'

SUBUNITS = {
    'SMARCE1': 'core (contacts)', 'SMARCB1': 'core (contacts)',
    'SMARCC1': 'core (contacts)', 'SMARCC2': 'core (contacts)',
    'SMARCD1': 'core (contacts)', 'SMARCD2': 'core (paralogue)',
    'DPF2': 'cBAF (contacts)', 'SMARCA4': 'ATPase', 'SMARCA2': 'ATPase',
    'ARID1B': 'paralogue', 'ACTB': 'core (no contact)',
    'ACTL6A': 'core (no contact)', 'PBRM1': 'PBAF', 'ARID2': 'PBAF',
    'BRD7': 'PBAF', 'PHF10': 'PBAF', 'BRD9': 'ncBAF', 'BICRA': 'ncBAF',
    'BICRAL': 'ncBAF', 'SS18': 'ncBAF',
}


def log(*a):
    print(*a, flush=True)


dep = pd.read_csv(os.path.join(BULK, 'depmap_CRISPRGeneEffect.csv'), index_col=0)
dep.columns = [c.split(' (')[0] for c in dep.columns]
dep = dep.loc[:, ~dep.columns.duplicated()]
mod = pd.read_csv(os.path.join(BULK, 'depmap_Model.csv'), low_memory=False)
idc = 'ModelID' if 'ModelID' in mod.columns else mod.columns[0]
tis = mod.set_index(idc)['OncotreeLineage']

# ---------------------------------------------------- A. the continuous axis
log('=' * 78)
log('CONTROL A -- continuous co-dependency, within tissue (must reproduce paper)')
log('=' * 78)
t = pd.Series({i: tis.get(i) for i in dep.index}).dropna()
rows = []
for g in SUBUNITS:
    if g not in dep.columns:
        continue
    rs, ns = [], 0
    for lineage, idx in t.groupby(t).groups.items():
        idx = [i for i in idx if i in dep.index]
        if len(idx) < 20:
            continue
        sub = dep.loc[idx, ['ARID1A', g]].dropna()
        if len(sub) < 20:
            continue
        r = stats.spearmanr(sub['ARID1A'], sub[g])[0]
        if np.isfinite(r):
            rs.append(r)
            ns += len(sub)
    if rs:
        rows.append(dict(gene=g, role=SUBUNITS[g], median_rho=float(np.median(rs)),
                         n_tissues=len(rs), n_lines=ns))
c = pd.DataFrame(rows).sort_values('median_rho', ascending=False)
log(f'{"subunit":10s}{"role":22s}{"median rho":>12s}{"tissues":>9s}')
for _, x in c.iterrows():
    log(f'{x.gene:10s}{x.role:22s}{x.median_rho:>+12.3f}{x.n_tissues:>9d}')

core = c[c.role.str.startswith(('core (contacts', 'cBAF'))].median_rho
other = c[c.role.isin(['PBAF', 'ncBAF'])].median_rho
u, p = stats.mannwhitneyu(core, other, alternative='two-sided')
log(f'\n  core/cBAF median {core.median():+.3f}  vs  PBAF/ncBAF {other.median():+.3f}'
    f'   Mann-Whitney P = {p:.4f}')
log('  -> the pipeline recovers the published separation' if p < 0.05
    else '  -> WARNING: pipeline does NOT recover the published separation')

# ------------------------------------------- B. a planted conditional effect
log('\n' + '=' * 78)
log('CONTROL B -- planted effect: SMARCA2 dependency in SMARCA4-low lines')
log('=' * 78)
log('SMARCA4-deficient cells are known to depend on SMARCA2; the pooled design')
log('must detect this if it can detect anything conditional.')


def zscore_within_tissue(frame):
    tt = pd.Series({i: tis.get(i) for i in frame.index})
    out = frame.copy().astype(float)
    for lineage, idx in tt.dropna().groupby(tt.dropna()).groups.items():
        idx = [i for i in idx if i in out.index]
        if len(idx) < 8:
            out.loc[idx] = np.nan
            continue
        blk = out.loc[idx]
        out.loc[idx] = (blk - blk.mean()) / blk.std(ddof=0)
    return out


hdr = pd.read_csv(os.path.join(G, 'CCLE_expression.csv'), nrows=0).columns.tolist()
col4 = [c for c in hdr if c.split(' (')[0] == 'SMARCA4'][0]
e4 = pd.read_csv(os.path.join(G, 'CCLE_expression.csv'),
                 usecols=[hdr[0], col4]).set_index(hdr[0])[col4]
lines = sorted(set(e4.dropna().index) & set(dep.index))
z = zscore_within_tissue(dep.loc[lines, ['SMARCA2']])
a = e4.loc[lines]
for q in (0.05, 0.10, 0.15):
    cut = a.quantile(q)
    lo = z.loc[[i for i in a[a <= cut].index if i in z.index], 'SMARCA2'].dropna()
    hi = z.loc[[i for i in a[a > cut].index if i in z.index], 'SMARCA2'].dropna()
    if len(lo) < 10:
        continue
    u, p = stats.mannwhitneyu(lo, hi, alternative='two-sided')
    log(f'  SMARCA4 bottom {int(q*100):2d}%  n={len(lo):3d}  '
        f'Δz = {lo.mean()-hi.mean():+.3f}   P = {p:.2e}'
        f'{"  <- detected" if p < 0.05 else "  <- NOT detected"}')
