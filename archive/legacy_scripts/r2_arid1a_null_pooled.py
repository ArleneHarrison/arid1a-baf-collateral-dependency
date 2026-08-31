# -*- coding: utf-8 -*-
"""Does the ARID1A-null state shift dependency toward PBAF/ncBAF?

The manuscript's genotype-conditional test splits 1,178 cell lines into five
per-tissue comparisons of 6-12 mutant lines each, then counts how many tissues
point the right way. That design throws away most of the panel and no comparison
survives correction, so it can neither confirm nor refute anything.

Two changes here:

  1. Define the null state by ARID1A PROTEIN rather than by somatic-mutation
     calls. Protein loss is the state that matters mechanistically and captures
     truncating mutations, deep deletions and silencing in one measurement.
  2. Standardise each gene's effect WITHIN tissue, then pool. Lineage is removed
     the same way the per-tissue design removes it, but the test uses every line
     at once instead of five small ones.

SMARCE1 is included: it is the top-ranked co-dependency in the manuscript and is
absent from the archived genotype table altogether.
"""
import io
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
BULK = '/data/s01011/ggn_project/02_targets/bulk_resources'
P = '/data/s01011/arid1a-nsclc/data/raw/sanger_proteomics'
G = '/data/s01011/00_Raw_Data/GDSC/GDSC_DrugSensitivity'
OUT = '/data/s01011/arid1a_module_generalisation'

SUBUNITS = {
    'SMARCE1': 'core (contacts ARID1A)', 'SMARCB1': 'core (contacts)',
    'SMARCC1': 'core (contacts)', 'SMARCC2': 'core (contacts)',
    'SMARCD1': 'core (contacts)', 'SMARCD2': 'core (paralogue)',
    'DPF2': 'cBAF (contacts)', 'SMARCA4': 'ATPase (drug target)',
    'SMARCA2': 'ATPase (drug target)', 'ARID1B': 'paralogue (synthetic lethal)',
    'ACTB': 'core (no contact)', 'ACTL6A': 'core (no contact)',
    'PBRM1': 'PBAF', 'ARID2': 'PBAF', 'BRD7': 'PBAF', 'PHF10': 'PBAF',
    'BRD9': 'ncBAF', 'BICRA': 'ncBAF', 'BICRAL': 'ncBAF', 'SS18': 'ncBAF',
}


def log(*a):
    print(*a, flush=True)


# ------------------------------------------------------------------- inputs
raw = pd.read_csv(os.path.join(P, 'proteomics_extracted',
                               'Protein_matrix_averaged_20250211.tsv'),
                  sep='\t', header=None, low_memory=False)
sym = raw.iloc[1, 2:].tolist()
dat = raw.iloc[3:].reset_index(drop=True)
prot = dat.iloc[:, 2:].apply(pd.to_numeric, errors='coerce')
prot.columns = sym
prot.index = dat.iloc[:, 1].values
prot = prot.loc[:, ~prot.columns.duplicated()]
ml = pd.read_csv(os.path.join(P, 'model_list.csv'), low_memory=False)
m2b = ml.dropna(subset=['BROAD_ID']).set_index('model_id')['BROAD_ID'].to_dict()
prot.index = [m2b.get(i) for i in prot.index]
prot = prot[[isinstance(i, str) and i.startswith('ACH-') for i in prot.index]]
prot = prot[~prot.index.duplicated()]

dep = pd.read_csv(os.path.join(BULK, 'depmap_CRISPRGeneEffect.csv'), index_col=0)
dep.columns = [c.split(' (')[0] for c in dep.columns]
dep = dep.loc[:, ~dep.columns.duplicated()]

mod = pd.read_csv(os.path.join(BULK, 'depmap_Model.csv'), low_memory=False)
idc = 'ModelID' if 'ModelID' in mod.columns else mod.columns[0]
tis = mod.set_index(idc)['OncotreeLineage']

# second, larger-n definition for a sensitivity check
hdr = pd.read_csv(os.path.join(G, 'CCLE_expression.csv'), nrows=0).columns.tolist()
acol = [c for c in hdr if c.split(' (')[0] == 'ARID1A'][0]
expr = pd.read_csv(os.path.join(G, 'CCLE_expression.csv'),
                   usecols=[hdr[0], acol]).set_index(hdr[0])[acol]

log(f'proteome {prot.shape[0]} | CRISPR {dep.shape[0]} | expression {len(expr)}')


def zscore_within_tissue(frame):
    """Standardise each gene within lineage, so pooling cannot import lineage."""
    t = pd.Series({i: tis.get(i) for i in frame.index})
    out = frame.copy().astype(float)
    for lineage, idx in t.dropna().groupby(t.dropna()).groups.items():
        idx = [i for i in idx if i in out.index]
        if len(idx) < 8:                      # too few to standardise reliably
            out.loc[idx] = np.nan
            continue
        blk = out.loc[idx]
        out.loc[idx] = (blk - blk.mean()) / blk.std(ddof=0)
    return out


def run(anchor_values, label, cuts=(0.05, 0.10, 0.15)):
    """Compare dependency of null lines against the rest, pooled after
    within-tissue standardisation."""
    lines = sorted(set(anchor_values.dropna().index) & set(dep.index))
    genes = [g for g in SUBUNITS if g in dep.columns]
    z = zscore_within_tissue(dep.loc[lines, genes])
    a = anchor_values.loc[lines]
    log(f'\n{"="*84}\n{label}: {len(lines)} lines, {len(genes)} subunits, '
        f'{z.notna().any(axis=1).sum()} usable after within-tissue standardisation\n{"="*84}')

    results = []
    for q in cuts:
        cut = a.quantile(q)
        null = a[a <= cut].index
        rest = a[a > cut].index
        rows = []
        for g in genes:
            x = z.loc[[i for i in null if i in z.index], g].dropna()
            y = z.loc[[i for i in rest if i in z.index], g].dropna()
            if len(x) < 10 or len(y) < 30:
                continue
            u, p = stats.mannwhitneyu(x, y, alternative='two-sided')
            rows.append(dict(cut=q, gene=g, role=SUBUNITS[g], n_null=len(x),
                             n_rest=len(y), z_null=x.mean(), z_rest=y.mean(),
                             delta=x.mean() - y.mean(), P=p))
        r = pd.DataFrame(rows)
        if r.empty:
            continue
        # BH within this cut
        r = r.sort_values('P')
        r['q'] = r.P * len(r) / np.arange(1, len(r) + 1)
        r['q'] = r['q'][::-1].cummin()[::-1].clip(upper=1)
        results.append(r)
        log(f'\n-- ARID1A-null = bottom {int(q*100)}%  (n_null = {r.n_null.iloc[0]})')
        log(f'{"subunit":10s}{"role":30s}{"Δz":>8s}{"P":>10s}{"q":>9s}')
        for _, x in r.sort_values('delta').iterrows():
            mark = ' *' if x.q < 0.05 else ''
            log(f'{x.gene:10s}{x.role:30s}{x.delta:>+8.3f}{x.P:>10.4f}{x.q:>9.3f}{mark}')
        log(f'   (Δz < 0 = null lines MORE dependent;  {int((r.q<0.05).sum())} of '
            f'{len(r)} pass BH)')
    return pd.concat(results) if results else pd.DataFrame()


prot_res = run(prot['ARID1A'], 'ARID1A-null defined by PROTEIN (primary)')
expr_res = run(expr, 'ARID1A-low defined by mRNA (larger n, weaker definition)')

os.makedirs(OUT, exist_ok=True)
if not prot_res.empty:
    prot_res.to_csv(os.path.join(OUT, 'r2_null_pooled_protein.csv'), index=False)
if not expr_res.empty:
    expr_res.to_csv(os.path.join(OUT, 'r2_null_pooled_mrna.csv'), index=False)
log(f'\nwritten to {OUT}')
