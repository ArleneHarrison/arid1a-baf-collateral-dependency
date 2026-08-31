# -*- coding: utf-8 -*-
"""Cross-complex validation of the collateral-dependency localization framework.

Step 0 (self-check): reproduce the published SWI/SNF anchor table (deep_anchor.csv)
        with this script. If the numbers do not match, nothing downstream is usable.
Step 1: apply the identical procedure to cohesin and to PRC2 -- two modular complexes
        with independently established, mutually exclusive accessory modules.

Module scoring, exactly as in the SWI/SNF analysis:
  * module score per cell line = mean gene effect over module members
  * the anchor gene is ALWAYS excluded from every module before scoring
  * Spearman correlation between anchor gene effect and module score
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

BASE = '/data/s01011/ggn_project/02_targets/bulk_resources'
MAT = os.path.join(BASE, 'depmap_CRISPRGeneEffect.csv')
MODEL = os.path.join(BASE, 'depmap_Model.csv')
OUT = '/data/s01011/arid1a-nsclc/results/tables'

# ----------------------------------------------------------------- complexes
COMPLEXES = {
    'SWI/SNF': {
        'modules': {
            'core':  ['SMARCB1', 'SMARCC1', 'SMARCC2', 'SMARCD1', 'SMARCD2',
                      'SMARCE1', 'ACTL6A', 'ACTB'],
            'cBAF':  ['ARID1A', 'ARID1B', 'DPF1', 'DPF2', 'DPF3', 'BCL7A',
                      'BCL7B', 'BCL7C', 'SS18', 'SS18L1'],
            'PBAF':  ['PBRM1', 'ARID2', 'BRD7', 'PHF10'],
            'ncBAF': ['BRD9', 'BICRA', 'BICRAL'],
        },
        # anchor -> module it is known to belong to ('none' = belongs to no sub-complex)
        'anchors': {'ARID1A': 'cBAF', 'ARID1B': 'cBAF', 'SMARCA4': 'core',
                    'SMARCB1': 'core', 'PBRM1': 'PBAF', 'BRD9': 'ncBAF',
                    'SMARCA2': 'none'},
    },

    # Cohesin. STAG1/STAG2 are mutually exclusive paralogues occupying the same
    # site on the ring -- structurally the same relationship as ARID1A/ARID1B in
    # cBAF, and likewise a documented synthetic-lethal pair.
    'Cohesin': {
        'modules': {
            'ring':      ['SMC1A', 'SMC3', 'RAD21'],
            'STAG':      ['STAG1', 'STAG2'],
            'loader':    ['NIPBL', 'MAU2'],
            'release':   ['WAPL', 'PDS5A', 'PDS5B'],
        },
        'anchors': {'STAG1': 'ring', 'STAG2': 'ring', 'RAD21': 'ring',
                    'SMC3': 'ring', 'NIPBL': 'loader', 'WAPL': 'release'},
    },

    # PRC2. PRC2.1 (MTF2/PHF1/PHF19/EPOP/PALI1) and PRC2.2 (AEBP2/JARID2) are
    # mutually exclusive accessory modules assembled on a shared catalytic core.
    'PRC2': {
        'modules': {
            'core':    ['EZH2', 'EED', 'SUZ12', 'RBBP4', 'RBBP7'],
            'PRC2.1':  ['MTF2', 'PHF1', 'PHF19', 'EPOP', 'LCOR'],
            'PRC2.2':  ['AEBP2', 'JARID2'],
        },
        'anchors': {'EZH2': 'core', 'SUZ12': 'core', 'EED': 'core',
                    'MTF2': 'PRC2.1', 'PHF19': 'PRC2.1',
                    'JARID2': 'PRC2.2', 'AEBP2': 'PRC2.2'},
    },
}

# ----------------------------------------------------------------- load data
print('loading matrix ...', flush=True)
df = pd.read_csv(MAT, index_col=0, low_memory=False)
# columns are "SYMBOL (ENTREZ)"
df.columns = [c.split(' (')[0] for c in df.columns]
df = df.loc[:, ~df.columns.duplicated()]
print(f'  {df.shape[0]} cell lines x {df.shape[1]} genes', flush=True)

lineage = None
if os.path.exists(MODEL):
    m = pd.read_csv(MODEL, low_memory=False)
    idcol = 'ModelID' if 'ModelID' in m.columns else m.columns[0]
    lincol = next((c for c in ('OncotreeLineage', 'lineage', 'Lineage')
                   if c in m.columns), None)
    if lincol:
        lineage = m.set_index(idcol)[lincol].reindex(df.index)
        print(f'  lineage annotation: {lineage.notna().sum()} lines, '
              f'{lineage.nunique()} lineages', flush=True)


def score(anchor, modules):
    """Spearman between anchor gene effect and each module score."""
    if anchor not in df.columns:
        return None
    a = df[anchor]
    out = {}
    for mod, members in modules.items():
        mem = [g for g in members if g in df.columns and g != anchor]
        if not mem:
            out[mod] = (np.nan, np.nan, 0, 0)
            continue
        s = df[mem].mean(axis=1)
        ok = a.notna() & s.notna()
        if ok.sum() < 30:
            out[mod] = (np.nan, np.nan, int(ok.sum()), len(mem))
            continue
        r, p = stats.spearmanr(a[ok], s[ok])
        out[mod] = (r, p, int(ok.sum()), len(mem))
    return out


# ----------------------------------------------------------------- run
rows = []
for cname, spec in COMPLEXES.items():
    print(f'\n===== {cname} =====', flush=True)
    mods = list(spec['modules'])
    header = f'{"anchor":10s}{"truth":9s}' + ''.join(f'{m:>11s}' for m in mods) + '   call   ok'
    print(header, flush=True)
    for anchor, truth in spec['anchors'].items():
        res = score(anchor, spec['modules'])
        if res is None:
            print(f'{anchor:10s} NOT IN MATRIX', flush=True)
            continue
        rs = {m: res[m][0] for m in mods}
        # assignment rule: module with the largest positive rho, provided it is
        # significant; otherwise 'none'
        best, bestr = None, -np.inf
        for m in mods:
            r, p = res[m][0], res[m][1]
            if not np.isnan(r) and r > bestr and p < 0.05 and r > 0.1:
                best, bestr = m, r
        call = best if best else 'none'
        hit = 'YES' if call == truth else ('YES' if (truth == 'none' and call == 'none') else 'no')
        print(f'{anchor:10s}{truth:9s}'
              + ''.join(f'{rs[m]:>11.3f}' if not np.isnan(rs[m]) else f'{"NA":>11s}'
                        for m in mods)
              + f'   {call:7s}{hit}', flush=True)
        for m in mods:
            r, p, n, k = res[m]
            rows.append(dict(complex=cname, anchor=anchor, known_module=truth,
                             module=m, rho=r, P=p, n_lines=n, n_members=k,
                             assigned=call, correct=(call == truth)))

res_df = pd.DataFrame(rows)
os.makedirs(OUT, exist_ok=True)
res_df.to_csv(os.path.join(OUT, 'cross_complex_anchor_validation.csv'), index=False)

# ----------------------------------------------------------------- summary
print('\n===== benchmark summary =====', flush=True)
per = res_df.drop_duplicates(['complex', 'anchor'])[['complex', 'anchor', 'correct']]
for cname in COMPLEXES:
    sub = per[per['complex'] == cname]
    print(f'  {cname:10s} {sub["correct"].sum()}/{len(sub)} anchors assigned to '
          f'their own module', flush=True)
print(f'  {"TOTAL":10s} {per["correct"].sum()}/{len(per)}', flush=True)

# ------------------------------------------- mutually exclusive paralogue pairs
print('\n===== mutually exclusive paralogue pairs =====', flush=True)
print('  (does a documented synthetic-lethal pair still show POSITIVE continuous', flush=True)
print('   co-dependency, as ARID1A~ARID1B does?)', flush=True)
for g1, g2, label in [('ARID1A', 'ARID1B', 'cBAF ARID paralogues'),
                      ('STAG1', 'STAG2', 'cohesin STAG paralogues'),
                      ('SMARCA4', 'SMARCA2', 'BAF ATPase paralogues')]:
    if g1 in df.columns and g2 in df.columns:
        ok = df[g1].notna() & df[g2].notna()
        r, p = stats.spearmanr(df.loc[ok, g1], df.loc[ok, g2])
        print(f'  {label:26s} {g1}~{g2}: rho={r:+.3f}  P={p:.2e}  n={ok.sum()}', flush=True)

print('\nwrote', os.path.join(OUT, 'cross_complex_anchor_validation.csv'), flush=True)
