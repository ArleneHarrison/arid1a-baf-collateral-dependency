# -*- coding: utf-8 -*-
"""Is the dependency a complex-assembly effect? A protein-vs-transcript contrast.

The manuscript's mechanism is that losing ARID1A destabilises the subunits it
contacts, so the residual complex becomes fragile. If that is right, the coupling
between ARID1A and those subunits should be visible at the PROTEIN level beyond
what is already present at the TRANSCRIPT level -- destabilisation is
post-transcriptional -- and the excess should be confined to the subunits ARID1A
actually touches.

Design: the same cell lines, the same partners, two layers.
    delta(g) = rho(ARID1A protein, g protein) - rho(ARID1A mRNA, g mRNA)
Prediction: delta > 0 for contacting subunits, ~0 for non-contacting ones, and
delta tracks interface size. A shared cell-line set is enforced so that sample
composition cannot produce the difference.

Proteomics: Sanger cell-line proteomics (the source already used in the paper).
Transcriptome: CCLE expression, matched through model_list BROAD_ID.
"""
import io
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
G = '/data/s01011/00_Raw_Data/GDSC/GDSC_DrugSensitivity'
P = '/data/s01011/arid1a-nsclc/data/raw/sanger_proteomics'
OUT = '/data/s01011/arid1a_module_generalisation'
os.makedirs(OUT, exist_ok=True)

ANCHOR = 'ARID1A'
# contact counts in 6LTJ, with SMARCC1 from 9A0K; 0 = does not contact ARID1A
CONTACT = {'SMARCC2': 754, 'SMARCA4': 557, 'SMARCD1': 553, 'SMARCB1': 300,
           'SMARCE1': 200, 'SMARCC1': 186, 'DPF2': 137,
           'ACTB': 0, 'ACTL6A': 0, 'PBRM1': 0, 'ARID2': 0, 'BRD7': 0,
           'PHF10': 0, 'BRD9': 0, 'BICRA': 0, 'SMARCA2': 0, 'ARID1B': 0}


def log(*a):
    print(*a, flush=True)


# ------------------------------------------------------------------ proteome
log('reading Sanger proteomics ...')
raw = pd.read_csv(os.path.join(P, 'proteomics_extracted',
                               'Protein_matrix_averaged_20250211.tsv'),
                  sep='\t', header=None, low_memory=False)
symbols = raw.iloc[1, 2:].tolist()            # row 1 carries gene symbols
data = raw.iloc[3:].reset_index(drop=True)    # rows 0-2 are headers
prot = data.iloc[:, 2:].apply(pd.to_numeric, errors='coerce')
prot.columns = symbols
prot.index = data.iloc[:, 1].values           # model_id
prot = prot.loc[:, ~prot.columns.duplicated()]
log(f'  proteome {prot.shape[0]} models x {prot.shape[1]} proteins')

# --------------------------------------------------------------- transcriptome
log('reading CCLE expression ...')
hdr = pd.read_csv(os.path.join(G, 'CCLE_expression.csv'), nrows=0).columns.tolist()
sym2col = {}
for c in hdr[1:]:
    sym2col.setdefault(c.split(' (')[0], c)
need = [ANCHOR] + list(CONTACT)
cols = [hdr[0]] + [sym2col[g] for g in need if g in sym2col]
expr = pd.read_csv(os.path.join(G, 'CCLE_expression.csv'), usecols=cols)
expr = expr.rename(columns={hdr[0]: 'DepMap_ID'}).set_index('DepMap_ID')
expr.columns = [c.split(' (')[0] for c in expr.columns]
log(f'  expression {expr.shape[0]} lines x {expr.shape[1]} genes of interest')

# ------------------------------------------------------------------- mapping
ml = pd.read_csv(os.path.join(P, 'model_list.csv'), low_memory=False)
m2b = ml.dropna(subset=['BROAD_ID']).set_index('model_id')['BROAD_ID'].to_dict()
prot.index = [m2b.get(i) for i in prot.index]
prot = prot[[isinstance(i, str) and i.startswith('ACH-') for i in prot.index]]
prot = prot[~prot.index.duplicated()]
log(f'  proteome models mapped to DepMap IDs: {prot.shape[0]}')

shared = sorted(set(prot.index) & set(expr.index))
log(f'  cell lines with BOTH proteome and transcriptome: {len(shared)}')
prot, expr = prot.loc[shared], expr.loc[shared]

# --------------------------------------------------------------------- test
rows = []
for g, nc in CONTACT.items():
    if g not in prot.columns or g not in expr.columns:
        continue
    # identical line set for both layers, so composition cannot drive delta
    ok = (prot[ANCHOR].notna() & prot[g].notna()
          & expr[ANCHOR].notna() & expr[g].notna())
    if ok.sum() < 60:
        continue
    rp = stats.spearmanr(prot.loc[ok, ANCHOR], prot.loc[ok, g])
    rm = stats.spearmanr(expr.loc[ok, ANCHOR], expr.loc[ok, g])
    rows.append(dict(gene=g, n_contacts=nc, n_lines=int(ok.sum()),
                     rho_protein=rp[0], P_protein=rp[1],
                     rho_mrna=rm[0], P_mrna=rm[1],
                     delta=rp[0] - rm[0],
                     contacting=nc > 0))

res = pd.DataFrame(rows).sort_values('n_contacts', ascending=False)
res.to_csv(os.path.join(OUT, 'protein_vs_mrna_delta.csv'), index=False)

log('\n' + '=' * 78)
log(f'{"gene":10s}{"contacts":>9s}{"n":>6s}{"rho_prot":>10s}{"rho_mRNA":>10s}{"delta":>9s}')
log('=' * 78)
for _, r in res.iterrows():
    mark = ' *' if r.contacting else ''
    log(f'{r.gene:10s}{r.n_contacts:>9.0f}{r.n_lines:>6.0f}'
        f'{r.rho_protein:>10.3f}{r.rho_mrna:>10.3f}{r.delta:>+9.3f}{mark}')

con = res[res.contacting].delta
non = res[~res.contacting].delta
summ = dict(
    n_lines_shared=len(shared),
    n_contacting=int(len(con)), n_noncontacting=int(len(non)),
    median_delta_contacting=float(con.median()),
    median_delta_noncontacting=float(non.median()),
    mannwhitney_two_sided_P=float(
        stats.mannwhitneyu(con, non, alternative='two-sided')[1]),
    spearman_contacts_vs_delta=float(
        stats.spearmanr(res.n_contacts, res.delta)[0]),
    spearman_contacts_vs_delta_P=float(
        stats.spearmanr(res.n_contacts, res.delta)[1]),
)
with open(os.path.join(OUT, 'protein_vs_mrna_summary.json'), 'w') as f:
    json.dump(summ, f, indent=2)
log('\n' + '=' * 78)
for k, v in summ.items():
    log(f'  {k:32s} {v}')
log('=' * 78)
