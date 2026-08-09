# -*- coding: utf-8 -*-
"""Re-audit BOTH structures after finding the fragment/copy error in 6C23.

Two bugs found in the PRC2 analysis:
  1. chains sharing a gene were assumed to be COPIES and collapsed with max();
     in 6C23 they are disjoint FRAGMENTS of one chain, so they must be SUMMED
  2. 'contact' was thresholded at > 0 atom pairs; one pair at 4.5 A is noise

Question now: does 6LTJ (the SWI/SNF structure used throughout the manuscript)
have the same fragment situation? The manuscript states SMARCC2 is a HOMODIMER
(chains N and O) and reports 754 = 269 + 485. If N and O are instead fragments of
a single SMARCC2 chain, that statement is wrong and must be corrected.
"""
import itertools
import numpy as np
import pandas as pd
from scipy import stats

SWI = '/data/s01011/arid1a-nsclc/data/prc2_struct/../../data/prc2_struct/6LTJ.pdb'
PRC2 = '/data/s01011/arid1a-nsclc/data/prc2_struct/6C23.pdb'
MAT = '/data/s01011/ggn_project/02_targets/bulk_resources/depmap_CRISPRGeneEffect.csv'
OUT = '/data/s01011/arid1a-nsclc/results/tables/prc2_refinement_CORRECTED.csv'


def load(pdb):
    d = {}
    for l in open(pdb):
        if l.startswith('ATOM') and l[76:78].strip() != 'H':
            c, rn = l[21], int(l[22:26])
            d.setdefault(c, {}).setdefault(rn, []).append(
                [float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return {c: {r: np.array(v) for r, v in x.items()} for c, x in d.items()}


def contacts(A, B, th=4.5):
    allB = np.vstack(list(B.values()))
    iface, npairs = set(), 0
    for rn, c in A.items():
        hit = np.sqrt(((c[:, None, :] - allB[None, :, :]) ** 2).sum(-1)) < th
        if hit.any():
            iface.add(rn)
            npairs += int(hit.sum())
    return iface, npairs


# ---------------------------------------------------------------- 6LTJ audit
import os
swi_path = None
for cand in ['/data/s01011/arid1a-nsclc/data/prc2_struct/6LTJ.pdb',
             '/data/s01011/arid1a-nsclc/data/structures/6LTJ.pdb']:
    if os.path.exists(cand):
        swi_path = cand
        break
if swi_path is None:
    os.system('cd /data/s01011/arid1a-nsclc/data/prc2_struct && '
              'wget -q https://files.rcsb.org/download/6LTJ.pdb -O 6LTJ.pdb')
    swi_path = '/data/s01011/arid1a-nsclc/data/prc2_struct/6LTJ.pdb'

st = load(swi_path)
print('=== 6LTJ: are SMARCC2 chains N and O copies or fragments? ===')
for c in ['N', 'O', 'L']:
    if c in st:
        res = sorted(st[c])
        print(f'  chain {c}: {len(res)} residues, range {res[0]}-{res[-1]}')
if 'N' in st and 'O' in st:
    sN, sO = set(st['N']), set(st['O'])
    ov = len(sN & sO)
    frac = ov / min(len(sN), len(sO))
    print(f'  overlap = {ov} residues ({frac:.0%} of the smaller chain)')
    print('  -> ' + ('SEPARATE COPIES (homodimer): summing is correct'
                     if frac > 0.5 else
                     'DISJOINT FRAGMENTS of one chain: summing is still correct, '
                     'but "homodimer" is the wrong description'))

print('\n=== 6LTJ: ARID1A contacts, per chain and summed ===')
PARTNERS = {'SMARCA4': ['I'], 'SMARCB1': ['M'], 'SMARCC2': ['N', 'O'],
            'SMARCD1': ['P'], 'SMARCE1': ['Q'], 'DPF2': ['R'],
            'ACTL6A': ['J'], 'ACTB': ['K']}
for g, chs in PARTNERS.items():
    tot_res, tot_pairs, per = set(), 0, []
    for c in chs:
        if c not in st:
            continue
        ifc, npr = contacts(st['L'], st[c])
        per.append(f'{c}:{npr}')
        tot_res |= ifc
        tot_pairs += npr
    print(f'  {g:9s} summed {tot_pairs:>5d} pairs, {len(tot_res):>3d} unique iface residues'
          f'   [{", ".join(per)}]')

# ---------------------------------------------------------------- PRC2 corrected
print('\n=== 6C23: corrected contacts (fragments SUMMED, threshold >= 3 iface residues) ===')
p = load(PRC2)
CH = {'SUZ12': ['A', 'M', 'Q'], 'EZH2': ['K', 'C'], 'EED': ['L'],
      'RBBP4': ['N'], 'AEBP2': ['P'], 'JARID2': ['E', 'O']}
cmap = {}
for anchor in ['EZH2', 'SUZ12']:
    cmap[anchor] = {}
    print(f'  anchor {anchor}:')
    for g, chs in CH.items():
        if g == anchor:
            continue
        tot_res, tot_pairs = set(), 0
        for ac in CH[anchor]:
            for pc in chs:
                if ac not in p or pc not in p:
                    continue
                ifc, npr = contacts(p[ac], p[pc])
                tot_res |= ifc
                tot_pairs += npr
        call = 'CONTACT' if len(tot_res) >= 3 else 'no contact'
        cmap[anchor][g] = (tot_pairs, len(tot_res), call)
        print(f'    {g:8s}{tot_pairs:>6d} pairs, {len(tot_res):>3d} iface residues   {call}')

# ---------------------------------------------------------------- module refinement
df = pd.read_csv(MAT, index_col=0, low_memory=False)
df.columns = [c.split(' (')[0] for c in df.columns]
df = df.loc[:, ~df.columns.duplicated()]
CORE = ['EZH2', 'EED', 'SUZ12', 'RBBP4', 'RBBP7']


def rho(anchor, members):
    mem = [g for g in members if g in df.columns and g != anchor]
    if not mem:
        return np.nan
    s = df[mem].mean(axis=1)
    a = df[anchor]
    ok = a.notna() & s.notna()
    return stats.spearmanr(a[ok], s[ok])[0]


print('\n=== module refinement with corrected contact calls ===')
rows = []
for anchor in ['EZH2', 'SUZ12']:
    curated = [g for g in CORE if g != anchor]
    contact_core = [g for g in curated if cmap[anchor].get(g, (0, 0, 'no'))[2] == 'CONTACT']
    print(f'\n  anchor {anchor}')
    r_cur = rho(anchor, curated)
    r_con = rho(anchor, contact_core)
    print(f'    curated              {curated}  rho = {r_cur:+.3f}')
    print(f'    contact-restricted   {contact_core}  rho = {r_con:+.3f}   '
          f'({r_con - r_cur:+.3f})')
    # exhaustive: best achievable subset
    best, bestr = None, -np.inf
    for k in range(1, len(curated) + 1):
        for sub in itertools.combinations(curated, k):
            r = rho(anchor, list(sub))
            if r == r and r > bestr:
                best, bestr = list(sub), r
    print(f'    empirical optimum    {best}  rho = {bestr:+.3f}   ({bestr - r_cur:+.3f})')
    agree = set(best) == set(contact_core)
    print(f'    contact rule == empirical optimum?  {"YES" if agree else "NO"}')
    rows += [dict(anchor=anchor, definition='curated', genes=';'.join(curated), rho=r_cur),
             dict(anchor=anchor, definition='contact_restricted',
                  genes=';'.join(contact_core), rho=r_con),
             dict(anchor=anchor, definition='empirical_optimum',
                  genes=';'.join(best), rho=bestr)]

pd.DataFrame(rows).to_csv(OUT, index=False)
print(f'\nwrote {OUT}')
