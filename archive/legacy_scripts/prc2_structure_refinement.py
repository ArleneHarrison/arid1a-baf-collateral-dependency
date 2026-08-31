# -*- coding: utf-8 -*-
"""Second complex: does structure-informed module refinement generalize to PRC2?

PRE-REGISTERED PREDICTION (stated before running):
  RBBP4/RBBP7 bind the N-terminal region of SUZ12, not EZH2 directly, and both are
  shared with other complexes (NuRD, CAF-1). By the same logic that made ACTB/ACTL6A
  dilution terms for the ARID1A anchor in SWI/SNF, restricting the PRC2 core module
  to the subunits that structurally contact EZH2 should RAISE the correlation.

Structure: PDB 6C23, human PRC2 with AEBP2 and JARID2 (Kasinath et al.).
  chains -- SUZ12: A,M,Q | EZH2: K,C | EED: L | RBBP4: N | AEBP2: P | JARID2: E,O
  RBBP7 is absent from the structure (paralogue of RBBP4), analogous to SMARCC1
  being absent from 6LTJ.

Method is identical to the SWI/SNF analysis: 4.5 A heavy-atom contact, anchor always
excluded from every module, Spearman against DepMap 24Q4 gene effect.
"""
import itertools
import numpy as np
import pandas as pd
from scipy import stats

PDB = '/data/s01011/arid1a-nsclc/data/prc2_struct/6C23.pdb'
MAT = '/data/s01011/ggn_project/02_targets/bulk_resources/depmap_CRISPRGeneEffect.csv'
OUT = '/data/s01011/arid1a-nsclc/results/tables/prc2_structure_refinement.csv'

CHAINS = {'SUZ12': ['A', 'M', 'Q'], 'EZH2': ['K', 'C'], 'EED': ['L'],
          'RBBP4': ['N'], 'AEBP2': ['P'], 'JARID2': ['E', 'O']}
# curated PRC2 core, as used in the benchmark section
CURATED_CORE = ['EZH2', 'EED', 'SUZ12', 'RBBP4', 'RBBP7']


def load(pdb):
    d = {}
    for l in open(pdb):
        if l.startswith('ATOM') and l[76:78].strip() != 'H':
            c = l[21]
            rn = int(l[22:26])
            d.setdefault(c, {}).setdefault(rn, []).append(
                [float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return {c: {r: np.array(v) for r, v in x.items()} for c, x in d.items()}


def contacts(A, B, th=4.5):
    """(interface residues of A, contacting atom PAIRS) -- same definition as 6LTJ work"""
    allB = np.vstack(list(B.values()))
    iface, npairs = set(), 0
    for rn, c in A.items():
        hit = np.sqrt(((c[:, None, :] - allB[None, :, :]) ** 2).sum(-1)) < th
        if hit.any():
            iface.add(rn)
            npairs += int(hit.sum())
    return iface, npairs


st = load(PDB)
print('chains present:', ''.join(sorted(st)))

print('\n=== step 1: which subunits does each anchor contact in 6C23? ===')
contact_map = {}
for anchor in ['EZH2', 'SUZ12']:
    print(f'\n  anchor = {anchor}')
    contact_map[anchor] = {}
    for partner, pch in CHAINS.items():
        if partner == anchor:
            continue
        best_pairs, best_res = 0, 0
        # take the maximum over all chain-copy combinations
        for ac, pc in itertools.product(CHAINS[anchor], pch):
            if ac not in st or pc not in st:
                continue
            ifc, npr = contacts(st[ac], st[pc])
            if npr > best_pairs:
                best_pairs, best_res = npr, len(ifc)
        contact_map[anchor][partner] = best_pairs
        flag = 'CONTACT' if best_pairs > 0 else 'no contact'
        print(f'    {partner:8s} {best_pairs:>5d} atom pairs, {best_res:>3d} iface residues   {flag}')

print('\n=== step 2: module refinement against DepMap 24Q4 ===')
df = pd.read_csv(MAT, index_col=0, low_memory=False)
df.columns = [c.split(' (')[0] for c in df.columns]
df = df.loc[:, ~df.columns.duplicated()]


def rho(anchor, members):
    mem = [g for g in members if g in df.columns and g != anchor]
    if not mem:
        return np.nan, np.nan, 0
    s = df[mem].mean(axis=1)
    a = df[anchor]
    ok = a.notna() & s.notna()
    r, p = stats.spearmanr(a[ok], s[ok])
    return r, p, len(mem)


rows = []
for anchor in ['EZH2', 'SUZ12']:
    cm = contact_map[anchor]
    curated = [g for g in CURATED_CORE if g != anchor]
    # structurally confirmed contacts among the curated core
    contact_core = [g for g in curated if cm.get(g, 0) > 0]
    # RBBP7 is not in the structure; treat by paralogy with RBBP4
    unresolved = [g for g in curated if g not in cm]

    r_cur, p_cur, k_cur = rho(anchor, curated)
    r_con, p_con, k_con = rho(anchor, contact_core)
    print(f'\n  anchor = {anchor}')
    print(f'    curated core            ({k_cur} genes: {", ".join(curated)})')
    print(f'      rho = {r_cur:+.3f}  P = {p_cur:.2e}')
    print(f'    contact-restricted core ({k_con} genes: {", ".join(contact_core)})')
    print(f'      rho = {r_con:+.3f}  P = {p_con:.2e}')
    print(f'    change = {r_con - r_cur:+.3f}   '
          f'{"PREDICTION HELD" if r_con > r_cur else "prediction NOT held"}')
    if unresolved:
        print(f'    (not resolved in structure, excluded from contact set: {", ".join(unresolved)})')
    rows.append(dict(anchor=anchor, definition='curated_core', genes=';'.join(curated),
                     k=k_cur, rho=r_cur, P=p_cur))
    rows.append(dict(anchor=anchor, definition='contact_restricted_core',
                     genes=';'.join(contact_core), k=k_con, rho=r_con, P=p_con))

    print(f'    individual non-contacting members:')
    for g in curated:
        if g not in contact_core and g in df.columns:
            a, b = df[anchor], df[g]
            ok = a.notna() & b.notna()
            r1, p1 = stats.spearmanr(a[ok], b[ok])
            tag = 'not in structure' if g in unresolved else 'no contact'
            print(f'      {g:8s} pairwise rho vs {anchor} = {r1:+.3f} (P = {p1:.2e})  [{tag}]')
            rows.append(dict(anchor=anchor, definition=f'single_{g}', genes=g, k=1,
                             rho=r1, P=p1))

    print(f'    leave-one-out from curated core:')
    for g in curated:
        others = [x for x in curated if x != g]
        r2, p2, _ = rho(anchor, others)
        cp = cm.get(g, None)
        tag = ('not in structure' if cp is None
               else ('contact' if cp > 0 else 'no contact'))
        print(f'      drop {g:8s} rho = {r2:+.3f}  delta = {r2 - r_cur:+.3f}   [{tag}]')
        rows.append(dict(anchor=anchor, definition=f'LOO_drop_{g}', genes=';'.join(others),
                         k=len(others), rho=r2, P=p2))

pd.DataFrame(rows).to_csv(OUT, index=False)
print(f'\nwrote {OUT}')
