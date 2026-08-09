# -*- coding: utf-8 -*-
"""ARID1A contact profile in 9A0K, using the manuscript's 6LTJ criterion.

9A0K is the deposited integrative model of the ENDOGENOUS human cBAF complex on
the nucleosome (3D-EM volumes + crosslinking-MS + comparative models). It is the
only PDB entry that places SMARCC1 in the assembled complex together with
ARID1A, so it supplies the data point the manuscript currently lists as missing.

Same criterion as 6LTJ: heavy-atom pairs within 4.5 A between the ARID1A chain
and each partner chain, plus the minimum heavy-atom distance.
"""
import io
import os
import sys

import numpy as np
from scipy.spatial import cKDTree

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
D = os.path.dirname(os.path.abspath(__file__))

CHAIN = {  # auth_asym_id -> gene, from the polymer_entity records
    'A': 'ARID1A', '1': 'SMARCC1', '2': 'SMARCC2', '4': 'SMARCA4',
    '6': 'ACTL6A', 'b': 'SMARCB1', 'B': 'ACTB', 'D': 'SMARCD1',
    'E': 'SMARCE1', 'P': 'DPF2',
}

xyz, chains, elems = [], [], []
with open(D + r'\9a0k.cif', encoding='utf-8', errors='replace') as fh:
    for line in fh:
        if not line.startswith(('ATOM ', 'HETATM ')):
            continue
        f = line.split()
        # columns follow the _atom_site loop read out of the file header
        elems.append(f[2])
        xyz.append((float(f[8]), float(f[9]), float(f[10])))
        chains.append(f[13])

xyz = np.asarray(xyz)
chains = np.asarray(chains)
elems = np.asarray(elems)
heavy = elems != 'H'
print(f'原子 {len(xyz)}，其中重原子 {heavy.sum()}')

a_mask = (chains == 'A') & heavy
A = xyz[a_mask]
print(f'ARID1A (链 A) 重原子 {len(A)}\n')
treeA = cKDTree(A)

# manuscript's 6LTJ values, for side-by-side comparison
LTJ = {'SMARCC2': 754, 'SMARCA4': 557, 'SMARCD1': 553, 'SMARCB1': 300,
       'SMARCE1': 200, 'DPF2': 137, 'ACTL6A': 0, 'ACTB': 0, 'SMARCC1': None}

print('=' * 76)
print('ARID1A 在 9A0K（内源 cBAF 整合模型）中的接触谱，判据 4.5 A 重原子对')
print('=' * 76)
print(f'{"伙伴":10s}{"接触原子对":>10s}{"最小距离(A)":>12s}   {"6LTJ 对应值":>12s}')
rows = []
for ch, gene in CHAIN.items():
    if gene == 'ARID1A':
        continue
    m = (chains == ch) & heavy
    if m.sum() == 0:
        continue
    P = xyz[m]
    pairs = treeA.query_ball_point(P, r=4.5)
    n = sum(len(p) for p in pairs)
    dmin = cKDTree(P).query(A)[0].min()
    ref = LTJ.get(gene)
    ref_s = '未解析' if ref is None else str(ref)
    rows.append((gene, n, dmin, ref))
    print(f'{gene:10s}{n:>10d}{dmin:>12.2f}   {ref_s:>12s}')

# rank agreement between the two structures, over subunits present in both
print()
both = [(g, n, r) for g, n, d, r in rows if r is not None]
from scipy import stats
if len(both) >= 4:
    g9 = [n for _, n, _ in both]
    g6 = [r for _, _, r in both]
    rho, p = stats.spearmanr(g9, g6)
    print(f'两结构接触量的秩相关（{len(both)} 个共有亚基）: Spearman rho = {rho:+.3f}, P = {p:.4f}')
    print(f'  共有亚基: {[g for g, _, _ in both]}')
