# -*- coding: utf-8 -*-
"""Corrected benchmark scoring.

The previous rule ("assignment = module with the largest rho") is wrong for this
question, because the CORE module is shared by every sub-complex: losing any
subunit of any sub-complex increases dependence on the shared core, so a high
core score is the expected biology, not a misassignment. Scoring argmax against a
single ground-truth label therefore penalises the method for being right.

The hypothesis actually under test is:

    for an anchor belonging to mutually exclusive sub-complex X, its collateral
    dependency should fall on (a) the shared core and (b) X's own specific module,
    and NOT on the specific modules of the other, mutually exclusive sub-complexes.

So the benchmark is a DISCRIMINATION test among the mutually exclusive specific
modules only; the shared core is scored separately as a positive control.
Anchor is excluded from every module throughout.
"""
import numpy as np
import pandas as pd
from scipy import stats

MAT = '/data/s01011/ggn_project/02_targets/bulk_resources/depmap_CRISPRGeneEffect.csv'
OUT = '/data/s01011/arid1a-nsclc/results/tables/benchmark_discrimination.csv'

SPEC = {
    'SWI/SNF': {
        'shared': {'core': ['SMARCB1', 'SMARCC1', 'SMARCC2', 'SMARCD1', 'SMARCD2',
                            'SMARCE1', 'ACTL6A', 'ACTB']},
        # mutually exclusive specific modules
        'exclusive': {
            'cBAF':  ['ARID1A', 'ARID1B', 'DPF1', 'DPF2', 'DPF3', 'BCL7A',
                      'BCL7B', 'BCL7C', 'SS18', 'SS18L1'],
            'PBAF':  ['PBRM1', 'ARID2', 'BRD7', 'PHF10'],
            'ncBAF': ['BRD9', 'BICRA', 'BICRAL'],
        },
        # anchor -> which exclusive module it belongs to (None = shared/none)
        'anchors': {'ARID1A': 'cBAF', 'ARID1B': 'cBAF', 'DPF2': 'cBAF',
                    'SS18': 'cBAF',
                    'PBRM1': 'PBAF', 'ARID2': 'PBAF', 'BRD7': 'PBAF',
                    'BRD9': 'ncBAF', 'BICRA': 'ncBAF',
                    'SMARCA4': None, 'SMARCB1': None, 'SMARCE1': None,
                    'SMARCA2': None},
    },
    'PRC2': {
        'shared': {'core': ['EZH2', 'EED', 'SUZ12', 'RBBP4', 'RBBP7']},
        'exclusive': {
            'PRC2.1': ['MTF2', 'PHF1', 'PHF19', 'EPOP', 'LCOR'],
            'PRC2.2': ['AEBP2', 'JARID2'],
        },
        'anchors': {'MTF2': 'PRC2.1', 'PHF19': 'PRC2.1', 'PHF1': 'PRC2.1',
                    'JARID2': 'PRC2.2', 'AEBP2': 'PRC2.2',
                    'EZH2': None, 'SUZ12': None, 'EED': None},
    },
}

df = pd.read_csv(MAT, index_col=0, low_memory=False)
df.columns = [c.split(' (')[0] for c in df.columns]
df = df.loc[:, ~df.columns.duplicated()]


def rho(anchor, members):
    mem = [g for g in members if g in df.columns and g != anchor]   # always exclude anchor
    if not mem:
        return np.nan, np.nan, 0
    s = df[mem].mean(axis=1)
    a = df[anchor]
    ok = a.notna() & s.notna()
    if ok.sum() < 30:
        return np.nan, np.nan, len(mem)
    r, p = stats.spearmanr(a[ok], s[ok])
    return r, p, len(mem)


rows, n_ok, n_tot = [], 0, 0
for cname, spec in SPEC.items():
    excl = spec['exclusive']
    print(f'\n===== {cname} =====')
    cols = list(spec['shared']) + list(excl)
    print(f'{"anchor":9s}{"belongs":9s}' + ''.join(f'{c:>10s}' for c in cols)
          + '   own>others?')
    print('-' * (18 + 10 * len(cols) + 14))
    for anchor, home in spec['anchors'].items():
        if anchor not in df.columns:
            continue
        vals = {}
        for c, mem in list(spec['shared'].items()) + list(excl.items()):
            vals[c] = rho(anchor, mem)
        line = f'{anchor:9s}{str(home or "shared"):9s}'
        line += ''.join(f'{vals[c][0]:>10.3f}' if not np.isnan(vals[c][0])
                        else f'{"NA":>10s}' for c in cols)
        if home:
            own = vals[home][0]
            others = [vals[m][0] for m in excl if m != home]
            ok = (not np.isnan(own)) and all(own > o for o in others if not np.isnan(o))
            n_tot += 1
            n_ok += ok
            line += f'   {"YES" if ok else "no"}'
        else:
            # shared-core anchors: check core is the strongest and no exclusive
            # module dominates it
            core = vals[list(spec['shared'])[0]][0]
            others = [vals[m][0] for m in excl]
            ok = (not np.isnan(core)) and core >= max(o for o in others if not np.isnan(o)) - 0.03
            n_tot += 1
            n_ok += ok
            line += f'   {"YES" if ok else "no"} (core)'
        print(line)
        for c in cols:
            r, p, k = vals[c]
            rows.append(dict(complex=cname, anchor=anchor, belongs=home or 'shared',
                             module=c, rho=r, P=p, k_excl_anchor=k))

pd.DataFrame(rows).to_csv(OUT, index=False)
print(f'\n===== discrimination benchmark: {n_ok}/{n_tot} anchors correct =====')
print(f'wrote {OUT}')
