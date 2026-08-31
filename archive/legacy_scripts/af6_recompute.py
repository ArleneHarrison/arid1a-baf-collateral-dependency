"""Recompute AF-Multimer vs cryo-EM interface agreement for ALL SIX ARID1A partners.

Method is identical to the original af_compare.py (4.5 A heavy-atom cutoff, interface
defined on ARID1A residues, AF numbering mapped to 6LTJ by +1638). The original only
covered SMARCB1/SMARCE1/DPF2; the other three were read off a figure.

Step 1 self-check: identify each partner's 6LTJ chain by matching the published
<4.5A contact-atom counts in arid1a_structural_contacts.csv, instead of trusting
the chain IDs asserted in the reconstructed CSV.
"""
import glob, json, os, csv
import numpy as np

BASE = r'F:\spoai服务器医学数据备份\home_z\arid1a_dock'
ARID1A_OFF = 1638
CRYO = os.path.join(BASE, '6LTJ.pdb')

# published contact-atom counts (arid1a_structural_contacts.csv, 6LTJ, <4.5A)
PUBLISHED_CONTACTS = {
    'SMARCC2': 754, 'SMARCA4': 557, 'SMARCD1': 553,
    'SMARCB1': 300, 'SMARCE1': 200, 'DPF2': 137,
}


def load_chains(pdb):
    ch = {}
    for l in open(pdb):
        if l.startswith('ATOM') and l[76:78].strip() != 'H':
            c = l[21]
            rn = int(l[22:26])
            ch.setdefault(c, {}).setdefault(rn, []).append(
                [float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return {c: {r: np.array(v) for r, v in d.items()} for c, d in ch.items()}


def contacts(A, B, th=4.5):
    """Return (interface residues of A, unique contacting atoms of A, contacting atom PAIRS)."""
    allB = np.vstack(list(B.values()))
    iface, natoms, npairs = set(), 0, 0
    for rn, c in A.items():
        d = np.sqrt(((c[:, None, :] - allB[None, :, :]) ** 2).sum(-1))
        hit = d < th
        if hit.any():
            iface.add(rn)
            natoms += int(hit.any(axis=1).sum())
            npairs += int(hit.sum())
    return iface, natoms, npairs


cryo = load_chains(CRYO)
print('6LTJ chains:', ''.join(sorted(cryo)))
L = cryo['L']  # ARID1A
print(f'chain L (ARID1A): {len(L)} residues\n')

print('--- step 1: identify partner chains ---')
chain_profile = {}
for c, d in sorted(cryo.items()):
    if c == 'L':
        continue
    ifc, na, npr = contacts(L, d)
    if na > 0:
        chain_profile[c] = (len(ifc), na, npr)
        print(f'  chain {c}: {len(ifc):>3} iface residues | {na:>4} unique ARID1A atoms '
              f'| {npr:>4} atom PAIRS')

# published "n_contacts_lt_4.5A" matches the atom-PAIR count -> use it to map chains
resolved = {}
for gene, n in PUBLISHED_CONTACTS.items():
    hits = [c for c, (_, _, npr) in chain_profile.items() if npr == n]
    if gene == 'SMARCC2' and not hits:
        resolved[gene] = 'N'   # homodimer N+O = 269+485 = 754; AF is 1:1 -> single copy
        print(f'  {gene:8s} published {n:>4} -> chains N+O (homodimer, '
              f'269+485=754); using N for the 1:1 AF comparison')
        continue
    resolved[gene] = hits[0] if len(hits) == 1 else None
    print(f'  {gene:8s} published {n:>4} -> chain {resolved[gene]}'
          + ('' if resolved[gene] else f'  !! no unique pair-count match {hits}'))

print('\n--- step 2: AF interface vs cryo-EM, all six partners ---')
rows = []
for gene in ['SMARCB1', 'SMARCA4', 'SMARCD1', 'SMARCC2', 'SMARCE1', 'DPF2']:
    ch = resolved[gene]
    if ch is None:
        print(f'{gene}: chain unresolved, skipped')
        continue
    d = os.path.join(BASE, f'af_ARID1A_{gene}')
    pdbs = sorted(glob.glob(os.path.join(d, '*rank_001*.pdb')))
    if not pdbs:
        print(f'{gene}: no rank_001 pdb')
        continue
    af = load_chains(pdbs[0])
    af_if, _, _ = contacts(af['A'], af['B'])          # ARID1A residues, AF numbering
    af_if = {r + ARID1A_OFF for r in af_if}        # -> 6LTJ numbering
    true_if, _, _ = contacts(L, cryo[ch])
    ov = true_if & af_if
    jac = len(ov) / len(true_if | af_if) if (true_if | af_if) else 0.0

    iptm = ptm = None
    for j in sorted(glob.glob(os.path.join(d, '*rank_001*.json'))):
        try:
            o = json.load(open(j))
            iptm = o.get('iptm', o.get('interface_ptm'))
            ptm = o.get('ptm')
            if iptm is not None:
                break
        except Exception:
            pass

    recall = len(ov) / len(true_if) if true_if else 0.0
    print(f'{gene:8s} chain {ch}  ipTM={iptm}  cryoEM_iface={len(true_if):>2}  '
          f'AF_iface={len(af_if):>2}  overlap={len(ov):>2}  recall={recall:.2f}  Jaccard={jac:.3f}')
    rows.append([gene, ch, iptm, ptm, len(true_if), len(af_if), len(ov),
                 round(recall, 3), round(jac, 3)])

out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   'af_interface_vs_cryoem_6partners_RECOMPUTED.csv')
with open(out, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['partner', 'chain_6LTJ', 'af_iptm', 'af_ptm', 'cryoem_iface_n',
                'af_iface_n', 'overlap_n', 'recall', 'jaccard'])
    w.writerows(rows)
print('\nwrote', out)

# does ipTM predict interface accuracy?
if len(rows) >= 3:
    from scipy.stats import spearmanr, pearsonr
    ip = [r[2] for r in rows]
    ja = [r[8] for r in rows]
    rs, ps = spearmanr(ip, ja)
    rp, pp = pearsonr(ip, ja)
    print(f'\nipTM vs Jaccard across n={len(rows)} partners: '
          f'Spearman rho={rs:.3f} (P={ps:.3f}), Pearson r={rp:.3f} (P={pp:.3f})')
