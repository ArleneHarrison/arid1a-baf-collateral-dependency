"""Harmonize existing Sanger proteins with pinned DepMap 24Q4 RNA/CNV.

Run on the user's server. No association or hypothesis tests in this stage.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

ALIASES = {'GLTSCR1': 'BICRA', 'GLTSCR1L': 'BICRAL'}


def digest(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8*1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def subset_matrix(path, genes):
    header = pd.read_csv(path, nrows=0).columns.tolist()
    mapping = {c: ALIASES.get(c.split(' (')[0], c.split(' (')[0]) for c in header[1:]}
    wanted = [c for c in header[1:] if mapping[c] in genes]
    if len({mapping[c] for c in wanted}) != len(wanted):
        raise ValueError(f'Ambiguous selected gene identifiers in {path.name}')
    f = pd.read_csv(path, usecols=[header[0]]+wanted, index_col=0)
    f = f.rename(columns=mapping)
    if not f.index.is_unique:
        raise ValueError(f'Duplicate ModelIDs in {path.name}')
    if not all(str(x).startswith('ACH-') for x in f.index):
        raise ValueError(f'Expected ModelIDs, found profile IDs or other rows in {path.name}')
    return f.apply(pd.to_numeric, errors='raise')


def main(args):
    raw, out = Path(args.raw), Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    spec = json.loads(Path(args.spec).read_text(encoding='utf-8'))
    genes = sorted(set([spec['anchor']] + [g for group in spec['modules'].values() for g in group]))
    expected = {
        'OmicsExpressionProteinCodingGenesTPMLogp1.csv': '71794802b750ce77c422dad0720a40af',
        'OmicsCNGene.csv': '4696565fa74dc805983bc278b4a41802',
        'depmap_Model.csv': '675210d17675f3517b0ce39a3c274f16',
        'depmap_CRISPRGeneEffect.csv': '6edf7ade09b9b34199210b559d4745d3',
    }
    provenance = []
    for name, checksum in expected.items():
        path = raw / name
        if not path.exists():
            if name == 'depmap_CRISPRGeneEffect.csv':
                print('DEPENDENCY_PENDING', flush=True)
                continue
            raise FileNotFoundError(path)
        actual = digest(path, 'md5')
        if actual != checksum:
            raise ValueError(f'Pinned 24Q4 checksum failed: {name}, retained untouched')
        provenance.append({'file':name, 'md5':actual, 'sha256':digest(path), 'bytes':path.stat().st_size})
        print('VERIFIED_INPUT', name, flush=True)
    path = raw / 'Protein_matrix_averaged_20250211.tsv'
    with path.open() as f:
        reader = csv.reader(f, delimiter='\t')
        uniprot, symbols, ids = next(reader), next(reader), next(reader)
    canon = [ALIASES.get(s, s) for s in symbols]
    keep = [i for i, s in enumerate(canon) if s in genes]
    if len(set(canon[i] for i in keep)) != len(keep):
        raise ValueError('Multiple measured protein entries map to one selected symbol; manual audit required')
    protein = pd.read_csv(path, sep='\t', skiprows=3, header=None, usecols=[1]+keep, index_col=0)
    protein.columns = [canon[i] for i in keep]
    protein = protein.apply(pd.to_numeric, errors='raise')
    if not protein.index.is_unique:
        raise ValueError('Duplicate protein source model IDs')
    original_protein_rows = len(protein)
    annotations = pd.read_csv(raw / 'model_list.csv', low_memory=False)
    mapping = annotations[['model_id','BROAD_ID']].dropna().drop_duplicates()
    if (mapping.groupby('model_id').BROAD_ID.nunique() > 1).any():
        raise ValueError('One Sanger model maps to multiple Broad IDs')
    mapping = mapping.set_index('model_id').BROAD_ID
    valid = protein.index.isin(mapping.index)
    missing_map = protein.index[~valid].tolist()
    protein = protein.loc[valid].sort_index()
    links = pd.DataFrame({'source_model_id':protein.index, 'sample_id':protein.index.map(mapping)})
    invalid = ~links.sample_id.astype(str).str.startswith('ACH-')
    if invalid.any():
        raise ValueError('Invalid Broad ModelID in Sanger mapping')
    links['duplicate_broad_id'] = links.sample_id.duplicated(keep='first')
    retained = links.loc[~links.duplicate_broad_id]
    protein = protein.loc[retained.source_model_id]
    protein.index = retained.sample_id.tolist()
    links.to_csv(out / 'r12_sanger_broad_mapping.csv', index=False)
    pd.DataFrame({'source_model_id':missing_map}).to_csv(out / 'r12_unmapped_sanger_models.csv', index=False)
    pd.DataFrame([{'uniprot':uniprot[i], 'source_symbol':symbols[i], 'gene_symbol':canon[i]}
                  for i in keep]).to_csv(out / 'r12_selected_protein_identifiers.csv', index=False)
    rna = subset_matrix(raw / 'OmicsExpressionProteinCodingGenesTPMLogp1.csv', genes)
    cnv = subset_matrix(raw / 'OmicsCNGene.csv', genes)
    if (cnv < 0).any().any():
        raise ValueError('Unexpected negative raw copy ratio; inspect scale before log2(x+1)')
    cnv = np.log2(cnv+1)
    metadata = pd.read_csv(raw / 'depmap_Model.csv', low_memory=False).set_index('ModelID')
    if not metadata.index.is_unique:
        raise ValueError('Duplicate model metadata')
    meta = pd.DataFrame({'patient_id':metadata.PatientID, 'lineage':metadata.OncotreeLineage,
                         'cell_line_name':metadata.CellLineName, 'disease':metadata.OncotreePrimaryDisease,
                         'cohort':'Sanger_DepMap24Q4', 'sample_type':'cancer_cell_line'})
    ids = sorted(set(protein.index) & set(rna.index) & set(cnv.index) & set(meta.index))
    for name, f in [('protein',protein),('rna',rna),('cnv',cnv),('metadata',meta)]:
        f.loc[ids].to_csv(out / f'{name}.csv', index_label='sample_id')
    if (raw / 'depmap_CRISPRGeneEffect.csv').exists():
        dependency = subset_matrix(raw / 'depmap_CRISPRGeneEffect.csv', genes)
        dependency.loc[dependency.index.intersection(ids)].to_csv(out / 'dependency.csv', index_label='sample_id')
    coverage = []
    for g in genes:
        row = {'gene':g}
        for name, f in [('protein',protein),('rna',rna),('cnv',cnv)]:
            row[name+'_n'] = int(f.loc[ids,g].notna().sum()) if g in f else 0
        coverage.append(row)
    pd.DataFrame(coverage).to_csv(out / 'r12_input_coverage.csv', index=False)
    meta.loc[ids].groupby(['lineage','disease'], dropna=False).size().rename('n').to_csv(out/'r12_lineage_disease_crosstab.csv')
    manifest = {'release':'DepMap Public 24Q4 v1', 'proteomics_release':'Sanger averaged 20250211',
                'original_protein_source_models':original_protein_rows,
                'mapped_protein_models':len(protein), 'all_omic_matched_models':len(ids),
                'duplicate_broad_mappings_excluded':int(links.duplicate_broad_id.sum()),
                'protein_scale':'Source averaged normalized values; not re-log-transformed or imputed',
                'rna_scale':'log2(TPM+1), unstranded model-level 24Q4 matrix as supplied',
                'cnv_scale':'log2(relative gene-level copy ratio+1), per official PortalOmicsCNGeneLog2 definition',
                'patient_selection':'Performed by analysis script after covariate filter, deterministic lexical ModelID',
                'alias_mapping':ALIASES, 'source_integrity':provenance,
                'protein_sha256':digest(path), 'mapping_sha256':digest(raw/'model_list.csv'),
                'spec_sha256':digest(args.spec)}
    (out / 'r12_input_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--raw', required=True)
    p.add_argument('--outdir', required=True)
    p.add_argument('--spec', required=True)
    main(p.parse_args())
