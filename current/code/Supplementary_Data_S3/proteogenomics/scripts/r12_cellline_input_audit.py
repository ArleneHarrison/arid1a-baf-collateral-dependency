"""Read-only audit of the existing server inputs. No outcome analyses."""
from pathlib import Path
import importlib
import json
import platform
import csv

BASE = Path('/data/s01011')
RAW = BASE / 'arid1a-nsclc/data/raw/sanger_proteomics'
BULK = BASE / 'ggn_project/02_targets/bulk_resources'
GENES = {'ARID1A', 'ARID1B', 'SMARCC1', 'SMARCC2', 'SMARCA4', 'SMARCA2',
         'SMARCD1', 'SMARCB1', 'SMARCE1', 'DPF2', 'ACTB', 'ACTL6A',
         'PBRM1', 'ARID2', 'BRD7', 'PHF10', 'BRD9', 'BICRA'}

print('PYTHON', platform.python_version())
for name in ['numpy', 'pandas', 'scipy', 'statsmodels', 'sklearn', 'requests', 'matplotlib']:
    try:
        module = importlib.import_module(name)
        print('PACKAGE', name, getattr(module, '__version__', 'unknown'))
    except ImportError:
        print('PACKAGE', name, 'NOT_INSTALLED')

for path in sorted((RAW / 'proteomics_extracted').glob('*')):
    print('PROTEOMICS_FILE', path.name, path.stat().st_size)

protein = RAW / 'proteomics_extracted/Protein_matrix_averaged_20250211.tsv'
with protein.open() as stream:
    reader = csv.reader(stream, delimiter='\t')
    rows = [next(reader) for _ in range(5)]
columns = [i for i, g in enumerate(rows[1]) if g in GENES]
print('PROTEOMICS_HEADER', json.dumps([[r[i] for i in [0, 1] + columns] for r in rows]))
print('AMBIGUOUS_ACTIN_SYMBOLS', [(rows[0][i], g) for i, g in enumerate(rows[1]) if 'ACTB' in g or 'BICRA' in g])

for path in [RAW / 'model_list.csv', BULK / 'depmap_Model.csv']:
    with path.open() as stream:
        reader = csv.DictReader(stream)
        first = next(reader)
        print('METADATA_FIELDS', path.name, list(first))
        print('METADATA_ID_EXAMPLE', path.name, {k: first[k] for k in ['model_id', 'BROAD_ID', 'ModelID', 'OncotreeLineage'] if k in first})

try:
    meta = json.loads((BULK / 'depmap_24q4_meta.json').read_text())
except ValueError:
    import requests
    print('EXISTING_META_INVALID_JSON_RETAINED_UNCHANGED')
    response = requests.get('https://api.figshare.com/v2/articles/27993248/versions/1', timeout=60)
    response.raise_for_status()
    meta = response.json()
print('RELEASE_META_TYPE', type(meta).__name__)
if isinstance(meta, dict):
    print('RELEASE_META_KEYS', list(meta))
    print('RELEASE_TITLE', meta.get('title'), meta.get('id'), meta.get('url'))
    for f in meta.get('files', []):
        name = f.get('name', '')
        if any(x in name.lower() for x in ['expression', 'copy', 'model.csv', 'crisprgeneeffect']):
            print('RELEASE_FILE', json.dumps(f))
else:
    print('RELEASE_META_PREVIEW', str(meta)[:1200])
