"""Server-only bounded staging gate and reproducible R12 cell-line run."""
from pathlib import Path
import os
import subprocess
import sys
import time

ROOT = Path('/root/autodl-tmp/arid1a_baf_r12_cellline')
RAW = ROOT / 'inputs'
expected = {'OmicsCNGene.csv':1394554236,
            'OmicsExpressionProteinCodingGenesTPMLogp1.csv':506628654,
            'Protein_matrix_averaged_20250211.tsv':50924586,
            'depmap_CRISPRGeneEffect.csv':428678699,
            'depmap_Model.csv':645696}
deadline = time.monotonic() + 3600
while True:
    waiting = [name for name, size in expected.items()
               if not (RAW/name).exists() or (RAW/name).stat().st_size != size]
    if not waiting and (RAW/'model_list.csv').exists():
        break
    if time.monotonic() > deadline:
        raise TimeoutError('Staging gate timeout; inputs retained, no inference run')
    print('WAITING_FOR_COMPLETE_TRANSFER', waiting, flush=True)
    time.sleep(30)

env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
subprocess.run([sys.executable, str(ROOT/'r12_protein_multiomics.py'), '--self-test'], check=True, env=env)
subprocess.run([sys.executable, str(ROOT/'r12_prepare_cellline.py'), '--raw', str(RAW),
                '--outdir', str(ROOT/'analysis_ready'), '--spec', str(ROOT/'r12_analysis_spec.json')], check=True, env=env)
subprocess.run([sys.executable, str(ROOT/'r12_protein_multiomics.py'), '--input', str(ROOT/'analysis_ready'),
                '--outdir', str(ROOT/'analysis_results'), '--spec', str(ROOT/'r12_analysis_spec.json'),
                '--cohort', 'Sanger_DepMap24Q4', '--categorical', 'lineage', '--deduplicate-patient'], check=True, env=env)
print('R12_CELLLINE_COMPLETED', flush=True)
