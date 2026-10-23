import os, sys, pickle
os.environ.setdefault('TF_CPP_MIN_LOG_LEVEL','3')
import colabfold.colabfold as cf
from pathlib import Path

NAME = sys.argv[1] if len(sys.argv)>1 else 'ARID1A_SMARCB1'
data = pickle.load(open(f'/home/z/arid1a_dock/msa_{NAME}.pkl','rb'))
seqA, seqB = data['seqA'], data['seqB']

# monkeypatch run_mmseqs2 to return prefetched MSA (fetched on the internet-connected local machine)
def fake_run_mmseqs2(x, prefix, *args, use_pairing=False, **kw):
    return data['paired'] if use_pairing else data['unpaired']
cf.run_mmseqs2 = fake_run_mmseqs2

# input csv (complex = chains joined by ':')
res_dir = Path(f'/home/z/arid1a_dock/af_{NAME}'); res_dir.mkdir(exist_ok=True)
csv = res_dir/'input.csv'
csv.write_text(f'id,sequence\n{NAME},{seqA}:{seqB}\n')

from colabfold.batch import get_queries, run
queries, is_complex = get_queries(str(csv))
print('is_complex=', is_complex, 'nqueries=', len(queries), flush=True)
run(queries=queries, result_dir=str(res_dir), num_models=3, is_complex=is_complex,
    model_type='auto', msa_mode='mmseqs2_uniref_env', pair_mode='unpaired_paired',
    num_recycles=3, rank_by='auto', num_seeds=1, keep_existing_results=False,
    data_dir='/home/z/colabfold_weights', host_url='https://api.colabfold.com', user_agent='arid1a/1.0')
print('AF_DONE', flush=True)
