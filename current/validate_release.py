"""Read-only audit of v1.1.1 assets and frozen numerical claims."""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
import pandas as pd
import fitz

W=Path(__file__).resolve().parent

def main():
    files=W/'metadata/current_file_manifest.json'
    if files.exists():
        for row in json.loads(files.read_text(encoding='utf-8')):
            p=W/row['file']
            assert p.is_file(),str(p)
            assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'],str(p)
    maintex=(W/'manuscript.tex').read_text(encoding='utf-8')
    site=maintex[maintex.index(r'\section{Introduction}'):]
    assert len(re.findall(r'\\begin\{figure\}',maintex))==4
    assert maintex.count(r'\begin{table}')==1
    for x in ['ARCHIVE_DOI','TODO','FIXME','PLACEHOLDER','sn-jnl','npj Systems','BMC Genomics']:
        assert x not in maintex,x
    assert '10.5281/zenodo.23587002' in maintex
    assert '+86 17625308115' in maintex
    assert '62572360' in maintex
    assert 'OpenAI Codex' in maintex
    si=(W/'supplementary_information.tex').read_text(encoding='utf-8')
    assert len(re.findall(r'\\section\*\{Supplementary Figure S\d+\}',si))==20
    for name in ['manuscript','supplementary_information']:
        d=fitz.open(W/f'{name}.pdf')
        text='\n'.join(p.get_text() for p in d)
        assert 'ARID1A' in text
        assert '\ufffd' not in text
        print(name,'pages',len(d))
    toc=pd.read_csv(W/'TOC_source_data.csv').set_index('cohort')
    assert toc['n'].to_dict()=={'LUAD':108,'LUSC':107}
    assert np.allclose(toc.estimate,[.250312,.170164])
    assert np.all(toc.q<.05)
    # When the full current/results tree is present, compare against frozen primary tests.
    matches=list((W/'results').rglob('r12_cptac_external_primary_tests.csv')) if (W/'results').exists() else []
    for p in matches:
        t=pd.read_csv(p)
        assert len(t)==4
        print('Primary external four-test family:',len(t))
        assert t.select_dtypes('number').notna().all().all()
        primary=t[t.test=='primary_own_adjusted'].set_index('cohort')
        for cohort,label,n,est,q in [('CPTAC_LUAD','LUAD',108,.2503120028199338,.0011998800119988001),('CPTAC_LUSC','LUSC',107,.17016399224122725,.0313968603139686)]:
            r=primary.loc[cohort]
            assert int(r.n_patients)==n
            assert np.isclose(r.estimate,est,atol=1e-12)
            assert np.isclose(r.q_external_four_tests,q,atol=1e-12)
            assert np.isclose(toc.loc[label,'estimate'],r.estimate,atol=1e-6)
        assert (t[t.test=='primary_raw_layer_difference'].q_external_four_tests>.05).all()
    for p in (W/'results').rglob('competitive_tests.csv') if (W/'results').exists() else []:
        t=pd.read_csv(p)
        assert len(t)==4
        assert np.allclose(t.BH_q_four_tests,1/20001)
        assert set(t.universe)=={7939,8649}
    for p in (W/'results').rglob('fixed_panel_KO_results.csv') if (W/'results').exists() else []:
        t=pd.read_csv(p);t=t[(t['query']=='ARID1A') & t.gene.isin(['SMARCC2','SMARCA4','SMARCD1','SMARCB1','SMARCE1','DPF2'])]
        assert len(t)==6 and (t.protein_log2FC<0).all()
        assert np.isclose(t.protein_log2FC.median(),-.4681511241411718)
    print('PASS: current files, fixed-panel outputs, data-availability DOI, contact details and declarations')

if __name__=='__main__':main()
