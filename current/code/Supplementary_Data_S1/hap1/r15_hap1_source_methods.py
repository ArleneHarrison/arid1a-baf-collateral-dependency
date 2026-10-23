"""Archive official open-access methods relevant to HAP1 interpretation."""
from pathlib import Path
import hashlib
import json
import urllib.request
import xml.etree.ElementTree as ET
import sys

out=Path(sys.argv[1]); out.mkdir(parents=True,exist_ok=True)
url='https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/PMC13281610/unicode'
content=urllib.request.urlopen(url,timeout=90).read()
(out/'source_manuscript_BioC.xml').write_bytes(content)
root=ET.fromstring(content)
terms=['technical replicates','independent gRNA library','log fold','dLFC','different sgRNAs','moderated t-test']
rows=[]
for passage in root.iter('passage'):
    texts=[t.text or '' for t in passage.findall('text')]
    for t in texts:
        if any(term.lower() in t.lower() for term in terms):
            rows.append({'section':[i.text for i in passage.findall('infon')],'text':t})
(out/'source_methods_selected.json').write_text(json.dumps({'url':url,'sha256':hashlib.sha256(content).hexdigest(),'passages':rows},indent=2),encoding='utf-8')
print(json.dumps(rows,indent=2))
