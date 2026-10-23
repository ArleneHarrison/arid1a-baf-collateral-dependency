from pathlib import Path
import requests,zipfile,io,json,hashlib
import pandas as pd
root=Path(__file__).resolve().parent/'ipsc';root.mkdir(exist_ok=True)
url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5583477/supplementaryFiles'
r=requests.get(url,timeout=90);r.raise_for_status()
z=zipfile.ZipFile(io.BytesIO(r.content));assert z.testzip() is None
files=[x for x in z.namelist() if Path(x).name in ['mmc7.xlsx','mmc1.pdf']]
for name in files:(root/Path(name).name).write_bytes(z.read(name))
manifest={'url':url,'article':'10.1016/j.celrep.2017.08.010','retrieved_files':[{ 'name':Path(n).name,'bytes':len(z.read(n)),'sha256':hashlib.sha256(z.read(n)).hexdigest()} for n in files]}
(root/'source_manifest.json').write_text(json.dumps(manifest,indent=2))
print(manifest)
book=pd.ExcelFile(root/'mmc7.xlsx');print(book.sheet_names)
for s in book.sheet_names:
    d=pd.read_excel(book,s,header=None,nrows=6);print(s,d.shape);print(d.to_string(index=False,header=False))
