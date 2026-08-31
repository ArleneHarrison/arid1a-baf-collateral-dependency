# -*- coding: utf-8 -*-
"""Pass 5, corrected abstract retrieval: parse <AbstractText> from the XML record.

The plain-text EFetch format puts the author-affiliation block where a naive
'first long paragraph' heuristic finds it, which is what went wrong on the first
attempt. XML is unambiguous.
"""
import io, json, os, re, ssl, time, urllib.request
import xml.etree.ElementTree as ET

D = (r'D:\桌面\_按项目整理\02_生信数据与医学AI项目\07_医学AI_深度学习_多项目孵化'
     r'\医学AI深度学习_导师课题项目群\范教授课题_NSCLC-SEER_分子对接多项目\纯生信'
     r'\ARID1A_BAF_最终项目')
ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE
UA = {'User-Agent': 'Mozilla/5.0'}
EUT = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'

raw = json.load(io.open(os.path.join(D, 'pass5_raw.json'), encoding='utf-8'))
entries = raw['entries']
pmids = [e['pmid'] for e in entries.values() if e['pmid']]

abstracts = {}
CH = 40
for i in range(0, len(pmids), CH):
    batch = ','.join(pmids[i:i + CH])
    url = f'{EUT}/efetch.fcgi?db=pubmed&id={batch}&retmode=xml'
    xml = None
    for t in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA),
                                        timeout=60, context=ctx) as r:
                xml = r.read()
            break
        except Exception:
            time.sleep(1.5 * (t + 1))
    if not xml:
        continue
    root = ET.fromstring(xml)
    for art in root.iter('PubmedArticle'):
        pid = art.findtext('.//PMID', '')
        segs = []
        for a in art.iter('AbstractText'):
            lbl = a.get('Label')
            txt = ''.join(a.itertext()).strip()
            segs.append(f'{lbl}: {txt}' if lbl else txt)
        abstracts[pid] = re.sub(r'\s+', ' ', ' '.join(segs)).strip()
    time.sleep(0.4)

got = sum(1 for v in abstracts.values() if v)
print(f'abstracts retrieved: {got}/{len(pmids)}')
missing = [k for k, e in entries.items() if e['pmid'] and not abstracts.get(e['pmid'])]
print(f'no abstract in PubMed (structural/data papers, editorials): {missing}')

pm2key = {e['pmid']: k for k, e in entries.items() if e['pmid']}
by_key = {pm2key[p]: a for p, a in abstracts.items() if p in pm2key}

out = io.open(os.path.join(D, 'PASS5_claim_source_review.md'), 'w', encoding='utf-8')
out.write('# Pass 5 - claim-source alignment review sheet\n\n')
out.write(f'{len(raw["claims"])} citation sites, {len(entries)} references, '
          f'{got} abstracts retrieved.\n\n---\n\n')
for i, c in enumerate(raw['claims'], 1):
    out.write(f'## Site {i}\n\n**Claim:** ...{c["claim"]}\n\n')
    for k in c['keys']:
        e = entries.get(k, {})
        ab = by_key.get(k, '')
        if not ab:
            ab = '[no abstract in PubMed - verify manually]'
        out.write(f'- **{k}** (PMID {e.get("pmid","-")}) {e.get("title","?")}\n')
        out.write(f'  > {ab[:850]}\n\n')
    out.write('\n')
out.close()
json.dump(by_key, io.open(os.path.join(D, 'pass5_abstracts.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=1)
print('wrote PASS5_claim_source_review.md')
