# -*- coding: utf-8 -*-
"""Pass 5: claim-source semantic alignment.

For every \\cite in main.tex, pull the surrounding claim sentence and the abstract
of each cited paper, so the pairing can be judged by reading rather than by keyword
overlap (the handbook is explicit that lexical overlap gives many false positives).
"""
import io, json, os, re, ssl, time, urllib.parse, urllib.request

D = (r'D:\桌面\_按项目整理\02_生信数据与医学AI项目\07_医学AI_深度学习_多项目孵化'
     r'\医学AI深度学习_导师课题项目群\范教授课题_NSCLC-SEER_分子对接多项目\纯生信'
     r'\ARID1A_BAF_最终项目\submission')
ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE
UA = {'User-Agent': 'Mozilla/5.0'}
EUT = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'


def fetch_abstract(pmid):
    url = (f'{EUT}/efetch.fcgi?db=pubmed&id={pmid}&rettype=abstract&retmode=text')
    for i in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA),
                                        timeout=30, context=ctx) as r:
                return r.read().decode('utf-8', 'replace')
        except Exception:
            time.sleep(1.0 * (i + 1))
    return ''


tex = io.open(os.path.join(D, 'main.tex'), encoding='utf-8').read()
bib = io.open(os.path.join(D, 'references.bib'), encoding='utf-8').read()

# key -> pmid, title
entries = {}
for blk in re.split(r'@article\{', bib)[1:]:
    key = blk.split(',')[0].strip()
    ti = re.search(r'title\s*=\s*\{\{(.+?)\}\}', blk, re.S)
    pm = re.search(r'PMID:\s*(\d+)', blk)
    entries[key] = dict(title=(ti.group(1).replace('\n', ' ').strip() if ti else ''),
                        pmid=(pm.group(1) if pm else ''))

# strip LaTeX comments, then find each cite with its sentence context
body = re.sub(r'(?<!\\)%.*', '', tex)
claims = []
for m in re.finditer(r'\\cite[a-z]*\{([^}]+)\}', body):
    keys = [k.strip() for k in m.group(1).split(',')]
    start = max(0, m.start() - 420)
    ctxt = body[start:m.end()]
    # keep the last sentence before the citation
    parts = re.split(r'(?<=[.!?])\s+(?=[A-Z\\])', ctxt)
    sent = ' '.join(parts[-2:]) if len(parts) > 1 else ctxt
    sent = re.sub(r'\\[a-zA-Z]+\*?(\[[^\]]*\])?(\{[^}]*\})?', ' ', sent)
    sent = re.sub(r'[{}$\\]', ' ', sent)
    sent = re.sub(r'\s+', ' ', sent).strip()
    claims.append(dict(keys=keys, claim=sent[-300:]))

print(f'{len(claims)} citation sites, {len(entries)} bib entries\n')

abstracts = {}
need = sorted({k for c in claims for k in c['keys']})
for k in need:
    pm = entries.get(k, {}).get('pmid', '')
    if pm:
        txt = fetch_abstract(pm)
        # keep only the abstract body
        blocks = [b.strip() for b in txt.split('\n\n') if len(b.strip()) > 150]
        abstracts[k] = (blocks[0] if blocks else txt)[:1400]
        time.sleep(0.36)
    else:
        abstracts[k] = '[not in PubMed - Crossref-only venue]'

out = io.open(os.path.join(D, '..', 'PASS5_claim_source_review.md'), 'w', encoding='utf-8')
out.write('# Pass 5 — claim-source alignment review sheet\n\n')
out.write(f'{len(claims)} citation sites across {len(entries)} references.\n\n')
for i, c in enumerate(claims, 1):
    out.write(f'## Site {i}\n\n**Claim:** ...{c["claim"]}\n\n')
    for k in c['keys']:
        e = entries.get(k, {})
        out.write(f'- **{k}** (PMID {e.get("pmid","-")}): {e.get("title","?")}\n')
        ab = abstracts.get(k, '')
        ab = re.sub(r'\s+', ' ', ab)[:700]
        out.write(f'  > {ab}\n\n')
    out.write('\n')
out.close()
json.dump(dict(claims=claims, abstracts=abstracts, entries=entries),
          io.open(os.path.join(D, '..', 'pass5_raw.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=1)
print('wrote PASS5_claim_source_review.md and pass5_raw.json')
