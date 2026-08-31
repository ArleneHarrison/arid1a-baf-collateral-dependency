# -*- coding: utf-8 -*-
"""Recover the 12 references whose remembered PMIDs pointed elsewhere.

Rule enforced here: resolve by TITLE, then require the returned title to match the
query before accepting. No PMID is ever taken on trust.
"""
import json, ssl, time, urllib.parse, urllib.request

ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE
UA = {'User-Agent': 'Mozilla/5.0'}
EUT = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'


def get(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA),
                                        timeout=30, context=ctx) as r:
                return json.loads(r.read().decode())
        except Exception:
            time.sleep(1.0 * (i + 1))
    return None


def esum(pmid):
    j = get(f'{EUT}/esummary.fcgi?db=pubmed&id={pmid}&retmode=json&version=2.0')
    try:
        o = j['result'][pmid]
    except Exception:
        return None
    if o.get('error') or not o.get('title'):
        return None
    doi = next((a['value'] for a in o.get('articleids', []) if a.get('idtype') == 'doi'), '')
    return dict(pmid=pmid, title=o.get('title', '').rstrip('.'),
                authors=[a['name'] for a in o.get('authors', []) if a.get('authtype') == 'Author'],
                year=(o.get('pubdate', '') or '')[:4], journal=o.get('source', ''),
                volume=o.get('volume', ''), issue=o.get('issue', ''),
                pages=o.get('pages', '') or o.get('elocationid', ''), doi=doi)


def norm(s):
    return ''.join(c.lower() for c in (s or '') if c.isalnum())


def resolve(title, must_contain):
    """search by title; accept only if the returned title genuinely matches"""
    for term in (f'{title}[Title]', title):
        j = get(f'{EUT}/esearch.fcgi?db=pubmed&term={urllib.parse.quote(term)}'
                f'&retmode=json&retmax=8&sort=relevance')
        for p in (j or {}).get('esearchresult', {}).get('idlist', []):
            s = esum(p)
            if not s:
                continue
            tl = s['title'].lower()
            if not all(w.lower() in tl for w in must_contain):
                continue
            a, b = norm(title), norm(s['title'])
            if a[:45] in b or b[:45] in a:
                return s
        time.sleep(0.3)
    return None


TARGETS = {
    'mashtalir2018':  ('Modular Organization and Assembly of SWI/SNF Family Chromatin Remodeling Complexes', ['swi/snf']),
    'he2020':         ('Structure of nucleosome-bound human BAF complex', ['baf', 'nucleosome']),
    'kasinath2018':   ('Structures of human PRC2 with its cofactors AEBP2 and JARID2', ['prc2']),
    'kasinath2021':   ('JARID2 and AEBP2 regulate PRC2 in the presence of H2AK119ub1 and other histone modifications', ['prc2']),
    'michel2018':     ('A non-canonical SWI/SNF complex is a synthetic lethal target in cancers driven by BAF complex perturbation', ['swi/snf']),
    'farnaby2019':    ('BAF complex vulnerabilities in cancer demonstrated via structure-based PROTAC design', ['baf']),
    'fedorov2015':    ('Selective targeting of the BRG/PB1 bromodomains impairs embryonic and trophoblast stem cell maintenance', ['bromodomain']),
    'centore2020':    ('Mammalian SWI/SNF Chromatin Remodeling Complexes: Emerging Mechanisms and Therapeutic Strategies', ['swi/snf']),
    'wanior2021':     ('Exploiting vulnerabilities of SWI/SNF chromatin remodelling complexes for cancer therapy', ['swi/snf']),
    'mathur2017':     ('ARID1A loss impairs enhancer-mediated gene regulation and drives colon cancer in mice', ['arid1a']),
    'trizzino2018':   ('The Tumor Suppressor ARID1A Controls Global Transcription via Pausing of RNA Polymerase II', ['arid1a']),
    'pacini2021':     ('Integrated cross-study datasets of genetic dependencies in cancer', ['dependenc']),
}

records = json.load(open('refs_final.json', encoding='utf-8'))
ok = fail = 0
for k, (title, must) in TARGETS.items():
    r = resolve(title, must)
    if r:
        records[k] = r
        ok += 1
        print(f'  OK  {k:15s} PMID {r["pmid"]:9s} {r["year"]}  {r["journal"][:24]:24s} {r["title"][:58]}')
    else:
        fail += 1
        print(f'  !!  {k:15s} UNRESOLVED -> {title[:65]}')
    time.sleep(0.45)


def esc(s):
    return (s or '').replace('&', r'\&').replace('%', r'\%').replace('_', r'\_')


def bibtex(k, r):
    au = ' and '.join(r['authors'][:12]) or 'Anonymous'
    if len(r['authors']) > 12:
        au += ' and others'
    f = [f'  author  = {{{esc(au)}}}', f'  title   = {{{{{esc(r["title"])}}}}}',
         f'  journal = {{{esc(r["journal"])}}}', f'  year    = {{{r["year"]}}}']
    for key, val in (('volume', r['volume']), ('number', r['issue']), ('pages', r['pages'])):
        if val:
            f.append(f'  {key:7s} = {{{val}}}')
    if r['doi']:
        f.append(f'  doi     = {{{r["doi"]}}}')
    if r['pmid']:
        f.append(f'  note    = {{PMID: {r["pmid"]}}}')
    return '@article{' + k + ',\n' + ',\n'.join(f) + '\n}\n'


with open('references.bib', 'w', encoding='utf-8') as fh:
    for k in sorted(records):
        fh.write(bibtex(k, records[k]) + '\n')
json.dump(records, open('refs_final.json', 'w', encoding='utf-8'), ensure_ascii=False, indent=1)

yrs = [int(r['year']) for r in records.values() if (r['year'] or '').isdigit()]
recent = sum(1 for y in yrs if y >= 2021)
print(f'\nrecovered {ok}, still unresolved {fail}')
print(f'=== FINAL: {len(records)} references, >=2021: {recent}/{len(yrs)} = '
      f'{100*recent/len(yrs):.1f}%  {"PASS" if recent/len(yrs)>=.5 else "FAIL"} ===')
