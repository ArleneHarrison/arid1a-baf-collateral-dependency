# -*- coding: utf-8 -*-
"""Pass 3: targeted replacement of the off-topic Pass-2 hits, resolution of the
three remaining unknowns, and journal-restricted search for target-journal papers.

Every hit is printed with its full title so semantic alignment can be judged
(handbook failure mode 2: a real paper that is about the wrong thing).
"""
import json, ssl, time, urllib.parse, urllib.request

ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE
UA = {'User-Agent': 'Mozilla/5.0'}
EUT = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'


def get(url, tries=3):
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
    doi = next((a['value'] for a in o.get('articleids', []) if a.get('idtype') == 'doi'), '')
    return dict(pmid=pmid, title=o.get('title', '').rstrip('.'),
                authors=[a['name'] for a in o.get('authors', [])][:3],
                year=(o.get('pubdate', '') or '')[:4], journal=o.get('source', ''),
                volume=o.get('volume', ''), issue=o.get('issue', ''),
                pages=o.get('pages', '') or o.get('elocationid', ''), doi=doi, source='pubmed')


def search(term, n=8):
    j = get(f'{EUT}/esearch.fcgi?db=pubmed&term={urllib.parse.quote(term)}'
            f'&retmode=json&retmax={n}&sort=relevance')
    return (j or {}).get('esearchresult', {}).get('idlist', [])


def show(tag, term, n=5):
    print(f'\n--- [{tag}]  {term}')
    out = []
    for p in search(term, n):
        s = esum(p)
        if not s:
            continue
        out.append(s)
        print(f'    {s["year"]}  {s["journal"][:26]:26s} PMID {s["pmid"]}')
        print(f'        {s["title"][:110]}')
        time.sleep(0.15)
    return out


QUERIES = [
    # replacements for off-topic Pass-2 hits
    ('swisnf_rev',   '(SWI/SNF OR BAF complex) AND (review[pt]) AND cancer AND 2022:2026[dp]'),
    ('arid1a_rev',   'ARID1A[Title] AND (cancer OR tumour OR tumor) AND review[pt] AND 2021:2026[dp]'),
    ('prc2_rev',     'PRC2[Title] AND (structure OR mechanism) AND review[pt] AND 2021:2026[dp]'),
    ('depmap_use',   '(Cancer Dependency Map OR DepMap)[Title/Abstract] AND 2022:2026[dp]'),
    ('crispr_rev',   'CRISPR screen[Title] AND (review[pt] OR method) AND cancer AND 2021:2026[dp]'),
    # the three still unresolved
    ('pbaf_struct',  'PBAF[Title] AND (nucleosome OR structure) AND (cryo-EM OR cryoEM OR structural)'),
    ('gtex',         'GTEx[Title] AND (genetic regulatory effects OR human tissues)'),
    # target journal: Briefings in Bioinformatics
    ('bib_1',        '"Brief Bioinform"[jour] AND (CRISPR OR dependency) AND 2022:2026[dp]'),
    ('bib_2',        '"Brief Bioinform"[jour] AND (protein complex OR co-essentiality OR module) AND 2022:2026[dp]'),
    ('bib_3',        '"Brief Bioinform"[jour] AND (cancer target OR synthetic lethality) AND 2022:2026[dp]'),
    ('bib_4',        '"Brief Bioinform"[jour] AND (AlphaFold OR protein structure prediction) AND 2022:2026[dp]'),
    # target journal alternative: iScience
    ('isci_1',       '"iScience"[jour] AND (chromatin remodeling OR SWI/SNF) AND 2022:2026[dp]'),
    # extra recent, on-topic
    ('smarcb1_2022', 'SMARCB1[Title] AND (rhabdoid OR loss) AND 2021:2026[dp]'),
    ('arid1a_immu',  'ARID1A[Title] AND lung AND 2021:2026[dp]'),
    ('coessential',  'co-essentiality[Title/Abstract] AND 2021:2026[dp]'),
    ('mr_method',    'Mendelian randomization[Title] AND (guideline OR best practice) AND 2021:2026[dp]'),
]

allres = {}
for tag, term in QUERIES:
    allres[tag] = show(tag, term)
    time.sleep(0.35)

json.dump(allres, open('refs_pass3_raw.json', 'w', encoding='utf-8'),
          ensure_ascii=False, indent=1)
print('\nwrote refs_pass3_raw.json  (titles above are for manual semantic screening)')
