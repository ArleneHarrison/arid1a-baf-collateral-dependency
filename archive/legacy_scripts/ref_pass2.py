# -*- coding: utf-8 -*-
"""Pass 2: repair bad matches, resolve the four not-found, and raise the
share of references from the last five years to >= 50%.

Handbook rules applied:
  - a Crossref hit whose venue is 'Faculty Opinions' is a post-publication review
    stub, never the primary paper -> reject and re-query
  - conference proceedings (NeurIPS) and stats journals are not in PubMed -> Crossref
  - to lift recency, ADD real recent papers tied to existing claims; never invent
"""
import json, ssl, time, urllib.parse, urllib.request

ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE
UA = {'User-Agent': 'Mozilla/5.0'}
EUT = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
BAD_VENUE = ('faculty opinions', 'f1000', 'post-publication')


def get(url, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA),
                                        timeout=30, context=ctx) as r:
                return json.loads(r.read().decode())
        except Exception:
            time.sleep(1.2 * (i + 1))
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


def pm_search(term, retmax=6):
    j = get(f'{EUT}/esearch.fcgi?db=pubmed&term={urllib.parse.quote(term)}&retmode=json&retmax={retmax}&sort=relevance')
    return (j or {}).get('esearchresult', {}).get('idlist', [])


def crossref(title, rows=5):
    u = ('https://api.crossref.org/works?query.bibliographic=' + urllib.parse.quote(title)
         + f'&rows={rows}&select=title,author,issued,container-title,DOI,volume,issue,page,type')
    j = get(u)
    out = []
    for it in (j or {}).get('message', {}).get('items', []):
        ct = (it.get('container-title') or [''])[0]
        if any(b in ct.lower() for b in BAD_VENUE):
            continue
        out.append(dict(pmid='', title=(it.get('title') or [''])[0],
                        authors=[f"{a.get('family','')} {(a.get('given','') or '')[:2]}".strip()
                                 for a in it.get('author', [])[:3]],
                        year=str((it.get('issued', {}).get('date-parts', [[None]])[0] or [''])[0]),
                        journal=ct, volume=it.get('volume', ''), issue=it.get('issue', ''),
                        pages=it.get('page', ''), doi=it.get('DOI', ''), source='crossref'))
    return out


def norm(s):
    return ''.join(c.lower() for c in (s or '') if c.isalnum())


def resolve(tag, title, prefer_crossref=False, want_year=None):
    """title-anchored resolution, PubMed first unless told otherwise"""
    cands = []
    if not prefer_crossref:
        for p in pm_search(f'{title}[Title]', 4) or pm_search(title, 4):
            s = esum(p)
            if s:
                cands.append(s)
    cands += crossref(title)
    tn = norm(title)
    for c in cands:
        cn = norm(c['title'])
        if not cn:
            continue
        if tn[:50] in cn or cn[:50] in tn:
            if want_year and (c['year'] or '').isdigit() and abs(int(c['year']) - want_year) > 2:
                continue
            c['tag'] = tag
            return c
    return None


REPAIR = [
    ('jumper2021',   'Highly accurate protein structure prediction with AlphaFold', False, 2021),
    ('evans2021',    'Protein complex prediction with AlphaFold-Multimer', True, None),
    ('eberhardt2021','AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings', False, 2021),
    ('karczewski2020','The mutational constraint spectrum quantified from variation in 141,456 humans', False, 2020),
    ('alpsoy2018',   'Glioma tumor suppressor candidate region gene 1 (GLTSCR1) and its paralog GLTSCR1-like form SWI/SNF chromatin remodeling subcomplexes', False, 2018),
    ('tcgaluad2014', 'Comprehensive molecular profiling of lung adenocarcinoma', False, 2014),
    # the four not found in Pass 1
    ('pbaf_struct',  'Structure of the human PBAF-nucleosome complex', False, None),
    ('lundberg2017', 'A Unified Approach to Interpreting Model Predictions', True, 2017),
    ('gtex2020',     'The GTEx Consortium atlas of genetic regulatory effects across human tissues', False, 2020),
    ('arid1a_enh',   'ARID1A loss impairs enhancer-mediated gene regulation and drives colon cancer in mice', False, 2017),
]

# recent papers (>=2021) tied to claims already in the manuscript
RECENT = [
    ('swisnf_rev24',  'SWI/SNF chromatin remodelling complexes in cancer therapeutic opportunities'),
    ('arid1a_rev23',  'ARID1A in cancer from mechanisms to therapeutic opportunities'),
    ('degrader2023',  'Discovery of a potent and selective SMARCA2 degrader'),
    ('fhd286_2024',   'FHD-286 a novel SMARCA2 SMARCA4 ATPase inhibitor'),
    ('smarca2_2023',  'Targeting SMARCA2 in SMARCA4-deficient cancers'),
    ('coess2023',     'Co-essentiality networks reveal functional gene modules cancer'),
    ('depmap2024',    'DepMap cancer dependency map integrative analysis 2024'),
    ('crispr_rev22',  'Genome-scale CRISPR screens for cancer target discovery'),
    ('paralog_2021',  'Systematic analysis of paralog synthetic lethality in cancer cells'),
    ('paralog_2021b', 'Paralog knockout profiling identifies DUSP4 and DUSP6 as a digenic dependency'),
    ('prc2_rev2022',  'PRC2 structure and function in development and disease'),
    ('af_bench2023',  'Benchmarking AlphaFold-Multimer for protein-protein interface prediction accuracy'),
    ('iptm_2023',     'Evaluation of AlphaFold-Multimer confidence metrics for protein complex prediction'),
    ('bib_crispr',    'Briefings in Bioinformatics CRISPR screen analysis computational methods review'),
    ('bib_complex',   'Briefings in Bioinformatics protein complex prediction machine learning'),
    ('bib_depend',    'Briefings in Bioinformatics cancer dependency prediction deep learning'),
    ('iscience_baf',  'iScience SWI/SNF chromatin remodeling cancer'),
    ('nsclc_arid1a',  'ARID1A mutation non-small cell lung cancer clinical outcome immunotherapy'),
    ('mr_guide2023',  'Guidelines for performing Mendelian randomization investigations update for summer 2023'),
    ('scfm_2024',     'Assessing the limits of zero-shot foundation models in single-cell biology'),
]

print('=== repairs ===')
fixed, still = [], []
for tag, title, pc, yr in REPAIR:
    r = resolve(tag, title, pc, yr)
    if r:
        fixed.append(r)
        print(f"  OK  [{tag:14s}] {(r['authors'] or ['?'])[0]:20s} {r['year']}  {r['journal'][:42]} ({r['source']})")
    else:
        still.append(dict(tag=tag, query=title))
        print(f"  !!  [{tag:14s}] still unresolved")
    time.sleep(0.4)

print('\n=== recent (>=2021) candidates, topic search ===')
recent_hits = []
for tag, q in RECENT:
    ids = pm_search(f'({q}) AND 2021:2026[dp]', 6)
    best = None
    for p in ids:
        s = esum(p)
        if s and s['year'].isdigit() and int(s['year']) >= 2021:
            best = s
            break
    if best:
        best['tag'] = tag
        recent_hits.append(best)
        print(f"  OK  [{tag:14s}] {(best['authors'] or ['?'])[0]:20s} {best['year']}  "
              f"{best['journal'][:34]:34s} PMID {best['pmid']}")
        print(f"          {best['title'][:95]}")
    else:
        print(f"  --  [{tag:14s}] no 2021+ hit for: {q[:60]}")
    time.sleep(0.4)

json.dump(dict(fixed=fixed, unresolved=still, recent=recent_hits),
          open('refs_pass2.json', 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
print(f'\nrepaired {len(fixed)}, still unresolved {len(still)}, recent found {len(recent_hits)}')
print('wrote refs_pass2.json')
