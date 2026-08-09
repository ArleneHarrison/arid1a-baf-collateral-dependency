# -*- coding: utf-8 -*-
"""Pass 4b: repair the mis-pointed PMIDs, then audit EVERY record's title against
the topic it is supposed to support (handbook failure modes 2 and 4).
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
    if o.get('error'):
        return None
    doi = next((a['value'] for a in o.get('articleids', []) if a.get('idtype') == 'doi'), '')
    return dict(pmid=pmid, title=o.get('title', '').rstrip('.'),
                authors=[a['name'] for a in o.get('authors', []) if a.get('authtype') == 'Author'],
                year=(o.get('pubdate', '') or '')[:4], journal=o.get('source', ''),
                volume=o.get('volume', ''), issue=o.get('issue', ''),
                pages=o.get('pages', '') or o.get('elocationid', ''), doi=doi)


def by_title(title):
    j = get(f'{EUT}/esearch.fcgi?db=pubmed&term={urllib.parse.quote(title + "[Title]")}'
            f'&retmode=json&retmax=3')
    ids = (j or {}).get('esearchresult', {}).get('idlist', [])
    for p in ids:
        s = esum(p)
        if s:
            return s
    return None


REPAIR = {
    'pacini2021':  'Integrated cross-study datasets of genetic dependencies in cancer',
    'boyle2018':   'High-resolution mapping of cancer cell networks using co-functional interactions',
    'jimenez2018': 'LightDock: a new multi-scale approach to protein-protein docking',
    'lundberg2020': 'From local explanations to global understanding with explainable AI for trees',
}

# key -> words that MUST appear in the real title (semantic anchor)
AUDIT = {
    'kadoch2013': ['swi/snf', 'malignancy'], 'mashtalir2018': ['swi/snf', 'assembly'],
    'centore2020': ['swi/snf'], 'mittal2020': ['swi/snf'], 'wanior2021': ['swi/snf'],
    'remodellers2024': ['remodell'], 'atpremod2024': ['remodell'],
    'helming2014': ['arid1b', 'arid1a'], 'hoffman2014': ['brm', 'brg1'],
    'michel2018': ['swi/snf', 'baf'], 'alpsoy2018': ['gltscr1'],
    'trizzino2018': ['arid1a'], 'mathur2017': ['arid1a'],
    'arid1atarget2021': ['arid1a'], 'arid1ahall2025': ['arid1a'], 'arid1aimmuno2024': ['arid1a'],
    'smarca4nsclc2023': ['smarca4'], 'swisnfnsclc2025': ['swi/snf'], 'brd9isci2026': ['brd9'],
    'he2020': ['baf', 'nucleosome'], 'yuan2022': ['pbaf'],
    'kasinath2018': ['prc2'], 'kasinath2021': ['prc2'], 'prc2struct2021': ['prc2'],
    'fedorov2015': ['bromodomain', 'brg'], 'farnaby2019': ['baf', 'protac'],
    'fhd286trial2025': ['fhd-286'], 'smd3236_2025': ['smarca2'], 'smarca4rev2025': ['smarca4'],
    'behan2019': ['crispr'], 'meyers2017': ['crispr'], 'tsherniak2017': ['dependency'],
    'pacini2021': ['dependenc'], 'depmap2025': ['dependency map'],
    'goncalves2022': ['proteomic'], 'tcgaluad2014': ['lung adenocarcinoma'],
    'kim2020': ['lung adenocarcinoma', 'single-cell'], 'corum2023': ['corum'],
    'pan2018': ['complex', 'fitness'], 'wainberg2021': ['co-essential'],
    'boyle2018': ['co-functional'], 'hassan2023': ['crispr'],
    'ito2021': ['paralog'], 'koferle2022': ['paralog'],
    'love2014': ['deseq2'], 'mayakonda2018': ['maftools'], 'xie2021': ['enrichr'],
    'traag2019': ['leiden'], 'cock2009': ['biopython'], 'jimenez2018': ['lightdock'],
    'jumper2021': ['alphafold'], 'mirdita2022': ['colabfold'], 'eastman2017': ['openmm'],
    'maier2015': ['ff14sb'], 'miller2012': ['mmpbsa'], 'roe2013': ['cpptraj'],
    'theodoris2023': ['transfer learning'], 'geneformerapp': ['geneformer'],
    'hemani2018': ['mr-base'], 'zhu2023af': ['alphafold-multimer'],
    'vosa2021': ['eqtl'], 'mckay2017': ['lung cancer'], 'gtex2020': ['gtex'],
    'karczewski2020': ['constraint'], 'bibinteractome': ['interactome'],
    'bibppi2023': ['protein-protein interaction'],
    'eberhardt2021': ['autodock vina'], 'evans2024': ['alphafold-multimer'],
    'benjamini1995': ['false discovery rate'], 'lundberg2020': ['explain'],
}

records = json.load(open('refs_final.json', encoding='utf-8'))

print('=== repairs ===')
for k, title in REPAIR.items():
    r = by_title(title)
    if r:
        records[k] = r
        print(f'  OK  {k:15s} PMID {r["pmid"]}  {r["year"]}  {r["journal"][:26]:26s} {r["title"][:60]}')
    else:
        print(f'  !!  {k:15s} still unresolved: {title}')
    time.sleep(0.4)
records.pop('lundberg2017', None)

print('\n=== semantic audit: does each title actually match its intended topic? ===')
bad = []
for k, r in sorted(records.items()):
    need = AUDIT.get(k)
    if not need:
        print(f'  ??  {k:16s} no audit anchor defined')
        continue
    t = (r['title'] or '').lower()
    if any(w in t for w in need):
        continue
    bad.append(k)
    print(f'  !!  {k:16s} expected {need}')
    print(f'        got: {r["title"][:95]}')
print(f'\n  mismatches: {len(bad)}  {bad if bad else "-- none --"}')

for k in bad:
    records.pop(k, None)
print(f'  removed {len(bad)} unverifiable entries')


def esc(s):
    return (s or '').replace('&', r'\&').replace('%', r'\%').replace('_', r'\_')


def bibtex(k, r):
    au = ' and '.join(r['authors'][:12]) or 'Anonymous'
    if len(r['authors']) > 12:
        au += ' and others'
    f = [f'  author  = {{{esc(au)}}}',
         f'  title   = {{{{{esc(r["title"])}}}}}',
         f'  journal = {{{esc(r["journal"])}}}',
         f'  year    = {{{r["year"]}}}']
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

yrs = [int(r['year']) for r in records.values() if (r['year'] or '').isdigit()]
recent = sum(1 for y in yrs if y >= 2021)
print(f'\n=== FINAL: {len(records)} references ===')
print(f'>= 2021: {recent}/{len(yrs)} = {100*recent/len(yrs):.1f}%  '
      f'{"PASS" if recent/len(yrs) >= .5 else "FAIL"}')
json.dump(records, open('refs_final.json', 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
print('wrote references.bib + refs_final.json')
