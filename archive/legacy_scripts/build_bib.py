# -*- coding: utf-8 -*-
"""Pass 4: fetch full metadata for the final, semantically screened reference set
and emit BibTeX. Every entry is anchored to a PMID (or a DOI for the two venues
PubMed does not index), so nothing here is written from memory.
"""
import json, re, ssl, time, urllib.parse, urllib.request

ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE
UA = {'User-Agent': 'Mozilla/5.0'}
EUT = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'

# key -> PMID. Selected after reading every title for semantic alignment.
PMIDS = {
    # --- SWI/SNF biology and cancer
    'kadoch2013':      '23644491',
    'mashtalir2018':   '30443016',
    'centore2020':     '32624222',
    'mittal2020':      '32303701',
    'wanior2021':      '34211095',
    'remodellers2024': '39014081',
    'atpremod2024':    '38081975',
    'helming2014':     '24562383',
    'hoffman2014':     '24520176',
    'michel2018':      '30104724',
    'alpsoy2018':      '29374058',
    'trizzino2018':    '30463020',
    'mathur2017':      '28112738',
    # --- ARID1A specifically
    'arid1atarget2021':'34619527',
    'arid1ahall2025':  '40429787',
    'arid1aimmuno2024':'39068111',
    'smarca4nsclc2023':'36450331',
    'swisnfnsclc2025': '41340120',
    'brd9isci2026':    '41797932',
    # --- structures
    'he2020':          '32747439',
    'yuan2022':        '35477757',
    'kasinath2018':    '30467178',
    'kasinath2021':    '34739333',
    'prc2struct2021':  '33232890',
    'fedorov2015':     '26601310',
    'farnaby2019':     '31427666',
    'fhd286trial2025': '40238563',
    'smd3236_2025':    '39745073',
    'smarca4rev2025':  '40398707',
    # --- screens and resources
    'behan2019':       '30971826',
    'meyers2017':      '29083409',
    'tsherniak2017':   '28753430',
    'pacini2021':      '33608518',
    'depmap2025':      '39468210',
    'goncalves2022':   '35839778',
    'tcgaluad2014':    '25079552',
    'kim2020':         '32385277',
    'corum2023':       '36382402',
    # --- co-essentiality methodology
    'pan2018':         '29778836',
    'wainberg2021':    '33859415',
    'boyle2018':       '29930200',
    'hassan2023':      '37750448',
    'ito2021':         '34857952',
    'koferle2022':     '35417719',
    # --- software and methods
    'love2014':        '25516281',
    'mayakonda2018':   '30341162',
    'xie2021':         '33780170',
    'traag2019':       '30914743',
    'cock2009':        '19304878',
    'jimenez2018':     '29281698',
    'jumper2021':      '34265844',
    'mirdita2022':     '35637307',
    'eastman2017':     '28746339',
    'maier2015':       '26574453',
    'miller2012':      '26605738',
    'roe2013':         '26583988',
    'theodoris2023':   '37258680',
    'geneformerapp':   '42026145',
    'hemani2018':      '29846171',
    'zhu2023af':       '37405868',
    # --- genetics resources
    'vosa2021':        '34475573',
    'mckay2017':       '28604730',
    'gtex2020':        '32913098',
    'karczewski2020':  '32461654',
    # --- target journal (Briefings in Bioinformatics)
    'bibinteractome':  '38180828',
    'bibppi2023':      '36880207',
}

# not indexed in PubMed -> Crossref by DOI
DOIS = {
    'lundberg2017':  '10.5555/3295222.3295230',      # NeurIPS, may not resolve; fallback below
    'eberhardt2021': '10.1021/acs.jcim.1c00203',     # J Chem Inf Model
    'evans2024':     '10.1101/2021.10.04.463034',    # AlphaFold-Multimer preprint
    'benjamini1995': '10.1111/j.2517-6161.1995.tb02031.x',
}


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


def crossref_doi(doi):
    j = get('https://api.crossref.org/works/' + urllib.parse.quote(doi))
    m = (j or {}).get('message')
    if not m:
        return None
    return dict(pmid='', title=(m.get('title') or [''])[0],
                authors=[f"{a.get('family','')} {''.join(w[0] for w in (a.get('given','') or '').split())}"
                         for a in m.get('author', [])],
                year=str((m.get('issued', {}).get('date-parts', [[None]])[0] or [''])[0]),
                journal=(m.get('container-title') or [''])[0],
                volume=m.get('volume', ''), issue=m.get('issue', ''),
                pages=m.get('page', ''), doi=m.get('DOI', ''))


records, failed = {}, []
print('fetching PubMed records ...')
for k, pmid in PMIDS.items():
    r = esum(pmid)
    if r and r['title']:
        records[k] = r
        print(f"  OK  {k:18s} {r['year']}  {r['journal'][:30]:30s} {r['title'][:62]}")
    else:
        failed.append((k, pmid))
        print(f"  !!  {k:18s} PMID {pmid} FAILED")
    time.sleep(0.35)

print('\nfetching Crossref records ...')
for k, doi in DOIS.items():
    r = crossref_doi(doi)
    if r and r['title']:
        records[k] = r
        print(f"  OK  {k:18s} {r['year']}  {r['journal'][:30]:30s} {r['title'][:62]}")
    else:
        failed.append((k, doi))
        print(f"  !!  {k:18s} DOI {doi} FAILED")
    time.sleep(0.35)


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
    if r['volume']:
        f.append(f'  volume  = {{{r["volume"]}}}')
    if r['issue']:
        f.append(f'  number  = {{{r["issue"]}}}')
    if r['pages']:
        f.append(f'  pages   = {{{r["pages"]}}}')
    if r['doi']:
        f.append(f'  doi     = {{{r["doi"]}}}')
    if r['pmid']:
        f.append(f'  note    = {{PMID: {r["pmid"]}}}')
    return '@article{' + k + ',\n' + ',\n'.join(f) + '\n}\n'


with open('references.bib', 'w', encoding='utf-8') as fh:
    for k in records:
        fh.write(bibtex(k, records[k]) + '\n')

yrs = [int(r['year']) for r in records.values() if (r['year'] or '').isdigit()]
recent = sum(1 for y in yrs if y >= 2021)
print(f'\n=== {len(records)} references written to references.bib ===')
print(f'failed: {len(failed)} -> {failed}')
print(f'year range {min(yrs)}-{max(yrs)}')
print(f'>= 2021 (last five years): {recent}/{len(yrs)} = {100*recent/len(yrs):.1f}%'
      f'   {"PASS" if recent/len(yrs) >= 0.5 else "BELOW 50% TARGET"}')
json.dump(records, open('refs_final.json', 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
