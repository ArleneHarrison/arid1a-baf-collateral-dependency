# -*- coding: utf-8 -*-
"""Pass 1: existence check for every candidate reference.

Per the reference-verification handbook:
  - PubMed esearch by title -> esummary for structured metadata (gold standard)
  - Crossref fallback for journals PubMed does not index (JASA, J Stat Softw, conferences)
  - anything found by neither is reported MISSING and must never be cited
"""
import json
import ssl
import sys
import time
import urllib.parse
import urllib.request

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE
UA = {'User-Agent': 'Mozilla/5.0 (reference verification)'}
EUT = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'


def get(url, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=30, context=ctx) as r:
                return json.loads(r.read().decode())
        except Exception:
            time.sleep(1.2 * (i + 1))
    return None


def pubmed_by_title(title):
    term = urllib.parse.quote(f'{title}[Title]')
    j = get(f'{EUT}/esearch.fcgi?db=pubmed&term={term}&retmode=json&retmax=3')
    ids = (j or {}).get('esearchresult', {}).get('idlist', [])
    if not ids:
        # relax: title words, not exact-title match
        term = urllib.parse.quote(title)
        j = get(f'{EUT}/esearch.fcgi?db=pubmed&term={term}&retmode=json&retmax=3')
        ids = (j or {}).get('esearchresult', {}).get('idlist', [])
    return ids


def esummary(pmid):
    j = get(f'{EUT}/esummary.fcgi?db=pubmed&id={pmid}&retmode=json&version=2.0')
    try:
        o = j['result'][pmid]
    except Exception:
        return None
    doi = ''
    for a in o.get('articleids', []):
        if a.get('idtype') == 'doi':
            doi = a.get('value', '')
    return dict(pmid=pmid, title=o.get('title', '').rstrip('.'),
                authors=[a['name'] for a in o.get('authors', [])][:3],
                year=(o.get('pubdate', '') or '')[:4],
                journal=o.get('source', ''), volume=o.get('volume', ''),
                issue=o.get('issue', ''), pages=o.get('pages', '') or o.get('elocationid', ''),
                doi=doi)


def crossref_by_title(title):
    u = ('https://api.crossref.org/works?query.bibliographic='
         + urllib.parse.quote(title) + '&rows=3&select=title,author,issued,container-title,DOI,volume,issue,page')
    j = get(u)
    items = (j or {}).get('message', {}).get('items', [])
    out = []
    for it in items:
        au = it.get('author', [])
        out.append(dict(pmid='', title=(it.get('title') or [''])[0],
                        authors=[f"{a.get('family','')} {a.get('given','')[:2]}".strip()
                                 for a in au[:3]],
                        year=str((it.get('issued', {}).get('date-parts', [[None]])[0] or [''])[0]),
                        journal=(it.get('container-title') or [''])[0],
                        volume=it.get('volume', ''), issue=it.get('issue', ''),
                        pages=it.get('page', ''), doi=it.get('DOI', '')))
    return out


def norm(s):
    return ''.join(c.lower() for c in (s or '') if c.isalnum())


CANDIDATES = [
    # tag, search title
    ('kadoch2013', 'Proteomic and bioinformatic analysis of mammalian SWI/SNF complexes identifies extensive roles in human malignancy'),
    ('mashtalir2018', 'Modular Organization and Assembly of SWI/SNF Family Chromatin Remodeling Complexes'),
    ('helming2014', 'ARID1B is a specific vulnerability in ARID1A-mutant cancers'),
    ('hoffman2014', 'Functional epigenetics approach identifies BRM/SMARCA2 as a critical synthetic lethal target in BRG1-deficient cancers'),
    ('michel2018', 'A non-canonical SWI/SNF complex is a synthetic lethal target in cancers driven by BAF complex perturbation'),
    ('alpsoy2018', 'Glioma tumor suppressor candidate region gene 1 and its paralog GLTSCR1-like form SWI/SNF chromatin remodeling subcomplexes'),
    ('he2020', 'Structure of nucleosome-bound human BAF complex'),
    ('yuan2022pbaf', 'Structural basis of human PBAF nucleosome recognition'),
    ('kasinath2018', 'Structures of human PRC2 with its cofactors AEBP2 and JARID2'),
    ('kasinath2021', 'JARID2 and AEBP2 regulate PRC2 in the presence of H2AK119ub1 and other histone modifications'),
    ('behan2019', 'Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens'),
    ('dempster2021', 'Chronos a cell population dynamics model of CRISPR experiments that improves inference of gene fitness'),
    ('meyers2017', 'Computational correction of copy number effect improves specificity of CRISPR-Cas9 essentiality screens in cancer cells'),
    ('tsherniak2017', 'Defining a Cancer Dependency Map'),
    ('pacini2021', 'Integrated cross-study datasets of genetic dependencies in cancer'),
    ('goncalves2022', 'Pan-cancer proteomic map of 949 human cell lines'),
    ('pan2018', 'Interrogation of Mammalian Protein Complex Structure Function and Membership Using Genome-Scale Fitness Screens'),
    ('wainberg2021', 'A genome-wide atlas of co-essential modules assigns function to uncharacterized genes'),
    ('boyle2018', 'High-resolution mapping of cancer cell networks using co-functional interactions'),
    ('love2014', 'Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2'),
    ('mayakonda2018', 'Maftools efficient and comprehensive analysis of somatic variants in cancer'),
    ('xie2021enrichr', 'Gene Set Knowledge Discovery with Enrichr'),
    ('traag2019', 'From Louvain to Leiden guaranteeing well-connected communities'),
    ('lundberg2017', 'A Unified Approach to Interpreting Model Predictions'),
    ('cock2009', 'Biopython freely available Python tools for computational molecular biology and bioinformatics'),
    ('jimenez2018', 'LightDock a new multi-scale approach to protein-protein docking'),
    ('roel2020lightdock', 'LightDock goes information-driven'),
    ('jumper2021', 'Highly accurate protein structure prediction with AlphaFold'),
    ('evans2021', 'Protein complex prediction with AlphaFold-Multimer'),
    ('mirdita2022', 'ColabFold making protein folding accessible to all'),
    ('eberhardt2021', 'AutoDock Vina 1.2.0 New Docking Methods Expanded Force Field and Python Bindings'),
    ('eastman2017', 'OpenMM 7 Rapid development of high performance algorithms for molecular dynamics'),
    ('maier2015', 'ff14SB Improving the Accuracy of Protein Side Chain and Backbone Parameters from ff99SB'),
    ('miller2012', 'MMPBSA.py An Efficient Program for End-State Free Energy Calculations'),
    ('roe2013', 'PTRAJ and CPPTRAJ Software for Processing and Analysis of Molecular Dynamics Trajectory Data'),
    ('theodoris2023', 'Transfer learning enables predictions in network biology'),
    ('hemani2018', 'The MR-Base platform supports systematic causal inference across the human phenome'),
    ('vosa2021', 'Large-scale cis- and trans-eQTL analyses identify thousands of genetic loci and polygenic scores that regulate blood gene expression'),
    ('mckay2017', 'Large-scale association analysis identifies new lung cancer susceptibility loci and heterogeneity in genetic susceptibility across histological subtypes'),
    ('gtex2020', 'The GTEx Consortium atlas of genetic regulatory effects across human tissues'),
    ('karczewski2020', 'The mutational constraint spectrum quantified from variation in 141456 humans'),
    ('kim2020', 'Single-cell RNA sequencing demonstrates the molecular and cellular reprogramming of metastatic lung adenocarcinoma'),
    ('tcgaluad2014', 'Comprehensive molecular profiling of lung adenocarcinoma'),
    ('tsitsiridis2023', 'CORUM the comprehensive resource of mammalian protein complexes 2022'),
    ('fedorov2015', 'Selective targeting of the BRG/PB1 bromodomains impairs embryonic and trophoblast stem cell maintenance'),
    ('farnaby2019', 'BAF complex vulnerabilities in cancer demonstrated via structure-based PROTAC design'),
    ('centore2020', 'Mammalian SWI/SNF Chromatin Remodeling Complexes Emerging Mechanisms and Therapeutic Strategies'),
    ('mittal2020', 'The SWI/SNF complex in cancer biology biomarkers and therapy'),
    ('wanior2021', 'Exploiting vulnerabilities of SWI/SNF chromatin remodelling complexes for cancer therapy'),
    ('reske2020', 'ARID1A-dependent maintenance of H3K4me1'),
    ('trizzino2018', 'The Tumor Suppressor ARID1A Controls Global Transcription via Pausing of RNA Polymerase II'),
]

print(f'Pass 1: existence check, {len(CANDIDATES)} candidates\n')
verified, needs_check = [], []
for tag, title in CANDIDATES:
    rec, src = None, ''
    for pmid in pubmed_by_title(title)[:3]:
        s = esummary(pmid)
        if s and (norm(title)[:55] in norm(s['title']) or norm(s['title'])[:55] in norm(title)):
            rec, src = s, 'pubmed'
            break
    if rec is None:
        for c in crossref_by_title(title):
            if norm(title)[:55] in norm(c['title']) or norm(c['title'])[:55] in norm(title):
                rec, src = c, 'crossref'
                break
    if rec:
        rec['tag'] = tag
        rec['source'] = src
        verified.append(rec)
        a = rec['authors'][0] if rec['authors'] else '?'
        print(f"  OK  [{tag:16s}] {a:22s} {rec['year']}  {rec['journal'][:38]:38s} ({src})")
    else:
        needs_check.append(dict(tag=tag, query=title))
        print(f"  !!  [{tag:16s}] NOT FOUND -> {title[:60]}")
    time.sleep(0.4)

print(f'\nverified {len(verified)} / {len(CANDIDATES)};  needs manual work: {len(needs_check)}')
yrs = [int(r['year']) for r in verified if (r['year'] or '').isdigit()]
recent = sum(1 for y in yrs if y >= 2021)
print(f'year range {min(yrs)}-{max(yrs)};  >=2021: {recent}/{len(yrs)} = {100*recent/len(yrs):.0f}%')
json.dump(dict(verified=verified, missing=needs_check),
          open('refs_pass1.json', 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
print('wrote refs_pass1.json')
