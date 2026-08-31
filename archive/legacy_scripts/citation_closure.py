# -*- coding: utf-8 -*-
"""Handbook Pass 1 closure check: every cited key defined, every defined key cited."""
import io, re, os

D = (r'D:\桌面\_按项目整理\02_生信数据与医学AI项目\07_医学AI_深度学习_多项目孵化'
     r'\医学AI深度学习_导师课题项目群\范教授课题_NSCLC-SEER_分子对接多项目\纯生信'
     r'\ARID1A_BAF_最终项目\submission')

tex = io.open(os.path.join(D, 'main.tex'), encoding='utf-8').read()
bib = io.open(os.path.join(D, 'references.bib'), encoding='utf-8').read()

cited = set()
for m in re.finditer(r'\\cite[a-z]*\{([^}]+)\}', tex):
    cited |= {k.strip() for k in m.group(1).split(',')}
defined = set(re.findall(r'@article\{([^,]+),', bib))

print(f'cited in text : {len(cited)}')
print(f'defined in bib: {len(defined)}')
print(f'missing  (cited, not in bib): {sorted(cited - defined) or "none"}')
print(f'orphan   (in bib, never cited): {sorted(defined - cited) or "none"}')

yrs = [int(y) for y in re.findall(r'year\s*=\s*\{(\d{4})\}', bib)]
rec = sum(1 for y in yrs if y >= 2021)
print(f'\nrecency: {rec}/{len(yrs)} = {100*rec/len(yrs):.1f}% published 2021 or later')

pm = len(re.findall(r'PMID:', bib))
doi = len(re.findall(r'doi\s*=', bib))
print(f'entries with DOI : {doi}/{len(yrs)}')
print(f'entries with PMID: {pm}/{len(yrs)}  (rest are Crossref-only venues)')
