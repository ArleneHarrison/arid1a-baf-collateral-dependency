# -*- coding: utf-8 -*-
"""Export the 15 main-text figures as separate numbered files for upload.

Briefings in Bioinformatics asks for figures supplied separately at 600 dpi for
line drawings and combination artwork. The vector PDFs already in
submission/figures are the primary deliverable; this script additionally
rasterises each at 600 dpi and renames both to Figure1..Figure15 in the order
they appear in the compiled manuscript, which is the order the submission system
expects.
"""
import os
import re
import shutil
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.abspath(os.path.join(HERE, '..', '..'))
SUB = os.path.join(PROJ, 'submission')
FIG = os.path.join(SUB, 'figures')
OUT = os.path.join(SUB, 'submission_package', 'figures')
os.makedirs(OUT, exist_ok=True)

# figure number -> source file, read from the compiled aux so the order can
# never drift from the manuscript
aux = open(os.path.join(SUB, 'main.aux'), encoding='utf-8').read()
numbers = {m.group(1): int(m.group(2))
           for m in re.finditer(r'\\newlabel\{(fig:[a-z0-9]+)\}\{\{(\d+)\}', aux)}

tex = open(os.path.join(SUB, 'main.tex'), encoding='utf-8').read()
blocks = re.findall(r'\\includegraphics\[[^\]]*\]\{([^}]+)\}.*?\\label\{(fig:[a-z0-9]+)\}',
                    tex, flags=re.S)

pairs = sorted(((numbers[lab], src) for src, lab in blocks if lab in numbers))
assert len(pairs) == len(blocks), 'a figure label is missing from main.aux'

print(f'exporting {len(pairs)} figures at 600 dpi\n')
for n, src in pairs:
    stem = src.replace('.pdf', '')
    pdf_in = os.path.join(FIG, stem + '.pdf')
    shutil.copy(pdf_in, os.path.join(OUT, f'Figure{n}.pdf'))
    subprocess.run(['pdftoppm', '-r', '600', '-png', '-singlefile',
                    pdf_in, os.path.join(OUT, f'Figure{n}')], check=True)
    kb = os.path.getsize(os.path.join(OUT, f'Figure{n}.png')) // 1024
    print(f'  Figure{n:<2d} <- {stem:22s} ({kb:,} KB at 600 dpi)')

print('\nwritten to', OUT)
