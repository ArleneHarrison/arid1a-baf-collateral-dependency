"""Download only pinned public 24Q4 inputs with complete-file integrity checks."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import requests

FILES = {
 'OmicsExpressionProteinCodingGenesTPMLogp1.csv': (51065489, 506628654, '71794802b750ce77c422dad0720a40af'),
 'OmicsCNGene.csv': (51065324, 1394554236, '4696565fa74dc805983bc278b4a41802'),
 'README.txt': (51065795, 43103, '54628c15a10b9ff6536db245cfed5231'),
}

def md5(path):
    h = hashlib.md5()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(8*1024*1024), b''):
            h.update(b)
    return h.hexdigest()

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--outdir', required=True)
    p.add_argument('--only', nargs='*', choices=list(FILES))
    args = p.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    manifest = []
    for name in args.only or FILES:
        fid, size, checksum = FILES[name]
        url = f'https://ndownloader.figshare.com/files/{fid}'
        dest = out / name
        if dest.exists():
            if dest.stat().st_size != size or md5(dest) != checksum:
                raise ValueError(f'Existing file failed integrity, retained untouched: {dest}')
            print('VERIFIED_EXISTING', name, flush=True)
        else:
            part = out / (name + '.part')
            if part.exists():
                part = out / (name + f'.{time.time_ns()}.part')
            print('DOWNLOADING', name, size, flush=True)
            with requests.get(url, stream=True, timeout=(20, 90)) as r:
                r.raise_for_status()
                h = hashlib.md5()
                written = 0
                tick = time.monotonic()
                with part.open('xb') as f:
                    for block in r.iter_content(2*1024*1024):
                        f.write(block)
                        h.update(block)
                        written += len(block)
                        if time.monotonic() - tick >= 25:
                            print('PROGRESS', name, written, '/', size, flush=True)
                            tick = time.monotonic()
            if written != size or h.hexdigest() != checksum:
                raise ValueError(f'Download failed integrity, partial file retained: {part}')
            part.rename(dest)
            print('VERIFIED_DOWNLOAD', name, flush=True)
        manifest.append(dict(file=name, url=url, bytes=size, md5=checksum,
                             release='DepMap Public 24Q4 v1',
                             doi='10.25452/figshare.plus.27993248.v1'))
    (out / 'r12_download_manifest.json').write_text(json.dumps(manifest, indent=2))

if __name__ == '__main__':
    main()
