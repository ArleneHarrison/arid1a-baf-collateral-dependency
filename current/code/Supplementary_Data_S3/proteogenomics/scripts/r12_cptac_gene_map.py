#!/usr/bin/env python3
"""Fetch the exact GENCODE 34 mapping used by the official cptac package."""
from pathlib import Path
import json
import pandas as pd
from r12_cptac_acquire_audit import download, session, write_json

root = Path("/root/autodl-tmp/arid1a_baf_r12_cptac/pancan_release")
raw = root / "raw_processed"
raw.mkdir(parents=True, exist_ok=True)
s = session()
record_file = raw / "cptac_package_zenodo_8394329.json"
record_info = download(s, "https://zenodo.org/api/records/8394329", record_file)
data = json.loads(record_file.read_text())
selected = [f for f in data["files"] if "bcm-luad-mapping-gencode.v34.basic.annotation-mapping.txt.gz" == f["key"]]
if len(selected) != 1:
    raise ValueError("Cannot uniquely find the cptac LUAD GENCODE 34 mapping")
sources = [record_info]
for f in selected:
    url = f["links"].get("self") or f["links"].get("content")
    target = raw / f["key"]
    info = download(s, url, target)
    sources.append({**info, "record_checksum": f["checksum"]})
    import hashlib
    if f["checksum"].startswith("md5:"):
        actual = hashlib.md5(target.read_bytes()).hexdigest()
        assert actual == f["checksum"].split(":", 1)[1], "Zenodo checksum mismatch"
    d = pd.read_csv(target, sep="\t", nrows=3)
    print(d.to_string(index=False), flush=True)
write_json(root / "gene_map_source_manifest.json", {"sources": sources,
    "package_download_source": "https://github.com/PayneLab/cptac/blob/master/cptac/tools/download_tools.py",
    "mapping_specification": "https://github.com/PayneLab/cptac/blob/master/cptac/cancers/bcm/bcmluad.py"})
