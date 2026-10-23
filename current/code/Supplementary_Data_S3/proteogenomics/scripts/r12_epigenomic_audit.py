#!/usr/bin/env python3
"""Bounded R12 public-data audit. Run on the user's Linux server only.

This phase downloads metadata and the small deposited processed files, never
SRA/FASTQ/BAM. It creates no biological result and performs no differential test.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import platform
import re
import tarfile
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

MAX_FILE_BYTES = 40 * 1024**2
UA = "ARID1A-BAF-public-reanalysis/12.0"


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_soft(data: bytes) -> list[dict]:
    records = []
    record = None
    for line in data.decode("utf-8", errors="replace").splitlines():
        if line.startswith("^") and " = " in line:
            typ, accession = line[1:].split(" = ", 1)
            record = {"record_type": typ, "accession": accession, "fields": {}}
            records.append(record)
        elif record is not None and line.startswith("!") and " = " in line:
            key, value = line[1:].split(" = ", 1)
            record["fields"].setdefault(key, []).append(value)
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default="/root/autodl-tmp/arid1a_baf_r12_epigenomic")
    args = parser.parse_args()
    if platform.system() != "Linux":
        raise SystemExit("Run on the authorized Linux server, not the workstation.")
    root = Path(args.root).resolve()
    if root != Path("/root/autodl-tmp/arid1a_baf_r12_epigenomic"):
        raise SystemExit("Unexpected output root; inspect the script before changing scope.")
    inputs, outputs = root / "inputs", root / "outputs"
    inputs.mkdir(parents=True, exist_ok=True)
    outputs.mkdir(parents=True, exist_ok=True)
    manifest = {"started_utc": datetime.now(timezone.utc).isoformat(), "host": platform.node(),
                "scope": "metadata and bounded processed-file audit only", "downloads": [], "errors": []}

    def fetch(url: str, name: str) -> bytes:
        path = inputs / name
        if path.exists():
            payload = path.read_bytes()
            cached = True
        else:
            payload = b""
            cached = False
            for attempt in range(3):
                try:
                    req = urllib.request.Request(url, headers={"User-Agent": UA})
                    with urllib.request.urlopen(req, timeout=35) as response:
                        if int(response.headers.get("Content-Length", 0)) > MAX_FILE_BYTES:
                            raise ValueError("File exceeds bounded audit size; no raw download started")
                        payload = response.read(MAX_FILE_BYTES + 1)
                    if len(payload) > MAX_FILE_BYTES:
                        raise ValueError("File exceeds bounded audit size")
                    if b"Checking your browser" in payload[:30000] or b"captcha" in payload[:3000].lower():
                        raise ValueError("Returned challenge HTML, not a data file")
                    path.write_bytes(payload)
                    break
                except Exception:
                    if attempt == 2:
                        raise
                    time.sleep(2 + attempt)
        manifest["downloads"].append({"url": url, "file": str(path), "bytes": len(payload),
                                      "sha256": sha256(payload), "cached": cached})
        print(f"AUDITED {name}: {len(payload)} bytes", flush=True)
        return payload

    all_records = {}
    for accession in ["GSE101966", "GSE164810"]:
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}&targ=all&form=text&view=full"
        try:
            payload = fetch(url, accession + "_family.soft")
            records = read_soft(payload)
            if not any(x["record_type"] == "SAMPLE" for x in records):
                raise ValueError("SOFT did not contain sample records")
            all_records[accession] = records
        except Exception as exc:
            manifest["errors"].append({"step": accession, "error": repr(exc)})

    inventory = {}
    processed = {
        "GSE101966_HCT116_FPKM.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE101nnn/GSE101966/suppl/GSE101966_HCT116_FPKM.txt.gz",
        "GSE101966_RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE101nnn/GSE101966/suppl/GSE101966_RAW.tar",
    }
    for name, url in processed.items():
        try:
            payload = fetch(url, name)
            if name.endswith(".gz"):
                lines = gzip.decompress(payload).decode("utf-8", errors="replace").splitlines()
                inventory[name] = {"line_count": len(lines), "header": lines[0], "preview": lines[1:4]}
            elif name.endswith(".tar"):
                members = []
                with tarfile.open(fileobj=io.BytesIO(payload), mode="r:") as archive:
                    for member in archive.getmembers():
                        if not member.isfile():
                            continue
                        stream = archive.extractfile(member)
                        assert stream is not None
                        content = stream.read()
                        if member.name.endswith(".gz"):
                            content = gzip.decompress(content)
                        lines = content.decode("utf-8", errors="replace").splitlines()
                        members.append({"name": member.name, "bytes": member.size,
                                        "header": lines[0] if lines else "", "preview": lines[1:3]})
                inventory[name] = members
        except Exception as exc:
            manifest["errors"].append({"step": name, "error": repr(exc)})

    try:
        url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=" + urllib.parse.quote("EXT_ID:28967863 AND SRC:MED") + "&format=json"
        search = json.loads(fetch(url, "kelso_europepmc_search.json"))
        hits = search.get("resultList", {}).get("result", [])
        pmcid = next((x.get("pmcid") for x in hits if x.get("pmcid")), None)
        if pmcid:
            xml = fetch(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML", "kelso_fulltext.xml")
            document = ET.fromstring(xml)
            paragraphs = [" ".join(x.itertext()) for x in document.iter("p")]
            keys = re.compile(r"replicat|clone|Horizon|shRNA|FPKM|count|download|RNA.seq", re.I)
            inventory["paper_methods_excerpt"] = [p for p in paragraphs if keys.search(p)]
            xlink = "{http://www.w3.org/1999/xlink}href"
            inventory["paper_supplements"] = [{"text": " ".join(el.itertext()), "href": el.attrib.get(xlink)}
                                              for el in document.iter() if el.attrib.get(xlink) and
                                              any(z in el.attrib.get(xlink, "").lower() for z in ["supp", "data", ".xls", ".csv", ".txt", ".zip"])]
    except Exception as exc:
        manifest["errors"].append({"step": "kelso_paper", "error": repr(exc)})

    rna_manifest = []
    for record in all_records.get("GSE101966", []):
        fields = record["fields"]
        if record["record_type"] != "SAMPLE" or fields.get("Sample_library_strategy") != ["RNA-Seq"]:
            continue
        title = " ".join(fields.get("Sample_title", []))
        arid1a_ko = int("ARID1A-/-" in title or "ARID1AKO" in title)
        arid1b_kd = int("ARID1BKD" in title)
        rna_manifest.append({"gsm": record["accession"], "title": title,
                             "arid1a_ko": arid1a_ko, "arid1b_kd": arid1b_kd,
                             "interaction": arid1a_ko * arid1b_kd,
                             "replicate_label": re.search(r"Rep(\d+)", title).group(1) if re.search(r"Rep(\d+)", title) else None,
                             "characteristics": fields.get("Sample_characteristics_ch1", []),
                             "source_name": fields.get("Sample_source_name_ch1", []),
                             "treatment": fields.get("Sample_treatment_protocol_ch1", []),
                             "growth": fields.get("Sample_growth_protocol_ch1", []),
                             "relations": fields.get("Sample_relation", [])})
    audit = {"metadata": all_records, "processed_inventory": inventory, "rna_sample_manifest": rna_manifest,
             "prespecified_question": "Does ARID1A knockout modify the effect of ARID1B knockdown?",
             "interaction_contrast": "(double perturbation - ARID1A KO) - (ARID1B KD - WT)",
             "formal_gate": "Raw comparable integer counts required for count-based DESeq2/edgeR; FPKM only permits explicitly exploratory log-linear analysis.",
             "biological_unit_gate": "Two labeled biological replicates per group do not establish two independent KO clones.",
             "integration_gate": "No individual-level paired integration unless the specimen correspondence is documented."}
    (outputs / "r12_epigenomic_data_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    manifest["completed_utc"] = datetime.now(timezone.utc).isoformat()
    (outputs / "r12_epigenomic_download_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"rna_samples": len(rna_manifest), "errors": manifest["errors"], "outputs": str(outputs)}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
