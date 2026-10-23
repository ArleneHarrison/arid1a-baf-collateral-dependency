"""Extract the GSE131907 sample-to-patient map from Supplementary Data 1.

The publisher workbook contains a malformed custom property that some
``openpyxl`` versions reject.  This reader intentionally uses only the OOXML
parts required for the first worksheet and records source provenance.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile

import pandas as pd


MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
NS = {"m": MAIN_NS}


def column_number(reference: str) -> int:
    letters = re.match(r"[A-Z]+", reference).group(0)
    number = 0
    for letter in letters:
        number = number * 26 + ord(letter) - ord("A") + 1
    return number - 1


def read_first_sheet(path: Path) -> list[list[str]]:
    with ZipFile(path) as archive:
        shared_xml = ET.fromstring(archive.read("xl/sharedStrings.xml"))
        shared = [
            "".join(node.text or "" for node in item.iter(f"{{{MAIN_NS}}}t"))
            for item in shared_xml.findall("m:si", NS)
        ]
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        rows: list[list[str]] = []
        for row in sheet.findall(".//m:sheetData/m:row", NS):
            values: dict[int, str] = {}
            for cell in row.findall("m:c", NS):
                value_node = cell.find("m:v", NS)
                value = "" if value_node is None else (value_node.text or "")
                if cell.attrib.get("t") == "s" and value:
                    value = shared[int(value)]
                values[column_number(cell.attrib["r"])] = value
            width = max(values, default=-1) + 1
            rows.append([values.get(i, "") for i in range(width)])
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    source = Path(args.input)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rows = read_first_sheet(source)
    header_index = next(i for i, row in enumerate(rows) if "Patient id" in row and "Samples" in row)
    headers = rows[header_index]
    records = []
    for row in rows[header_index + 1 :]:
        padded = row + [""] * (len(headers) - len(row))
        record = dict(zip(headers, padded))
        if record.get("Patient id") and record.get("Samples"):
            records.append(record)
    frame = pd.DataFrame(records)
    mapping = frame[["Samples", "Patient id", "Tissue origins"]].rename(
        columns={"Samples": "Sample", "Patient id": "patient_id", "Tissue origins": "tissue_origin"}
    )
    if mapping.Sample.duplicated().any():
        raise ValueError("Supplementary Data 1 contains duplicate sample names")
    if mapping.patient_id.nunique() != 44 or len(mapping) != 58:
        raise ValueError(
            f"Unexpected mapping dimensions: {len(mapping)} samples, "
            f"{mapping.patient_id.nunique()} patients"
        )
    mapping.to_csv(out / "gse131907_sample_patient_map.csv", index=False)
    sha256 = hashlib.sha256(source.read_bytes()).hexdigest()
    provenance = {
        "source": str(source.resolve()),
        "source_sha256": sha256,
        "source_description": "Supplementary Data 1 of Kim et al., Nature Communications 2020",
        "doi": "10.1038/s41467-020-16164-1",
        "n_samples": int(len(mapping)),
        "n_patients": int(mapping.patient_id.nunique()),
    }
    (out / "gse131907_sample_patient_map_provenance.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8"
    )
    print(mapping.to_string(index=False))
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()
