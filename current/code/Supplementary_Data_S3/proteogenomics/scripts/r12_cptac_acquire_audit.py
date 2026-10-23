#!/usr/bin/env python3
"""Acquire published CPTAC lung processed matrices and audit sample linkage.

No association tests, imputation, rescaling, or outcome-driven filtering.
Run on the user server. Sources are the original-cohort LinkedOmics portals,
not a mixture of reprocessed pan-cancer and original-release matrices.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import socket
from datetime import datetime, timezone
from html import unescape
from pathlib import Path
from urllib.parse import urljoin, urlparse

import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BAF = ["ARID1A", "ARID1B", "SMARCA4", "SMARCA2", "SMARCB1", "SMARCC1",
       "SMARCC2", "SMARCD1", "SMARCD2", "SMARCD3", "SMARCE1", "ACTL6A",
       "ACTL6B", "ACTB", "DPF1", "DPF2", "DPF3", "ARID2", "PBRM1", "BRD7",
       "PHF10", "BRD9", "BICRA", "BICRAL", "SS18", "SS18L1"]
PORTALS = {"LUAD": "https://linkedomics.org/data_download/CPTAC-LUAD/",
           "LUSC": "https://linkedomics.org/data_download/CPTAC-LSCC/"}
DOIS = {"LUAD": "10.1016/j.cell.2020.06.013", "LUSC": "10.1016/j.cell.2021.07.016"}


def now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False), encoding="utf-8")


def session():
    s = requests.Session()
    s.headers["User-Agent"] = "ARID1A-BAF-CPTAC-reanalysis/1.0 research-processed-data"
    retry = Retry(total=2, connect=2, read=2, backoff_factor=2,
                  status_forcelist=[429, 500, 502, 503, 504])
    s.mount("https://", HTTPAdapter(max_retries=retry))
    return s


def download(s, url, path, cap_mb=150):
    sidecar = path.with_suffix(path.suffix + ".source.json")
    if path.exists() and sidecar.exists():
        previous = json.loads(sidecar.read_text(encoding="utf-8"))
        if previous.get("url") == url and sha256(path) == previous.get("sha256"):
            return previous
        raise ValueError(f"Existing file fails provenance check: {path}")
    print(f"DOWNLOAD {url}", flush=True)
    with s.get(url, stream=True, timeout=(20, 90)) as r:
        r.raise_for_status()
        cap = cap_mb * 1024 * 1024
        if int(r.headers.get("Content-Length", 0)) > cap:
            raise ValueError(f"Processed file exceeds {cap_mb} MB safety cap")
        tmp = path.with_suffix(path.suffix + ".part")
        size = 0
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1024 * 512):
                size += len(chunk)
                if size > cap:
                    raise ValueError("Downloaded response exceeded file size cap")
                fh.write(chunk)
        tmp.replace(path)
        info = {"url": url, "resolved_url": r.url, "downloaded_utc": now(),
                "last_modified": r.headers.get("Last-Modified"),
                "etag": r.headers.get("ETag"), "bytes": size,
                "sha256": sha256(path), "file": str(path)}
        write_json(sidecar, info)
    return info


def catalog(html, base):
    result = []
    for row in re.findall(r"<tr\b[^>]*>(.*?)</tr>", html, re.I | re.S):
        cells = re.findall(r"<td\b[^>]*>(.*?)</td>", row, re.I | re.S)
        if len(cells) < 7:
            continue
        clean = [re.sub(r"\s+", " ", unescape(re.sub(r"<[^>]*>", " ", c))).strip()
                 for c in cells]
        links = re.findall(r"href\s*=\s*['\"]([^'\"]+)", row, re.I)
        if not links:
            continue
        result.append({"label": clean[0], "url": urljoin(base, links[0]),
                       "description": clean[2], "level": clean[3], "unit": clean[4],
                       "expected_samples": clean[5], "expected_features": clean[6]})
    return result


def selected(row):
    label = row["label"]
    if label.startswith("RNAseq ("):
        return "rna"
    if label.startswith("Proteome ("):
        return "protein"
    if label.startswith("SCNV (") and "log ratio" in row["unit"].lower():
        return "cnv"
    if label.startswith("Clinical") or label.startswith("Molecular phenotypes"):
        return "metadata"
    return None


def tissue_type(row):
    label = row["label"].lower()
    if "normal" in label:
        return "normal"
    if "tumor" in label or label == "clinical":
        return "tumor"
    raise ValueError(f"Cannot determine sample type from portal label: {row}")


def patient_id(raw, tissue):
    """Only documented CPTAC patient syntax; no arbitrary split or fuzzy joins."""
    value = str(raw).strip()
    # LinkedOmics uses R's make.names notation (dots for hyphens; X before digits).
    match = re.fullmatch(r"(C3[NL])[.-](\d{5})(?:[._-]([NT]))?", value)
    if match:
        stem = f"{match[1]}-{match[2]}"
        suffix = match[3]
    else:
        match = re.fullmatch(r"X?(\d{2}LU\d{3})(?:[._-]([NT]))?", value)
        if not match:
            raise ValueError(f"Unrecognized patient ID (manual audit needed): {raw}")
        stem, suffix = match[1], match[2]
    if suffix is not None and suffix != {"tumor": "T", "normal": "N"}[tissue]:
        raise ValueError(f"Sample suffix conflicts with source tissue type: {raw}")
    return stem


def read_matrix(path, row, cohort, modality, tissue, out):
    # Inspect exact header before pandas can automatically mangle duplicate names.
    with open(path, encoding="utf-8-sig", newline="") as fh:
        header = next(csv.reader(fh, delimiter="\t"))
    if len(set(header[1:])) != len(header[1:]):
        raise ValueError(f"Duplicate sample columns in {path}")
    raw = pd.read_csv(path, sep="\t", index_col=0, na_values=["NA", "NaN", "", "null"])
    if raw.shape[1] != int(row["expected_samples"]):
        raise ValueError(f"Sample count disagrees with portal: {path}: {raw.shape}")
    raw.index = raw.index.astype(str).str.strip()
    matrix = raw.apply(pd.to_numeric, errors="coerce")
    # A text feature annotation column must never be interpreted as an assay sample.
    newly_missing = matrix.isna() & raw.notna()
    if newly_missing.to_numpy().any():
        raise ValueError(f"Unexpected nonnumeric measurement in {path}")
    duplicate_rows = matrix.loc[matrix.index.duplicated(keep=False)].copy()
    if not duplicate_rows.empty:
        duplicate_rows.to_csv(out / f"{cohort}_{modality}_{tissue}_duplicate_gene_rows.tsv", sep="\t")
        # Keep deterministic summary without choosing the row that gives a desired effect.
        matrix = matrix.groupby(level=0, sort=False).median()
    ids = [patient_id(x, tissue) for x in matrix.columns]
    if len(set(ids)) != len(ids):
        raise ValueError(f"Canonical IDs collide in {path}")
    mapping = pd.DataFrame({"cohort": cohort, "modality": modality,
                            "sample_type": tissue, "raw_sample_id": matrix.columns,
                            "patient_id": ids,
                            "sample_id": [f"{cohort}:{p}:{tissue}" for p in ids],
                            "source_url": row["url"]})
    matrix.columns = mapping["sample_id"]
    matrix.index.name = "gene_symbol"
    output = out / f"{cohort}_{modality}_{tissue}_gene_symbol.tsv.gz"
    matrix.to_csv(output, sep="\t", compression="gzip")
    coverage = []
    for gene in BAF:
        measured = matrix.loc[gene].notna() if gene in matrix.index else pd.Series(False, index=matrix.columns)
        coverage.append({"cohort": cohort, "modality": modality, "sample_type": tissue,
                         "gene_symbol": gene, "present": gene in matrix.index,
                         "n_nonmissing": int(measured.sum()), "n_total": len(measured),
                         "fraction_nonmissing": float(measured.mean())})
    audit = {"cohort": cohort, "modality": modality, "sample_type": tissue,
             "n_samples": matrix.shape[1], "n_gene_symbols": matrix.shape[0],
             "n_original_rows": raw.shape[0], "duplicate_rows": duplicate_rows.shape[0],
             "units": row["unit"], "source_url": row["url"],
             "output": str(output), "sha256": sha256(output),
             "transformation": "none; retain the portal's already-normalized numerical values",
             "duplicate_rule": "median across rows with the same exact symbol; exported audit",
             "missing_value_policy": "preserved; no imputation"}
    return matrix, mapping, coverage, audit


def read_metadata(path, row, cohort, tissue, out):
    raw = pd.read_csv(path, sep="\t", dtype=str, na_values=["NA", "NaN", ""])
    raw = raw.rename(columns={raw.columns[0]: "raw_sample_id"})
    # TSI stores a data-type declaration in its first data row.
    is_type = raw.iloc[:, 1:].apply(lambda col: col.isin(["BIN", "CAT", "CON", "ORD"])).all(axis=1)
    raw = raw.loc[~is_type].copy()
    raw["patient_id"] = raw["raw_sample_id"].map(lambda x: patient_id(x, tissue))
    raw["cohort"], raw["sample_type"] = cohort, tissue
    raw["sample_id"] = raw["patient_id"].map(lambda p: f"{cohort}:{p}:{tissue}")
    if raw["sample_id"].duplicated().any():
        raise ValueError(f"Duplicate patient metadata records in {path}")
    raw["source_url"] = row["url"]
    label = "molecular" if row["label"].startswith("Molecular") else "clinical"
    raw.to_csv(out / f"{cohort}_{label}_{tissue}_metadata.tsv", sep="\t", index=False)
    return raw


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("/root/autodl-tmp/arid1a_baf_r12_cptac"))
    ap.add_argument("--cohort", nargs="+", choices=list(PORTALS), default=list(PORTALS))
    args = ap.parse_args()
    root = args.root.resolve()
    source_dir, out = root / "raw_processed", root / "audit"
    source_dir.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    s = session()
    manifest = {"created_utc": now(), "hostname": socket.gethostname(), "sources": [],
                "cohort_releases": DOIS, "scope": "processed data acquisition and matching only"}
    coverages, mappings, assays, flow, catalog_rows, metadata_tables = [], [], [], [], [], []
    matrices = {}
    for cohort in args.cohort:
        portal = PORTALS[cohort]
        html_path = source_dir / f"{cohort}_portal.html"
        manifest["sources"].append(download(s, portal, html_path))
        rows = catalog(html_path.read_text(encoding="utf-8"), portal)
        if not rows:
            raise ValueError("No datasets found in portal HTML")
        for row in rows:
            catalog_rows.append({"cohort": cohort, **row})
            modality = selected(row)
            if modality is None:
                continue
            tissue = tissue_type(row)
            path = source_dir / Path(urlparse(row["url"]).path).name
            provenance = download(s, row["url"], path)
            manifest["sources"].append({"cohort": cohort, "modality": modality,
                                         "sample_type": tissue, **provenance})
            write_json(out / "r12_cptac_source_manifest.json", manifest)
            if modality == "metadata":
                metadata_tables.append(read_metadata(path, row, cohort, tissue, out))
            else:
                mat, mapping, coverage, audit = read_matrix(path, row, cohort, modality, tissue, out)
                matrices[(cohort, modality, tissue)] = mat
                mappings.append(mapping)
                coverages.extend(coverage)
                assays.append(audit)
        for tissue in ["tumor", "normal"]:
            sets = {m: set(matrices[(cohort, m, tissue)].columns)
                    for m in ["protein", "rna", "cnv"] if (cohort, m, tissue) in matrices}
            for m, ids in sets.items():
                flow.append({"cohort": cohort, "sample_type": tissue, "step": m,
                             "n_samples": len(ids)})
            matched = sets.get("protein", set()) & sets.get("rna", set())
            all_ids = set.union(*sets.values()) if sets else set()
            sample_rows = []
            for sid in sorted(all_ids):
                sample_rows.append({"sample_id": sid, "cohort": cohort,
                                    "patient_id": sid.split(":")[1], "sample_type": tissue,
                                    **{f"has_{m}": sid in sets.get(m, set()) for m in ["protein", "rna", "cnv"]},
                                    "rna_protein_matched": sid in matched,
                                    "rna_protein_cnv_matched": sid in matched & sets.get("cnv", set())})
            sample_df = pd.DataFrame(sample_rows)
            sample_df.to_csv(out / f"{cohort}_{tissue}_sample_match.tsv", sep="\t", index=False)
            for step, ids in [("protein_RNA", matched), ("protein_RNA_CNV", matched & sets.get("cnv", set()))]:
                flow.append({"cohort": cohort, "sample_type": tissue, "step": step, "n_samples": len(ids)})
            # Select identical sample columns only; genes are left in their own assay universe.
            for m in sets:
                sel = sorted(matched & sets.get("cnv", set())) if tissue == "tumor" else sorted(matched)
                matrices[(cohort, m, tissue)].reindex(columns=sel).to_csv(
                    out / f"{cohort}_{m}_{tissue}_matched.tsv.gz", sep="\t", compression="gzip")
    pd.concat(mappings, ignore_index=True).to_csv(out / "r12_cptac_sample_mapping.tsv", sep="\t", index=False)
    pd.DataFrame(coverages).to_csv(out / "r12_cptac_baf_coverage.tsv", sep="\t", index=False)
    pd.DataFrame(flow).to_csv(out / "r12_cptac_sample_flow.tsv", sep="\t", index=False)
    pd.DataFrame(catalog_rows).to_csv(out / "r12_cptac_source_catalog.tsv", sep="\t", index=False)
    pd.DataFrame(assays).to_csv(out / "r12_cptac_assay_audit.tsv", sep="\t", index=False)
    fields = []
    for frame in metadata_tables:
        for col in frame.columns:
            fields.append({"cohort": frame["cohort"].iloc[0], "sample_type": frame["sample_type"].iloc[0],
                           "field": col, "n_records": len(frame), "n_nonmissing": int(frame[col].notna().sum()),
                           "candidate_purity": bool(re.search("purity|cellularity", col, re.I)),
                           "candidate_batch": bool(re.search("batch|plex|channel|instrument", col, re.I))})
    pd.DataFrame(fields).to_csv(out / "r12_cptac_metadata_coverage.tsv", sep="\t", index=False)
    manifest["finished_utc"] = now()
    manifest["results"] = assays
    write_json(out / "r12_cptac_source_manifest.json", manifest)
    print(pd.DataFrame(flow).to_string(index=False), flush=True)
    print("FINISHED acquisition/matching; no association analysis performed", flush=True)


if __name__ == "__main__":
    main()
