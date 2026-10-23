#!/usr/bin/env python3
"""Acquire the separately versioned pan-cancer CPTAC release, then inspect schemas."""
import argparse
import json
from pathlib import Path
from urllib.parse import urlparse
import pandas as pd
from r12_cptac_acquire_audit import download, session, catalog, write_json, now


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=Path("/root/autodl-tmp/arid1a_baf_r12_cptac/pancan_release"))
    p.add_argument("--cohort", nargs="+", default=["LUAD", "LSCC"])
    args = p.parse_args()
    raw = args.root / "raw_processed"
    raw.mkdir(parents=True, exist_ok=True)
    s = session()
    manifest = {"started_utc": now(), "release": "CPTAC pan-cancer; linked version, distinct from original cohort portals",
                "sources": []}
    datasets = []
    labels = {"RNAseq (Gene level, Tumor)": ("rna", "tumor"),
              "RNAseq (Gene level, Normal)": ("rna", "normal"),
              "Proteomics (Tumor)": ("protein", "tumor"),
              "Proteomics (Normal)": ("protein", "normal"),
              "CNV (log2 ratio)": ("cnv", "tumor"),
              "Meta table": ("metadata", "tumor"),
              "Phenotype": ("phenotype", "tumor")}
    for cohort in args.cohort:
        portal = f"https://linkedomics.org/data_download/CPTAC-pancan-{cohort}/"
        page = raw / f"{cohort}_pancan_portal.html"
        manifest["sources"].append(download(s, portal, page))
        rows = catalog(page.read_text(encoding="utf-8"), portal)
        for row in rows:
            if row["label"] not in labels:
                continue
            modality, tissue = labels[row["label"]]
            fn = raw / Path(urlparse(row["url"]).path).name
            info = download(s, row["url"], fn)
            entry = {"cohort": "LUSC" if cohort == "LSCC" else cohort,
                     "modality": modality, "sample_type": tissue, **row, **info}
            manifest["sources"].append(entry)
            datasets.append(entry)
            write_json(args.root / "source_manifest.json", manifest)
            d = pd.read_csv(fn, sep="\t", nrows=2)
            print(json.dumps({"file": fn.name, "columns": d.columns.tolist()[:12],
                              "first_records": d.iloc[:, :5].to_dict(orient="records")}), flush=True)
    pd.DataFrame(datasets).to_csv(args.root / "source_catalog.tsv", sep="\t", index=False)
    manifest["finished_utc"] = now()
    write_json(args.root / "source_manifest.json", manifest)


if __name__ == "__main__":
    main()
