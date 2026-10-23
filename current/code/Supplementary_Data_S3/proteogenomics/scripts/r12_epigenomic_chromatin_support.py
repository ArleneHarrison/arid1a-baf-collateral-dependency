#!/usr/bin/env python3
"""Server-side descriptive promoter peak support and metadata consistency audit.

Per-sample peak calls are thresholded and use different intervals. They are NOT
a shared-region count matrix, and this program never assigns differential
accessibility p-values or equates peak non-detection with chromatin closure.
"""
from __future__ import annotations
import gzip
import hashlib
import io
import json
import platform
import re
import tarfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/root/autodl-tmp/arid1a_baf_r12_epigenomic")


def overlap_any(peaks, promoters):
    answer = np.zeros(len(promoters), dtype=bool)
    for chromosome, idx in promoters.groupby("chromosome", sort=False).groups.items():
        idx = np.asarray(list(idx))
        intervals = peaks[peaks.chromosome == chromosome].sort_values("start")
        if intervals.empty:
            continue
        starts = intervals.start.to_numpy(dtype=np.int64)
        maximum_ends = np.maximum.accumulate(intervals.end.to_numpy(dtype=np.int64))
        q = promoters.loc[idx]
        j = np.searchsorted(starts, q.promoter_end.to_numpy(), side="left") - 1
        valid = j >= 0
        answer[idx[valid]] = maximum_ends[j[valid]] > q.promoter_start.to_numpy()[valid]
    return answer


def main():
    if platform.system() != "Linux":
        raise SystemExit("Run only on the authorized Linux server")
    outputs = ROOT / "outputs"
    audit = json.loads((outputs / "r12_epigenomic_data_audit.json").read_text())
    samples = {r["accession"]: r["fields"] for r in audit["metadata"]["GSE101966"] if r["record_type"] == "SAMPLE"}
    rna = pd.read_csv(outputs / "r12_rna_primary_interaction.csv")
    rna = rna.sort_values("gene").reset_index(drop=True)
    rna["annotated_TSS"] = np.where(rna.strand == "+", rna.start, rna.end - 1)
    rna["promoter_start"] = np.maximum(0, rna.annotated_TSS - np.where(rna.strand == "+", 2000, 500))
    rna["promoter_end"] = rna.annotated_TSS + np.where(rna.strand == "+", 500, 2000)
    support = rna[["gene", "transcript_id", "chromosome", "strand", "annotated_TSS", "promoter_start", "promoter_end"]].copy()
    inventory = []
    groups = defaultdict(list)
    with tarfile.open(ROOT / "inputs" / "GSE101966_RAW.tar", "r:") as archive:
        for member in archive.getmembers():
            if not member.isfile():
                continue
            match = re.match(r"(GSM\d+)_", Path(member.name).name)
            if not match:
                continue
            gsm = match.group(1)
            fields = samples[gsm]
            title = fields["Sample_title"][0]
            assay = next((mark for mark in ["ATAC", "H3K27ac", "H3K27me3", "H3K4me3", "H3K4me"] if mark in title), None)
            if assay is None:
                continue
            a = int("ARID1A-/-" in title or "ARID1AKO" in title)
            b = int("ARID1BKD" in title)
            group = {(0,0): "WT", (0,1): "B_KD", (1,0): "A_KO", (1,1): "double"}[(a,b)]
            stream = archive.extractfile(member)
            compressed = stream.read()
            data = gzip.decompress(compressed).decode("utf-8", errors="strict")
            table = pd.read_csv(io.StringIO(data), sep="\t", comment="#", header=None,
                                usecols=[0,1,2,3,5], names=["peak_id", "chromosome", "start", "end", "normalized_tag_count"])
            table["start"] = pd.to_numeric(table.start, errors="coerce")
            table["end"] = pd.to_numeric(table.end, errors="coerce")
            table = table.dropna(subset=["chromosome", "start", "end"])
            table[["start", "end"]] = table[["start", "end"]].astype(np.int64)
            assert (table.end > table.start).all()
            called = overlap_any(table, rna)
            column = f"{assay}_{group}_{gsm}_peak_detected"
            support[column] = called.astype(int)
            groups[(assay, group)].append(column)
            declared = re.search(r"# total peaks = (\d+)", data)
            if declared:
                assert int(declared.group(1)) == len(table), (member.name, declared.group(1), len(table))
            total_tags = re.search(r"# Total tags = ([\d.]+)", data)
            inventory.append({"gsm": gsm, "title": title, "assay": assay, "group": group,
                               "archive_member": member.name, "sha256_compressed_member": hashlib.sha256(compressed).hexdigest(),
                               "called_regions": len(table), "median_region_width": float(np.median(table.end - table.start)),
                               "promoters_with_any_called_region": int(called.sum()),
                               "reported_total_tags": float(total_tags.group(1)) if total_tags else None,
                               "data_type": "per-sample called regions with normalized tags; not common-region counts",
                               "formal_differential_testing": "not_possible_from_these_thresholded_files_alone"})
    for (assay, group), columns in groups.items():
        assert len(columns) == 2, (assay, group, columns)
        support[f"{assay}_{group}_replicates_detected"] = support[columns].sum(axis=1)
        support[f"{assay}_{group}_both_replicates_detected"] = (support[columns].sum(axis=1) == 2).astype(int)
    support["interpretation"] = "Physical overlap with annotated promoter; peak absence is not quantitative loss; no distal target assignment"
    support.to_csv(outputs / "r12_promoter_peak_support.csv", index=False)
    pd.DataFrame(inventory).to_csv(outputs / "r12_chromatin_file_audit.csv", index=False)
    joined = rna.merge(support.drop(columns=["chromosome", "strand", "annotated_TSS", "promoter_start", "promoter_end"]),
                        on=["gene", "transcript_id"], how="left", validate="one_to_one")
    joined.to_csv(outputs / "r12_gene_molecular_evidence.csv", index=False)
    frozen = ROOT / "inputs" / "r12_analysis_spec.json"
    if frozen.exists():
        modules = json.loads(frozen.read_text())["modules"]
        membership = pd.DataFrame([{"module": name, "gene": gene} for name, genes in modules.items() for gene in genes])
        # Verified one-to-one historic HGNC symbols in this 2017 deposit.
        # Retain both names so missing aliases cannot masquerade as unmeasured.
        legacy = {"BICRA": "GLTSCR1", "BICRAL": "GLTSCR1L"}
        membership["source_gene"] = membership.gene.map(lambda g: g if g in set(joined.gene) else legacy.get(g, g))
        module_result = membership.merge(joined.rename(columns={"gene": "source_gene"}), on="source_gene", how="left", validate="many_to_one")
        raw_expression = pd.read_csv(ROOT / "inputs" / "GSE101966_HCT116_FPKM.txt.gz", sep="\t")
        deposited_symbols = set(raw_expression["Annotation/Divergence"].str.split("|").str[0])
        module_result["response_availability"] = np.where(module_result.interaction_log2_effect.notna(), "tested",
            np.where(module_result.source_gene.isin(deposited_symbols), "present_but_below_prespecified_expression_filter", "not_in_deposit"))
        module_result["interpretation_module"] = "External transcriptional-response check; expression response is not cell dependency"
        module_result.to_csv(outputs / "r12_prespecified_module_gene_response.csv", index=False)
    known_panel = ["MET", "TGFA", "TGFBR2", "PLAU", "EGFR", "ARID1A", "ARID1B"]
    joined[joined.gene.isin(known_panel)].to_csv(outputs / "r12_literature_anchor_gene_response.csv", index=False)

    # The same Ab-number occurs under contradictory genotype labels across
    # assays. It is a filename clue, NOT proof of a shared biological specimen.
    mouse = []
    token_rows = defaultdict(list)
    for record in audit["metadata"].get("GSE164810", []):
        if record["record_type"] != "SAMPLE":
            continue
        fields = record["fields"]
        filenames = [value for key, values in fields.items() if key.startswith("Sample_supplementary_file") for value in values]
        chars = fields.get("Sample_characteristics_ch1", [])
        genotype = next((x.split(":",1)[1].strip() for x in chars if x.lower().startswith("genotype:")), "unknown")
        tokens = sorted(set(re.findall(r"Ab\d+", " ".join(filenames), flags=re.I)))
        row = {"gsm": record["accession"], "title": fields.get("Sample_title", [""])[0],
                "genotype_as_deposited": genotype, "assay": fields.get("Sample_library_strategy", [""])[0],
                "files": ";".join(filenames), "filename_Ab_tokens": ";".join(tokens),
                "SRA_BioSample_links": ";".join(fields.get("Sample_relation", [])),
                "sample_level_pairing": "not_established"}
        mouse.append(row)
        for token in tokens:
            token_rows[token.lower()].append(row)
    flags = []
    for token, rows in token_rows.items():
        genotypes = sorted(set(x["genotype_as_deposited"] for x in rows))
        if len(genotypes) > 1:
            flags.append({"filename_token": token, "genotypes": genotypes,
                          "gsms": [x["gsm"] for x in rows], "status": "unresolved_filename_metadata_inconsistency",
                          "action": "Do not relabel; do not assume matched-mouse multiomics; exclude confirmatory paired analysis pending source resolution"})
    pd.DataFrame(mouse).to_csv(outputs / "r12_gse164810_sample_metadata.csv", index=False)
    conclusion = {"checked_utc": datetime.now(timezone.utc).isoformat(), "host": platform.node(),
        "GSE101966_chromatin_files": len(inventory),
        "formal_chromatin_interaction_status": "Not estimable from independently called per-sample peaks; common-region read counts absent",
        "scope_of_peak_support": "Descriptive same-promoter overlap only; no inferred enhancer-gene targets or differential p values",
        "verified_module_aliases": {"BICRA": "GLTSCR1", "BICRAL": "GLTSCR1L"},
        "alias_primary_sources": ["https://www.ncbi.nlm.nih.gov/gene/29998", "https://www.ncbi.nlm.nih.gov/gene/23506"],
        "GSE164810_genotype_filename_flags": flags,
        "GSE164810_design": "RNA 2 KP + 2 KPA; ATAC 2 KP + 2 KPA; KP ARID1A CUT&Tag n=1 plus IgG n=1",
        "GSE164810_decision": "No relabeling or mouse differential analysis in this phase; public deposit does not establish individual pairing",
        "GSE164810_remaining_gate": "Source publication supplement/submitter mapping must resolve filename inconsistency before matched integration"}
    (outputs / "r12_chromatin_and_mouse_audit_summary.json").write_text(json.dumps(conclusion, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(conclusion, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
