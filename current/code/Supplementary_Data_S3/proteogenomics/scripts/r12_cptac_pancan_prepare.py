#!/usr/bin/env python3
"""Freeze analysis-ready gene matrices and patient metadata; no association tests."""
import argparse
import csv
import re
from pathlib import Path
import numpy as np
import pandas as pd
from r12_cptac_acquire_audit import BAF, patient_id, sha256, write_json, now


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("/root/autodl-tmp/arid1a_baf_r12_cptac/pancan_release"))
    args = parser.parse_args()
    raw = args.root / "raw_processed"
    out = args.root / "audit"
    out.mkdir(exist_ok=True)
    cat = pd.read_csv(args.root / "source_catalog.tsv", sep="\t")
    gm = pd.read_csv(raw / "bcm-luad-mapping-gencode.v34.basic.annotation-mapping.txt.gz", sep="\t")
    gm = gm[["gene", "gene_name"]].dropna().drop_duplicates()
    ambiguous = gm[gm.gene.duplicated(keep=False)]
    ambiguous.to_csv(out / "gene_mapping_ambiguous_ids.tsv", sep="\t", index=False)
    gm = gm.loc[~gm.gene.isin(ambiguous.gene)].set_index("gene")["gene_name"]
    primary = ["ARID1A", "SMARCC2", "SMARCA4", "SMARCD1", "SMARCB1", "SMARCE1", "DPF2",
               "ACTL6A", "PBRM1", "ARID2", "BRD7", "PHF10", "BRD9"]
    all_genes = list(dict.fromkeys(BAF + ["BCL7A", "BCL7B", "BCL7C"]))
    assay_audit, coverage, mappings, flows, metadata_coverage, outputs = [], [], [], [], [], []
    exported_entry_checks = []
    matrices = {}
    source_snapshots = {}
    for _, rec in cat.loc[cat.modality.isin(["rna", "protein", "cnv"])].iterrows():
        file = Path(rec["file"])
        with file.open(encoding="utf-8-sig", newline="") as fh:
            head = next(csv.reader(fh, delimiter="\t"))
        assert len(head[1:]) == len(set(head[1:])), f"Duplicate original columns: {file}"
        d = pd.read_csv(file, sep="\t", index_col=0)
        assert d.columns.size == int(rec.expected_samples), f"Portal sample count mismatch: {file}"
        assert d.index.is_unique, f"Repeated Ensembl measurement rows: {file}"
        vals = d.apply(pd.to_numeric, errors="raise")
        if np.isinf(vals.to_numpy()).any():
            raise ValueError(f"Infinite input measurements: {file}")
        key = (rec.cohort, rec.modality, rec.sample_type)
        ids = [patient_id(x, rec.sample_type) for x in vals.columns]
        assert len(ids) == len(set(ids)), f"Patient canonicalization collision: {file}"
        sid = [f"{rec.cohort}:{p}:{rec.sample_type}" for p in ids]
        map_df = pd.DataFrame({"cohort": rec.cohort, "release": "pancan_data_freeze_v1.2",
                               "modality": rec.modality, "sample_type": rec.sample_type,
                               "raw_sample_id": vals.columns, "patient_id": ids,
                               "sample_id": sid, "source_url": rec.url})
        mappings.append(map_df)
        vals.columns = sid
        annot = pd.DataFrame({"ensembl_versioned_id": vals.index,
                              "gene_symbol": vals.index.map(gm)})
        annot["mapped"] = annot.gene_symbol.notna()
        source_primary = vals.loc[annot.gene_symbol.isin(primary).to_numpy()].copy()
        source_primary.index = annot.loc[annot.gene_symbol.isin(primary), "gene_symbol"].to_numpy()
        assert source_primary.index.is_unique, "Primary gene mapping is ambiguous"
        source_snapshots[key] = source_primary
        annot.to_csv(out / ("_".join(key) + "_gene_mapping.tsv.gz"), sep="\t", index=False, compression="gzip")
        mapped = vals.loc[annot["mapped"].to_numpy()].copy()
        mapped.index = annot.loc[annot.mapped, "gene_symbol"].to_numpy()
        duplicates = mapped.index[mapped.index.duplicated(keep=False)].unique().tolist()
        if duplicates:
            mapped.loc[duplicates].to_csv(out / ("_".join(key) + "_duplicate_symbols.tsv.gz"), sep="\t", compression="gzip")
        assert not set(primary) & set(duplicates), "Primary symbol maps to multiple gene measurements"
        mapped = mapped.groupby(level=0, sort=False).median()
        mapped.index.name = "gene_symbol"
        matrices[key] = mapped
        dest = out / ("_".join(key) + "_gene_symbol.tsv.gz")
        mapped.to_csv(dest, sep="\t", compression="gzip")
        outputs.append({"file": str(dest), "sha256": sha256(dest)})
        assay_audit.append({"cohort": rec.cohort, "modality": rec.modality, "sample_type": rec.sample_type,
                            "n_samples": mapped.shape[1], "n_ensembl_rows": vals.shape[0],
                            "n_symbol_mapped_rows": int(annot.mapped.sum()), "n_final_gene_symbols": len(mapped),
                            "n_unmapped_ensembl_rows": int((~annot.mapped).sum()),
                            "n_duplicate_symbols_collapsed": len(duplicates),
                            "unit": rec.unit, "source_url": rec.url,
                            "processing": "Exact versioned GENCODE34 mapping; duplicate symbols median; no imputation or numerical rescaling"})
        for g in all_genes:
            n = int(mapped.loc[g].notna().sum()) if g in mapped.index else 0
            coverage.append({"cohort": rec.cohort, "modality": rec.modality, "sample_type": rec.sample_type,
                             "gene_symbol": g, "present": g in mapped.index, "n_nonmissing": n,
                             "n_total": len(ids), "fraction_nonmissing": n / len(ids)})
    for cohort in ["LUAD", "LUSC"]:
        metadata = None
        for _, rec in cat.loc[cat.cohort.eq(cohort) & cat.modality.isin(["metadata", "phenotype"])].iterrows():
            m = pd.read_csv(rec["file"], sep="\t", dtype=str)
            idcol = m.columns[0]
            m = m.loc[m[idcol].ne("data_type")].copy()
            ids = [patient_id(x, "tumor") for x in m[idcol]]
            assert len(ids) == len(set(ids)), "Metadata canonical IDs collide"
            m.index = [f"{cohort}:{p}:tumor" for p in ids]
            m.index.name = "sample_id"
            m.columns = [rec.modality + "__" + c for c in m.columns]
            for col in m.columns:
                metadata_coverage.append({"cohort": cohort, "field": col, "n_records": len(m),
                                          "n_nonmissing": int(m[col].notna().sum()),
                                          "candidate_purity": bool(re.search("purity|cellularity", col, re.I)),
                                          "candidate_batch": bool(re.search("batch|plex|tmt|channel|instrument", col, re.I))})
            metadata = m if metadata is None else metadata.join(m, how="outer", validate="one_to_one")
        metadata["cohort"] = cohort
        metadata["patient_id"] = metadata.index.str.split(":").str[1]
        metadata["sample_type"] = "tumor"
        for source, target in [("phenotype__WES_purity", "purity_wes"), ("phenotype__WGS_purity", "purity_wgs")]:
            metadata[target] = pd.to_numeric(metadata[source], errors="coerce") if source in metadata else np.nan
        metadata["purity_primary"] = metadata["purity_wes"]
        metadata["purity_primary_source"] = "WES_purity, same pan-cancer freeze; not imputed from WGS"
        metadata.to_csv(out / f"{cohort}_tumor_metadata_full.tsv", sep="\t")
        for tissue in ["tumor", "normal"]:
            sets = {m: set(matrices[(cohort, m, tissue)].columns)
                    for m in ["protein", "rna", "cnv"] if (cohort, m, tissue) in matrices}
            union = set.union(*sets.values())
            sample_df = pd.DataFrame(index=sorted(union))
            sample_df.index.name = "sample_id"
            sample_df["cohort"] = cohort
            sample_df["patient_id"] = sample_df.index.str.split(":").str[1]
            sample_df["sample_type"] = tissue
            for mod in ["protein", "rna", "cnv"]:
                sample_df[f"has_{mod}"] = sample_df.index.isin(sets.get(mod, set()))
            both = sets["protein"] & sets["rna"]
            joint = both & sets.get("cnv", set())
            sample_df["rna_protein_matched"] = sample_df.index.isin(both)
            sample_df["rna_protein_cnv_matched"] = sample_df.index.isin(joint)
            sample_df.to_csv(out / f"{cohort}_{tissue}_sample_match.tsv", sep="\t")
            for mod in sets:
                flows.append({"cohort": cohort, "sample_type": tissue, "step": mod, "n_samples": len(sets[mod])})
            for label, samples in [("protein_RNA", both), ("protein_RNA_CNV", joint)]:
                flows.append({"cohort": cohort, "sample_type": tissue, "step": label, "n_samples": len(samples)})
            analysis_ids = sorted(joint if tissue == "tumor" else both)
            for mod in sets:
                matrices[(cohort, mod, tissue)].reindex(columns=analysis_ids).to_csv(
                    out / f"{cohort}_{mod}_{tissue}_matched.tsv.gz", sep="\t", compression="gzip")
            if tissue == "tumor":
                md = metadata.reindex(analysis_ids)
                assert md["patient_id"].notna().all(), "Missing clinical mapping in matched tumor assay"
                md.to_csv(out / f"{cohort}_tumor_metadata_matched.tsv", sep="\t")
                flows.append({"cohort": cohort, "sample_type": tissue, "step": "protein_RNA_CNV_WES_purity",
                              "n_samples": int(md.purity_primary.notna().sum())})
                # Fixed-panel feasibility per anchor-partner pair; no correlations.
                pmat, rmat, cmat = [matrices[(cohort, m, tissue)].reindex(columns=analysis_ids)
                                    for m in ["protein", "rna", "cnv"]]
                pair_counts = []
                for gene in all_genes:
                    mask = pd.Series(True, index=analysis_ids)
                    for g in ["ARID1A", gene]:
                        for mat in [pmat, rmat, cmat]:
                            mask &= mat.loc[g].notna() if g in mat.index else False
                    pair_counts.append({"cohort": cohort, "gene_symbol": gene,
                                        "n_anchor_partner_RNA_CNV_protein_complete": int(mask.sum()),
                                        "n_also_WES_purity_complete": int((mask & md.purity_primary.notna()).sum()),
                                        "primary_fixed_panel": gene in primary})
                pd.DataFrame(pair_counts).to_csv(out / f"{cohort}_anchor_partner_pair_coverage.tsv", sep="\t", index=False)
                ready = args.root / "analysis_ready" / cohort
                ready.mkdir(parents=True, exist_ok=True)
                for mod, mat in [("protein", pmat), ("rna", rmat), ("cnv", cmat)]:
                    export = mat.reindex(index=all_genes, columns=analysis_ids).T
                    export.index.name = "sample_id"
                    file = ready / f"{mod}.csv"
                    export.to_csv(file)
                    check = pd.read_csv(file, index_col=0)
                    assert check.index.tolist() == analysis_ids
                    assert check.columns.tolist() == all_genes
                    np.testing.assert_allclose(check.to_numpy(), export.to_numpy(), rtol=1e-12, atol=1e-12, equal_nan=True)
                    for sid in [analysis_ids[0], analysis_ids[len(analysis_ids) // 2], analysis_ids[-1]]:
                        for gene in primary:
                            if gene not in mat.index:
                                continue
                            exported_entry_checks.append({"cohort": cohort, "modality": mod,
                                "sample_id": sid, "gene_symbol": gene,
                                "source_downloaded_value": source_snapshots[(cohort, mod, tissue)].loc[gene, sid],
                                "prepared_matrix_value": mat.loc[gene, sid],
                                "roundtrip_export_value": check.loc[sid, gene],
                                "equal_within_1e_12": bool(np.isclose(source_snapshots[(cohort, mod, tissue)].loc[gene, sid], check.loc[sid, gene], rtol=1e-12, atol=1e-12, equal_nan=True))})
                    outputs.append({"file": str(file), "sha256": sha256(file)})
                meta_ready = md[["patient_id", "sample_type", "cohort", "purity_primary", "purity_wgs"]].rename(columns={"purity_primary": "purity"})
                assert meta_ready.patient_id.is_unique
                assert meta_ready.sample_type.eq("tumor").all()
                meta_ready.to_csv(ready / "metadata.csv")
                write_json(ready / "input_manifest.json", {"cohort": cohort,
                    "release": "CPTAC pan-cancer data_freeze_v1.2_reorganized",
                    "source_portal": "https://linkedomics.org/data_download/CPTAC-pancan-" + ("LSCC" if cohort == "LUSC" else cohort) + "/",
                    "sample_orientation": "Rows are matched tumor samples, columns are exact mapped gene symbols",
                    "n_tumors": len(meta_ready), "n_WES_purity_nonmissing": int(meta_ready.purity.notna().sum()),
                    "protein_scale": "Log2 reference intensity normalized protein abundance (portal: log2 MS1 intensity); kept as supplied",
                    "rna_scale": "Upper quartile 1500 normalized log2 RSEM gene expression; kept as supplied",
                    "cnv_scale": "WES gene-level log2 copy-number ratio; kept as supplied",
                    "purity_source": "phenotype table WES_purity; no WGS substitution or imputation",
                    "batch": "Not available in supplied meta/phenotype tables; no batch column fabricated",
                    "gene_mapping": "Exact versioned Ensembl to GENCODE34 symbols; median for duplicate symbols, primary genes have no duplicate symbols",
                    "sample_mapping": "Exact patient key with cohort and tumor type; no normal tissue copied into tumor data",
                    "absent_gene_handling": "Unmeasured prespecified genes are all-NaN columns, not zero-filled",
                    "metadata_sha256": sha256(ready / "metadata.csv"),
                    "assay_sha256": {m: sha256(ready / f"{m}.csv") for m in ["protein", "rna", "cnv"]}})
    pd.concat(mappings, ignore_index=True).to_csv(out / "sample_mapping.tsv", sep="\t", index=False)
    pd.DataFrame(assay_audit).to_csv(out / "assay_audit.tsv", sep="\t", index=False)
    pd.DataFrame(coverage).to_csv(out / "baf_coverage.tsv", sep="\t", index=False)
    pd.DataFrame(metadata_coverage).to_csv(out / "metadata_coverage.tsv", sep="\t", index=False)
    pd.DataFrame(flows).to_csv(out / "sample_flow.tsv", sep="\t", index=False)
    pd.DataFrame(exported_entry_checks).to_csv(out / "export_entryvalue_spotchecks.tsv", sep="\t", index=False)
    write_json(out / "preparation_manifest.json", {"finished_utc": now(), "input_release": "data_freeze_v1.2_reorganized",
        "sample_key": "cohort:canonical_patient_id:source_tissue_type; no fuzzy joins",
        "gene_annotation": "Exact versioned Ensembl gene identifiers using official cptac GENCODE34 mapping",
        "scope": "Acquisition, annotation and coverage only. No associations or fitted models.",
        "unavailable_covariates": "No TMT/plex batch column found in downloaded metadata; do not claim batch correction",
        "normal_tissue": "Retained separately; somatic CNV and tumor purity not copied into normal samples",
        "outputs": outputs})
    print(pd.DataFrame(flows).to_string(index=False))
    cov = pd.DataFrame(coverage)
    print(cov.loc[cov.sample_type.eq("tumor")].pivot(index="gene_symbol", columns=["cohort", "modality"], values="n_nonmissing").to_string())


if __name__ == "__main__":
    main()
