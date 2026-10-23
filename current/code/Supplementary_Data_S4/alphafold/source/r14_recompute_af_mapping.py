#!/usr/bin/env python3
"""Repair three gapped ARID1A residue mappings; reuse the fixed six archived rank-1 models.

Run only on the user's server. No AlphaFold invocation, new pose selection,
experimental-contact modification, or manuscript edit is performed.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy
from scipy.stats import spearmanr

AA = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(), "ARNDCQEGHILKMFPSTWYV"))
ORDER = ["SMARCB1", "SMARCA4", "SMARCD1", "SMARCC2", "SMARCE1", "DPF2"]
CHAINS = dict(zip(ORDER, "MIPNQR"))
GAPPED = {"SMARCA4", "SMARCD1", "SMARCC2"}
SEGMENTS = [(1639,1746), (1802,1862), (1954,2025), (2047,2210), (2225,2285)]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_chains(pdb):
    """Preserve the legacy ATOM/non-H parser exactly for distance calculations."""
    ch = {}
    for line in Path(pdb).read_text().splitlines():
        if line.startswith("ATOM") and line[76:78].strip() != "H":
            chain, rn = line[21], int(line[22:26])
            ch.setdefault(chain, {}).setdefault(rn, []).append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return {c: {r: np.array(v) for r, v in d.items()} for c, d in ch.items()}


def ca_sequence(pdb, chain):
    out = {}
    for line in Path(pdb).read_text().splitlines():
        if line.startswith("ATOM") and line[21] == chain and line[12:16].strip() == "CA":
            rn = int(line[22:26])
            assert not line[26].strip(), "Insertion codes require an explicit mapping audit"
            amino = AA[line[17:20].strip()]
            if rn in out:
                assert out[rn] == amino
            out[rn] = amino
    return out


def contacts(a, b, threshold=4.5):
    """Identical residue, unique-atom and atom-pair definitions to af6_recompute.py."""
    all_b = np.vstack(list(b.values()))
    interface, natoms, npairs = set(), 0, 0
    for rn, coord in a.items():
        distances = np.sqrt(((coord[:, None, :] - all_b[None, :, :]) ** 2).sum(-1))
        hit = distances < threshold
        if hit.any():
            interface.add(rn)
            natoms += int(hit.any(axis=1).sum())
            npairs += int(hit.sum())
    return interface, natoms, npairs


def agreement(predicted, reference):
    overlap = predicted & reference
    union = predicted | reference
    return {"overlap_n": len(overlap), "recall": len(overlap) / len(reference),
            "jaccard": len(overlap) / len(union) if union else 0.0}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=False)
    inp = args.root / "inputs"
    cryo_path = inp / "6LTJ.pdb"
    old_path = inp / "af_interface_vs_cryoem_6partners_RECOMPUTED.csv"
    old = {r["partner"]: r for r in read_csv(old_path)}
    assert list(old) == ORDER
    reference_chains = load_chains(cryo_path)
    reference_seq = ca_sequence(cryo_path, "L")
    reference_residues = list(reference_seq)
    gapped_reference_sequence = "".join(reference_seq.values())
    continuous_sequence = read_csv(inp / "SMARCB1/input.csv")[0]["sequence"].split(":")[0]
    assert len(reference_residues) == 466 and len(continuous_sequence) == 647
    # Two independent index constructions: observed ordered PDB CA identifiers and piecewise source segments.
    piecewise_residues = [r for lo, hi in SEGMENTS for r in range(lo, hi + 1)]
    assert reference_residues == piecewise_residues
    assert "".join(continuous_sequence[r - 1639] for r in reference_residues) == gapped_reference_sequence
    rows, comparisons, mapping_rows, checks, sets = [], [], [], [], {}
    input_paths = [cryo_path, old_path]
    for gene in ORDER:
        folder = inp / gene
        csv_path, config_path = folder / "input.csv", folder / "config.json"
        seq_a, seq_b = read_csv(csv_path)[0]["sequence"].split(":")
        config = json.loads(config_path.read_text())
        assert config["model_type"] == "alphafold2_multimer_v3" and config["use_templates"] is False
        pdbs = list(folder.glob("*unrelaxed_rank_001*.pdb"))
        score_files = list(folder.glob("*scores_rank_001*.json"))
        assert len(pdbs) == len(score_files) == 1
        pdb, score_path = pdbs[0], score_files[0]
        score = json.loads(score_path.read_text())
        assert score["iptm"] == float(old[gene]["af_iptm"])
        assert score["ptm"] == float(old[gene]["af_ptm"])
        af_chains = load_chains(pdb)
        af_seq_a, af_seq_b = ca_sequence(pdb, "A"), ca_sequence(pdb, "B")
        assert list(af_seq_a) == list(range(1, len(seq_a) + 1))
        assert "".join(af_seq_a.values()) == seq_a
        assert "".join(af_seq_b.values()) == seq_b
        if gene in GAPPED:
            assert seq_a == gapped_reference_sequence
            assert continuous_sequence.find(seq_a) == -1
            mapping = dict(zip(af_seq_a, reference_residues))
            independent_mapping = dict(enumerate(piecewise_residues, start=1))
        else:
            assert seq_a == continuous_sequence
            mapping = {r: r + 1638 for r in af_seq_a}
            independent_mapping = dict(enumerate(range(1639, 2286), start=1))
        assert mapping == independent_mapping
        assert len(mapping) == len(set(mapping.values()))
        for af_rn, source_rn in mapping.items():
            assert af_seq_a[af_rn] == continuous_sequence[source_rn - 1639]
            mapping_rows.append({"partner": gene, "af_chain_A_residue": af_rn, "amino_acid": af_seq_a[af_rn],
                                 "legacy_plus1638": af_rn + 1638, "mapped_6LTJ_residue": source_rn,
                                 "mapping_changed": source_rn != af_rn + 1638})
        raw_interface, pred_atoms, pred_pairs = contacts(af_chains["A"], af_chains["B"])
        reference_interface, ref_atoms, ref_pairs = contacts(reference_chains["L"], reference_chains[CHAINS[gene]])
        assert old[gene]["chain_6LTJ"] == CHAINS[gene]
        legacy_interface = {r + 1638 for r in raw_interface}
        mapped_interface = {mapping[r] for r in raw_interface}
        old_agreement = agreement(legacy_interface, reference_interface)
        new_agreement = agreement(mapped_interface, reference_interface)
        assert len(reference_interface) == int(old[gene]["cryoem_iface_n"])
        assert len(raw_interface) == int(old[gene]["af_iface_n"])
        assert old_agreement["overlap_n"] == int(old[gene]["overlap_n"])
        assert round(old_agreement["jaccard"], 3) == float(old[gene]["jaccard"])
        assert round(old_agreement["recall"], 3) == float(old[gene]["recall"])
        if gene not in GAPPED:
            assert new_agreement == old_agreement
        row = {"partner": gene, "chain_6LTJ": CHAINS[gene], "af_iptm": score["iptm"], "af_ptm": score["ptm"],
               "cryoem_iface_n": len(reference_interface), "af_iface_n": len(raw_interface),
               "overlap_n": new_agreement["overlap_n"], "recall": round(new_agreement["recall"], 3),
               "jaccard": round(new_agreement["jaccard"], 3)}
        rows.append(row)
        comparisons.append({**row, "arid1a_input_length": len(seq_a), "partner_input_length": len(seq_b),
                            "old_overlap_n": old_agreement["overlap_n"], "old_recall_exact": old_agreement["recall"],
                            "new_recall_exact": new_agreement["recall"], "old_jaccard_exact": old_agreement["jaccard"],
                            "new_jaccard_exact": new_agreement["jaccard"], "old_jaccard_display": float(old[gene]["jaccard"]),
                            "mapping_changed": gene in GAPPED, "fixed_rank1_filename": pdb.name})
        checks.append({"partner":gene, "legacy_csv_fully_reproduced": True, "input_matches_AF_chain_A":True,
                       "input_matches_AF_chain_B":True, "mapping_matches_independent_piecewise_construction":True,
                       "mapped_residue_amino_acids_match_647aa_reference":True,
                       "experimental_chain_preserved":True, "experimental_interface_residues_preserved":True,
                       "prediction_atom_pairs_unchanged":pred_pairs, "experimental_atom_pairs_unchanged":ref_pairs,
                       "predicted_interface_n_preserved":len(raw_interface), "reference_interface_n_preserved":len(reference_interface)})
        sets[gene] = {"predicted_interface_AF_indices":sorted(raw_interface),
                      "legacy_predicted_interface_6LTJ_indices":sorted(legacy_interface),
                      "corrected_predicted_interface_6LTJ_indices":sorted(mapped_interface),
                      "unchanged_reference_interface_6LTJ_indices":sorted(reference_interface)}
        input_paths += [csv_path, config_path, pdb, score_path]
    correlations = {}
    iptm = [r["af_iptm"] for r in rows]
    for name, jac in [("legacy_display_3dp", [float(old[g]["jaccard"]) for g in ORDER]),
                      ("corrected_display_3dp", [r["jaccard"] for r in rows]),
                      ("legacy_exact", [r["old_jaccard_exact"] for r in comparisons]),
                      ("corrected_exact", [r["new_jaccard_exact"] for r in comparisons])]:
        statistic, pvalue = spearmanr(iptm, jac)
        correlations[name] = {"n":6, "spearman_rho":float(statistic), "two_sided_scipy_asymptotic_p":float(pvalue),
                              "jaccard_values_in_fixed_order":jac}
    write_csv(args.outdir / "af_interface_vs_cryoem_6partners_R14.csv", rows)
    write_csv(args.outdir / "r14_af_old_new_comparison.csv", comparisons)
    write_csv(args.outdir / "r14_af_residue_mapping.csv", mapping_rows)
    write_csv(args.outdir / "r14_af_consistency_checks.csv", checks)
    (args.outdir / "r14_af_interface_sets.json").write_text(json.dumps(sets, indent=2))
    (args.outdir / "r14_af_correlations.json").write_text(json.dumps(correlations, indent=2))
    manifest = {"created_utc":datetime.now(timezone.utc).isoformat(), "host":platform.node(),
                "execution":"user_5090_server; only existing fixed rank-1 interface bookkeeping and the same six-partner Spearman association",
                "model_rerun":False, "new_pose_selection":False, "experimental_contacts_modified":False,
                "metric":"ARID1A-side residue interface; legacy ATOM/non-H parser and strict distance <4.5 Angstrom; Jaccard intersection/union",
                "SMARCC2_reference":"Existing single copy N retained, no N+O pooling",
                "all_legacy_table_entries_reproduced":True, "all_mapping_checks_passed":True,
                "changes_limited_to_three_gapped_ARID1A_index_maps":True,
                "partners_fixed_order":ORDER, "software":{"numpy":np.__version__,"scipy":scipy.__version__},
                "script_sha256":sha(__file__),
                "input_sha256":{str(p.relative_to(args.root)):sha(p) for p in input_paths},
                "output_sha256":{p.name:sha(p) for p in sorted(args.outdir.iterdir())}}
    (args.outdir / "r14_af_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps({"comparison":[{k:r[k] for k in ("partner","old_overlap_n","overlap_n","old_jaccard_exact","new_jaccard_exact")} for r in comparisons],
                      "correlations":correlations,"all_integrity_checks_passed":True},indent=2))


if __name__ == "__main__":
    main()
