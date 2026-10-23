"""Reproduce only final six-partner interface metrics from fixed source inputs.

The parser, contact definition and residue mapping logic use the archived
implementation; this entry point has no dependency on historical metric tables.
It is provided for reproducibility and was not executed during data curation.
"""
import argparse
import csv
import json
from pathlib import Path

from scipy.stats import spearmanr
from r14_recompute_af_mapping import (
    ORDER, CHAINS, GAPPED, SEGMENTS, read_csv, load_chains, ca_sequence, contacts,
    agreement, sha, write_csv,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=False)
    inputs = args.root / "inputs"
    reference = load_chains(inputs / "6LTJ.pdb")
    reference_sequence = ca_sequence(inputs / "6LTJ.pdb", "L")
    reference_residues = list(reference_sequence)
    assert reference_residues == [r for lo, hi in SEGMENTS for r in range(lo, hi + 1)]
    continuous = read_csv(inputs / "SMARCB1/input.csv")[0]["sequence"].split(":")[0]
    assert len(continuous) == 647 and len(reference_residues) == 466
    assert "".join(continuous[r - 1639] for r in reference_residues) == "".join(reference_sequence.values())
    rows, crosswalk, exact, sources = [], [], [], {"inputs/6LTJ.pdb":sha(inputs / "6LTJ.pdb")}
    for gene in ORDER:
        folder = inputs / gene
        seq_a, seq_b = read_csv(folder / "input.csv")[0]["sequence"].split(":")
        config = json.loads((folder / "config.json").read_text())
        assert config["model_type"] == "alphafold2_multimer_v3" and config["use_templates"] is False
        pdbs, scores = list(folder.glob("*unrelaxed_rank_001*.pdb")), list(folder.glob("*scores_rank_001*.json"))
        assert len(pdbs) == len(scores) == 1
        pdb, score_path = pdbs[0], scores[0]
        score = json.loads(score_path.read_text())
        chains = load_chains(pdb)
        aa, bb = ca_sequence(pdb, "A"), ca_sequence(pdb, "B")
        assert list(aa) == list(range(1, len(seq_a)+1))
        assert "".join(aa.values()) == seq_a and "".join(bb.values()) == seq_b
        if gene in GAPPED:
            assert seq_a == "".join(reference_sequence.values()) and continuous.find(seq_a) == -1
            mapping = dict(zip(aa, reference_residues))
        else:
            assert seq_a == continuous
            mapping = {r:r+1638 for r in aa}
        assert len(mapping) == len(set(mapping.values()))
        for index, residue in mapping.items():
            assert aa[index] == continuous[residue-1639]
            crosswalk.append({"partner":gene,"af_chain_A_residue":index,"amino_acid":aa[index],"mapped_6LTJ_residue":residue})
        predicted, _, _ = contacts(chains["A"], chains["B"])
        experimental, _, _ = contacts(reference["L"], reference[CHAINS[gene]])
        overlap = agreement({mapping[r] for r in predicted}, experimental)
        exact.append(overlap["jaccard"])
        rows.append({"partner":gene,"chain_6LTJ":CHAINS[gene],"af_iptm":score["iptm"],"af_ptm":score["ptm"],
                     "cryoem_iface_n":len(experimental),"af_iface_n":len(predicted),"overlap_n":overlap["overlap_n"],
                     "recall":round(overlap["recall"],3),"jaccard":round(overlap["jaccard"],3)})
        for path in (folder/"input.csv",folder/"config.json",pdb,score_path):
            sources[path.relative_to(args.root).as_posix()] = sha(path)
    write_csv(args.outdir/"interface_metrics.csv",rows)
    write_csv(args.outdir/"residue_crosswalk.csv",crosswalk)
    correlations = {}
    for label, values in (("exact",exact),("display_3dp",[r["jaccard"] for r in rows])):
        rho, p = spearmanr([r["af_iptm"] for r in rows],values)
        correlations[label] = {"n":6,"spearman_rho":float(rho),"two_sided_scipy_asymptotic_p":float(p),"jaccard_values_in_fixed_order":values}
    (args.outdir/"confidence_interface_correlations.json").write_text(json.dumps(correlations,indent=2))
    (args.outdir/"input_hashes.json").write_text(json.dumps(sources,indent=2))


if __name__ == "__main__":
    main()
