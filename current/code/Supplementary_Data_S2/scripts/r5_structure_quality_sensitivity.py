# -*- coding: utf-8 -*-
"""Validate the local structural basis of the ARID1A contact analysis.

Runs on the server against deposited mmCIF coordinates and wwPDB validation XML.
It reports global/chain map-model quality, residue-level Q scores at the ARID1A
interface, contact calls over multiple distance cutoffs, a bottom-decile-Q
coordinate sensitivity analysis, and genetic group-comparison sensitivities.

The output is descriptive/robustness evidence.  A coordinate model remains an
interpretation of the cryo-EM map; no result here is called a mechanism.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser
from scipy import stats
from scipy.spatial import cKDTree


CHAIN_MAP = {
    "6LTJ": {
        "ARID1A": ["L"], "SMARCA4": ["I"], "ACTL6A": ["J"],
        "ACTB": ["K"], "SMARCB1": ["M"], "SMARCC2": ["N", "O"],
        "SMARCD1": ["P"], "SMARCE1": ["Q"], "DPF2": ["R"],
    },
    "6LTH": {
        "ARID1A": ["L"], "SMARCA4": ["I"], "SMARCB1": ["M"],
        "SMARCC2": ["N", "O"], "SMARCD1": ["P"],
        "SMARCE1": ["Q"], "DPF2": ["R"],
    },
}

CONTACT_GENES = ["SMARCC2", "SMARCA4", "SMARCD1", "SMARCB1", "SMARCE1", "DPF2"]
NONCONTACT_GENES = ["ACTB", "ACTL6A", "PBRM1", "ARID2", "BRD7", "PHF10", "BRD9", "BICRA"]
INTEGRATIVE_CONTACT = "SMARCC1"
CUTOFFS = (4.0, 4.5, 5.0)
SEED = 20260901


def load_validation(path):
    root = ET.fromstring(gzip.open(path, "rb").read())
    entry = root.find("Entry").attrib
    chains = {e.attrib.get("chain"): e.attrib for e in root.findall("ModelledEntityInstance")}
    residues = {}
    for e in root.findall("ModelledSubgroup"):
        a = e.attrib
        try:
            key = (a["chain"], int(a["resnum"]), a.get("icode", " ").strip())
        except (KeyError, ValueError):
            continue
        residues[key] = {
            "q_score": float(a["Q_score"]) if "Q_score" in a else np.nan,
            "residue_inclusion": float(a["residue_inclusion"])
            if "residue_inclusion" in a else np.nan,
            "rama": a.get("rama"), "rota": a.get("rota"),
        }
    return entry, chains, residues


def residue_key(chain_id, residue):
    return (chain_id, int(residue.id[1]), str(residue.id[2]).strip())


def atoms_for(structure, chains, qmap=None, q_floor_by_chain=None):
    coords, keys = [], []
    for chain_id in chains:
        chain = structure[0][chain_id]
        for residue in chain:
            if residue.id[0] != " ":
                continue
            key = residue_key(chain_id, residue)
            if q_floor_by_chain is not None:
                q = qmap.get(key, {}).get("q_score", np.nan)
                floor = q_floor_by_chain.get(chain_id, -np.inf)
                if not np.isfinite(q) or q < floor:
                    continue
            for atom in residue:
                if (atom.element or "").upper() == "H":
                    continue
                coords.append(atom.coord)
                keys.append(key)
    return np.asarray(coords, dtype=float), keys


def contact_metrics(a_xyz, a_keys, b_xyz, b_keys, cutoff):
    if len(a_xyz) == 0 or len(b_xyz) == 0:
        return dict(pair_count=0, min_distance=np.nan,
                    a_residues=set(), b_residues=set())
    ta, tb = cKDTree(a_xyz), cKDTree(b_xyz)
    neighbours = ta.query_ball_tree(tb, cutoff)
    pair_count = sum(len(v) for v in neighbours)
    a_res, b_res = set(), set()
    for i, js in enumerate(neighbours):
        if js:
            a_res.add(a_keys[i])
            b_res.update(b_keys[j] for j in js)
    min_distance = float(ta.query(b_xyz, k=1)[0].min())
    return dict(pair_count=pair_count, min_distance=min_distance,
                a_residues=a_res, b_residues=b_res)


def chain_q_floors(qmap, chains, quantile=0.10):
    out = {}
    for chain in chains:
        vals = [v["q_score"] for (c, _, _), v in qmap.items()
                if c == chain and np.isfinite(v["q_score"])]
        out[chain] = float(np.quantile(vals, quantile)) if vals else -np.inf
    return out


def bootstrap_median_difference(x, y, n_boot=20000):
    rng = np.random.default_rng(SEED)
    vals = np.empty(n_boot)
    for i in range(n_boot):
        xb = x[rng.integers(0, len(x), len(x))]
        yb = y[rng.integers(0, len(y), len(y))]
        vals[i] = np.median(xb) - np.median(yb)
    return [float(v) for v in np.quantile(vals, [0.025, 0.975])]


def cliffs_delta(x, y):
    return float((np.greater.outer(x, y).sum() - np.less.outer(x, y).sum()) /
                 (len(x) * len(y)))


def group_test(df, cohort, contact, noncontact, label):
    col = f"{cohort}_rho"
    x = df.set_index("gene").loc[contact, col].astype(float).to_numpy()
    y = df.set_index("gene").loc[noncontact, col].astype(float).to_numpy()
    u = stats.mannwhitneyu(x, y, alternative="two-sided", method="exact")
    ci_low, ci_high = bootstrap_median_difference(x, y)
    return {
        "analysis": label, "cohort": cohort,
        "n_contacting": len(x), "n_noncontacting": len(y),
        "median_contacting": float(np.median(x)),
        "median_noncontacting": float(np.median(y)),
        "median_difference": float(np.median(x) - np.median(y)),
        "bootstrap_95ci_low": ci_low,
        "bootstrap_95ci_high": ci_high,
        "cliffs_delta": cliffs_delta(x, y),
        "mannwhitney_U": float(u.statistic), "exact_two_sided_P": float(u.pvalue),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--rho-csv", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    parser = MMCIFParser(QUIET=True, auth_chains=True, auth_residues=True)
    all_contact_rows, global_rows, chain_rows, residue_rows = [], [], [], []
    interface_by_structure = {}
    qmaps = {}

    for pdbid in ("6LTJ", "6LTH"):
        cif = os.path.join(args.input_dir, f"{pdbid}.cif")
        xml = os.path.join(args.input_dir, f"{pdbid.lower()}_validation.xml.gz")
        structure = parser.get_structure(pdbid, cif)
        entry, chain_meta, qmap = load_validation(xml)
        qmaps[pdbid] = qmap
        global_rows.append({
            "pdbid": pdbid,
            "resolution_A": float(entry["PDB-resolution"]),
            "emdb_id": entry.get("emdb_id"),
            "global_Q_score": float(entry.get("Q_score", "nan")),
            "global_Q_absolute_percentile": float(entry.get("absolute-percentile-qrelative", "nan")),
            "global_Q_same_resolution_percentile": float(entry.get("relative-percentile-qrelative", "nan")),
            "atom_inclusion_backbone": float(entry.get("atom_inclusion_backbone", "nan")),
            "atom_inclusion_all_atoms": float(entry.get("atom_inclusion_all_atoms", "nan")),
            "clashscore": float(entry.get("clashscore", "nan")),
            "rama_outliers_percent": float(entry.get("percent-rama-outliers", "nan")),
            "rotamer_outliers_percent": float(entry.get("percent-rota-outliers", "nan")),
        })
        for gene, chains in CHAIN_MAP[pdbid].items():
            for chain in chains:
                meta = chain_meta.get(chain, {})
                vals = [v["q_score"] for (c, _, _), v in qmap.items()
                        if c == chain and np.isfinite(v["q_score"])]
                inc = [v["residue_inclusion"] for (c, _, _), v in qmap.items()
                       if c == chain and np.isfinite(v["residue_inclusion"])]
                chain_rows.append({
                    "pdbid": pdbid, "gene": gene, "chain": chain,
                    "n_residues_with_Q": len(vals),
                    "median_residue_Q": float(np.median(vals)) if vals else np.nan,
                    "mean_residue_Q": float(np.mean(vals)) if vals else np.nan,
                    "mean_residue_inclusion": float(np.mean(inc)) if inc else np.nan,
                    "deposited_chain_Q": float(meta.get("Q_score", "nan")),
                    "deposited_chain_average_residue_inclusion":
                        float(meta.get("average_residue_inclusion", "nan")),
                })

        all_chains = sorted({c for cs in CHAIN_MAP[pdbid].values() for c in cs})
        floors = chain_q_floors(qmap, all_chains, 0.10)
        a_xyz, a_keys = atoms_for(structure, CHAIN_MAP[pdbid]["ARID1A"])
        aq_xyz, aq_keys = atoms_for(
            structure, CHAIN_MAP[pdbid]["ARID1A"], qmap, floors
        )
        interface_by_structure[pdbid] = set()
        for gene, chains in CHAIN_MAP[pdbid].items():
            if gene == "ARID1A":
                continue
            b_xyz, b_keys = atoms_for(structure, chains)
            bq_xyz, bq_keys = atoms_for(structure, chains, qmap, floors)
            for cutoff in CUTOFFS:
                m = contact_metrics(a_xyz, a_keys, b_xyz, b_keys, cutoff)
                mq = contact_metrics(aq_xyz, aq_keys, bq_xyz, bq_keys, cutoff)
                if cutoff == 4.5 and gene in CONTACT_GENES:
                    interface_by_structure[pdbid].update(m["a_residues"])
                all_contact_rows.append({
                    "pdbid": pdbid, "gene": gene, "cutoff_A": cutoff,
                    "pair_count": m["pair_count"],
                    "min_distance_A": m["min_distance"],
                    "arid1a_interface_residues": len(m["a_residues"]),
                    "partner_interface_residues": len(m["b_residues"]),
                    "pair_count_after_bottom10pct_Q_filter": mq["pair_count"],
                    "contact_call": int(m["pair_count"] > 0),
                    "contact_call_after_Q_filter": int(mq["pair_count"] > 0),
                })

        arid_keys = [k for k in qmap if k[0] == "L" and np.isfinite(qmap[k]["q_score"])]
        interface = interface_by_structure[pdbid]
        for key in arid_keys:
            residue_rows.append({
                "pdbid": pdbid, "chain": key[0], "resnum": key[1],
                "icode": key[2], "is_interface_4.5A": int(key in interface),
                "q_score": qmap[key]["q_score"],
                "residue_inclusion": qmap[key]["residue_inclusion"],
            })

    contact_df = pd.DataFrame(all_contact_rows)
    contact_df.to_csv(os.path.join(args.out, "structure_contact_cutoff_sensitivity.csv"), index=False)
    pd.DataFrame(global_rows).to_csv(os.path.join(args.out, "structure_global_validation.csv"), index=False)
    pd.DataFrame(chain_rows).to_csv(os.path.join(args.out, "structure_chain_validation.csv"), index=False)
    residue_df = pd.DataFrame(residue_rows)
    residue_df.to_csv(os.path.join(args.out, "arid1a_residue_Q_validation.csv"), index=False)

    q_summary = []
    for pdbid, d in residue_df.groupby("pdbid"):
        x = d[d["is_interface_4.5A"] == 1].q_score.to_numpy()
        y = d[d["is_interface_4.5A"] == 0].q_score.to_numpy()
        u = stats.mannwhitneyu(x, y, alternative="two-sided")
        q_summary.append({
            "pdbid": pdbid, "n_interface": len(x), "n_noninterface": len(y),
            "median_Q_interface": float(np.median(x)),
            "median_Q_noninterface": float(np.median(y)),
            "mannwhitney_two_sided_P": float(u.pvalue),
        })
    pd.DataFrame(q_summary).to_csv(os.path.join(args.out, "arid1a_interface_Q_summary.csv"), index=False)

    # Paired Q-score comparison of the same ARID1A residues between deposits.
    q1 = residue_df[residue_df.pdbid == "6LTJ"].set_index(["chain", "resnum", "icode"])
    q2 = residue_df[residue_df.pdbid == "6LTH"].set_index(["chain", "resnum", "icode"])
    shared = q1.join(q2, lsuffix="_6LTJ", rsuffix="_6LTH", how="inner")
    shared.to_csv(os.path.join(args.out, "arid1a_Q_6LTJ_vs_6LTH.csv"))
    paired_q = stats.wilcoxon(shared.q_score_6LTH, shared.q_score_6LTJ,
                              alternative="greater", zero_method="wilcox")

    rho = pd.read_csv(args.rho_csv)
    tests = []
    for cohort in ("lung", "pancancer"):
        tests.append(group_test(
            rho, cohort, CONTACT_GENES + [INTEGRATIVE_CONTACT], NONCONTACT_GENES,
            "original_binary_partition_7v8"
        ))
        tests.append(group_test(
            rho, cohort, CONTACT_GENES, NONCONTACT_GENES,
            "experimental_6LTJ_6LTH_contacts_only_6v8"
        ))
        tests.append(group_test(
            rho, cohort, CONTACT_GENES, ["ACTB", "ACTL6A"],
            "subunits_physically_present_in_6LTJ_only_6v2"
        ))
        # Leave-one-gene sensitivity of the original 7-v-8 partition.
        for gene in CONTACT_GENES + [INTEGRATIVE_CONTACT] + NONCONTACT_GENES:
            c = [g for g in CONTACT_GENES + [INTEGRATIVE_CONTACT] if g != gene]
            n = [g for g in NONCONTACT_GENES if g != gene]
            tests.append(group_test(rho, cohort, c, n, f"leave_one_out:{gene}"))
    pd.DataFrame(tests).to_csv(os.path.join(args.out, "structural_group_sensitivity.csv"), index=False)

    summary = {
        "6LTJ_interface_residue_count": len(interface_by_structure["6LTJ"]),
        "6LTH_interface_residue_count": len(interface_by_structure["6LTH"]),
        "shared_ARID1A_residues_with_Q": int(len(shared)),
        "median_Q_6LTJ_shared": float(shared.q_score_6LTJ.median()),
        "median_Q_6LTH_shared": float(shared.q_score_6LTH.median()),
        "paired_Q_wilcoxon_one_sided_P": float(paired_q.pvalue),
        "all_six_contacts_persist_4_0_to_5_0A": bool(
            contact_df[(contact_df.gene.isin(CONTACT_GENES))]
            .groupby(["pdbid", "gene"]).contact_call.min().eq(1).all()
        ),
        "all_six_contacts_persist_after_bottom10pct_Q_filter": bool(
            contact_df[(contact_df.gene.isin(CONTACT_GENES))]
            .groupby(["pdbid", "gene"]).contact_call_after_Q_filter.min().eq(1).all()
        ),
        "random_seed": SEED,
    }
    with open(os.path.join(args.out, "structure_validation_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
