"""Frozen R15 paired strong-null sensitivity analysis; run only on seeta.

Uses existing curated eligible complexes without selecting on new outcomes.
Same PatientID outer splits and symmetric anchor/drop selection for both arms.
No claim that a stable dropped member is biologically an incorrect member.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import numpy as np
import pandas as pd
import scipy
from scipy import stats

SPEC = {
    "version": "R15-strong-null-v1",
    "seed": 20260906,
    "n_splits": 10,
    "null_replicates_per_kind": 5,
    "min_complete_rows_per_half": 150,
    "strong_dependency_threshold": -1.0,
    "training_matching_features": ["mean", "sd_ddof0", "fraction_observed_below_minus1", "missing_fraction"],
    "matching_scale": "population SD across train-eligible background genes; zero scale replaced by 1",
    "matching_caliper_each_feature_sd": 0.5,
    "matching_nearest_candidates": 30,
    "matching_sampling": "greedy targets in ascending candidate-count order; random inverse-distance weighted draw among nearest available 30; no gene reuse within a null module",
    "matching_failure": "record failure; never relax caliper or resample a failed replicate",
    "pool": "all genes with >=150 observed training rows and training SD>1e-8, excluding current real-complex genes; no test statistics used",
    "split": "10 global half splits of PatientID groups; missing PatientID gets unique ModelID group; sorted group list permuted by fixed seed; same split for every complex and null",
    "missingness": "for each real-null pair separately, restrict both arms to identical complete-case rows within each frozen half; >=150 per half; no imputation; count every exclusion",
    "selection": "maximize train Spearman(anchor, mean other members), then maximize train leave-one-member-out gain; anchor always excluded; lexical gene ordering resolves ties",
    "heldout": "carry train anchor/drop into test; test best-drop search only assesses identity concordance for that fixed anchor, not test anchor selection",
    "primary_eligible": "same complex set with all 10 splits x 5 replicates x both null kinds successful; retain all failed attempts in attrition table",
    "endpoints": ["drop_identity_concordance", "heldout_selected_gain"],
    "aggregation": "pair real-minus-null; average 5 null replicates then 10 splits within complex; average complex differences within overlap family; equal-weight families primary",
    "overlap_families": "connected components of full frozen eligible catalog at gene-set Jaccard>=0.5, before observing outcomes",
    "inference": "two-sided sign-flip P with 20000 draws and +1 correction; mean-difference percentile bootstrap 10000 draws; BH across 2 endpoints x 2 null kinds x 2 units; complex-level sensitivity does not assume biological independence",
    "limitations": ["no genome-wide network-degree matching", "no novel data or independent cell-screen resource", "no tissue-stratified splits", "matching only marginal training profiles, not covariance", "overlap-family grouping may leave dependence across families", "train-dependent null members can change across splits", "five null replicates cannot provide high-resolution per-complex empirical P", "complete-case common-support subset may limit generalizability", "stability is not biological membership correctness"],
}
STATE = {}


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def stable_seed(label):
    return int.from_bytes(hashlib.sha256(f'{SPEC["seed"]}|{label}'.encode()).digest()[:8], "little") % (2**32 - 1)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")


def build_catalog(columns, manifest, cp_path, u2s_path):
    old = pd.read_csv(manifest)
    if old.groupby("source_complex").k_total.nunique().max() != 1:
        raise ValueError("R5 module sizes conflict")
    frozen = old.drop_duplicates("source_complex").set_index("source_complex")
    if len(frozen) != 277:
        raise ValueError(f"Expected 277 original eligible complexes, got {len(frozen)}")
    up = pd.read_csv(u2s_path, sep="\t")
    u2s = dict(zip(up.Entry, up["Gene Names (primary)"]))
    cp = pd.read_csv(cp_path, sep="\t", low_memory=False)
    mol = next(c for c in cp.columns if "Identifiers" in c)
    out = []
    for _, row in cp.iterrows():
        cid = str(row.iloc[0])
        if cid not in frozen.index:
            continue
        accessions = {a.split("(")[0].split("-")[0] for a in str(row[mol]).split("|") if a}
        genes = sorted({u2s[a] for a in accessions if a in u2s and isinstance(u2s[a], str)} & set(columns))
        if len(genes) != int(frozen.loc[cid, "k_total"]):
            raise ValueError(f"Frozen size mismatch for {cid}: {len(genes)}")
        out.append({"source_complex": cid, "source_name": str(frozen.loc[cid, "source_name"]), "k_total": len(genes), "genes": ";".join(genes)})
    out = pd.DataFrame(out).sort_values("source_complex").reset_index(drop=True)
    if len(out) != 277 or out.source_complex.duplicated().any():
        raise ValueError("Incomplete/duplicated frozen catalog")
    sets = [set(g.split(";")) for g in out.genes]
    parent = list(range(len(out)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    overlaps = []
    for i in range(len(sets)):
        for j in range(i):
            jac = len(sets[i] & sets[j]) / len(sets[i] | sets[j])
            if jac >= 0.5:
                parent[root(i)] = root(j)
                overlaps.append({"complex_a": out.iloc[i].source_complex, "complex_b": out.iloc[j].source_complex, "jaccard": jac})
    roots = [root(i) for i in range(len(out))]
    labels = {r: f"family_{n:03d}" for n, r in enumerate(sorted(set(roots)), 1)}
    out["overlap_family"] = [labels[r] for r in roots]
    return out, pd.DataFrame(overlaps)


def prepare(args):
    folder = Path(args.out)
    folder.mkdir(parents=True, exist_ok=True)
    if (folder / "r15_frozen_spec.json").exists():
        raise RuntimeError("Frozen specification exists; refusing overwrite")
    write_json(folder / "r15_frozen_spec.json", SPEC)
    inputs = {"dependency": args.dep, "models": args.models, "cp": str(folder / "cp_frozen.tsv"), "u2s": str(folder / "up2sym_frozen.tsv"), "r5_manifest": str(folder / "r5_null_manifest_frozen.csv"), "script": str(Path(__file__).resolve())}
    write_json(folder / "r15_input_hashes.json", {k: {"path": p, "sha256": sha(p), "bytes": Path(p).stat().st_size} for k, p in inputs.items()})
    dep = pd.read_csv(args.dep, index_col=0, low_memory=False)
    dep.columns = [str(c).split(" (")[0] for c in dep.columns]
    dup_genes = int(dep.columns.duplicated().sum())
    dep = dep.loc[:, ~dep.columns.duplicated()]
    if dep.index.duplicated().any():
        raise ValueError("Duplicate ModelID rows must be resolved explicitly")
    meta = pd.read_csv(args.models, low_memory=False).set_index("ModelID")
    if meta.index.duplicated().any():
        raise ValueError("Duplicate ModelID metadata")
    info = meta.reindex(dep.index)
    patient = info.PatientID.fillna("").astype(str)
    missing_patient = patient.eq("")
    groups = np.where(missing_patient, "MODEL:" + dep.index.astype(str), "PATIENT:" + patient)
    arr = dep.to_numpy(dtype=np.float64)
    np.save(folder / "dependency.npy", arr)
    write_json(folder / "genes.json", dep.columns.tolist())
    catalog, overlaps = build_catalog(dep.columns, inputs["r5_manifest"], inputs["cp"], inputs["u2s"])
    catalog.to_csv(folder / "r15_frozen_complexes.csv", index=False)
    overlaps.to_csv(folder / "r15_overlap_edges.csv", index=False)
    base = pd.DataFrame({"row": np.arange(len(dep)), "ModelID": dep.index, "PatientID_group": groups, "missing_PatientID": missing_patient.to_numpy(), "lineage": info.OncotreeLineage.fillna("Unknown").to_numpy()})
    uniques = np.unique(groups)
    split_rows, split_qc, features = [], [], []
    for split in range(SPEC["n_splits"]):
        rng = np.random.default_rng(stable_seed(f"outer:{split}"))
        shuffled = rng.permutation(uniques)
        train_groups = set(shuffled[:len(shuffled) // 2])
        train = np.array([g in train_groups for g in groups])
        b = base.copy()
        b["split"] = split
        b["partition"] = np.where(train, "train", "test")
        split_rows.append(b)
        x = arr[train]
        valid = np.isfinite(x)
        n = valid.sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = np.nanmean(x, axis=0)
            sd = np.nanstd(x, axis=0)
            frac = ((x < -1) & valid).sum(axis=0) / n
            missing = 1 - n / len(x)
        feat = np.column_stack([mean, sd, frac, missing, n])
        features.append(feat)
        pool = (n >= 150) & (sd > 1e-8) & np.isfinite(feat[:, :4]).all(axis=1)
        split_qc.append({"split": split, "seed": stable_seed(f"outer:{split}"), "n_train": int(train.sum()), "n_test": int((~train).sum()), "train_patients": len(train_groups), "test_patients": len(uniques) - len(train_groups), "patient_overlap": len(set(groups[train]) & set(groups[~train])), "eligible_background_genes": int(pool.sum())})
    pd.concat(split_rows).to_csv(folder / "r15_patient_splits.csv", index=False)
    pd.DataFrame(split_qc).to_csv(folder / "r15_split_qc.csv", index=False)
    np.save(folder / "training_features.npy", np.stack(features))
    preflight = {"n_models": len(dep), "n_genes": dep.shape[1], "n_unique_patient_groups": len(uniques), "n_missing_patient_ids": int(missing_patient.sum()), "duplicate_gene_symbols_first_retained": dup_genes, "eligible_complexes": len(catalog), "overlap_families": int(catalog.overlap_family.nunique()), "matrix_nonfinite_fraction": float((~np.isfinite(arr)).mean()), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__, "python": platform.python_version(), "spec_sha256": sha(folder / "r15_frozen_spec.json"), "catalog_sha256": sha(folder / "r15_frozen_complexes.csv")}
    write_json(folder / "r15_preflight.json", preflight)
    print(json.dumps(preflight, indent=2), flush=True)


def init_worker(folder):
    folder = Path(folder)
    STATE["arr"] = np.load(folder / "dependency.npy", mmap_mode="r")
    STATE["features"] = np.load(folder / "training_features.npy", mmap_mode="r")
    STATE["genes"] = json.loads((folder / "genes.json").read_text())
    STATE["gene_index"] = {g: i for i, g in enumerate(STATE["genes"])}
    splits = pd.read_csv(folder / "r15_patient_splits.csv")
    STATE["splits"] = {s: (d.loc[d.partition.eq("train"), "row"].to_numpy(), d.loc[d.partition.eq("test"), "row"].to_numpy()) for s, d in splits.groupby("split")}


def rankcorr_columns(x, y):
    if y.ndim == 1:
        y = y[:, None]
    if x.ndim == 1:
        x = x[:, None]
    rx = stats.rankdata(x, axis=0)
    ry = stats.rankdata(y, axis=0)
    rx -= rx.mean(axis=0)
    ry -= ry.mean(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return (rx * ry).sum(axis=0) / np.sqrt((rx * rx).sum(axis=0) * (ry * ry).sum(axis=0))


def evaluate(train, test, genes):
    k = train.shape[1]
    means = (train.sum(axis=1, keepdims=True) - train) / (k - 1)
    rhos = rankcorr_columns(train, means)
    if not np.isfinite(rhos).any():
        return None
    ai = int(np.nanargmax(rhos))
    members = np.array([j for j in range(k) if j != ai])
    x = train[:, members]
    loomeans = (x.sum(axis=1, keepdims=True) - x) / (k - 2)
    gains = rankcorr_columns(train[:, ai], loomeans) - rhos[ai]
    if not np.isfinite(gains).any():
        return None
    di_local = int(np.nanargmax(gains))
    di = int(members[di_local])
    y = test[:, members]
    test_base = float(rankcorr_columns(test[:, ai], y.mean(axis=1))[0])
    test_loo = (y.sum(axis=1, keepdims=True) - y) / (k - 2)
    test_gains = rankcorr_columns(test[:, ai], test_loo) - test_base
    if not np.isfinite(test_gains).all():
        return None
    bi = int(np.argmax(test_gains))
    return {"anchor": genes[ai], "selected_drop": genes[di], "test_best_drop": genes[int(members[bi])], "concordance": int(bi == di_local), "train_anchor_rho": float(rhos[ai]), "test_anchor_rho": test_base, "train_selected_gain": float(gains[di_local]), "heldout_selected_gain": float(test_gains[di_local])}


def match_null(real_idx, split, kind, rep, cid):
    feat = np.array(STATE["features"][split])
    eligible = (feat[:, 4] >= 150) & (feat[:, 1] > 1e-8) & np.isfinite(feat[:, :4]).all(axis=1)
    pool = np.flatnonzero(eligible)
    scale = np.std(feat[pool, :4], axis=0)
    scale[scale < 1e-12] = 1
    pool = np.setdiff1d(pool, real_idx)
    rng = np.random.default_rng(stable_seed(f"{cid}:{split}:{kind}:{rep}"))
    if kind == "size":
        picks = rng.choice(pool, len(real_idx), replace=False)
        return sorted(picks.tolist(), key=lambda j: STATE["genes"][j]), [], "ok"
    if not eligible[real_idx].all():
        return None, [], "real_target_not_training_eligible"
    candidates = {}
    for ri in real_idx:
        delta = (feat[pool, :4] - feat[ri, :4]) / scale
        valid = np.max(np.abs(delta), axis=1) <= SPEC["matching_caliper_each_feature_sd"]
        inds = pool[valid]
        distances = np.sum(delta[valid] ** 2, axis=1)
        order = np.lexsort((inds, distances))
        candidates[ri] = (inds[order], distances[order])
    target_order = sorted(real_idx, key=lambda ri: (len(candidates[ri][0]), STATE["genes"][ri]))
    used, balances = set(), []
    for ri in target_order:
        ci, ds = candidates[ri]
        keep = np.array([g not in used for g in ci], dtype=bool)
        ci, ds = ci[keep][:30], ds[keep][:30]
        if len(ci) == 0:
            return None, balances, f"no_caliper_candidate:{STATE['genes'][ri]}"
        weights = 1 / (np.sqrt(ds) + 0.05)
        ni = int(rng.choice(ci, p=weights / weights.sum()))
        used.add(ni)
        row = {"source_complex": cid, "split": split, "null_rep": rep, "real_gene": STATE["genes"][ri], "matched_gene": STATE["genes"][ni], "n_caliper_candidates": len(candidates[ri][0])}
        for fi, fn in enumerate(SPEC["training_matching_features"]):
            row[f"real_{fn}"] = float(feat[ri, fi])
            row[f"matched_{fn}"] = float(feat[ni, fi])
            row[f"standardized_delta_{fn}"] = float((feat[ni, fi] - feat[ri, fi]) / scale[fi])
        balances.append(row)
    return sorted(used, key=lambda j: STATE["genes"][j]), balances, "ok"


def analyze_complex(record):
    cid = record["source_complex"]
    real_genes = record["genes"].split(";")
    real_idx = [STATE["gene_index"][g] for g in real_genes]
    arr = STATE["arr"]
    rows, matching, manifests = [], [], []
    for split in range(10):
        tr, te = STATE["splits"][split]
        for kind in ("strong", "size"):
            for rep in range(1, 6):
                base = {"source_complex": cid, "source_name": record["source_name"], "overlap_family": record["overlap_family"], "split": split, "kind": kind, "null_rep": rep, "k_total": len(real_idx), "train_initial": len(tr), "test_initial": len(te)}
                ni, balance, status = match_null(real_idx, split, kind, rep, cid)
                matching.extend(balance)
                if ni is None:
                    rows.append(dict(base, status=status))
                    continue
                null_genes = [STATE["genes"][j] for j in ni]
                manifests.append(dict(base, genes=";".join(null_genes), seed=stable_seed(f"{cid}:{split}:{kind}:{rep}")))
                combined = arr[:, real_idx + ni]
                complete = np.isfinite(combined).all(axis=1)
                trc, tec = tr[complete[tr]], te[complete[te]]
                base.update(train_complete=len(trc), test_complete=len(tec), train_excluded=len(tr)-len(trc), test_excluded=len(te)-len(tec))
                if len(trc) < 150 or len(tec) < 150:
                    rows.append(dict(base, status="insufficient_paired_complete_cases"))
                    continue
                real = evaluate(arr[np.ix_(trc, real_idx)], arr[np.ix_(tec, real_idx)], real_genes)
                null = evaluate(arr[np.ix_(trc, ni)], arr[np.ix_(tec, ni)], null_genes)
                if real is None or null is None:
                    rows.append(dict(base, status="undefined_rank_correlation"))
                    continue
                result = dict(base, status="ok")
                result.update({f"real_{k}": v for k, v in real.items()})
                result.update({f"null_{k}": v for k, v in null.items()})
                result["concordance_difference"] = real["concordance"] - null["concordance"]
                result["heldout_gain_difference"] = real["heldout_selected_gain"] - null["heldout_selected_gain"]
                rows.append(result)
    return rows, matching, manifests


def inference(values, label):
    values = np.asarray(values, dtype=float)
    n = len(values)
    if not n:
        return {"n": 0, "mean_difference": None, "median_difference": None, "ci_low": None, "ci_high": None, "p_two_sided": None}
    rng = np.random.default_rng(stable_seed(label))
    boot = np.mean(values[rng.integers(0, n, size=(10000, n))], axis=1)
    signs = rng.choice([-1, 1], size=(20000, n))
    perm = (signs * values).mean(axis=1)
    estimate = float(values.mean())
    p = (1 + int((np.abs(perm) >= abs(estimate) - 1e-15).sum())) / 20001
    lo, hi = np.quantile(boot, [0.025, 0.975])
    return {"n": n, "mean_difference": estimate, "median_difference": float(np.median(values)), "ci_low": float(lo), "ci_high": float(hi), "p_two_sided": float(p)}


def summarize(folder):
    folder = Path(folder)
    results = pd.read_csv(folder / "r15_pair_results.csv")
    ok = results.loc[results.status.eq("ok")].copy()
    counts = ok.groupby(["source_complex", "kind"]).size().unstack(fill_value=0)
    for k in ("strong", "size"):
        if k not in counts:
            counts[k] = 0
    keep = counts.index[counts.strong.eq(50) & counts["size"].eq(50)]
    ok["primary_complete_complex"] = ok.source_complex.isin(keep)
    numeric = ["concordance_difference", "heldout_gain_difference", "real_concordance", "null_concordance", "real_heldout_selected_gain", "null_heldout_selected_gain"]
    # Each null kind has identical numbers of reps/splits in the primary set.
    primary = ok.loc[ok.primary_complete_complex].groupby(["source_complex", "overlap_family", "kind"])[numeric].mean().reset_index()
    primary.to_csv(folder / "r15_primary_complex_results.csv", index=False)
    families = primary.groupby(["overlap_family", "kind"])[numeric].mean().reset_index()
    families.to_csv(folder / "r15_primary_family_results.csv", index=False)
    rows = []
    for unit, table in (("family_primary", families), ("complex_sensitivity", primary)):
        for kind in ("strong", "size"):
            for endpoint in ("concordance_difference", "heldout_gain_difference"):
                r = inference(table.loc[table.kind.eq(kind), endpoint], f"infer:{unit}:{kind}:{endpoint}")
                rows.append(dict(unit=unit, kind=kind, endpoint=endpoint, **r))
    tests = pd.DataFrame(rows)
    valid = tests.p_two_sided.notna()
    tests.loc[valid, "BH_q_8_tests"] = stats.false_discovery_control(tests.loc[valid, "p_two_sided"].to_numpy(), method="bh")
    tests.to_csv(folder / "r15_primary_inference.csv", index=False)
    results.groupby(["kind", "status"]).size().rename("n_attempts").reset_index().to_csv(folder / "r15_attrition_summary.csv", index=False)
    results.loc[~results.status.eq("ok")].to_csv(folder / "r15_failed_attempts.csv", index=False)
    # Descriptive only: all available attempts, never substitutes primary results.
    ok.groupby(["source_complex", "overlap_family", "kind"])[numeric].agg(["mean", "count"]).to_csv(folder / "r15_available_case_descriptive.csv")
    balance = pd.read_csv(folder / "r15_gene_matching_balance.csv")
    bcols = [c for c in balance if c.startswith("standardized_delta_")]
    balance_summary = []
    for c in bcols:
        b = balance[c].to_numpy()
        balance_summary.append({"feature": c, "n_gene_matches": len(b), "mean_standardized_delta": float(b.mean()), "mean_abs_standardized_delta": float(abs(b).mean()), "p95_abs_standardized_delta": float(np.quantile(abs(b), .95)), "max_abs_standardized_delta": float(abs(b).max()), "n_caliper_violations": int((abs(b) > .5 + 1e-10).sum())})
    pd.DataFrame(balance_summary).to_csv(folder / "r15_matching_balance_summary.csv", index=False)
    summary = {"eligible_complexes": int(results.source_complex.nunique()), "attempts": len(results), "successful_attempts": len(ok), "primary_complete_complexes": len(keep), "primary_families": int(primary.overlap_family.nunique()), "full_frozen_family_count": int(pd.read_csv(folder / "r15_frozen_complexes.csv").overlap_family.nunique()), "primary_set_rule": SPEC["primary_eligible"], "matching_caliper_violations": sum(x["n_caliper_violations"] for x in balance_summary), "spec_sha256": sha(folder / "r15_frozen_spec.json"), "script_sha256": sha(Path(__file__)), "limitations": SPEC["limitations"]}
    write_json(folder / "r15_summary.json", summary)
    print(json.dumps(summary, indent=2), flush=True)
    print(tests.to_string(index=False), flush=True)


def run(args):
    folder = Path(args.out)
    if json.loads((folder / "r15_frozen_spec.json").read_text()) != SPEC:
        raise RuntimeError("Specification differs from frozen version")
    frozen_hash = json.loads((folder / "r15_input_hashes.json").read_text())["script"]["sha256"]
    if frozen_hash != sha(Path(__file__)):
        raise RuntimeError("Code changed after preregistered freeze")
    if (folder / "r15_pair_results.csv").exists():
        raise RuntimeError("Results exist; refusing overwrite")
    catalog = pd.read_csv(folder / "r15_frozen_complexes.csv")
    records = catalog.to_dict("records")
    start = time.monotonic()
    rows, matching, manifests = [], [], []
    with ProcessPoolExecutor(max_workers=min(args.workers, 8), initializer=init_worker, initargs=(str(folder),)) as pool:
        for i, (a, b, c) in enumerate(pool.map(analyze_complex, records, chunksize=1), 1):
            rows.extend(a)
            matching.extend(b)
            manifests.extend(c)
            if i % 10 == 0 or i == len(records):
                elapsed = time.monotonic() - start
                print(f"completed {i}/{len(records)} complexes; {elapsed:.1f}s; ETA {elapsed/i*(len(records)-i):.1f}s", flush=True)
    pd.DataFrame(rows).to_csv(folder / "r15_pair_results.csv", index=False)
    pd.DataFrame(matching).to_csv(folder / "r15_gene_matching_balance.csv", index=False)
    pd.DataFrame(manifests).to_csv(folder / "r15_null_gene_manifest.csv", index=False)
    summarize(folder)
    write_json(folder / "r15_output_hashes.json", {p.name: {"sha256": sha(p), "bytes": p.stat().st_size} for p in sorted(folder.glob("r15_*")) if p.is_file() and p.name != "r15_output_hashes.json"})
    print(f"FINISHED in {time.monotonic()-start:.1f}s", flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--stage", choices=["prepare", "run"], required=True)
    p.add_argument("--out", default="/root/autodl-tmp/arid1a_baf_r15/strong_null")
    p.add_argument("--dep", default="/root/autodl-tmp/arid1a_baf_r12_cellline/inputs/depmap_CRISPRGeneEffect.csv")
    p.add_argument("--models", default="/root/autodl-tmp/arid1a_baf_r12_cellline/inputs/depmap_Model.csv")
    p.add_argument("--workers", type=int, default=8)
    args = p.parse_args()
    (prepare if args.stage == "prepare" else run)(args)


if __name__ == "__main__":
    main()
