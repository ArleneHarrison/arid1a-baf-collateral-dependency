# -*- coding: utf-8 -*-
"""Nested, selection-symmetric validation of complex-module refinement.

The earlier 227-complex survey selected the highest-correlating anchor for each
real complex, whereas its size-matched random null used an arbitrary anchor.
This script removes that asymmetry.  In every outer split, BOTH a curated
complex and each size-matched random module:

  1. select the best anchor using discovery-half cell lines only;
  2. nominate the member whose removal most improves the discovery correlation;
  3. carry that anchor and nomination into the held-out half;
  4. record whether the same member is the held-out best drop and whether its
     held-out gain is positive.

The primary unit is the curated complex.  Each curated complex is paired with
two independently sampled, size-matched random modules.  Inference compares a
complex's validation concordance with the mean concordance of its two matched
null modules.  No result from this script is confirmatory; it is a robustness
analysis of the module-refinement procedure.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from scipy import stats


DEP: pd.DataFrame | None = None
N_SPLITS = 30
MIN_HALF = 150
BASE_SEED = 20260901


def stable_seed(label: str) -> int:
    digest = hashlib.sha256(f"{BASE_SEED}|{label}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "little") % (2**32 - 1)


def spearman_fast(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 3 or np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return np.nan
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    return float(np.corrcoef(rx, ry)[0, 1])


def select_anchor_and_drop(arr: np.ndarray, genes: list[str], idx: np.ndarray):
    """Select best anchor, then best leave-one-out drop, on one data half."""
    block = arr[idx]
    best_anchor_idx = None
    best_anchor_rho = -np.inf
    for ai in range(block.shape[1]):
        module_idx = [j for j in range(block.shape[1]) if j != ai]
        rho = spearman_fast(block[:, ai], block[:, module_idx].mean(axis=1))
        if np.isfinite(rho) and rho > best_anchor_rho:
            best_anchor_idx = ai
            best_anchor_rho = rho
    if best_anchor_idx is None:
        return None

    module_idx = [j for j in range(block.shape[1]) if j != best_anchor_idx]
    anchor = block[:, best_anchor_idx]
    base = spearman_fast(anchor, block[:, module_idx].mean(axis=1))
    best_drop_idx = None
    best_gain = -np.inf
    for di in module_idx:
        rest = [j for j in module_idx if j != di]
        rho = spearman_fast(anchor, block[:, rest].mean(axis=1))
        gain = rho - base
        if np.isfinite(gain) and gain > best_gain:
            best_drop_idx = di
            best_gain = gain
    if best_drop_idx is None:
        return None
    return {
        "anchor_idx": best_anchor_idx,
        "anchor": genes[best_anchor_idx],
        "anchor_rho": best_anchor_rho,
        "drop_idx": best_drop_idx,
        "drop": genes[best_drop_idx],
        "gain": best_gain,
    }


def validate_selected(arr: np.ndarray, genes: list[str], idx: np.ndarray,
                      anchor_idx: int, selected_drop_idx: int):
    """Evaluate the discovery-selected anchor/drop on held-out lines."""
    block = arr[idx]
    module_idx = [j for j in range(block.shape[1]) if j != anchor_idx]
    anchor = block[:, anchor_idx]
    base = spearman_fast(anchor, block[:, module_idx].mean(axis=1))
    gains = {}
    for di in module_idx:
        rest = [j for j in module_idx if j != di]
        rho = spearman_fast(anchor, block[:, rest].mean(axis=1))
        gains[di] = rho - base
    best_drop_idx = max(gains, key=gains.get)
    return {
        "anchor_rho": base,
        "best_drop_idx": best_drop_idx,
        "best_drop": genes[best_drop_idx],
        "selected_gain": gains[selected_drop_idx],
        "best_gain": gains[best_drop_idx],
    }


def analyze_one(task):
    global DEP
    module_id, source_complex, kind, null_rep, genes = task
    frame = DEP[genes].dropna(axis=0, how="any")
    if len(frame) < 2 * MIN_HALF:
        return [], None
    arr = frame.to_numpy(dtype=float)
    rng = np.random.default_rng(stable_seed(module_id))
    rows = []
    for split in range(N_SPLITS):
        perm = rng.permutation(len(arr))
        cut = len(perm) // 2
        discovery, validation = perm[:cut], perm[cut:]
        if len(discovery) < MIN_HALF or len(validation) < MIN_HALF:
            continue
        selected = select_anchor_and_drop(arr, genes, discovery)
        if selected is None:
            continue
        held = validate_selected(
            arr, genes, validation,
            selected["anchor_idx"], selected["drop_idx"]
        )
        rows.append({
            "module_id": module_id,
            "source_complex": source_complex,
            "kind": kind,
            "null_rep": null_rep,
            "k_total": len(genes),
            "n_complete": len(frame),
            "split": split,
            "selected_anchor": selected["anchor"],
            "selected_drop": selected["drop"],
            "validation_best_drop": held["best_drop"],
            "same_drop": int(selected["drop_idx"] == held["best_drop_idx"]),
            "train_anchor_rho": selected["anchor_rho"],
            "validation_anchor_rho": held["anchor_rho"],
            "train_selected_gain": selected["gain"],
            "validation_selected_gain": held["selected_gain"],
            "validation_gain_positive": int(held["selected_gain"] > 0),
        })
    if not rows:
        return [], None
    d = pd.DataFrame(rows)
    summary = {
        "module_id": module_id,
        "source_complex": source_complex,
        "kind": kind,
        "null_rep": null_rep,
        "k_total": len(genes),
        "n_complete": len(frame),
        "n_splits": len(d),
        "concordance": float(d.same_drop.mean()),
        "positive_gain_fraction": float(d.validation_gain_positive.mean()),
        "median_validation_gain": float(d.validation_selected_gain.median()),
        "median_train_anchor_rho": float(d.train_anchor_rho.median()),
        "median_validation_anchor_rho": float(d.validation_anchor_rho.median()),
    }
    return rows, summary


def cliff_delta(x: np.ndarray, y: np.ndarray) -> float:
    return float((np.greater.outer(x, y).sum() - np.less.outer(x, y).sum()) /
                 (len(x) * len(y)))


def bootstrap_paired_median_diff(diff: np.ndarray, n_boot=20000):
    rng = np.random.default_rng(BASE_SEED)
    n = len(diff)
    vals = np.empty(n_boot)
    for i in range(n_boot):
        vals[i] = np.median(diff[rng.integers(0, n, size=n)])
    return [float(v) for v in np.quantile(vals, [0.025, 0.975])]


def build_complexes(dep_columns: set[str], cp_path: str, u2s_path: str):
    u2s_df = pd.read_csv(u2s_path, sep="\t")
    u2s = dict(zip(u2s_df["Entry"], u2s_df["Gene Names (primary)"]))
    cp = pd.read_csv(cp_path, sep="\t", low_memory=False)
    mol_col = next(c for c in cp.columns if "Identifiers" in c)
    ac_col = cp.columns[0]
    name_col = "Recommended name"
    records = []
    for _, row in cp.iterrows():
        raw = row[mol_col]
        if not isinstance(raw, str):
            continue
        accessions = {a.split("(")[0].split("-")[0] for a in raw.split("|") if a}
        genes = sorted({u2s[a] for a in accessions
                        if a in u2s and isinstance(u2s.get(a), str)} & dep_columns)
        # Match the earlier resampling range: anchor + 5--20 scored members.
        if 6 <= len(genes) <= 21:
            records.append((str(row[ac_col]), str(row[name_col]), genes))
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bulk", required=True,
                        help="Directory containing depmap_CRISPRGeneEffect.csv")
    parser.add_argument("--cp", required=True,
                        help="Complex Portal complex-membership TSV")
    parser.add_argument("--u2s", required=True,
                        help="UniProt-accession-to-gene-symbol TSV")
    parser.add_argument("--out", required=True)
    parser.add_argument("--splits", type=int, default=30)
    parser.add_argument("--null-reps", type=int, default=2)
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()

    global DEP, N_SPLITS
    N_SPLITS = args.splits
    os.makedirs(args.out, exist_ok=True)
    dep_path = os.path.join(args.bulk, "depmap_CRISPRGeneEffect.csv")
    DEP = pd.read_csv(dep_path, index_col=0, low_memory=False)
    DEP.columns = [c.split(" (")[0] for c in DEP.columns]
    DEP = DEP.loc[:, ~DEP.columns.duplicated()]

    complexes = build_complexes(set(DEP.columns), args.cp, args.u2s)
    all_complex_members = set(g for _, _, genes in complexes for g in genes)
    pool = [g for g in DEP.columns if DEP[g].notna().sum() >= 2 * MIN_HALF]
    tasks = []
    null_manifest = []
    for ac, name, genes in complexes:
        tasks.append((f"real:{ac}", ac, "real", 0, genes))
        for rep in range(1, args.null_reps + 1):
            rng = np.random.default_rng(stable_seed(f"null:{ac}:{rep}"))
            # Sampling from the full well-covered genome is conservative: a
            # chance draw may include genuine co-complex genes.
            pick = sorted(rng.choice(pool, size=len(genes), replace=False).tolist())
            tasks.append((f"null:{ac}:{rep}", ac, "null", rep, pick))
            null_manifest.append({
                "source_complex": ac, "source_name": name, "null_rep": rep,
                "k_total": len(genes), "genes": ";".join(pick),
                "n_overlap_any_curated_member": len(set(pick) & all_complex_members),
            })

    print(f"DepMap {DEP.shape}; {len(complexes)} curated complexes; "
          f"{len(tasks)} total modules; {N_SPLITS} nested splits", flush=True)
    split_rows, summaries = [], []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futures = {ex.submit(analyze_one, task): task[0] for task in tasks}
        for i, fut in enumerate(as_completed(futures), 1):
            rows, summary = fut.result()
            split_rows.extend(rows)
            if summary is not None:
                summaries.append(summary)
            if i % 50 == 0 or i == len(futures):
                print(f"  completed {i}/{len(futures)}", flush=True)

    split_df = pd.DataFrame(split_rows)
    module_df = pd.DataFrame(summaries)
    split_df.to_csv(os.path.join(args.out, "nested_split_results.csv"), index=False)
    module_df.to_csv(os.path.join(args.out, "nested_module_summary.csv"), index=False)
    pd.DataFrame(null_manifest).to_csv(
        os.path.join(args.out, "nested_null_manifest.csv"), index=False
    )

    real = module_df[module_df.kind == "real"].set_index("source_complex")
    null = (module_df[module_df.kind == "null"]
            .groupby("source_complex")
            .agg(null_concordance=("concordance", "mean"),
                 null_positive_gain_fraction=("positive_gain_fraction", "mean"),
                 null_median_validation_gain=("median_validation_gain", "mean")))
    paired = real.join(null, how="inner")
    paired["concordance_difference"] = (
        paired.concordance - paired.null_concordance
    )
    paired["positive_gain_fraction_difference"] = (
        paired.positive_gain_fraction - paired.null_positive_gain_fraction
    )
    paired.to_csv(os.path.join(args.out, "nested_paired_comparison.csv"))

    diff = paired.concordance_difference.to_numpy()
    wil = stats.wilcoxon(diff, alternative="greater", zero_method="wilcox")
    real_vals = paired.concordance.to_numpy()
    null_vals = paired.null_concordance.to_numpy()
    summary = {
        "design": "nested discovery-half anchor/drop selection with held-out validation",
        "n_curated_complexes_available": len(complexes),
        "n_paired_complexes_analyzed": int(len(paired)),
        "n_splits_per_module": N_SPLITS,
        "n_size_matched_null_modules_per_complex": args.null_reps,
        "real_median_concordance": float(np.median(real_vals)),
        "matched_null_median_concordance": float(np.median(null_vals)),
        "paired_median_difference": float(np.median(diff)),
        "paired_median_difference_bootstrap_95ci": bootstrap_paired_median_diff(diff),
        "paired_wilcoxon_statistic": float(wil.statistic),
        "paired_wilcoxon_one_sided_P": float(wil.pvalue),
        "cliffs_delta_real_vs_paired_null": cliff_delta(real_vals, null_vals),
        "real_median_positive_gain_fraction": float(
            paired.positive_gain_fraction.median()),
        "matched_null_median_positive_gain_fraction": float(
            paired.null_positive_gain_fraction.median()),
        "random_seed": BASE_SEED,
    }
    with open(os.path.join(args.out, "nested_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
