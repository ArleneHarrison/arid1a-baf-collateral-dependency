"""Independent scalar-Spearman replay and integrity checks for frozen R15 nulls."""
import hashlib
import json
import os
from pathlib import Path

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path("/root/autodl-tmp/arid1a_baf_r15/strong_null")


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(8 << 20), b""):
            h.update(b)
    return h.hexdigest()


def independent(train, test, genes):
    correlations = []
    for i in range(len(genes)):
        others = [j for j in range(len(genes)) if j != i]
        correlations.append(spearmanr(train[:, i], train[:, others].mean(axis=1)).statistic)
    ai = int(np.nanargmax(correlations))
    others = [j for j in range(len(genes)) if j != ai]
    gains = []
    for di in others:
        rest = [j for j in others if j != di]
        gains.append(spearmanr(train[:, ai], train[:, rest].mean(axis=1)).statistic-correlations[ai])
    di_local = int(np.nanargmax(gains))
    base_test = spearmanr(test[:, ai], test[:, others].mean(axis=1)).statistic
    held = []
    for di in others:
        rest = [j for j in others if j != di]
        held.append(spearmanr(test[:, ai], test[:, rest].mean(axis=1)).statistic-base_test)
    bi = int(np.nanargmax(held))
    return {"anchor": genes[ai], "selected_drop": genes[others[di_local]], "test_best_drop": genes[others[bi]], "concordance": int(di_local==bi), "train_anchor_rho": correlations[ai], "test_anchor_rho": base_test, "train_selected_gain": gains[di_local], "heldout_selected_gain": held[di_local]}


def main():
    results = pd.read_csv(ROOT/"r15_pair_results.csv")
    nulls = pd.read_csv(ROOT/"r15_null_gene_manifest.csv")
    catalog = pd.read_csv(ROOT/"r15_frozen_complexes.csv").set_index("source_complex")
    splits = pd.read_csv(ROOT/"r15_patient_splits.csv")
    matrix = np.load(ROOT/"dependency.npy", mmap_mode="r")
    genes = json.loads((ROOT/"genes.json").read_text())
    gi = {g:i for i,g in enumerate(genes)}
    checks = {}
    checks["patient_group_overlap_zero"] = all(not (set(d.loc[d.partition.eq("train"), "PatientID_group"]) & set(d.loc[d.partition.eq("test"), "PatientID_group"])) for _,d in splits.groupby("split"))
    checks["ten_splits_each_model_once"] = len(splits)==10*len(matrix) and not splits.duplicated(["split","ModelID"]).any()
    checks["all_277_times_10_times_2_times_5_attempted"] = len(results)==27700 and not results.duplicated(["source_complex","split","kind","null_rep"]).any()
    checks["no_null_duplicate_or_real_overlap"] = all(len(set(r.genes.split(";")))==r.k_total and not (set(r.genes.split(";")) & set(catalog.loc[r.source_complex,"genes"].split(";"))) for r in nulls.itertuples())
    balance = pd.read_csv(ROOT/"r15_gene_matching_balance.csv")
    cols = [c for c in balance if c.startswith("standardized_delta_")]
    checks["strong_caliper_all_four_features"] = bool((balance[cols].abs()<=.5+1e-10).all().all())
    pairs = results.loc[results.status.eq("ok")].sort_values(["source_complex","split","kind","null_rep"])
    # Deterministic systematic sample spans the full ordered catalog and both arms.
    take = pairs.iloc[np.linspace(0,len(pairs)-1,min(20,len(pairs))).astype(int)]
    nulls = nulls.set_index(["source_complex","split","kind","null_rep"])
    replay = []
    for r in take.itertuples():
        rg = catalog.loc[r.source_complex,"genes"].split(";")
        ng = nulls.loc[(r.source_complex,r.split,r.kind,r.null_rep),"genes"].split(";")
        ri,ni = [gi[g] for g in rg],[gi[g] for g in ng]
        block = matrix[:,ri+ni]
        valid = np.isfinite(block).all(axis=1)
        sd = splits.loc[splits.split.eq(r.split)]
        tr = sd.loc[sd.partition.eq("train"),"row"].to_numpy()
        te = sd.loc[sd.partition.eq("test"),"row"].to_numpy()
        tr,te = tr[valid[tr]],te[valid[te]]
        for prefix,idx,names in (("real",ri,rg),("null",ni,ng)):
            calc = independent(matrix[np.ix_(tr,idx)],matrix[np.ix_(te,idx)],names)
            for key,value in calc.items():
                reported = getattr(r,prefix+"_"+key)
                equal = (value==reported) if isinstance(value,str) else bool(np.isclose(value,reported,atol=1e-12,rtol=1e-12))
                replay.append({"source_complex":r.source_complex,"split":r.split,"kind":r.kind,"null_rep":r.null_rep,"arm":prefix,"metric":key,"reported":reported,"independent":value,"matches":equal})
        checks[f"paired_row_counts:{r.source_complex}:{r.split}:{r.kind}:{r.null_rep}"] = len(tr)==r.train_complete and len(te)==r.test_complete
    replay = pd.DataFrame(replay)
    replay.to_csv(ROOT/"r15_independent_replay.csv",index=False)
    checks["independent_scalar_spearman_all_matches"] = bool(replay.matches.all())
    hashes = json.loads((ROOT/"r15_input_hashes.json").read_text())
    checks["all_original_input_and_code_hashes_unchanged"] = all(sha(v["path"])==v["sha256"] for v in hashes.values())
    primary = pd.read_csv(ROOT/"r15_primary_complex_results.csv")
    full = set(results.loc[results.status.eq("ok")].groupby("source_complex").size().loc[lambda x:x.eq(100)].index)
    checks["primary_set_exact_fully_successful_catalog"] = set(primary.source_complex)==full
    audit = {"status":"PASS" if all(checks.values()) else "FAIL", "checks":checks,"replayed_pairs":len(take),"independent_metric_checks":len(replay),"audit_script_sha256":sha(Path(__file__))}
    (ROOT/"r15_independent_audit.json").write_text(json.dumps(audit,indent=2),encoding="utf-8")
    manifest = {p.name:{"sha256":sha(p),"bytes":p.stat().st_size} for p in sorted(ROOT.iterdir()) if p.is_file() and p.suffix not in (".npy",) and p.name!="r15_release_hashes.json"}
    (ROOT/"r15_release_hashes.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    print(json.dumps(audit,indent=2),flush=True)
    if audit["status"]!="PASS":
        raise RuntimeError("Independent audit failed")


if __name__=="__main__":
    main()
