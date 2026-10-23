#!/usr/bin/env python3
"""Run frozen CPTAC purity-primary and base-sensitivity analyses on the server."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
import argparse

ROOT = Path("/root/autodl-tmp/arid1a_baf_r12_cptac")
OUT = ROOT / "analysis_results"
PYTHON = "/root/autodl-tmp/arid1a_baf_virtual_ko_v2/.venv/bin/python"


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def stamp():
    return datetime.now(timezone.utc).isoformat()


def main():
    global OUT
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', default='analysis_results')
    args = parser.parse_args()
    OUT = ROOT / args.outdir
    if any((OUT / name).exists() for name in ['frozen', 'LUAD_purity', 'LUSC_purity', 'LUAD_base', 'LUSC_base']):
        raise SystemExit(f"Refusing to overwrite existing formal results: {OUT}")
    OUT.mkdir(exist_ok=True)
    frozen = OUT / "frozen"
    frozen.mkdir()
    for name in ["r12_protein_multiomics.py", "r12_analysis_spec.json"]:
        shutil.copy2(ROOT / name, frozen / name)
    env = os.environ.copy()
    for key in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
        env[key] = "1"
    script, spec = frozen / "r12_protein_multiomics.py", frozen / "r12_analysis_spec.json"
    assert json.loads(spec.read_text())["bootstrap_replicates"] == 10000
    with (OUT / "self_test.log").open("w") as log:
        result = subprocess.run([PYTHON, "-u", str(script), "--self-test"],
                                stdout=log, stderr=subprocess.STDOUT, env=env)
    if result.returncode != 0:
        raise SystemExit("Frozen script self-test failed")
    print("FROZEN_SELF_TEST_PASS", flush=True)
    report = {"started_utc": stamp(), "server_python": PYTHON,
              "script_sha256": sha(script), "spec_sha256": sha(spec),
              "process_limit": 2, "threads_per_process": 1,
              "design": "Purity model is primary; RNA/CNV-only base is predeclared sensitivity",
              "bootstrap_replicates": 10000,
              "multiplicity": "Final BH across four external primary tests will be applied by root",
              "runs": []}

    def run_one(cohort, model):
        start = time.monotonic()
        args = [PYTHON, "-u", str(script), "--input", str(ROOT / "pancan_release/analysis_ready" / cohort),
                "--outdir", str(OUT / f"{cohort}_{model}"), "--spec", str(spec),
                "--cohort", f"CPTAC_{cohort}", "--deduplicate-patient"]
        if model == "purity":
            args += ["--numeric", "purity"]
        print(f"START {cohort} {model}", flush=True)
        with (OUT / f"{cohort}_{model}.log").open("w") as log:
            proc = subprocess.run(args, stdout=log, stderr=subprocess.STDOUT, env=env)
        item = {"cohort": cohort, "model": model, "exit_code": proc.returncode,
                "elapsed_seconds": time.monotonic() - start, "command": args,
                "completed_utc": stamp()}
        print(f"DONE {cohort} {model} code={proc.returncode}", flush=True)
        return item

    for model in ["purity", "base"]:
        with ThreadPoolExecutor(max_workers=2) as executor:
            jobs = [executor.submit(run_one, c, model) for c in ["LUAD", "LUSC"]]
            results = [j.result() for j in jobs]
        report["runs"].extend(results)
        (OUT / "run_report.json").write_text(json.dumps(report, indent=2))
        if any(r["exit_code"] != 0 for r in results):
            raise SystemExit("Formal analysis failed; retain logs without replacing results")
    report["finished_utc"] = stamp()
    report["artifact_sha256"] = {str(p.relative_to(OUT)): sha(p) for p in OUT.rglob("*")
                                  if p.is_file() and p.name != "run_report.json"}
    (OUT / "run_report.json").write_text(json.dumps(report, indent=2))
    print("ALL_FOUR_FORMAL_ANALYSES_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
