"""Summarise ARID1A mutation frequency and variant classes from TCGA MAF files."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from common import ensure_parent


def labelled(values: list[str]) -> dict[str, str]:
    result = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"Expected LABEL=VALUE, received {value!r}")
        label, item = value.split("=", 1)
        result[label] = item
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maf", action="append", required=True, help="COHORT=path/to/file.maf[.gz]")
    parser.add_argument("--denominator", action="append", required=True, help="COHORT=number_sequenced")
    parser.add_argument("--gene", default="ARID1A")
    parser.add_argument("--frequency-output", required=True)
    parser.add_argument("--spectrum-output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    mafs = labelled(args.maf)
    denominators = {key: int(value) for key, value in labelled(args.denominator).items()}
    frequency_rows = []
    spectrum_frames = []
    for cohort, path in mafs.items():
        maf = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
        required = {"Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification"}
        missing = required - set(maf.columns)
        if missing:
            raise ValueError(f"{path} lacks columns: {sorted(missing)}")
        variants = maf.loc[maf.Hugo_Symbol.astype(str).str.casefold() == args.gene.casefold()].copy()
        n_mutant = variants.Tumor_Sample_Barcode.astype(str).nunique()
        denominator = denominators[cohort]
        frequency_rows.append({
            "cohort": cohort, "n_mut_patients": n_mutant,
            "n_total": denominator, "freq_pct": 100 * n_mutant / denominator,
        })
        spectrum_frames.append(
            variants.groupby("Variant_Classification", as_index=False).size().rename(columns={"size": "N"}).assign(cohort=cohort)
        )
    frequency = pd.DataFrame(frequency_rows)
    if len(frequency) > 1:
        n_mutant = frequency.n_mut_patients.sum()
        denominator = frequency.n_total.sum()
        frequency.loc[len(frequency)] = ["NSCLC", n_mutant, denominator, 100 * n_mutant / denominator]
    frequency.to_csv(ensure_parent(args.frequency_output), index=False)
    pd.concat(spectrum_frames, ignore_index=True).to_csv(ensure_parent(args.spectrum_output), index=False)
    print(frequency.to_string(index=False))


if __name__ == "__main__":
    main()

