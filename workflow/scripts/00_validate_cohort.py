#!/usr/bin/env python3
"""Audit every frozen figure-specific cohort before a pipeline run."""
from __future__ import annotations
import hashlib, json, os, sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "config"))
from settings import QC1_INPUT_CSV, QC2_INPUT_CSV, QC3_INPUT_CSV, QC4_INPUT_CSV, OUTPUT_ROOT

OUTDIR = Path(os.environ.get("MICROBIOME_OUTPUT", OUTPUT_ROOT)) / "00_cohort"
COHORTS = {"qc1": QC1_INPUT_CSV, "qc2": QC2_INPUT_CSV, "qc3": QC3_INPUT_CSV, "qc4": QC4_INPUT_CSV}

def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    records = []
    for name, path in COHORTS.items():
        df = pd.read_csv(path, low_memory=False)
        if df.sample_id.duplicated().any(): raise ValueError(f"{name}: duplicate sample IDs")
        records.append({"cohort": name, "input": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                        "samples": len(df), "unique_sample_ids": int(df.sample_id.nunique()),
                        "taxa_features": int(len(df.columns) - 3), "age_nonmissing": int(df.age.notna().sum()),
                        "weight_nonmissing": int(df.weight.notna().sum()), "age_ge_1": int((df.age >= 1).sum()),
                        "age_ge_1_weight_nonmissing": int(((df.age >= 1) & df.weight.notna()).sum())})
        df[["sample_id", "age", "weight"]].to_csv(OUTDIR / f"{name}_sample_manifest.csv", index=False)
    pd.DataFrame(records).to_csv(OUTDIR / "cohort_audit.csv", index=False)
    (OUTDIR / "cohort_audit.json").write_text(json.dumps(records, indent=2) + "\n")
    print(pd.DataFrame(records).to_string(index=False))

if __name__ == "__main__": main()
