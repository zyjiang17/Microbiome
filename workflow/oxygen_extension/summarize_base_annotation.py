#!/usr/bin/env python3
"""Summarize coverage of the evidence-based base oxygen annotation."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent
queue = pd.read_csv(ROOT / "outputs/species_annotation_review_queue.csv")
ann = pd.read_csv(ROOT / "base_oxygen_annotation_v1.csv")
if ann.species.duplicated().any():
    raise ValueError("Duplicate species in base annotation")
merged = queue.merge(ann, on="species", how="left", validate="one_to_one")
merged["oxygen_phenotype"] = merged.oxygen_phenotype.fillna("unreviewed")
merged["analysis_group"] = merged.analysis_group.fillna("excluded")
merged["evidence_level"] = merged.evidence_level.fillna("none")
merged.to_csv(ROOT / "outputs/base_oxygen_annotation_all_species.csv", index=False)

total = merged.mean_original_abundance.sum()
summary = (merged.groupby(["analysis_group", "oxygen_phenotype"], dropna=False)
           .agg(species_features=("species", "size"),
                mean_original_abundance=("mean_original_abundance", "sum"))
           .reset_index())
summary["share_of_direct_species_abundance"] = summary.mean_original_abundance / total
summary.to_csv(ROOT / "outputs/base_oxygen_annotation_coverage.csv", index=False)

reviewed = merged.evidence_level.ne("none")
strict = merged.evidence_level.eq("high") & merged.analysis_group.isin(["oxygen_tolerant", "anaerobic"])
broad = merged.evidence_level.isin(["high", "moderate", "low"]) & merged.analysis_group.isin(["oxygen_tolerant", "anaerobic"])
print(f"Reviewed features: {reviewed.sum()}/{len(merged)}; abundance coverage={merged.loc[reviewed, 'mean_original_abundance'].sum()/total:.3%}")
print(f"Strict high-confidence assigned coverage={merged.loc[strict, 'mean_original_abundance'].sum()/total:.3%}")
print(f"Broad assigned coverage={merged.loc[broad, 'mean_original_abundance'].sum()/total:.3%}")
print(summary.to_string(index=False))
