#!/usr/bin/env python3
"""Species-level Shannon diversity by age and weight group.

Age diversity uses QC3 (the 1,043-sample cluster/noise-excluded cohort).
Weight diversity uses QC4 samples with recorded weight (n=618). Shannon is
calculated from direct, non-strain species abundances with no prevalence filter;
each sample is re-normalized across its classified species before calculation.
"""
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT / "config"))
from settings import QC3_INPUT_CSV, QC4_INPUT_CSV, OUTPUT_ROOT

AGE_INPUT = Path(os.environ.get("MICROBIOME_QC3_INPUT", QC3_INPUT_CSV))
WEIGHT_INPUT = Path(os.environ.get("MICROBIOME_QC4_INPUT", QC4_INPUT_CSV))
OUTDIR = Path(os.environ.get("MICROBIOME_OUTPUT", OUTPUT_ROOT)) / "08_figure5_shannon"


def species_relative_abundance(df: pd.DataFrame) -> pd.DataFrame:
    """Return direct species abundances, normalized within each sample."""
    cols = [
        col for col in df.columns
        if "|s__" in col and "|t__" not in col and not any(tag in col for tag in ("GGB", "FGB", "OFGB", "CFGB"))
    ]
    if not cols:
        raise ValueError("No direct non-strain species columns were found")
    values = df[cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    return values.div(values.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)


def shannon(values: pd.DataFrame) -> pd.Series:
    positive = values.where(values > 0, 1.0)
    return -(values * np.log(positive)).sum(axis=1)


def group_plot(
    data: pd.DataFrame,
    *,
    group_column: str,
    order: list[str],
    display_labels: list[str],
    x_label: str,
    output_stem: Path,
) -> None:
    present = [group for group in order if (data[group_column].astype(str) == group).any()]
    values = [data.loc[data[group_column].astype(str) == group, "shannon"].to_numpy() for group in present]
    fig, ax = plt.subplots(figsize=(6, 5))
    violin = ax.violinplot(values, showmeans=False, showmedians=False, showextrema=False)
    for body in violin["bodies"]:
        body.set_facecolor("#377BA8")
        body.set_edgecolor("#333333")
        body.set_linewidth(1.15)
        body.set_alpha(0.96)
    rng = np.random.default_rng(42)
    for index, value in enumerate(values, start=1):
        ax.scatter(index + rng.uniform(-0.13, 0.13, len(value)), value,
                   s=10, alpha=0.32, color="#152C3D", linewidths=0, zorder=3)
    ax.boxplot(
        values, widths=0.09, patch_artist=True, showfliers=False,
        boxprops={"facecolor": "#242424", "edgecolor": "#242424", "alpha": 0.86},
        whiskerprops={"color": "#242424", "linewidth": 1.15},
        capprops={"color": "#242424", "linewidth": 1.15},
        medianprops={"color": "white", "linewidth": 1.45},
    )
    ax.set_xticks(range(1, len(present) + 1))
    ax.set_xticklabels(display_labels, fontsize=9)
    ax.set(xlabel=x_label, ylabel="Shannon Diversity")
    ax.xaxis.label.set_size(11)
    ax.yaxis.label.set_size(11)
    ax.tick_params(axis="y", labelsize=9)
    fig.tight_layout()
    fig.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    age_df = pd.read_csv(AGE_INPUT)
    age_df = age_df[age_df["age"].notna()].copy()
    age_df["age_group"] = pd.cut(
        age_df["age"], bins=[0, 2, 7, np.inf], labels=["Puppy (0–2)", "Adult (>2–7)", "Mature (>7)"], include_lowest=True
    )
    age_df["shannon"] = shannon(species_relative_abundance(age_df)).to_numpy()
    age_table = age_df[["sample_id", "age", "age_group", "shannon"]].copy()
    age_table.to_csv(OUTDIR / "figure5_species_shannon_age_qc3.csv", index=False)
    group_plot(
        age_table, group_column="age_group", order=["Puppy (0–2)", "Adult (>2–7)", "Mature (>7)"],
        display_labels=["Puppy", "Adult", "Mature"],
        x_label="Dog Age Group",
        output_stem=OUTDIR / "figure5_species_shannon_by_age_group",
    )

    weight_df = pd.read_csv(WEIGHT_INPUT)
    weight_df = weight_df[weight_df["weight"].notna()].copy()
    weight_df["weight_group"] = pd.cut(
        weight_df["weight"], bins=[0, 20, 60, np.inf], labels=["Small (<20 lb)", "Medium (20–<60 lb)", "Large (≥60 lb)"],
        right=False, include_lowest=True,
    )
    weight_df["shannon"] = shannon(species_relative_abundance(weight_df)).to_numpy()
    weight_table = weight_df[["sample_id", "age", "weight", "weight_group", "shannon"]].copy()
    weight_table.to_csv(OUTDIR / "figure5_species_shannon_weight_qc4.csv", index=False)
    group_plot(
        weight_table, group_column="weight_group", order=["Small (<20 lb)", "Medium (20–<60 lb)", "Large (≥60 lb)"],
        display_labels=["Small", "Medium", "Large"],
        x_label="Dog Weight Group",
        output_stem=OUTDIR / "figure5_species_shannon_by_weight_group",
    )

    rows = []
    for trait, data, variable, cohort in [
        ("age", age_table, "age", "QC3"),
        ("weight", weight_table, "weight", "QC4"),
    ]:
        rho, p = spearmanr(data[variable], data["shannon"])
        rows.append({"trait": trait, "cohort": cohort, "n": len(data), "spearman_rho": rho, "spearman_p": p,
                     "diversity_definition": "direct species relative abundance, re-normalized within sample; no prevalence filter"})
    pd.DataFrame(rows).to_csv(OUTDIR / "figure5_species_shannon_statistics.csv", index=False)


if __name__ == "__main__":
    main()
