#!/usr/bin/env python3
"""Sample-specific oxygen figures for the canine oral microbiome manuscript."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-dog-oxygen")
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "inputs/qc3_duplicates_cluster1_and_minus1_removed.csv"
DEFAULT_EFFECTS = ROOT / "oxygen_extension/analysis_results/annotated_species_age_effects_qc3.csv"
DEFAULT_OXYGEN = ROOT / "oxygen_extension/analysis_results/oxygen_group_abundance_qc3_expanded.csv"
DEFAULT_OUTPUT = ROOT / "oxygen_extension/analysis_results/figures"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def age_window_quantiles(
    x: np.ndarray,
    y: np.ndarray,
    half_width_years: float = 0.75,
    min_samples: int = 20,
):
    grid = np.arange(float(x.min()), float(x.max()) + 0.05, 0.1)
    centers, medians, lower, upper = [], [], [], []
    for center in grid:
        values = y[np.abs(x - center) <= half_width_years]
        if len(values) < min_samples:
            continue
        centers.append(center)
        medians.append(float(np.median(values)))
        lower.append(float(np.quantile(values, 0.25)))
        upper.append(float(np.quantile(values, 0.75)))
    return np.asarray(centers), np.asarray(medians), np.asarray(lower), np.asarray(upper)


def plot_balance(oxygen_rows: list[dict[str, str]], output: Path) -> None:
    age = np.asarray([float(row["age"]) for row in oxygen_rows])
    oxygen_tolerant = np.asarray([float(row["Oxygen-tolerant"]) for row in oxygen_rows])
    anaerobic = np.asarray([float(row["Anaerobic"]) for row in oxygen_rows])
    classified_total = oxygen_tolerant + anaerobic
    keep = classified_total > 0
    age = age[keep]
    oxygen_fraction = 100.0 * oxygen_tolerant[keep] / classified_total[keep]
    cx, median, lower, upper = age_window_quantiles(age, oxygen_fraction)

    fig, ax = plt.subplots(figsize=(8.2, 5.8))
    ax.scatter(age, oxygen_fraction, s=18, alpha=0.28, color="#565f6b", edgecolor="none", rasterized=True)
    ax.fill_between(cx, lower, upper, color="#7aa6d8", alpha=0.25, linewidth=0)
    ax.plot(cx, median, color="#2f6fb0", linewidth=2.4)
    ax.axhline(50, color="#4d4d4d", linewidth=1, linestyle="--")
    ax.set_xlim(left=0)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Age (years)")
    ax.set_ylabel("Oxygen-tolerant fraction of classified oxygen groups (%)")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    stem = output / "sample_level_oxygen_tolerant_fraction_vs_age"
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_species_heatmap(input_path: Path, effect_rows: list[dict[str, str]], output: Path) -> None:
    oxygen_tolerant = [row for row in effect_rows if row["analysis_group"] == "oxygen_tolerant"]
    oxygen_tolerant.sort(key=lambda row: float(row["mean_abundance"]), reverse=True)
    selected = oxygen_tolerant[:12]
    selected_species = [row["species"] for row in selected]

    with input_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("Input abundance table has no header")
        clade_by_species = {}
        for column in reader.fieldnames:
            if "|s__" in column and "|t__" not in column:
                clade_by_species[column.rsplit("|s__", 1)[1]] = column
        missing = [species for species in selected_species if species not in clade_by_species]
        if missing:
            raise ValueError(f"Selected species absent from QC3 input: {missing}")

        sample_ids, age, matrix_rows = [], [], []
        for row in reader:
            sample_ids.append(row["sample_id"])
            age.append(float(row["age"]))
            matrix_rows.append([float(row[clade_by_species[species]]) for species in selected_species])

    age_array = np.asarray(age)
    matrix = np.asarray(matrix_rows, dtype=float).T
    order = np.argsort(age_array, kind="stable")
    sorted_age = age_array[order]
    log_matrix = np.log10(matrix[:, order] + 1e-5)

    phenotype_labels = {
        "Obligate aerobe": "obligate aerobe",
        "Aerobe": "aerobe",
        "Facultative anaerobe": "facultative",
        "Microaerophile": "microaerophile",
        "Variable oxygen-tolerant": "variable",
    }
    row_labels = [
        f"{row['display_species']} ({phenotype_labels.get(row['phenotype'], row['phenotype'].lower())})"
        for row in selected
    ]

    fig = plt.figure(figsize=(15.5, 7.0))
    grid = fig.add_gridspec(2, 1, height_ratios=[0.36, 6], hspace=0.04)
    age_ax = fig.add_subplot(grid[0])
    heat_ax = fig.add_subplot(grid[1], sharex=age_ax)

    age_image = age_ax.imshow(
        sorted_age[np.newaxis, :],
        aspect="auto",
        cmap="viridis",
        interpolation="nearest",
        extent=[0, len(sorted_age), 0, 1],
    )
    age_ax.set_yticks([])
    age_ax.set_xticks([])
    age_bar = fig.colorbar(age_image, ax=age_ax, orientation="horizontal", fraction=0.8, pad=0.25, aspect=45)
    age_bar.set_label("Age (years)")

    image = heat_ax.imshow(
        log_matrix,
        aspect="auto",
        cmap="magma",
        vmin=-5,
        vmax=0,
        interpolation="nearest",
        extent=[0, len(sorted_age), len(selected_species), 0],
    )
    heat_ax.set_yticks(np.arange(len(selected_species)) + 0.5)
    heat_ax.set_yticklabels(row_labels)
    heat_ax.set_xlabel("Individual dogs ordered by age")
    heat_ax.set_xticks([])

    for boundary_age in (2, 7):
        boundary = int(np.searchsorted(sorted_age, boundary_age, side="right"))
        heat_ax.axvline(boundary, color="white", linewidth=1.1, alpha=0.9)
        age_ax.axvline(boundary, color="white", linewidth=1.1, alpha=0.9)

    abundance_bar = fig.colorbar(image, ax=heat_ax, fraction=0.018, pad=0.015)
    abundance_bar.set_label(r"log$_{10}$(relative abundance + 10$^{-5}$)")
    fig.tight_layout()
    stem = output / "sample_level_oxygen_tolerant_species_heatmap"
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--effects", type=Path, default=DEFAULT_EFFECTS)
    parser.add_argument("--oxygen", type=Path, default=DEFAULT_OXYGEN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    oxygen_rows = read_rows(args.oxygen)
    effect_rows = read_rows(args.effects)
    plot_balance(oxygen_rows, args.output)
    plot_species_heatmap(args.input, effect_rows, args.output)


if __name__ == "__main__":
    main()
