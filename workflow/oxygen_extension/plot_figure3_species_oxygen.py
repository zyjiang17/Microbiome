#!/usr/bin/env python3
"""Create revised Figure 3 and the species-resolution coverage audit."""

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
DEFAULT_TAXON_UMAP = (
    ROOT
    / "oxygen_extension/analysis_results/06_figure3_bacterial_umap/bacterial_order_umap_coordinates.csv"
)
DEFAULT_OXYGEN = ROOT / "oxygen_extension/analysis_results/oxygen_age_group_composition_expanded.csv"
DEFAULT_SAMPLE_OXYGEN = ROOT / "oxygen_extension/analysis_results/oxygen_group_abundance_qc3_expanded.csv"
DEFAULT_OUTPUT = ROOT / "oxygen_extension/analysis_results/figures"

AGE_GROUPS = ["Puppy", "Adult", "Mature"]
AGE_DISPLAY = ["Puppy\n(0–2)", "Adult\n(>2–7)", "Mature\n(>7)"]

def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def short_age_group(value: str) -> str:
    if value.startswith("Puppy"):
        return "Puppy"
    if value.startswith("Adult"):
        return "Adult"
    if value.startswith("Mature"):
        return "Mature"
    raise ValueError(f"Unexpected age group: {value}")


def prepare_species_composition(rows: list[dict[str, str]]) -> list[dict[str, float | str]]:
    by_group = {short_age_group(row["age_group"]): row for row in rows}
    if set(by_group) != set(AGE_GROUPS):
        raise ValueError(f"Expected age groups {AGE_GROUPS}, found {sorted(by_group)}")

    result: list[dict[str, float | str]] = []
    oxygen_columns = [
        "Obligate aerobe",
        "Aerobe",
        "Facultative anaerobe",
        "Microaerophile",
        "Variable oxygen-tolerant",
    ]
    for group in AGE_GROUPS:
        row = by_group[group]
        oxygen_tolerant = sum(float(row[column]) for column in oxygen_columns)
        anaerobic = float(row["Anaerobe"])
        unknown = float(row["Unassigned direct species"])
        unresolved = float(row["Unresolved below species"])
        species_total = oxygen_tolerant + anaerobic + unknown
        if abs((species_total + unresolved) - 100.0) > 0.05:
            raise ValueError(f"Age-group components do not sum to 100 for {group}")
        result.append(
            {
                "Age group": group,
                "Oxygen-tolerant species": 100.0 * oxygen_tolerant / species_total,
                "Anaerobic species": 100.0 * anaerobic / species_total,
                "Unknown oxygen phenotype": 100.0 * unknown / species_total,
                "Species-resolved coverage": species_total,
                "Not resolved to species": unresolved,
            }
        )
    return result


def write_summary(rows: list[dict[str, float | str]], path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def taxon_umap_panel(fig, ax, rows: list[dict[str, str]]) -> None:
    x = np.asarray([float(row["UMAP1"]) for row in rows])
    y = np.asarray([float(row["UMAP2"]) for row in rows])
    rho = np.asarray([float(row["spearman_r"]) for row in rows])
    abundance = np.asarray([float(row["mean_abundance"]) for row in rows])
    size = 30 + 330 * np.sqrt(abundance / abundance.max())
    points = ax.scatter(
        x, y, c=rho, cmap="RdBu_r", s=size, alpha=0.9,
        edgecolor="black", linewidth=0.4, vmin=-0.7, vmax=0.7,
    )
    center_x, center_y = float(np.median(x)), float(np.median(y))
    manual_offsets = {
        "Xanthomonadales": (-10, -10),
        "Propionibacteriales": (10, -10),
        "Tissierellales": (8, -8),
        "Endomicrobiales": (-8, -8),
    }
    for i, (xi, yi, row) in enumerate(zip(x, y, rows)):
        x_offset = 9 if xi <= center_x else -9
        y_offset = 7 if yi <= center_y else -7
        if i % 2:
            y_offset *= -1
        x_offset, y_offset = manual_offsets.get(row["order"], (x_offset, y_offset))
        ax.annotate(
            row["order"],
            xy=(xi, yi),
            xytext=(x_offset, y_offset),
            textcoords="offset points",
            fontsize=6.4,
            ha="left" if x_offset > 0 else "right",
            va="bottom" if y_offset > 0 else "top",
            arrowprops=dict(arrowstyle="-", color="gray", lw=0.45, shrinkA=4, shrinkB=4),
        )
    colorbar = fig.colorbar(points, ax=ax, pad=0.02, fraction=0.05)
    colorbar.set_label("Spearman correlation with host age")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.margins(x=0.08, y=0.08)
    ax.text(0.01, 0.99, "A", transform=ax.transAxes, va="top", fontweight="bold")


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


def sample_level_panel(ax, rows: list[dict[str, str]], panel_label: str = "B") -> None:
    age = np.asarray([float(row["age"]) for row in rows])
    oxygen_tolerant = np.asarray([float(row["Oxygen-tolerant"]) for row in rows])
    anaerobic = np.asarray([float(row["Anaerobic"]) for row in rows])
    classified_total = oxygen_tolerant + anaerobic
    keep = classified_total > 0
    age = age[keep]
    fraction = 100.0 * oxygen_tolerant[keep] / classified_total[keep]
    cx, median, lower, upper = age_window_quantiles(age, fraction)

    ax.scatter(age, fraction, s=18, alpha=0.28, color="#565f6b", edgecolor="none", rasterized=True)
    ax.fill_between(cx, lower, upper, color="#7aa6d8", alpha=0.25, linewidth=0)
    ax.plot(cx, median, color="#2f6fb0", linewidth=2.4, label="Age-window median")
    ax.axhline(50, color="#4d4d4d", linewidth=1, linestyle="--", label="Equal fractions")
    ax.set_xlim(left=0)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Age (years)")
    ax.set_ylabel("Oxygen-tolerant fraction of classified oxygen groups (%)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(0.01, 0.99, panel_label, transform=ax.transAxes, va="top", fontweight="bold")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--taxon-umap", type=Path, default=DEFAULT_TAXON_UMAP)
    parser.add_argument("--oxygen", type=Path, default=DEFAULT_OXYGEN)
    parser.add_argument("--sample-oxygen", type=Path, default=DEFAULT_SAMPLE_OXYGEN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    taxon_umap_rows = read_rows(args.taxon_umap)
    oxygen_rows = prepare_species_composition(read_rows(args.oxygen))
    sample_oxygen_rows = read_rows(args.sample_oxygen)
    write_summary(oxygen_rows, args.output / "figure3b_species_oxygen_composition_data.csv")

    # Standalone revised Figure 3B: one point per dog.
    fig, ax = plt.subplots(figsize=(8.3, 6.2))
    sample_level_panel(ax, sample_oxygen_rows)
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    stem = args.output / "figure3b_sample_level_oxygen_fraction"
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    # Combined Figure 3: taxon UMAP (A) and one point per dog (B).
    fig, axes = plt.subplots(1, 2, figsize=(17.5, 7.2), gridspec_kw={"width_ratios": [1.2, 1]})
    taxon_umap_panel(fig, axes[0], taxon_umap_rows)
    sample_level_panel(axes[1], sample_oxygen_rows)
    axes[1].legend(frameon=False, loc="upper right")
    fig.tight_layout()
    stem = args.output / "figure3_taxon_umap_and_sample_oxygen_fraction"
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    # Separate coverage audit; no physiology categories are mixed into this plot.
    coverage = [float(row["Species-resolved coverage"]) for row in oxygen_rows]
    fig, ax = plt.subplots(figsize=(7.4, 5.4))
    bars = ax.bar(AGE_DISPLAY, coverage, width=0.65, color="#6f7782")
    ax.set_ylim(0, 100)
    ax.set_ylabel("Total abundance resolved to species (%)")
    ax.spines[["top", "right"]].set_visible(False)
    for bar, value in zip(bars, coverage):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 1.0, f"{value:.1f}%", ha="center")
    fig.tight_layout()
    stem = args.output / "supp_species_resolution_coverage_by_age_group"
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
