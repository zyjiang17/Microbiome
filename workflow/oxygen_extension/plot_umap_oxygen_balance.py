#!/usr/bin/env python3
"""Color the frozen dog-level UMAP by age and oxygen balance.

The UMAP coordinates are never refit here. The left panel reproduces the
existing continuous-age coloring for all QC2 dogs. The right panel joins the
QC3 expanded oxygen annotation by sample ID and colors matched dogs by the
oxygen-tolerant:anaerobe log10 ratio. Samples without a matched oxygen result
are omitted from the oxygen panel.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-dog-oxygen")
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_COORDS = ROOT / "analysis_outputs/02_figure2/dog_umap_coordinates.csv"
DEFAULT_OXYGEN = ROOT / "oxygen_extension/analysis_results/oxygen_group_abundance_qc3_expanded.csv"
DEFAULT_OUTPUT = ROOT / "oxygen_extension/analysis_results/figures"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coordinates", type=Path, default=DEFAULT_COORDS)
    parser.add_argument("--oxygen", type=Path, default=DEFAULT_OXYGEN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    coordinates = read_csv(args.coordinates)
    oxygen_rows = read_csv(args.oxygen)
    if not coordinates or not oxygen_rows:
        raise ValueError("Coordinate and oxygen tables must both be non-empty")

    oxygen_by_sample = {row["sample_id"]: float(row["Log-ratio"]) for row in oxygen_rows}
    if len(oxygen_by_sample) != len(oxygen_rows):
        raise ValueError("Oxygen table contains duplicate sample IDs")

    sample_ids = [row["sample_id"] for row in coordinates]
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("UMAP coordinate table contains duplicate sample IDs")

    x = [float(row["UMAP1"]) for row in coordinates]
    y = [float(row["UMAP2"]) for row in coordinates]
    age = [float(row["age"]) for row in coordinates]
    matched = [i for i, sample_id in enumerate(sample_ids) if sample_id in oxygen_by_sample]
    unmatched = [i for i, sample_id in enumerate(sample_ids) if sample_id not in oxygen_by_sample]
    log_ratio = [oxygen_by_sample[sample_ids[i]] for i in matched]

    missing_coordinates = sorted(set(oxygen_by_sample).difference(sample_ids))
    if missing_coordinates:
        raise ValueError(
            f"{len(missing_coordinates)} oxygen samples lack frozen UMAP coordinates; "
            f"first missing ID: {missing_coordinates[0]}"
        )

    ratio_min, ratio_max = min(log_ratio), max(log_ratio)
    if not ratio_min < 0 < ratio_max:
        raise ValueError("The oxygen log-ratio must span zero for a centered diverging scale")

    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.2), sharex=True, sharey=True)
    for ax in axes:
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.spines[["top", "right"]].set_visible(False)

    age_points = axes[0].scatter(
        x,
        y,
        c=age,
        cmap="viridis",
        s=24,
        alpha=0.88,
        edgecolor="black",
        linewidth=0.25,
        rasterized=True,
    )
    age_bar = fig.colorbar(age_points, ax=axes[0], fraction=0.046, pad=0.03)
    age_bar.set_label("Age (years)")
    axes[0].text(0.02, 0.98, "A", transform=axes[0].transAxes, va="top", fontweight="bold")
    ratio_points = axes[1].scatter(
        [x[i] for i in matched],
        [y[i] for i in matched],
        c=log_ratio,
        cmap="coolwarm_r",
        norm=TwoSlopeNorm(vmin=ratio_min, vcenter=0.0, vmax=ratio_max),
        s=24,
        alpha=0.9,
        edgecolor="black",
        linewidth=0.25,
        rasterized=True,
    )
    ratio_bar = fig.colorbar(ratio_points, ax=axes[1], fraction=0.046, pad=0.03)
    ratio_bar.set_label(r"log$_{10}$(oxygen-tolerant / anaerobic)")
    axes[1].text(0.02, 0.98, "B", transform=axes[1].transAxes, va="top", fontweight="bold")

    fig.tight_layout()
    args.output.mkdir(parents=True, exist_ok=True)
    stem = args.output / "dog_umap_age_and_oxygen_balance"
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    audit = {
        "coordinate_file": str(args.coordinates),
        "oxygen_file": str(args.oxygen),
        "umap_was_refit": False,
        "qc2_coordinate_samples": len(coordinates),
        "qc3_oxygen_matched_samples": len(matched),
        "qc2_only_gray_samples": len(unmatched),
        "oxygen_log_ratio_min": ratio_min,
        "oxygen_log_ratio_max": ratio_max,
        "output_png": str(stem.with_suffix(".png")),
        "output_pdf": str(stem.with_suffix(".pdf")),
    }
    (args.output / "dog_umap_age_and_oxygen_balance_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
