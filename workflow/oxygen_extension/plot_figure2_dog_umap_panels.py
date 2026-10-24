#!/usr/bin/env python3
"""Assemble revised Figure 2 from frozen dog-level UMAP coordinates.

Panels are: A, age; B, oxygen balance; C, pseudotime; and D, pseudotime
against age.  This script never refits UMAP or pseudotime.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-dog-figure2")
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.nonparametric.smoothers_lowess import lowess


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_COORDS = ROOT / "analysis_outputs/02_figure2/dog_umap_coordinates.csv"
DEFAULT_PSEUDOTIME = ROOT / "analysis_outputs/02_figure2/pseudotime_main_component.csv"
DEFAULT_OXYGEN = ROOT / "oxygen_extension/analysis_results/oxygen_group_abundance_qc3_expanded.csv"
DEFAULT_OUTPUT = ROOT / "final_figures"
OXYGEN_VMIN, OXYGEN_VMAX = -1.0, 6.0
OXYGEN_CMAP = LinearSegmentedColormap.from_list(
    "oxygen_balance_asymmetric",
    [(0.00, "#8c2d04"), (1 / 7, "#f2f2f2"), (0.34, "#9ecae1"),
     (0.67, "#3182bd"), (1.00, "#08306b")],
)


def label_panel(ax, label: str) -> None:
    ax.text(0.01, 0.99, label, transform=ax.transAxes, va="top", ha="left",
            fontsize=13, fontweight="bold")


def style_umap(ax) -> None:
    ax.set(xlabel="UMAP 1", ylabel="UMAP 2")
    ax.spines[["top", "right"]].set_visible(False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coordinates", type=Path, default=DEFAULT_COORDS)
    parser.add_argument("--pseudotime", type=Path, default=DEFAULT_PSEUDOTIME)
    parser.add_argument("--oxygen", type=Path, default=DEFAULT_OXYGEN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    coords = pd.read_csv(args.coordinates)
    pseudo = pd.read_csv(args.pseudotime)
    oxygen = pd.read_csv(args.oxygen)[["sample_id", "Log-ratio"]]
    if coords.sample_id.duplicated().any() or pseudo.sample_id.duplicated().any() or oxygen.sample_id.duplicated().any():
        raise ValueError("Sample IDs must be unique in every Figure 2 input")
    data = coords.merge(oxygen, on="sample_id", how="left").merge(
        pseudo[["sample_id", "pseudotime"]], on="sample_id", how="left"
    )

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 11.5), sharex="col", sharey="row")
    ax = axes[0, 0]
    points = ax.scatter(data.UMAP1, data.UMAP2, c=data.age, cmap="viridis", s=23,
                        alpha=.88, edgecolor="black", linewidth=.22, rasterized=True)
    fig.colorbar(points, ax=ax, fraction=.046, pad=.025, label="Age (years)")
    style_umap(ax); label_panel(ax, "A")

    ax = axes[0, 1]
    unmatched = data["Log-ratio"].isna()
    ax.scatter(data.loc[unmatched, "UMAP1"], data.loc[unmatched, "UMAP2"], s=18,
               color="#d4d4d4", alpha=.65, edgecolor="none", rasterized=True)
    matched = ~unmatched
    points = ax.scatter(data.loc[matched, "UMAP1"], data.loc[matched, "UMAP2"],
                        c=data.loc[matched, "Log-ratio"], cmap=OXYGEN_CMAP,
                        norm=Normalize(OXYGEN_VMIN, OXYGEN_VMAX, clip=True), s=23,
                        alpha=.9, edgecolor="black", linewidth=.22, rasterized=True)
    fig.colorbar(points, ax=ax, fraction=.046, pad=.025,
                 label=r"log$_{10}$(oxygen-tolerant / anaerobic)")
    style_umap(ax); label_panel(ax, "B")

    ax = axes[1, 0]
    missing_pseudo = data.pseudotime.isna()
    ax.scatter(data.loc[missing_pseudo, "UMAP1"], data.loc[missing_pseudo, "UMAP2"],
               s=18, color="#d4d4d4", alpha=.65, edgecolor="none", rasterized=True)
    keep = ~missing_pseudo
    points = ax.scatter(data.loc[keep, "UMAP1"], data.loc[keep, "UMAP2"],
                        c=data.loc[keep, "pseudotime"], cmap="magma", vmin=0, vmax=1,
                        s=23, alpha=.9, edgecolor="black", linewidth=.22, rasterized=True)
    fig.colorbar(points, ax=ax, fraction=.046, pad=.025, label="Pseudotime")
    style_umap(ax); label_panel(ax, "C")

    ax = axes[1, 1]
    trend_data = data.loc[keep, ["age", "pseudotime"]].dropna().sort_values("age")
    curve = lowess(trend_data.pseudotime, trend_data.age, frac=.4, it=0)
    rho, p_value = spearmanr(trend_data.age, trend_data.pseudotime)
    ax.scatter(trend_data.age, trend_data.pseudotime, s=15, alpha=.32,
               color="#5e6472", edgecolor="none", rasterized=True)
    ax.plot(curve[:, 0], curve[:, 1], color="#d55e00", linewidth=2.3)
    ax.set(xlabel="Age (years)", ylabel="Pseudotime (0–1)", ylim=(-.03, 1.03))
    ax.spines[["top", "right"]].set_visible(False)
    p_text = "p < 0.001" if p_value < .001 else f"p = {p_value:.3f}"
    ax.text(.97, .04, f"Spearman ρ = {rho:.3f}\n{p_text}", transform=ax.transAxes,
            ha="right", va="bottom")
    label_panel(ax, "D")

    fig.tight_layout()
    args.output.mkdir(parents=True, exist_ok=True)
    stem = args.output / "figure2_dog_umap_age_oxygen_pseudotime"
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    audit = {
        "coordinates": str(args.coordinates), "pseudotime": str(args.pseudotime),
        "oxygen": str(args.oxygen), "umap_was_refit": False,
        "coordinate_samples": int(len(coords)), "oxygen_matched_samples": int(matched.sum()),
        "pseudotime_samples": int(keep.sum()), "oxygen_scale": [OXYGEN_VMIN, OXYGEN_VMAX],
        "oxygen_zero_color_position": 1 / 7,
    }
    (args.output / "figure2_dog_umap_age_oxygen_pseudotime_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
