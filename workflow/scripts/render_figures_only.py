#!/usr/bin/env python3
"""Render publication figures from saved intermediate data only.

This script never fits UMAP or Elastic Net models and never invokes MetaPhlAn.
It intentionally consumes the coordinate, abundance, and held-out prediction
tables created by the analytical workflow.  The source directory must contain
the folders listed in ``render_inputs/README.md``.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.stats import pearsonr, spearmanr
from statsmodels.nonparametric.smoothers_lowess import lowess


ORDER_COLORS = {
    "Moraxellales": "#6783cf", "Neisseriales": "#2e5dc2", "Flavobacteriales": "#c5d3ef",
    "Burkholderiales": "#297f5a", "Pasteurellales": "#64ba85", "Enterobacterales": "#a9ce72",
    "Bacteroidales": "#ed994b", "Campylobacterales": "#dc3b29", "Desulfovibrionales": "#b81d15",
    "Other": "#e6e6e6",
}


def save(fig: plt.Figure, stem: Path) -> None:
    """Write review-friendly PNG and archival vector PDF versions."""
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def require(source: Path, name: str) -> Path:
    path = source / name
    if not path.exists():
        raise FileNotFoundError(f"Missing render input: {path}")
    return path


def format_p(p: float) -> str:
    """Format a calculated p-value consistently across title-free figures."""
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}"


def render_figure1(source: Path, output: Path) -> None:
    """Recreate only the visual clustering layout from saved order abundances."""
    order = pd.read_csv(require(source, "01_figure1/figure1_order_abundance_qc2_duplicates_1pct.csv"), index_col=0)
    manifest = pd.read_csv(require(source, "01_figure1/qc2_sample_manifest.csv"))
    age_group = pd.cut(manifest.age, bins=[0, 2, 7, np.inf], labels=["Puppy", "Adult", "Mature"], include_lowest=True)
    age_palette = {"Puppy": "#ACDFFF", "Adult": "#FFEB56", "Mature": "#07D4A1"}
    aerobic = {"Neisseriales", "Pasteurellales", "Lactobacillales", "Bacillales", "Micrococcales",
               "Corynebacteriales", "Burkholderiales", "Flavobacteriales", "Mycoplasmatales"}
    anaerobic = {"Bacteroidales", "Campylobacterales", "Fusobacteriales", "Spirochaetales",
                 "Erysipelotrichales", "Clostridiales", "Eubacteriales", "Desulfovibrionales",
                 "Synergistales", "Tissierellales", "Veillonellales"}
    category = lambda x: "aerobic" if x in aerobic else "anaerobic" if x in anaerobic else "other"
    mat = np.log10(order + 1e-5).T
    row_linkage = linkage(mat.rank(axis=1), method="average", metric="correlation")
    col_linkage = linkage(mat.rank(axis=1).T, method="average", metric="correlation")
    row_order, col_order = leaves_list(row_linkage), leaves_list(col_linkage)
    ordered = mat.iloc[row_order, col_order]
    row_colors = ["#1F77B4" if category(t) == "aerobic" else "#D62728" if category(t) == "anaerobic" else "#D0D2D3" for t in mat.index]
    sample_color = dict(zip(manifest.sample_id.astype(str), age_group.astype(str).map(age_palette)))
    col_colors = [sample_color.get(str(x), "#D3D3D3") for x in mat.columns]
    ordered_row_colors = [row_colors[i] for i in row_order]
    ordered_col_colors = [col_colors[i] for i in col_order]
    fig = plt.figure(figsize=(18, 11), dpi=300)
    ax_cbar = fig.add_axes([.02, .78, .04, .15]); ax_cd = fig.add_axes([.20, .755, .70, .16])
    ax_cc = fig.add_axes([.20, .718, .70, .03]); ax_rd = fig.add_axes([.015, .10, .16, .61])
    ax_rc = fig.add_axes([.17, .10, .03, .61]); ax = fig.add_axes([.20, .10, .70, .61])
    for dend_ax, data, kwargs in [(ax_cd, col_linkage, {}), (ax_rd, row_linkage, {"orientation": "left"})]:
        dendrogram(data, ax=dend_ax, no_labels=True, color_threshold=0, above_threshold_color="#666666",
                   link_color_func=lambda _: "#666666", **kwargs)
        dend_ax.set_xticks([]); dend_ax.set_yticks([])
        for spine in dend_ax.spines.values(): spine.set_visible(False)
    col_lookup = {x: i for i, x in enumerate(dict.fromkeys(ordered_col_colors))}
    row_lookup = {x: i for i, x in enumerate(dict.fromkeys(ordered_row_colors))}
    ax_cc.imshow([[col_lookup[x] for x in ordered_col_colors]], aspect="auto", cmap=ListedColormap(list(col_lookup)))
    ax_rc.imshow([[row_lookup[x]] for x in ordered_row_colors], aspect="auto", cmap=ListedColormap(list(row_lookup)))
    for axis in (ax_cc, ax_rc):
        axis.set_xticks([]); axis.set_yticks([])
        for spine in axis.spines.values(): spine.set_visible(False)
    image = ax.imshow(ordered, aspect="auto", cmap="RdYlBu_r", vmin=-5, vmax=0, interpolation="nearest")
    fig.colorbar(image, cax=ax_cbar, label="log10(relative abundance)")
    ax.set(xlabel=f"Samples (n={ordered.shape[1]:,})", ylabel="")
    ax.set_xticks([]); ax.set_yticks(range(len(ordered.index))); ax.set_yticklabels(ordered.index)
    ax.yaxis.tick_right(); ax.tick_params(axis="y", labelright=True, labelleft=False, right=True, left=False)
    for tick in ax.get_yticklabels(): tick.set_color("#1F77B4" if category(tick.get_text()) == "aerobic" else "#D62728" if category(tick.get_text()) == "anaerobic" else "black")
    plt.setp(ax.get_yticklabels(), fontsize=12)
    fig.legend(handles=[Patch(facecolor=age_palette["Puppy"], label="Puppy (0–2)"),
                        Patch(facecolor=age_palette["Adult"], label="Adult (2–7)"),
                        Patch(facecolor=age_palette["Mature"], label="Mature (7+)"),
                        Patch(facecolor="#1F77B4", label="Aerobic / Facultative"),
                        Patch(facecolor="#D62728", label="Anaerobic"), Patch(facecolor="#B0BEC5", label="Other")],
               loc="upper left", bbox_to_anchor=(.09, .93), frameon=False, fontsize=10)
    save(fig, output / "figure1_order_clustermap")


def render_figure2(source: Path, output: Path) -> None:
    coords = pd.read_csv(require(source, "02_figure2/dog_umap_coordinates.csv"))
    fig, ax = plt.subplots(figsize=(7, 6))
    points = ax.scatter(coords.UMAP1, coords.UMAP2, c=coords.age, cmap="viridis", s=28,
                        alpha=.85, edgecolor="black", linewidth=.25)
    ax.set(xlabel="UMAP 1", ylabel="UMAP 2")
    fig.colorbar(points, ax=ax, label="Dog age (years)")
    save(fig, output / "figure2_dog_sample_umap")

    pseudo = pd.read_csv(require(source, "02_figure2/pseudotime_main_component.csv"))
    rho, p = spearmanr(pseudo.age, pseudo.pseudotime)
    curve = lowess(pseudo.pseudotime, pseudo.age, frac=.4)
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(pseudo.age, pseudo.pseudotime, alpha=.5, edgecolors="none", s=10)
    ax.plot(curve[:, 0], curve[:, 1], color="#d55e00", linewidth=2)
    ax.set(xlabel="Age (years)", ylabel="Pseudotime (0–1)")
    ax.text(.96, .03, f"Spearman ρ = {rho:.3f}\n{format_p(p)}", transform=ax.transAxes,
            va="bottom", ha="right")
    save(fig, output / "figure2_pseudotime_vs_age")


def render_figure3(source: Path, output: Path) -> None:
    means = pd.read_csv(require(source, "03_figure3/figure3_order_mean_composition_by_age_group_within_orders.csv"), index_col=0)
    ax = means.plot(kind="bar", stacked=True, figsize=(10, 6), width=.8,
                    color=[ORDER_COLORS.get(c, "#d9d9d9") for c in means.columns])
    ax.set_xlabel("")
    ax.set_ylabel("Mean relative abundance among classified orders")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    ax.figure.tight_layout()
    save(ax.figure, output / "figure3_order_composition_by_age_group")


def render_figure4(source: Path, output: Path) -> None:
    data = pd.read_csv(require(source, "04_figure4/figure4_genus_relative_abundance_qc3.csv"))
    for genus in ("Porphyromonas", "Conchiformibius"):
        col = f"g__{genus}"
        rho, p = spearmanr(data.age, data[col])
        curve = lowess(data[col], data.age, frac=.65, it=0)
        fig, ax = plt.subplots(figsize=(7, 5.5))
        ax.scatter(data.age, data[col], alpha=.55, s=20)
        ax.plot(curve[:, 0], curve[:, 1], color="#d55e00", linewidth=2)
        ax.set_xlim(left=0)
        ax.set_ylim(-.05, 1)
        ax.set(xlabel="Age (years)", ylabel="MetaPhlAn relative abundance")
        ax.text(.96, .96, f"Spearman ρ = {rho:.3f}\n{format_p(p)}", transform=ax.transAxes,
                va="top", ha="right")
        save(fig, output / f"figure4_age_and_{genus.lower()}_abundance")


def render_figure5(source: Path, output: Path) -> None:
    outcomes = pd.read_csv(require(source, "05_figure5/figure5_cv_outcomes_no_1pct.csv"))
    for trait, data in outcomes.groupby("trait", sort=False):
        actual, predicted = data.actual.to_numpy(), data.predicted_cv.to_numpy()
        r, p = pearsonr(actual, predicted)
        slope, intercept = np.polyfit(actual, predicted, 1)
        line = np.linspace(min(actual.min(), predicted.min()), max(actual.max(), predicted.max()), 100)
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.scatter(actual, predicted, s=20, alpha=.65)
        ax.plot(line, line, "--", label="Ideal: y = x")
        ax.plot(line, intercept + slope * line, label="Trend")
        ax.set(xlabel=f"Actual {trait.lower()}", ylabel=f"Predicted {trait.lower()}")
        ax.text(.05, .95, f"Nested CV\nPearson r = {r:.3f}\n{format_p(p)}", transform=ax.transAxes,
                va="top", bbox={"facecolor": "white", "edgecolor": "none", "alpha": .88, "pad": 2})
        ax.legend(loc="lower right", frameon=True, framealpha=.88)
        fig.tight_layout()
        save(fig, output / f"figure5_{trait.lower()}_enet_prediction")

    corr = pd.read_csv(require(source, "05_figure5/figure5_taxon_age_weight_spearman_qc4.csv"))
    labels = pd.read_csv(require(source, "05_figure5/figure5_taxon_age_weight_spearman_labeled_qc4.csv"))
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(corr.r_age, corr.r_weight, s=12, alpha=.5, color="lightgray")
    palette = {"Age up / Weight down": "#C53A33", "Age down / Weight up": "#3B75AF"}
    texts, target_x, target_y = [], [], []
    for group, color in palette.items():
        subset = labels[labels.label_group == group]
        ax.scatter(subset.r_age, subset.r_weight, s=35, color=color, label=group)
        for _, row in subset.iterrows():
            texts.append(ax.text(row.r_age, row.r_weight, row.label, fontsize=8, color=color))
            target_x.append(row.r_age)
            target_y.append(row.r_weight)
    adjust_text(
        texts, x=target_x, y=target_y, target_x=target_x, target_y=target_y, ax=ax,
        expand_points=(2.2, 2.2), expand_text=(1.8, 1.8), force_text=1.5, force_points=1.0, lim=500,
        arrowprops=dict(arrowstyle="-", color="0.35", lw=.7, shrinkA=4, shrinkB=4),
    )
    ax.axhline(0, color="black", linewidth=.5)
    ax.axvline(0, color="black", linewidth=.5)
    ax.set(xlabel="Spearman r (taxon vs age)", ylabel="Spearman r (taxon vs weight)")
    ax.margins(x=.15, y=.15)
    ax.legend(frameon=False)
    fig.tight_layout()
    save(fig, output / "figure5_taxon_age_weight_spearman")


def render_figure6(source: Path, output: Path) -> None:
    data = pd.read_csv(require(source, "06_figure3_bacterial_umap/bacterial_order_umap_coordinates.csv"))
    size = 30 + 330 * np.sqrt(data.mean_abundance / data.mean_abundance.max())
    fig, ax = plt.subplots(figsize=(9.5, 7.5))
    points = ax.scatter(data.UMAP1, data.UMAP2, c=data.spearman_r, cmap="RdBu_r", s=size,
                        alpha=.9, edgecolor="black", linewidth=.4, vmin=-.7, vmax=.7)
    texts = [ax.text(row.UMAP1, row.UMAP2, row.order, fontsize=7) for _, row in data.iterrows()]
    adjust_text(texts, objects=points, expand=(1.5, 1.6), force_static=(.8, 1.3), force_text=(.5, .7),
                arrowprops=dict(arrowstyle="-", color="gray", lw=.5, shrinkA=5))
    fig.colorbar(points, ax=ax, label="Spearman correlation with host age")
    ax.set(xlabel="UMAP 1", ylabel="UMAP 2")
    save(fig, output / "figure6_bacterial_order_umap")


def render_shannon_diversity(source: Path, output: Path) -> None:
    """Render Shannon group plots without recomputing diversity values."""
    specs = [
        ("08_figure5_shannon/figure5_species_shannon_age_qc3.csv", "age_group",
         ["Puppy (0–2)", "Adult (>2–7)", "Mature (>7)"],
         ["Puppy", "Adult", "Mature"], "Dog Age Group", "age", "figure5_species_shannon_by_age_group"),
        ("08_figure5_shannon/figure5_species_shannon_weight_qc4.csv", "weight_group",
         ["Small (<20 lb)", "Medium (20–<60 lb)", "Large (≥60 lb)"],
         ["Small", "Medium", "Large"], "Dog Weight Group", "weight", "figure5_species_shannon_by_weight_group"),
    ]
    for filename, group_col, order, display_labels, x_label, statistic_column, stem in specs:
        data = pd.read_csv(require(source, filename))
        present = [group for group in order if (data[group_col].astype(str) == group).any()]
        values = [data.loc[data[group_col].astype(str) == group, "shannon"].to_numpy() for group in present]
        fig, ax = plt.subplots(figsize=(6, 5))
        violin = ax.violinplot(values, showmeans=False, showmedians=False, showextrema=False)
        for body in violin["bodies"]:
            body.set_facecolor("#377BA8"); body.set_edgecolor("#333333"); body.set_linewidth(1.15); body.set_alpha(.96)
        rng = np.random.default_rng(42)
        for index, value in enumerate(values, start=1):
            ax.scatter(index + rng.uniform(-.13, .13, len(value)), value,
                       s=10, alpha=.32, color="#152C3D", linewidths=0, zorder=3)
        ax.boxplot(values, widths=.09, patch_artist=True, showfliers=False,
                   boxprops={"facecolor": "#242424", "edgecolor": "#242424", "alpha": .86},
                   whiskerprops={"color": "#242424", "linewidth": 1.15},
                   capprops={"color": "#242424", "linewidth": 1.15},
                   medianprops={"color": "white", "linewidth": 1.45})
        ax.set_xticks(range(1, len(present) + 1)); ax.set_xticklabels(display_labels, fontsize=9)
        ax.set(xlabel=x_label, ylabel="Shannon Diversity")
        ax.xaxis.label.set_size(11); ax.yaxis.label.set_size(11); ax.tick_params(axis="y", labelsize=9)
        rho, p = spearmanr(data[statistic_column], data.shannon)
        ax.text(.03, .97, f"Spearman ρ = {rho:.3f}\n{format_p(p)}", transform=ax.transAxes,
                va="top", ha="left", fontsize=9,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": .80, "pad": 1.5})
        fig.tight_layout()
        save(fig, output / stem)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path, help="Saved render-input directory")
    parser.add_argument("--output", required=True, type=Path, help="Directory for freshly rendered figures")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    render_figure1(args.source, args.output)
    render_figure2(args.source, args.output)
    render_figure3(args.source, args.output)
    render_figure4(args.source, args.output)
    render_figure5(args.source, args.output)
    render_figure6(args.source, args.output)
    render_shannon_diversity(args.source, args.output)


if __name__ == "__main__":
    main()
