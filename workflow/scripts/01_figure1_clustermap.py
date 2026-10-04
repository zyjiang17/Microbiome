#!/usr/bin/env python3
"""Figure 1 clustermap: QC2 input with oxygen-preference order annotations."""

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
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage


ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT / "config"))
from settings import QC2_INPUT_CSV, OUTPUT_ROOT
INPUT_FILE = Path(os.environ.get("MICROBIOME_QC2_INPUT", QC2_INPUT_CSV))
OUTDIR = Path(os.environ.get("MICROBIOME_OUTPUT", OUTPUT_ROOT)) / "01_figure1"
META_COLS = {"sample_id", "age", "weight"}

AGE_PALETTE = {"Puppy": "#ACDFFF", "Adult": "#FFEB56", "Mature": "#07D4A1"}
AerobicFacultative = {
    "Neisseriales",
    "Pasteurellales",
    "Lactobacillales",
    "Bacillales",
    "Micrococcales",
    "Corynebacteriales",
    "Burkholderiales",
    "Flavobacteriales",
    "Mycoplasmatales",
}
ANAEROBIC_ORDERS = {
    "Bacteroidales",
    "Campylobacterales",
    "Fusobacteriales",
    "Spirochaetales",
    "Erysipelotrichales",
    "Clostridiales",
    "Eubacteriales",
    "Desulfovibrionales",
    "Synergistales",
    "Tissierellales",
    "Veillonellales",
}


def taxa_matrix(df: pd.DataFrame) -> pd.DataFrame:
    taxa_cols = [c for c in df.columns if c not in META_COLS]
    return df[taxa_cols].apply(pd.to_numeric, errors="coerce").fillna(0)


def extract_taxon_level(clade: str, prefix: str) -> str:
    for part in str(clade).split("|"):
        if part.startswith(prefix):
            return part.replace(prefix, "")
    return "Unclassified"


def collapse_taxa(X: pd.DataFrame, prefix: str) -> pd.DataFrame:
    level_map = {col: extract_taxon_level(col, prefix) for col in X.columns}
    out = X.T.groupby(level_map).sum().T
    if "Unclassified" in out.columns:
        out = out.drop(columns=["Unclassified"])
    return out.loc[:, ~out.columns.str.startswith(("OFG", "OGF", "OGFB", "GGB"))]


def prevalence_filter(X: pd.DataFrame, min_prev: float = 0.01) -> pd.DataFrame:
    return X.loc[:, (X > 0).mean(axis=0) >= min_prev]


def age_groups(df: pd.DataFrame) -> pd.Series:
    return pd.cut(
        df["age"],
        bins=[0, 2, 7, np.inf],
        labels=["Puppy", "Adult", "Mature"],
        include_lowest=True,
    )


def order_category(order: str) -> str:
    if order in AerobicFacultative:
        return "aerobic"
    if order in ANAEROBIC_ORDERS:
        return "anaerobic"
    return "other"


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT_FILE)
    df = df[df["age"].notna()].copy()
    df.index = df["sample_id"].astype(str)

    X = taxa_matrix(df)
    X.index = df.index
    order = prevalence_filter(collapse_taxa(X, "o__"), min_prev=0.01)
    order.to_csv(OUTDIR / "figure1_order_abundance_qc2_duplicates_1pct.csv")

    plot_mat = np.log10(order + 1e-5).T
    row_colors = [
        "#1F77B4" if order_category(t) == "aerobic" else "#D62728" if order_category(t) == "anaerobic" else "#D0D2D3"
        for t in plot_mat.index
    ]

    meta = df[["sample_id", "age"]].copy()
    meta["age_group"] = age_groups(meta)
    meta = meta.drop_duplicates("sample_id")
    sample_to_color = dict(zip(meta["sample_id"].astype(str), meta["age_group"].map(AGE_PALETTE)))
    col_colors = [sample_to_color.get(str(sample), "#D3D3D3") for sample in plot_mat.columns]

    row_linkage = linkage(plot_mat.rank(axis=1), method="average", metric="correlation")
    col_linkage = linkage(plot_mat.rank(axis=1).T, method="average", metric="correlation")
    row_order = leaves_list(row_linkage)
    col_order = leaves_list(col_linkage)
    ordered = plot_mat.iloc[row_order, col_order]
    ordered_row_colors = [row_colors[i] for i in row_order]
    ordered_col_colors = [col_colors[i] for i in col_order]

    fig = plt.figure(figsize=(18, 11), dpi=300)
    ax_cbar = fig.add_axes([0.02, 0.78, 0.04, 0.15])
    ax_col_dend = fig.add_axes([0.20, 0.755, 0.70, 0.16])
    ax_col_colors = fig.add_axes([0.20, 0.718, 0.70, 0.030])
    ax_row_dend = fig.add_axes([0.015, 0.10, 0.16, 0.61])
    ax_row_colors = fig.add_axes([0.17, 0.10, 0.03, 0.61])
    ax_heatmap = fig.add_axes([0.20, 0.10, 0.70, 0.61])

    dendrogram(
        col_linkage,
        ax=ax_col_dend,
        no_labels=True,
        color_threshold=0,
        above_threshold_color="#666666",
        link_color_func=lambda _: "#666666",
    )
    for line in ax_col_dend.collections:
        line.set_linewidth(0.5)

    dendrogram(
        row_linkage,
        ax=ax_row_dend,
        orientation="left",
        no_labels=True,
        color_threshold=0,
        above_threshold_color="#666666",
        link_color_func=lambda _: "#666666",
    )
    for ax in [ax_col_dend, ax_row_dend]:
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

    col_lookup = {color: i for i, color in enumerate(dict.fromkeys(ordered_col_colors))}
    row_lookup = {color: i for i, color in enumerate(dict.fromkeys(ordered_row_colors))}
    ax_col_colors.imshow([[col_lookup[c] for c in ordered_col_colors]], aspect="auto", cmap=ListedColormap(list(col_lookup)))
    ax_row_colors.imshow([[row_lookup[c]] for c in ordered_row_colors], aspect="auto", cmap=ListedColormap(list(row_lookup)))
    for ax in [ax_col_colors, ax_row_colors]:
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

    im = ax_heatmap.imshow(ordered, aspect="auto", cmap="RdYlBu_r", vmin=-5, vmax=0, interpolation="nearest")
    fig.colorbar(im, cax=ax_cbar, label="log10(relative abundance)")
    ax_heatmap.set_xlabel(f"Samples (n={ordered.shape[1]:,})")
    ax_heatmap.set_ylabel("")
    ax_heatmap.set_xticks([])
    ax_heatmap.set_yticks(range(len(ordered.index)))
    ax_heatmap.set_yticklabels(ordered.index)
    ax_heatmap.yaxis.tick_right()
    ax_heatmap.tick_params(axis="y", labelright=True, labelleft=False, right=True, left=False)
    for tick in ax_heatmap.get_yticklabels():
        category = order_category(tick.get_text())
        tick.set_color("#1F77B4" if category == "aerobic" else "#D62728" if category == "anaerobic" else "black")
    plt.setp(ax_heatmap.get_yticklabels(), fontsize=12)

    handles = [
        Patch(facecolor=AGE_PALETTE["Puppy"], label="Puppy (0-2)"),
        Patch(facecolor=AGE_PALETTE["Adult"], label="Adult (2-7)"),
        Patch(facecolor=AGE_PALETTE["Mature"], label="Mature (7+)"),
        Patch(facecolor="#1F77B4", label="Aerobic / Facultative\n(health-associated)"),
        Patch(facecolor="#D62728", label="Anaerobic\n(disease-associated)"),
        Patch(facecolor="#B0BEC5", label="Other / non-core"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.09, 0.93), frameon=False, fontsize=10)
    fig.savefig(OUTDIR / "figure1_order_clustermap_qc2_duplicates_1pct.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
