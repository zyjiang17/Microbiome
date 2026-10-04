#!/usr/bin/env python3
"""Figure 3 composition among classified orders only.

This version excludes features lacking an order-level `o__` annotation, collapses
the remaining features to order level, then renormalizes within each sample.
"""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT / "config"))
from settings import QC2_INPUT_CSV, OUTPUT_ROOT
INPUT_FILE = Path(os.environ.get("MICROBIOME_QC2_INPUT", QC2_INPUT_CSV))
OUTDIR = Path(os.environ.get("MICROBIOME_OUTPUT", OUTPUT_ROOT)) / "03_figure3"
META_COLS = {"sample_id", "age", "weight"}

ORDER_COLORS = {
    "Moraxellales": "#6783cf",
    "Neisseriales": "#2e5dc2",
    "Flavobacteriales": "#c5d3ef",
    "Burkholderiales": "#297f5a",
    "Pasteurellales": "#64ba85",
    "Enterobacterales": "#a9ce72",
    "Bacteroidales": "#ed994b",
    "Campylobacterales": "#dc3b29",
    "Desulfovibrionales": "#b81d15",
    "Other": "#e6e6e6",
}
STACK_ORDER = [
    "Other",
    "Neisseriales",
    "Moraxellales",
    "Flavobacteriales",
    "Burkholderiales",
    "Pasteurellales",
    "Enterobacterales",
    "Bacteroidales",
    "Campylobacterales",
    "Desulfovibrionales",
]


def taxa_matrix(df: pd.DataFrame) -> pd.DataFrame:
    taxa_cols = [c for c in df.columns if c not in META_COLS]
    return df[taxa_cols].apply(pd.to_numeric, errors="coerce").fillna(0)


def extract_order(clade: str) -> str | None:
    for part in str(clade).split("|"):
        if part.startswith("o__"):
            return part.replace("o__", "")
    return None


def collapse_to_classified_orders(X: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    order_map = {col: extract_order(col) for col in X.columns}
    classified_cols = [col for col, order in order_map.items() if order is not None]
    excluded_cols = [col for col, order in order_map.items() if order is None]

    collapsed = X[classified_cols].T.groupby({col: order_map[col] for col in classified_cols}).sum().T
    excluded = pd.DataFrame(
        {
            "taxon_feature": excluded_cols,
            "nonzero_samples": [(X[col] > 0).sum() for col in excluded_cols],
            "total_abundance": [X[col].sum() for col in excluded_cols],
        }
    ).sort_values("total_abundance", ascending=False)
    return collapsed, excluded


def age_groups(df: pd.DataFrame) -> pd.Series:
    return pd.cut(
        df["age"],
        bins=[0, 2, 7, float("inf")],
        labels=["Puppy", "Adult", "Mature"],
        include_lowest=True,
    )


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT_FILE)
    df = df[df["age"].notna()].copy()
    df["age_group"] = age_groups(df)

    classified_orders, excluded = collapse_to_classified_orders(taxa_matrix(df))
    excluded.to_csv(OUTDIR / "figure3_within_orders_excluded_no_order_annotation.csv", index=False)

    order_rel = classified_orders.div(classified_orders.sum(axis=1).replace(0, 1), axis=0)
    named_orders = [order for order in STACK_ORDER if order != "Other"]

    plot_df = pd.DataFrame(index=order_rel.index)
    for order in named_orders:
        plot_df[order] = order_rel[order] if order in order_rel.columns else 0.0
    other_cols = [col for col in order_rel.columns if col not in named_orders]
    plot_df["Other"] = order_rel[other_cols].sum(axis=1) if other_cols else 0.0
    plot_df = plot_df[STACK_ORDER]
    plot_df["age_group"] = df["age_group"].values

    means = plot_df.groupby("age_group", observed=True).mean()
    means.to_csv(OUTDIR / "figure3_order_mean_composition_by_age_group_within_orders.csv")

    colors = [ORDER_COLORS[col] for col in means.columns]
    ax = means.plot(kind="bar", stacked=True, figsize=(10, 6), width=0.8, color=colors)
    ax.set_xlabel("")
    ax.set_ylabel("Mean relative abundance among classified orders")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    plt.tight_layout()
    plt.savefig(OUTDIR / "figure3_order_composition_by_age_group_within_orders.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    main()
