#!/usr/bin/env python3
"""Runnable QC2 bacterial-order UMAP from the collaborator notebook."""
from __future__ import annotations
import os, sys
from pathlib import Path
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap
from adjustText import adjust_text
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "config"))
from settings import QC2_INPUT_CSV, OUTPUT_ROOT
INPUT_FILE = Path(os.environ.get("MICROBIOME_QC2_INPUT", QC2_INPUT_CSV))
OUTDIR = Path(os.environ.get("MICROBIOME_OUTPUT", OUTPUT_ROOT)) / "06_figure3_bacterial_umap"

def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT_FILE)
    features = df.drop(columns=["sample_id", "age", "weight"])
    cols = [c for c in features if "|o__" in c and "|f__" not in c and "|o__OFGB" not in c]
    orders = features[cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    orders = orders.div(orders.sum(axis=1).replace(0, 1), axis=0)
    x = np.log10(orders + 1e-5)
    coords = umap.UMAP(metric="correlation", n_neighbors=12, min_dist=.25, random_state=42).fit_transform(x.T)
    labels = [c.split("|o__")[-1] for c in cols]
    summary = pd.DataFrame({"order": labels, "UMAP1": coords[:, 0], "UMAP2": coords[:, 1],
                            "spearman_r": [spearmanr(orders[c], df.age).statistic for c in cols],
                            "mean_abundance": [orders[c].mean() for c in cols]})
    summary.to_csv(OUTDIR / "bacterial_order_umap_statistics.csv", index=False)
    # Preserve the embedding so subsequent publication-format renders do not
    # need to repeat UMAP fitting.
    summary.to_csv(OUTDIR / "bacterial_order_umap_coordinates.csv", index=False)
    size = 30 + 330 * np.sqrt(summary.mean_abundance / summary.mean_abundance.max())
    fig, ax = plt.subplots(figsize=(9.5, 7.5))
    points = ax.scatter(summary.UMAP1, summary.UMAP2, c=summary.spearman_r, cmap="RdBu_r", s=size,
                        alpha=.9, edgecolor="black", linewidth=.4, vmin=-.7, vmax=.7)
    texts = [ax.text(row.UMAP1, row.UMAP2, row.order, fontsize=7) for _, row in summary.iterrows()]
    adjust_text(texts, objects=points, expand=(1.5, 1.6), force_static=(.8, 1.3), force_text=(.5, .7),
                arrowprops=dict(arrowstyle="-", color="gray", lw=.5, shrinkA=5))
    fig.colorbar(points, ax=ax, label="Spearman Correlation Coefficient (with host age)")
    ax.set(xlabel="UMAP1", ylabel="UMAP2")
    fig.savefig(OUTDIR / "bacterial_order_umap.pdf", bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__": main()
