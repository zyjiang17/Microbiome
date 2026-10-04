#!/usr/bin/env python3
"""QC3 genus-specific trends using original MetaPhlAn relative abundance."""
from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.nonparametric.smoothers_lowess import lowess

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "config"))
from settings import QC3_INPUT_CSV, OUTPUT_ROOT

INPUT_FILE = Path(os.environ.get("MICROBIOME_QC3_INPUT", QC3_INPUT_CSV))
OUTDIR = Path(os.environ.get("MICROBIOME_OUTPUT", OUTPUT_ROOT)) / "04_figure4"


def genus_columns(df: pd.DataFrame) -> list[str]:
    return [c for c in df if "|g__" in c and "|s__" not in c and not c.split("|")[-1].startswith("g__GGB")]


def plot_trend(data: pd.DataFrame, genus: str, target: Path) -> tuple[float, float]:
    fitted = lowess(data[genus], data.age, frac=.65, it=0)
    rho, p = spearmanr(data.age, data[genus])
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.scatter(data.age, data[genus], alpha=.35, edgecolors="none", s=8)
    ax.plot(fitted[:, 0], fitted[:, 1], color="red", linewidth=1.5)
    ax.set_xlim(left=0)
    ax.set_ylim(-.05, 1); ax.set(xlabel="Age (years)", ylabel="MetaPhlAn relative abundance")
    fig.savefig(target, bbox_inches="tight")
    plt.close(fig)
    return float(rho), float(p)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT_FILE, low_memory=False)
    cols = genus_columns(df)
    names = {}
    for genus in ("Conchiformibius", "Porphyromonas"):
        match = [c for c in cols if c.endswith(f"g__{genus}")]
        if len(match) != 1: raise ValueError(f"Expected one genus column for {genus}; found {len(match)}")
        names[genus] = match[0]
    data = df[["sample_id", "age", "weight"]].copy()
    for genus, column in names.items():
        data[f"g__{genus}"] = pd.to_numeric(df[column], errors="coerce").fillna(0)
    data.to_csv(OUTDIR / "figure4_genus_relative_abundance_qc3.csv", index=False)
    stats = []
    for genus in ("Porphyromonas", "Conchiformibius"):
        label = f"g__{genus}"
        rho, p = plot_trend(data, label, OUTDIR / f"figure4_age_and_{genus.lower()}_abundance.pdf")
        stats.append({"taxon": label, "n": len(data), "spearman_r": rho, "spearman_p": p,
                      "x_axis": "linear age in years", "smoother": "LOWESS frac=0.65, it=0",
                      "abundance_definition": "MetaPhlAn relative abundance"})
    pd.DataFrame(stats).to_csv(OUTDIR / "figure4_age_trend_statistics.csv", index=False)


if __name__ == "__main__":
    main()
