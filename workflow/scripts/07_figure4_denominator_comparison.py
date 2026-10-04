#!/usr/bin/env python3
"""Compare two abundance denominators for the prespecified Figure 4 genera."""
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
from settings import QC2_INPUT_CSV, OUTPUT_ROOT

INPUT_FILE = Path(os.environ.get("MICROBIOME_QC2_INPUT", QC2_INPUT_CSV))
OUTDIR = Path(os.environ.get("MICROBIOME_OUTPUT", OUTPUT_ROOT)) / "07_denominator_comparison"
GENERA = ("Porphyromonas", "Conchiformibius")


def genus_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df if "|g__" in c and "|s__" not in c and not c.split("|")[-1].startswith("g__GGB")]


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT_FILE, low_memory=False)
    cols = genus_cols(df)
    genus_total = df[cols].apply(pd.to_numeric, errors="coerce").fillna(0).sum(axis=1).replace(0, np.nan)
    records, stats = [], []
    for genus in GENERA:
        match = [c for c in cols if c.endswith(f"g__{genus}")]
        if len(match) != 1:
            raise ValueError(f"Expected exactly one column for {genus}; found {len(match)}")
        col = match[0]
        raw = pd.to_numeric(df[col], errors="coerce").fillna(0)
        within_genus = (raw / genus_total).fillna(0)
        for definition, y in [("MetaPhlAn relative abundance", raw), ("Proportion among genus-level abundances", within_genus)]:
            rho, p = spearmanr(df.age, y)
            stats.append({"genus": genus, "definition": definition, "n": len(df), "spearman_r": rho, "spearman_p": p,
                          "age_axis": "linear years", "smoother": "LOWESS frac=0.65, it=0"})
            records.append(pd.DataFrame({"sample_id": df.sample_id, "age": df.age, "genus": genus,
                                         "definition": definition, "value": y}))
    long = pd.concat(records, ignore_index=True)
    long.to_csv(OUTDIR / "figure4_denominator_comparison_data.csv", index=False)
    pd.DataFrame(stats).to_csv(OUTDIR / "figure4_denominator_comparison_statistics.csv", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(14, 9), dpi=300, sharex=True)
    definitions = ["MetaPhlAn relative abundance", "Proportion among genus-level abundances"]
    for i, genus in enumerate(GENERA):
        for j, definition in enumerate(definitions):
            ax = axes[i, j]
            d = long[(long.genus == genus) & (long.definition == definition)]
            fit = lowess(d.value, d.age, frac=.65, it=0, return_sorted=True)
            summary = next(x for x in stats if x["genus"] == genus and x["definition"] == definition)
            ax.scatter(d.age, d.value, alpha=.35, edgecolors="none", s=10)
            ax.plot(fit[:, 0], fit[:, 1], color="red", linewidth=1.8)
            ax.set(xlabel="Age (years)", ylabel=definition, title=f"{genus}: r = {summary['spearman_r']:.3f}")
            ax.set_xlim(left=0)
    fig.tight_layout()
    fig.savefig(OUTDIR / "figure4_denominator_comparison_linear_age.png", bbox_inches="tight")
    fig.savefig(OUTDIR / "figure4_denominator_comparison_linear_age.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
