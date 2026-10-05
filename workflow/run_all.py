#!/usr/bin/env python3
"""Run the figure-specific cohorts in manuscript workflow order."""
from __future__ import annotations
import subprocess, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent
STEPS = ["00_validate_cohort.py", "01_figure1_clustermap.py", "02_figure2_umap_pseudotime.py", "03_figure3_composition.py", "04_figure4_taxa_diversity.py", "05_figure5_enet.py", "06_figure3_bacterial_umap.py", "08_figure5_shannon_diversity.py"]
for step in STEPS:
    print(f"\n==> {step}", flush=True)
    subprocess.run([sys.executable, str(ROOT / "scripts" / step)], check=True)
