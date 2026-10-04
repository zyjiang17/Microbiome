#!/usr/bin/env python3
"""Faithful runnable version of the collaborator's QC2 dog UMAP notebook."""
from __future__ import annotations

import hashlib
import json
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
import umap
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, shortest_path
from scipy.stats import spearmanr
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
from statsmodels.nonparametric.smoothers_lowess import lowess

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "config"))
from settings import QC2_INPUT_CSV, OUTPUT_ROOT

INPUT_FILE = Path(os.environ.get("MICROBIOME_QC2_INPUT", QC2_INPUT_CSV))
OUTDIR = Path(os.environ.get("MICROBIOME_OUTPUT", OUTPUT_ROOT)) / "02_figure2"
N_NEIGHBORS, MIN_DIST, SEED, PSEUDOTIME_K = 12, 0.25, 42, 10


def order_matrix(df: pd.DataFrame) -> pd.DataFrame:
    features = df.drop(columns=["sample_id", "age", "weight"])
    cols = [c for c in features if "|o__" in c and "|f__" not in c and "|o__OFGB" not in c]
    x = features[cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    return x.div(x.sum(axis=1).replace(0, 1), axis=0)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT_FILE)
    if df.sample_id.duplicated().any():
        raise ValueError("QC2 must already contain unique sample_id values")
    metadata = df[["sample_id", "age", "weight"]].copy()
    matrix = order_matrix(df)
    embedding = umap.UMAP(metric="correlation", n_neighbors=N_NEIGHBORS,
                          min_dist=MIN_DIST, random_state=SEED).fit_transform(np.log10(matrix + 1e-5))
    result = metadata.copy()
    result[["UMAP1", "UMAP2"]] = embedding

    fig, ax = plt.subplots(figsize=(9.5, 7.5))
    points = ax.scatter(result.UMAP1, result.UMAP2, c=result.age, cmap="viridis", s=35,
                        alpha=.9, edgecolor="black", linewidth=.4)
    ax.set(xlabel="UMAP1", ylabel="UMAP2")
    fig.colorbar(points, ax=ax, label="Age (years)")
    fig.savefig(OUTDIR / "dog_sample_umap.pdf", bbox_inches="tight")
    plt.close(fig)

    dbscan = DBSCAN(eps=.55, min_samples=10).fit_predict(embedding)
    result["dbscan_cluster"] = dbscan
    result.to_csv(OUTDIR / "dog_umap_coordinates.csv", index=False)
    result.loc[result.dbscan_cluster.isin([1, -1]), ["sample_id", "dbscan_cluster"]].to_csv(
        OUTDIR / "cluster1_and_minus1_samples.csv", index=False)

    nn = NearestNeighbors(n_neighbors=PSEUDOTIME_K + 1, metric="euclidean").fit(embedding)
    distances, indices = nn.kneighbors(embedding)
    starts, ends, weights = [], [], []
    for i in range(len(embedding)):
        for distance, neighbor in zip(distances[i][1:], indices[i][1:]):
            starts.append(i); ends.append(neighbor); weights.append(float(distance))
    graph = csr_matrix((weights, (starts, ends)), shape=(len(embedding), len(embedding)))
    graph = graph.maximum(graph.T)
    component_count, labels = connected_components(graph, directed=False)
    counts = np.bincount(labels)
    main_component = int(np.argmax(counts))
    root = int(np.argmin(result.age.to_numpy()))
    if labels[root] != main_component:
        raise RuntimeError("Youngest sample is not in the main pseudotime component")
    keep = labels == main_component
    trimmed = graph[keep][:, keep]
    root_trimmed = int(np.cumsum(keep)[root] - 1)
    distance = shortest_path(trimmed, directed=False, indices=root_trimmed)
    pseudo = distance / distance.max()
    age = result.loc[keep, "age"].to_numpy()
    rho, p = spearmanr(age, pseudo)
    curve = lowess(pseudo, age, frac=.4)
    pd.DataFrame({"sample_id": result.loc[keep, "sample_id"], "age": age, "pseudotime": pseudo}).to_csv(
        OUTDIR / "pseudotime_main_component.csv", index=False)
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(age, pseudo, alpha=.5, edgecolors="none", s=10)
    ax.plot(curve[:, 0], curve[:, 1], color="red", linewidth=1)
    ax.set(xlabel="Age (years)", ylabel="Pseudotime (0-1)")
    ax.text(.96, .03, f"ρ = {rho:.3f}\np < 0.001", transform=ax.transAxes, va="bottom", ha="right")
    fig.savefig(OUTDIR / "pseudotime_vs_age.pdf", bbox_inches="tight")
    plt.close(fig)
    audit = {"input": str(INPUT_FILE), "sha256": hashlib.sha256(INPUT_FILE.read_bytes()).hexdigest(),
             "samples": len(df), "order_features": matrix.shape[1], "metric": "correlation",
             "n_neighbors": N_NEIGHBORS, "min_dist": MIN_DIST, "dbscan_eps": .55,
             "dbscan_min_samples": 10, "pseudotime_knn": PSEUDOTIME_K,
             "connected_components": int(component_count), "main_component_n": int(counts[main_component]),
             "excluded_component_n": int(len(df) - counts[main_component]), "spearman_rho": float(rho), "spearman_p": float(p)}
    (OUTDIR / "umap_pseudotime_audit.json").write_text(json.dumps(audit, indent=2) + "\n")


if __name__ == "__main__":
    main()
