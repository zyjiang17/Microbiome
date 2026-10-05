#!/usr/bin/env python3
"""Figure 5 age/weight analyses: QC4, no prevalence cutoff."""

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
from scipy.stats import pearsonr, spearmanr
from adjustText import adjust_text
from sklearn.linear_model import ElasticNet, LinearRegression
from sklearn.metrics import mean_absolute_error, median_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT / "config"))
from settings import QC4_INPUT_CSV, OUTPUT_ROOT
INPUT_FILE = Path(os.environ.get("MICROBIOME_QC4_INPUT", QC4_INPUT_CSV))
OUTDIR = Path(os.environ.get("MICROBIOME_OUTPUT", OUTPUT_ROOT)) / "05_figure5"
META_COLS = {"sample_id", "age", "weight"}
# Broad grid retained from the earlier workflow; selected independently in
# each outer-fold training set through the inner GridSearchCV.
PARAM_GRID = {
    "enet__alpha": [1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1, 10, 100],
    "enet__l1_ratio": [0.2, 0.4, 0.6, 0.8],
}


def format_p(p: float) -> str:
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}"


def taxa_matrix(df: pd.DataFrame) -> pd.DataFrame:
    taxa_cols = [c for c in df.columns if c not in META_COLS]
    return df[taxa_cols].apply(pd.to_numeric, errors="coerce").fillna(0)


def feature_spearman_r(X: pd.DataFrame, y: pd.Series) -> pd.DataFrame:
    y = pd.to_numeric(y, errors="coerce")
    rows = []
    for col in X.columns:
        x = pd.to_numeric(X[col], errors="coerce")
        mask = x.notna() & y.notna()
        if mask.sum() < 3 or x[mask].nunique() <= 1:
            continue
        rows.append((col, spearmanr(x[mask], y[mask])[0]))
    return pd.DataFrame(rows, columns=["clade", "r"])


def make_label(clade: str) -> str | float:
    parts = str(clade).split("|")
    if any(part.startswith("t__") for part in parts):
        return np.nan
    return parts[-1]


def make_pipeline(random_state: int = 42) -> Pipeline:
    return Pipeline(
        [
            ("scale", StandardScaler()),
            ("enet", ElasticNet(max_iter=5000, tol=1e-3, random_state=random_state, selection="random")),
        ]
    )


def fold_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    mse = mean_squared_error(y_true, y_pred)
    out = {
        "mae": mean_absolute_error(y_true, y_pred),
        "median_absolute_error": median_absolute_error(y_true, y_pred),
        "r2": r2_score(y_true, y_pred),
        "mse": mse,
        "rmse": float(np.sqrt(mse)),
        "actual_min": float(np.min(y_true)),
        "actual_max": float(np.max(y_true)),
        "predicted_min": float(np.min(y_pred)),
        "predicted_max": float(np.max(y_pred)),
        "actual_range": float(np.max(y_true) - np.min(y_true)),
        "predicted_range": float(np.max(y_pred) - np.min(y_pred)),
        "bias_pred_minus_actual": float(np.mean(y_pred - y_true)),
    }
    out["predicted_actual_range_ratio"] = out["predicted_range"] / out["actual_range"] if out["actual_range"] else np.nan
    if len(y_true) >= 3 and np.std(y_true) > 0 and np.std(y_pred) > 0:
        r, p = pearsonr(y_true, y_pred)
        out.update({"r": r, "p": p})
    else:
        out.update({"r": np.nan, "p": np.nan})
    return out


def plot_taxon_age_weight_correlation(X_both: pd.DataFrame, df_both: pd.DataFrame) -> None:
    age_corr = feature_spearman_r(X_both, df_both["age"]).rename(columns={"r": "r_age"})
    weight_corr = feature_spearman_r(X_both, df_both["weight"]).rename(columns={"r": "r_weight"})
    plot_df = age_corr.merge(weight_corr, on="clade")

    plot_df = plot_df.loc[~plot_df["clade"].str.contains(r"\|t__", regex=True, na=False)].copy()
    plot_df = plot_df.loc[~plot_df["clade"].str.contains(r"GGB|FGB|OFGB|CFGB", regex=True, na=False)].copy()
    plot_df["label"] = plot_df["clade"].apply(make_label)

    age_pos_weight_neg = (plot_df["r_age"] > 0) & (plot_df["r_weight"] < 0)
    age_neg_weight_pos = (plot_df["r_age"] < 0) & (plot_df["r_weight"] > 0)
    plot_df["discord_up_age"] = np.where(age_pos_weight_neg, plot_df["r_age"].abs() + plot_df["r_weight"].abs(), np.nan)
    plot_df["discord_up_weight"] = np.where(age_neg_weight_pos, plot_df["r_age"].abs() + plot_df["r_weight"].abs(), np.nan)

    top_age_up = plot_df.nlargest(5, "discord_up_age").copy()
    top_weight_up = plot_df.nlargest(5, "discord_up_weight").copy()

    fig, ax = plt.subplots(figsize=(8, 8), dpi=150)
    ax.scatter(plot_df["r_age"], plot_df["r_weight"], s=12, alpha=0.5, color="lightgray")
    ax.scatter(top_age_up["r_age"], top_age_up["r_weight"], s=35, color="#C53A33", label="Age up / Weight down")
    ax.scatter(top_weight_up["r_age"], top_weight_up["r_weight"], s=35, color="#3B75AF", label="Age down / Weight up")

    texts = []
    target_x = []
    target_y = []
    bacteroidales_rank = 0
    for _, row in top_age_up.iterrows():
        point_x = row["r_age"]
        point_y = row["r_weight"]
        text_x = point_x
        text_y = point_y
        ha = "left"
        if "Bacteroidales" in row["clade"] and point_x > 0 and point_y < 0:
            text_x += 0.012 + 0.006 * bacteroidales_rank
            text_y += 0.018 + 0.018 * bacteroidales_rank
            ha = "right" if bacteroidales_rank % 2 else "left"
            bacteroidales_rank += 1
        texts.append(ax.text(text_x, text_y, row["label"], fontsize=10, color="#C53A33", ha=ha))
        target_x.append(point_x)
        target_y.append(point_y)

    for _, row in top_weight_up.iterrows():
        point_x = row["r_age"]
        point_y = row["r_weight"]
        texts.append(ax.text(point_x, point_y, row["label"], fontsize=10, color="#3B75AF"))
        target_x.append(point_x)
        target_y.append(point_y)

    adjust_text(
        texts,
        x=target_x,
        y=target_y,
        target_x=target_x,
        target_y=target_y,
        ax=ax,
        expand_points=(2.5, 2.5),
        expand_text=(2.0, 2.0),
        force_text=2.0,
        force_points=1.2,
        lim=500,
        arrowprops=dict(arrowstyle="-", color="0.35", lw=0.7, shrinkA=4, shrinkB=4),
    )

    ax.axhline(0, color="black", linewidth=0.5)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_xlabel("Spearman r (taxon vs age)")
    ax.set_ylabel("Spearman r (taxon vs weight)")
    ax.legend(frameon=False)
    plt.tight_layout()
    fig.savefig(OUTDIR / "figure5_taxon_age_weight_spearman_qc4.png", dpi=300)
    plt.close(fig)

    plot_df.to_csv(OUTDIR / "figure5_taxon_age_weight_spearman_qc4.csv", index=False)
    pd.concat(
        [
            top_age_up.assign(label_group="Age up / Weight down"),
            top_weight_up.assign(label_group="Age down / Weight up"),
        ],
        ignore_index=True,
    ).to_csv(OUTDIR / "figure5_taxon_age_weight_spearman_labeled_qc4.csv", index=False)


def run_nested_cv(X: pd.DataFrame, y: pd.Series, sample_ids: pd.Series, trait: str) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    y_arr = pd.to_numeric(y, errors="coerce").to_numpy()
    sample_ids = sample_ids.astype(str).to_numpy()
    outer_cv = KFold(n_splits=10, shuffle=True, random_state=42)
    rows = []
    fold_rows = []

    for fold, (train_idx, test_idx) in enumerate(outer_cv.split(X), start=1):
        inner_cv = KFold(n_splits=5, shuffle=True, random_state=42 + fold)
        grid = GridSearchCV(
            make_pipeline(random_state=42 + fold),
            PARAM_GRID,
            cv=inner_cv,
            scoring="neg_mean_squared_error",
            n_jobs=1,
        )
        grid.fit(X.iloc[train_idx], y_arr[train_idx])
        pred = grid.predict(X.iloc[test_idx])
        best = grid.best_params_

        for idx, predicted in zip(test_idx, pred):
            rows.append(
                {
                    "trait": trait,
                    "fold": fold,
                    "sample_id": sample_ids[idx],
                    "actual": y_arr[idx],
                    "predicted_cv": predicted,
                    "residual": y_arr[idx] - predicted,
                    "best_alpha": best["enet__alpha"],
                    "best_l1_ratio": best["enet__l1_ratio"],
                }
            )

        metrics = fold_metrics(y_arr[test_idx], pred)
        fold_rows.append(
            {
                "trait": trait,
                "fold": fold,
                "train_n": len(train_idx),
                "test_n": len(test_idx),
                "best_alpha": best["enet__alpha"],
                "best_l1_ratio": best["enet__l1_ratio"],
                **metrics,
            }
        )

    outcomes = pd.DataFrame(rows).sort_values(["fold", "sample_id"])
    fold_summary = pd.DataFrame(fold_rows)
    overall = fold_metrics(outcomes["actual"].to_numpy(), outcomes["predicted_cv"].to_numpy())
    summary = {
        "trait": trait,
        "n": len(outcomes),
        "n_features": X.shape[1],
        "outer_cv_folds": 10,
        "inner_cv_folds": 5,
        "alpha_grid": ";".join(map(str, PARAM_GRID["enet__alpha"])),
        "l1_ratio_grid": ";".join(map(str, PARAM_GRID["enet__l1_ratio"])),
        **overall,
    }
    return summary, outcomes, fold_summary


def plot_predictions(outcomes: pd.DataFrame, trait: str) -> None:
    y_true = outcomes["actual"].to_numpy()
    y_pred = outcomes["predicted_cv"].to_numpy()
    metrics = fold_metrics(y_true, y_pred)
    lr = LinearRegression().fit(y_true.reshape(-1, 1), y_pred)
    line_x = np.linspace(min(y_true.min(), y_pred.min()), max(y_true.max(), y_pred.max()), 100)

    fig, ax = plt.subplots(figsize=(6, 5), dpi=300)
    ax.scatter(y_true, y_pred, s=20, alpha=0.65)
    ax.plot(line_x, line_x, "--", label="ideal: y=x")
    ax.plot(line_x, lr.intercept_ + lr.coef_[0] * line_x, label="trend")
    ax.set_xlabel(f"Actual {trait.lower()}")
    ax.set_ylabel(f"Predicted {trait.lower()}")
    ax.text(
        0.05,
        0.95,
        f"Nested CV\nPearson r = {metrics['r']:.3f}\n{format_p(metrics['p'])}\n"
        f"R² = {metrics['r2']:.3f}\nRMSE = {metrics['rmse']:.2f}",
        transform=ax.transAxes,
        va="top",
    )
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUTDIR / f"figure5_{trait.lower()}_enet_prediction_no_1pct.png")
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT_FILE)
    X = np.log10(taxa_matrix(df) + 1e-6)

    age_df = df[df["age"].notna() & (df["age"] >= 1)].copy()
    weight_df = df[df["age"].notna() & df["weight"].notna() & (df["age"] >= 1)].copy()
    plot_taxon_age_weight_correlation(X.loc[weight_df.index], weight_df)

    results = []
    all_outcomes = []
    all_fold_summaries = []
    for trait, trait_df, y_col in [("Age", age_df, "age"), ("Weight", weight_df, "weight")]:
        summary, outcomes, fold_summary = run_nested_cv(
            X.loc[trait_df.index],
            trait_df[y_col],
            trait_df["sample_id"],
            trait,
        )
        results.append(summary)
        all_outcomes.append(outcomes)
        all_fold_summaries.append(fold_summary)
        outcomes.to_csv(OUTDIR / f"figure5_{trait.lower()}_cv_outcomes_no_1pct.csv", index=False)
        plot_predictions(outcomes, trait)

    pd.concat(all_outcomes, ignore_index=True).to_csv(OUTDIR / "figure5_cv_outcomes_no_1pct.csv", index=False)
    pd.concat(all_fold_summaries, ignore_index=True).to_csv(OUTDIR / "figure5_fold_metrics_no_1pct.csv", index=False)
    pd.DataFrame(results).to_csv(OUTDIR / "figure5_enet_summary_results_no_1pct.csv", index=False)


if __name__ == "__main__":
    main()
