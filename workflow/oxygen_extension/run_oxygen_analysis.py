#!/usr/bin/env python3
"""Run evidence-tiered oxygen analyses and create manuscript-style figures.

All primary calculations use original MetaPhlAn relative-abundance units.
Annotation-only normalization is never used. QC3 is the primary cohort and
QC4 (age >= 1 year) is the sensitivity cohort.
"""
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
from scipy.stats import spearmanr, rankdata
import statsmodels.api as sm
from statsmodels.nonparametric.smoothers_lowess import lowess
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = HERE / "analysis_results"
FIG = OUT / "figures"
OUT.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

INPUTS = {
    "qc3": ROOT / "inputs/qc3_duplicates_cluster1_and_minus1_removed.csv",
    "qc4": ROOT / "inputs/qc4_duplicates_cluster_1_minus_1_removed_age_ge_1.csv",
}

COLORS = {
    "Obligate aerobe": "#2155A4",
    "Aerobe": "#3F83CC",
    "Facultative anaerobe": "#62B0C9",
    "Microaerophile": "#78C7A3",
    "Variable oxygen-tolerant": "#B7D977",
    "Anaerobe": "#D6604D",
    "Unassigned direct species": "#B8B8B8",
    "Unresolved below species": "#E8E8E8",
    "Oxygen-tolerant": "#2E78B7",
    "Anaerobic": "#C94A3F",
    "Log-ratio": "#684A9E",
}

DISPLAY = {
    "obligate_aerobe": "Obligate aerobe",
    "aerobe": "Aerobe",
    "facultative_anaerobe": "Facultative anaerobe",
    "microaerophile": "Microaerophile",
    "variable_microaerophile_facultative": "Variable oxygen-tolerant",
    "variable_aerobe_facultative_microaerophile": "Variable oxygen-tolerant",
    "anaerobe": "Anaerobe",
}


def load_annotation() -> pd.DataFrame:
    queue = pd.read_csv(HERE / "outputs/species_annotation_review_queue.csv")[["clade", "species", "mean_original_abundance"]]
    ann = pd.read_csv(HERE / "base_oxygen_annotation_v1.csv")
    out = queue.merge(ann, on="species", how="left", validate="one_to_one")
    out["evidence_level"] = out.evidence_level.fillna("none")
    out["analysis_group"] = out.analysis_group.fillna("excluded")
    out["oxygen_phenotype"] = out.oxygen_phenotype.fillna("unreviewed")
    out["annotation_level"] = out.annotation_level.fillna("unreviewed")
    return out


def annotation_mask(ann: pd.DataFrame, tier: str) -> pd.Series:
    if tier == "strict":
        return ann.evidence_level.eq("high") & ann.annotation_level.eq("species")
    if tier == "expanded":
        return ann.evidence_level.isin(["high", "moderate"])
    if tier == "broad":
        return ann.evidence_level.isin(["high", "moderate", "low"])
    raise ValueError(tier)


def group_abundance(df: pd.DataFrame, ann: pd.DataFrame, tier: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    mask = annotation_mask(ann, tier)
    active = ann.loc[mask].copy()
    values = df[active.clade].apply(pd.to_numeric, errors="raise")
    groups = {}
    for phenotype, part in active.groupby("oxygen_phenotype"):
        label = DISPLAY.get(phenotype, phenotype)
        abundance = values[part.clade].sum(axis=1)
        groups[label] = groups.get(label, 0) + abundance
    result = pd.DataFrame(groups, index=df.index)
    for expected in list(COLORS)[:6]:
        if expected not in result:
            result[expected] = 0.0
    result["Oxygen-tolerant"] = sum(
        result[name] for name in ["Obligate aerobe", "Aerobe", "Facultative anaerobe", "Microaerophile", "Variable oxygen-tolerant"]
    )
    result["Anaerobic"] = result["Anaerobe"]
    result["Log-ratio"] = np.log10((result["Oxygen-tolerant"] + 1e-6) / (result["Anaerobic"] + 1e-6))
    return result, active


def association_rows(df: pd.DataFrame, grouped: pd.DataFrame, cohort: str, tier: str) -> list[dict]:
    rows = []
    for outcome in grouped.columns:
        rho, p = spearmanr(df.age, grouped[outcome])
        row = {"cohort": cohort, "tier": tier, "outcome": outcome, "n": len(df), "spearman_rho": rho, "spearman_p": p}
        complete = df.age.notna() & df.weight.notna() & grouped[outcome].notna()
        y = rankdata(grouped.loc[complete, outcome])
        age = rankdata(df.loc[complete, "age"])
        weight = rankdata(df.loc[complete, "weight"])
        fit = sm.OLS(y, sm.add_constant(np.column_stack([age, weight]))).fit(cov_type="HC3")
        row.update({"weight_adjusted_n": int(complete.sum()), "rank_age_beta_weight_adjusted": fit.params[1],
                    "rank_age_p_weight_adjusted": fit.pvalues[1]})
        rows.append(row)
    return rows


def smooth_xy(age: pd.Series, y: pd.Series, frac: float = 0.42) -> np.ndarray:
    return lowess(y.to_numpy(), age.to_numpy(), frac=frac, it=0, return_sorted=True)


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIG / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_multicategory_trends(df: pd.DataFrame, grouped: pd.DataFrame, tier: str) -> None:
    categories = ["Obligate aerobe", "Aerobe", "Facultative anaerobe", "Microaerophile", "Variable oxygen-tolerant", "Anaerobe"]
    fig, ax = plt.subplots(figsize=(8.8, 5.4))
    for category in categories:
        smooth = smooth_xy(df.age, grouped[category] * 100)
        ax.plot(smooth[:, 0], np.maximum(smooth[:, 1], 0), lw=2.2, color=COLORS[category], label=category)
    ax.set(xlabel="Age (years)", ylabel="Original MetaPhlAn relative abundance (%)", xlim=(0, 17))
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=.2)
    ax.legend(frameon=False, ncol=3, fontsize=8.5, loc="upper center", bbox_to_anchor=(.5, 1.16))
    fig.tight_layout(rect=[0, 0, 1, .92])
    save_figure(fig, f"oxygen_multicategory_age_trends_{tier}")


def plot_binary_sensitivity(df: pd.DataFrame, grouped_by_tier: dict[str, pd.DataFrame]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.2))
    tier_style = {"strict": ("Strict", "-"), "expanded": ("Expanded", "--"), "broad": ("Broad", ":")}
    for ax, outcome in zip(axes, ["Oxygen-tolerant", "Anaerobic", "Log-ratio"]):
        for tier, grouped in grouped_by_tier.items():
            y = grouped[outcome] if outcome == "Log-ratio" else grouped[outcome] * 100
            smooth = smooth_xy(df.age, y)
            plotted = smooth[:, 1] if outcome == "Log-ratio" else np.maximum(smooth[:, 1], 0)
            ax.plot(smooth[:, 0], plotted, lw=2.2, ls=tier_style[tier][1], label=tier_style[tier][0],
                    color=COLORS[outcome])
        ax.set_xlim(0, 17)
        ax.set_xlabel("Age (years)")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=.2)
        ax.set_title(outcome)
    axes[0].set_ylabel("Original relative abundance (%)")
    axes[1].set_ylabel("Original relative abundance (%)")
    axes[2].set_ylabel("log10(oxygen-tolerant / anaerobic)")
    axes[0].legend(frameon=False)
    fig.tight_layout()
    save_figure(fig, "oxygen_binary_annotation_sensitivity")


def plot_age_group_composition(df: pd.DataFrame, grouped: pd.DataFrame, ann: pd.DataFrame, tier: str) -> None:
    direct_total = df[ann.clade].sum(axis=1)
    shown = ["Obligate aerobe", "Aerobe", "Facultative anaerobe", "Microaerophile", "Variable oxygen-tolerant", "Anaerobe"]
    plot = grouped[shown].copy()
    plot["Unassigned direct species"] = (direct_total - plot.sum(axis=1)).clip(lower=0)
    plot["Unresolved below species"] = (1 - direct_total).clip(lower=0)
    plot["age_group"] = pd.cut(df.age, [0, 2, 7, np.inf], labels=["Puppy (0–2)", "Adult (>2–7)", "Mature (>7)"], include_lowest=True)
    means = plot.groupby("age_group", observed=True).mean() * 100
    counts = plot.groupby("age_group", observed=True).size()
    means.to_csv(OUT / f"oxygen_age_group_composition_{tier}.csv")
    fig, ax = plt.subplots(figsize=(9.2, 5.7))
    bottom = np.zeros(len(means))
    for category in shown + ["Unassigned direct species", "Unresolved below species"]:
        labels = [f"{group}\n(n={counts.loc[group]})" for group in means.index]
        ax.bar(labels, means[category], bottom=bottom, color=COLORS[category], width=.68, label=category)
        bottom += means[category].to_numpy()
    ax.set(ylabel="Mean original relative abundance (%)", ylim=(0, 100), xlabel="")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, fontsize=9)
    fig.tight_layout()
    save_figure(fig, f"oxygen_age_group_composition_{tier}")


def plot_species_effects(df: pd.DataFrame, ann: pd.DataFrame) -> pd.DataFrame:
    rows = []
    active = ann[annotation_mask(ann, "broad") & ann.analysis_group.isin(["oxygen_tolerant", "anaerobic"])].copy()
    for row in active.itertuples():
        rho, p = spearmanr(df.age, df[row.clade])
        rows.append({"species": row.species, "display_species": row.species.replace("_", " "), "phenotype": DISPLAY.get(row.oxygen_phenotype, row.oxygen_phenotype),
                     "analysis_group": row.analysis_group, "evidence_level": row.evidence_level,
                     "mean_abundance": df[row.clade].mean(), "prevalence": (df[row.clade] > 0).mean(), "rho": rho, "p": p})
    effects = pd.DataFrame(rows)
    effects["q_bh_within_displayed_base"] = multipletests(effects.p, method="fdr_bh")[1]
    effects.to_csv(OUT / "annotated_species_age_effects_qc3.csv", index=False)
    top = effects.nlargest(22, "mean_abundance").sort_values("rho")
    fig, ax = plt.subplots(figsize=(8.4, 7.4))
    colors = [COLORS[DISPLAY.get(p, p)] for p in top.phenotype]
    ax.barh(top.display_species.str.replace("sp ", "sp. ", regex=False), top.rho, color=colors)
    ax.axvline(0, color="#444444", lw=.8)
    ax.set(xlabel="Spearman correlation with age", ylabel="")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", alpha=.2)
    fig.tight_layout()
    save_figure(fig, "dominant_annotated_species_age_correlations")
    return effects


def plot_neisseria(df: pd.DataFrame, ann: pd.DataFrame) -> None:
    neis = ann[ann.species.str.startswith("Neisseria_")].copy()
    neis["mean"] = [df[c].mean() for c in neis.clade]
    neis = neis.nlargest(6, "mean")
    fig, axes = plt.subplots(2, 3, figsize=(12.2, 7.2), sharex=True)
    rng = np.random.default_rng(42)
    for ax, row in zip(axes.flat, neis.itertuples()):
        y = df[row.clade] * 100
        ax.scatter(df.age + rng.normal(0, .025, len(df)), y, s=7, alpha=.18, color="#254E70", linewidths=0)
        smooth = smooth_xy(df.age, y, frac=.5)
        ax.plot(smooth[:, 0], np.maximum(smooth[:, 1], 0), color="#C94A3F", lw=2)
        rho, p = spearmanr(df.age, y)
        ax.set_title(row.species.replace("_", " "), fontsize=10, fontstyle="italic")
        ax.text(.97, .95, f"rho={rho:.2f}\np={p:.2g}", transform=ax.transAxes, ha="right", va="top", fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xlim(0, 17)
        ax.set_xlabel("Age (years)")
        ax.set_ylabel("Original abundance (%)")
    fig.tight_layout()
    save_figure(fig, "neisseria_species_age_trends")


def robustness_table(df: pd.DataFrame, ann: pd.DataFrame) -> pd.DataFrame:
    active = ann[annotation_mask(ann, "broad") & ann.analysis_group.eq("oxygen_tolerant")].copy()
    base = df[active.clade].sum(axis=1)
    scenarios = {"All broad oxygen-tolerant taxa": base}
    neis = active[active.species.str.startswith("Neisseria_")].clade
    neisseriales = active[active.clade.str.contains("|o__Neisseriales|", regex=False)].clade
    conch = active[active.species.eq("Conchiformibius_steedae")].clade
    scenarios["Without Neisseria"] = base - df[neis].sum(axis=1)
    scenarios["Without Neisseriales"] = base - df[neisseriales].sum(axis=1)
    scenarios["Without Conchiformibius steedae"] = base - df[conch].sum(axis=1)
    for row in active.nlargest(8, "mean_original_abundance").itertuples():
        label = f"Without {row.species.replace('_', ' ')}"
        if label not in scenarios:
            scenarios[label] = base - df[row.clade]
    rows = []
    for scenario, values in scenarios.items():
        rho, p = spearmanr(df.age, values)
        rows.append({"scenario": scenario, "n": len(df), "spearman_rho": rho, "p": p, "mean_original_abundance": values.mean()})
    out = pd.DataFrame(rows)
    out["q_bh"] = multipletests(out.p, method="fdr_bh")[1]
    out.to_csv(OUT / "oxygen_tolerant_leave_out_sensitivity_qc3.csv", index=False)
    plot = out.sort_values("spearman_rho")
    fig, ax = plt.subplots(figsize=(8.4, 5.9))
    ax.barh(plot.scenario, plot.spearman_rho, color="#3F83CC")
    ax.axvline(0, color="#444444", lw=.8)
    ax.set(xlabel="Spearman correlation with age", ylabel="")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", alpha=.2)
    fig.tight_layout()
    save_figure(fig, "oxygen_tolerant_leave_out_sensitivity")
    return out


def main() -> None:
    ann = load_annotation()
    all_stats = []
    primary_grouped = {}
    for cohort, path in INPUTS.items():
        df = pd.read_csv(path)
        if df.sample_id.duplicated().any():
            raise ValueError(f"{cohort} contains duplicate sample IDs")
        for tier in ["strict", "expanded", "broad"]:
            grouped, active = group_abundance(df, ann, tier)
            grouped.insert(0, "sample_id", df.sample_id.values)
            grouped.insert(1, "age", df.age.values)
            grouped.insert(2, "weight", df.weight.values)
            grouped.to_csv(OUT / f"oxygen_group_abundance_{cohort}_{tier}.csv", index=False)
            numerical = grouped.drop(columns=["sample_id", "age", "weight"])
            all_stats.extend(association_rows(df, numerical, cohort, tier))
            if cohort == "qc3":
                primary_grouped[tier] = numerical

    stats = pd.DataFrame(all_stats)
    for (cohort, tier), indexes in stats.groupby(["cohort", "tier"]).groups.items():
        idx = list(indexes)
        stats.loc[idx, "spearman_q_bh_within_tier"] = multipletests(stats.loc[idx, "spearman_p"], method="fdr_bh")[1]
        stats.loc[idx, "weight_adjusted_q_bh_within_tier"] = multipletests(stats.loc[idx, "rank_age_p_weight_adjusted"], method="fdr_bh")[1]
    stats.to_csv(OUT / "oxygen_group_age_associations.csv", index=False)

    primary = pd.read_csv(INPUTS["qc3"])
    plot_multicategory_trends(primary, primary_grouped["strict"], "strict")
    plot_multicategory_trends(primary, primary_grouped["expanded"], "expanded")
    plot_binary_sensitivity(primary, primary_grouped)
    plot_age_group_composition(primary, primary_grouped["expanded"], ann, "expanded")
    plot_species_effects(primary, ann)
    plot_neisseria(primary, ann)
    robustness_table(primary, ann)
    print(stats[stats.outcome.isin(["Oxygen-tolerant", "Anaerobic", "Log-ratio"])].to_string(index=False))


if __name__ == "__main__":
    main()
