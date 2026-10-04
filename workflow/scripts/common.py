"""Shared, documented transformations used by all figure scripts."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "config"))
from settings import INPUT_CSV, OUTPUT_ROOT

META = {"sample_id", "age", "weight"}

def load_data(path: Path = INPUT_CSV) -> pd.DataFrame:
    data = pd.read_csv(path)
    required = {"sample_id", "age", "weight"}
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"Input is missing required columns: {sorted(missing)}")
    if data["sample_id"].duplicated().any():
        raise ValueError("Input contains duplicate sample_id values. Freeze a deduplicated cohort first.")
    return data

def taxa(data: pd.DataFrame) -> pd.DataFrame:
    return data[[c for c in data.columns if c not in META]].apply(pd.to_numeric, errors="coerce").fillna(0)

def age_group(age: pd.Series) -> pd.Categorical:
    return pd.cut(age, [0, 2, 7, np.inf], labels=["Puppy (0-2)", "Adult (>2-7)", "Mature (>7)"], include_lowest=True)

def level_name(clade: str, prefix: str) -> str:
    for piece in str(clade).split("|"):
        if piece.startswith(prefix):
            return piece[len(prefix):]
    return "Unclassified"

def collapse(X: pd.DataFrame, prefix: str, drop_unclassified: bool = True) -> pd.DataFrame:
    groups = {column: level_name(column, prefix) for column in X.columns}
    out = X.T.groupby(groups).sum().T
    if drop_unclassified and "Unclassified" in out:
        out = out.drop(columns="Unclassified")
    return out.loc[:, ~out.columns.str.startswith(("OFG", "OGF", "OGFB", "GGB", "FGB", "CFGB"))]

def prevalence(X: pd.DataFrame, cutoff: float) -> pd.DataFrame:
    return X.loc[:, (X > 0).mean(axis=0) >= cutoff]

def outdir(name: str) -> Path:
    path = OUTPUT_ROOT / name
    path.mkdir(parents=True, exist_ok=True)
    return path
