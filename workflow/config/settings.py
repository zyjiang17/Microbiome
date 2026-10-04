"""Single configuration point for the reproducible analysis release."""
from pathlib import Path
import os

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = PIPELINE_ROOT.parent
INPUT_DIR = Path(os.environ.get("MICROBIOME_INPUT_DIR", PIPELINE_ROOT / "inputs"))

# These names are the frozen collaborator-uploaded cohorts.  Analyses select
# their input explicitly; do not substitute a global "final" table.
QC1_INPUT_CSV = INPUT_DIR / "qc1_duplicates_removed_only.csv"
QC2_INPUT_CSV = INPUT_DIR / "qc2_duplicates_cluster1_and_minus1_removed.csv"
QC3_INPUT_CSV = INPUT_DIR / "qc3_duplicates_removed_prevalence_1pct.csv"
QC4_INPUT_CSV = INPUT_DIR / "qc4_duplicates_cluster_1_minus_1_removed_age_ge_1.csv"

# Backward-compatible default for helpers that have not yet been given a
# figure-specific input. New scripts should import one of the constants above.
INPUT_CSV = Path(os.environ.get("MICROBIOME_INPUT", QC2_INPUT_CSV))
UMAP_INPUT_CSV = Path(os.environ.get("MICROBIOME_UMAP_INPUT", QC3_INPUT_CSV))
OUTPUT_ROOT = Path(os.environ.get("MICROBIOME_OUTPUT", PIPELINE_ROOT / "results"))
RANDOM_SEED = 42
MIN_PREVALENCE = 0.01
AGE_WEIGHT_MIN_AGE = 1.0
