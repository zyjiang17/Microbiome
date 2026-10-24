#!/usr/bin/env bash
#SBATCH --job-name=dog_humann_age
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=humann_%j.out
#SBATCH --error=humann_%j.err
# Submit from this directory with: sbatch submit_humann.slurm
set -euo pipefail

PROJECT_ROOT="${HUMANN_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
PROJECT_ROOT="$(cd "$PROJECT_ROOT" && pwd)"
INPUT_DIR="${HUMANN_INPUT_DIR:-$PROJECT_ROOT/../humann_inputs}"
OUTPUT_ROOT="${HUMANN_OUTPUT_ROOT:-$PROJECT_ROOT/cluster_runs}"
PIPELINE_ROOT="${MICROBIOME_PIPELINE_DIR:-/u/home/z/zyjiang/canine_oral_microbiome_repro/pipeline}"
RUN_DIR="$OUTPUT_ROOT/slurm_${SLURM_JOB_ID:-manual}"
mkdir -p "$RUN_DIR"/results "$RUN_DIR"/figures

module load python/3.10.20
VENV="${HUMANN_VENV:-$PIPELINE_ROOT/.venv_preview_20261002}"
if [[ ! -f "$VENV/bin/activate" ]]; then
  echo "ERROR: venv not found at $VENV. Set HUMANN_VENV to an environment with requirements_humann.txt installed." >&2
  exit 1
fi
source "$VENV/bin/activate"
export MPLBACKEND=Agg
export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib_humann_${SLURM_JOB_ID:-manual}"
export XDG_CACHE_HOME="${TMPDIR:-/tmp}/cache_humann_${SLURM_JOB_ID:-manual}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export HUMANN_INPUT_DIR="$INPUT_DIR"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"

for f in qc1_pathabundance_cpm_n_reads_ge_2_mil.tsv qc1_pathcoverage_n_reads_ge_2_mil.tsv; do
 [[ -s "$INPUT_DIR/$f" ]] || { echo "ERROR: missing $INPUT_DIR/$f" >&2; exit 2; }
done
python3 -c 'import pandas,numpy,scipy,statsmodels,sklearn,matplotlib; print("Python libraries OK")'
python3 "$PROJECT_ROOT/code/preflight_humann.py" --input-dir "$INPUT_DIR" --out "$RUN_DIR/preflight.json"
cp "$PROJECT_ROOT/code/"*.py "$RUN_DIR/"
cp "$PROJECT_ROOT/code/export_official_metacyc_mapping.R" "$RUN_DIR/"
sha256sum "$INPUT_DIR"/qc1_*.tsv > "$RUN_DIR/input_sha256.txt"
cd "$RUN_DIR"
python3 -u run_humann_analysis.py \
 --abundance "$INPUT_DIR/qc1_pathabundance_cpm_n_reads_ge_2_mil.tsv" \
 --coverage "$INPUT_DIR/qc1_pathcoverage_n_reads_ge_2_mil.tsv" \
 --output "$RUN_DIR"
python3 -u run_focus_panels.py
python3 -u run_weight_analysis.py
# Compare selected headline values to the supplied exploratory reference snapshot.
python3 "$PROJECT_ROOT/code/compare_reference.py" \
 --reference "$PROJECT_ROOT/reference_results/analysis_summary.json" \
 --actual "$RUN_DIR/results/analysis_summary.json" \
 --output "$RUN_DIR/reference_comparison.json"
if [[ -d "$PIPELINE_ROOT/inputs" ]]; then
  python3 "$PROJECT_ROOT/code/check_taxonomy_overlap.py" \
    --input-dir "$INPUT_DIR" --taxa-dir "$PIPELINE_ROOT/inputs" --out-dir "$RUN_DIR/results" \
    || echo "WARNING: optional taxonomic cohort comparison failed; main analysis completed."
fi
echo "COMPLETED HUMAnN pathway workflow; figures/results: $RUN_DIR"
