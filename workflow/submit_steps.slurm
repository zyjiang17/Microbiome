#!/usr/bin/env bash
# Submit selected steps, for example:
#   sbatch submit_steps.slurm 00 01 03 05
#   sbatch submit_steps.slurm 02 06
#SBATCH --job-name=dog_microbiome_steps
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/dog_microbiome_steps_%j.out
#SBATCH --error=logs/dog_microbiome_steps_%j.err

set -euo pipefail

PIPELINE_DIR="/u/home/z/zyjiang/canine_oral_microbiome_repro/pipeline"
cd "$PIPELINE_DIR"
mkdir -p logs results

module load python/3.10.20
source "$PIPELINE_DIR/.venv_preview_20261002/bin/activate"
export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib"
export XDG_CACHE_HOME="${TMPDIR:-/tmp}/cache"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
export MICROBIOME_INPUT_DIR="${MICROBIOME_INPUT_DIR:-$PIPELINE_DIR/inputs}"
export MICROBIOME_OUTPUT="${MICROBIOME_OUTPUT:-$PIPELINE_DIR/results/slurm_${SLURM_JOB_ID}}"

if (( $# == 0 )); then
  echo "Provide one or more step numbers: 00 01 02 03 04 05 06 07" >&2
  exit 2
fi

python3 -c 'import numpy, pandas, scipy, sklearn, matplotlib, umap'
for step in "$@"; do
  case "$step" in
    00) script="00_validate_cohort.py" ;;
    01) script="01_figure1_clustermap.py" ;;
    02) script="02_figure2_umap_pseudotime.py" ;;
    03) script="03_figure3_composition.py" ;;
    04) script="04_figure4_taxa_diversity.py" ;;
    05) script="05_figure5_enet.py" ;;
    06) script="06_figure3_bacterial_umap.py" ;;
    07) script="07_figure4_denominator_comparison.py" ;;
    *) echo "Unknown step: $step" >&2; exit 2 ;;
  esac
  echo "==> Running $script"
  python3 "scripts/$script"
done
