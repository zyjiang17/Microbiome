#!/usr/bin/env bash
# Render existing figure data only; no MetaPhlAn, UMAP fitting, or ENet fitting.
# Usage: sbatch submit_render_only.slurm [render-input-directory]
# Example: sbatch submit_render_only.slurm render_inputs/current
#SBATCH --job-name=dog_microbiome_render
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/dog_microbiome_render_%j.out
#SBATCH --error=logs/dog_microbiome_render_%j.err

set -euo pipefail
PIPELINE_DIR="/u/home/z/zyjiang/canine_oral_microbiome_repro/pipeline"
cd "$PIPELINE_DIR"
module load python/3.10.20
source "$PIPELINE_DIR/.venv_preview_20261002/bin/activate"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib"
export XDG_CACHE_HOME="${TMPDIR:-/tmp}/cache"

SOURCE_DIR="${1:-render_inputs/current}"
OUTPUT_DIR="results/render_${SLURM_JOB_ID}"
python3 scripts/render_figures_only.py --source "$SOURCE_DIR" --output "$OUTPUT_DIR"
