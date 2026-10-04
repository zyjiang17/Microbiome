#!/usr/bin/env bash
# Submit with: sbatch submit_all.slurm
#SBATCH --job-name=dog_microbiome_all
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/dog_microbiome_all_%j.out
#SBATCH --error=logs/dog_microbiome_all_%j.err

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

python3 -c 'import numpy, pandas, scipy, sklearn, matplotlib, umap'
python3 run_all.py
