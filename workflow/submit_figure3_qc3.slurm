#!/usr/bin/env bash
# Refit the bacterial-order UMAP on QC3 and rebuild the combined Figure 3.
# Submit with: sbatch submit_figure3_qc3.slurm
# Uses the same Python environment as the latest workflow on GitHub main.
#SBATCH --job-name=dog_figure3_qc3
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/dog_figure3_qc3_%j.out
#SBATCH --error=logs/dog_figure3_qc3_%j.err

set -euo pipefail

PIPELINE_DIR="/u/home/z/zyjiang/canine_oral_microbiome_repro/pipeline"
cd "$PIPELINE_DIR"
mkdir -p logs oxygen_extension/analysis_results

module load python/3.10.20
source "$PIPELINE_DIR/.venv_preview_20261002/bin/activate"
python3 -c 'import sys; assert sys.version_info[:2] == (3, 10), sys.version'

export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib"
export XDG_CACHE_HOME="${TMPDIR:-/tmp}/cache"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"

export MICROBIOME_INPUT_DIR="${MICROBIOME_INPUT_DIR:-$PIPELINE_DIR/inputs}"
export MICROBIOME_OUTPUT="$PIPELINE_DIR/oxygen_extension/analysis_results"

python3 -c 'import numpy, pandas, scipy, sklearn, matplotlib, umap'
python3 scripts/06_figure3_bacterial_umap.py
python3 oxygen_extension/plot_figure3_species_oxygen.py

echo "Figure 3 QC3 outputs:"
ls -lh oxygen_extension/analysis_results/figures/figure3_taxon_umap_and_sample_oxygen_fraction.*
