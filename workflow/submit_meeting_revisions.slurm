#!/usr/bin/env bash
# Rebuild the manuscript figures revised at the October meeting.
# Submit from the Slurm-preview login host with:
#   sbatch submit_meeting_revisions.slurm
#SBATCH --job-name=dog_meeting_figs
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=logs/dog_meeting_figs_%j.out
#SBATCH --error=logs/dog_meeting_figs_%j.err

set -euo pipefail

PIPELINE_DIR="/u/home/z/zyjiang/canine_oral_microbiome_repro/pipeline"
cd "$PIPELINE_DIR"
mkdir -p logs final_figures oxygen_extension/analysis_results

module load python/3.10.20
source "$PIPELINE_DIR/.venv_preview_20261002/bin/activate"
python3 -c 'import sys; assert sys.version_info[:2] == (3, 10), sys.version'
python3 -c 'import numpy,pandas,scipy,statsmodels,sklearn,matplotlib,umap,adjustText'

export MPLBACKEND=Agg
export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib_meeting_${SLURM_JOB_ID}"
export XDG_CACHE_HOME="${TMPDIR:-/tmp}/cache_meeting_${SLURM_JOB_ID}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"

# Figure 1: same clustering, corrected age-strip alignment, no oxygen annotation.
python3 scripts/render_figures_only.py --source render_inputs/current --output final_figures

# Figure 2: frozen dog coordinates shown by age, oxygen, pseudotime, then age vs pseudotime.
python3 oxygen_extension/plot_figure2_dog_umap_panels.py

# Figure 3A is the only UMAP refit in this job and is explicitly fit on QC3.
export MICROBIOME_INPUT_DIR="${MICROBIOME_INPUT_DIR:-$PIPELINE_DIR/inputs}"
export MICROBIOME_OUTPUT="$PIPELINE_DIR/oxygen_extension/analysis_results"
python3 scripts/06_figure3_bacterial_umap.py
python3 oxygen_extension/plot_figure3_species_oxygen.py

echo "Meeting-revision outputs:"
ls -lh final_figures/figure1_order_clustermap.* \
  final_figures/figure2_dog_umap_age_oxygen_pseudotime.* \
  final_figures/figure4_age_associated_taxa_and_diversity.* \
  final_figures/figure5_age_enet_prediction.* \
  final_figures/figure6_weight_enet_prediction.* \
  final_figures/supplementary_taxon_age_weight_spearman.* \
  oxygen_extension/analysis_results/figures/figure3_bacterial_umap_and_oxygen_age_trends.*
