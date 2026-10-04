# Canine oral microbiome workflow

Reproducible analysis and figure-rendering workflow for the canine oral microbiome manuscript. The workflow keeps the analytical cohort, saved intermediate tables, Slurm submission scripts, final title-free figures, and the code that generated them together in one repository.

## Repository contents

| Path | Contents |
| --- | --- |
| `inputs/` | Four curated taxonomic abundance matrices used by the workflow. |
| `scripts/` | Python code for cohort validation, figure analyses, and figure-only rendering. |
| `config/settings.py` | Cohort-to-analysis mapping and output settings. |
| `submit_all.slurm` | Full Hoffman2 workflow. |
| `submit_steps.slurm` | Selected analytical steps only. |
| `submit_render_only.slurm` | Redraw figures from saved outputs without refitting models. |
| `render_inputs/current/` | Saved tables used by the render-only workflow. |
| `analysis_outputs/` | Current analytical outputs supporting the manuscript figures. |
| `final_figures/` | Title-free manuscript figures in PNG (300 dpi) and PDF formats. |
| `provenance/collaborator_github/` | Original collaborator code retained for provenance. |

Raw FASTQ files and MetaPhlAn input/output files are deliberately not included in this GitHub package because of their size and controlled data-upload workflow. They should be archived separately in the approved data repository.

## Cohorts

The pipeline uses a figure-specific cohort rather than a single matrix for all analyses. The mapping is defined in `config/settings.py`.

| Dataset | Samples | Purpose |
| --- | ---: | --- |
| QC1 | 1,125 | Duplicate removal only; retained for provenance. |
| QC2 | 1,125 | Duplicate removal plus 1% prevalence filtering; used for clustermap and UMAP analyses. |
| QC3 | 1,043 | Duplicate removal plus exclusion of the designated UMAP island and DBSCAN noise; used for composition and genus-age trend figures. |
| QC4 | 675 | QC3 exclusions with age ≥1 year; used for age prediction. Weight prediction uses the 618 samples with recorded weight. |

## Analytical workflow

1. `00_validate_cohort.py` — validates input cohorts and writes sample manifests.
2. `01_figure1_clustermap.py` — order-level clustermap using QC2.
3. `02_figure2_umap_pseudotime.py` — dog-level UMAP and pseudotime using QC2.
4. `03_figure3_composition.py` — order composition by age group using QC3.
5. `04_figure4_taxa_diversity.py` — raw MetaPhlAn relative-abundance trends for *Porphyromonas* and *Conchiformibius* using QC3.
6. `05_figure5_enet.py` — nested cross-validated Elastic Net models using QC4.
7. `06_figure3_bacterial_umap.py` — bacterial-order UMAP using QC2.

The Elastic Net workflow applies log10(abundance + 1e-6), standardization within the model pipeline, 10 outer folds, and 5 inner folds for parameter selection. The archived search grid is alpha = 1e-5 to 100 and L1 ratio = 0.2, 0.4, 0.6, or 0.8.

## Hoffman2 setup

The Slurm scripts are configured for the Hoffman2 Slurm-preview environment with Python 3.10.20. From the repository root:

```bash
module load python/3.10.20
python3 -m venv .venv_preview
source .venv_preview/bin/activate
pip install -r requirements.txt
```

Update the virtual-environment path in the `submit_*.slurm` scripts if using a different environment name or location.

## Run the analysis

```bash
sbatch submit_all.slurm
```

Run selected steps, for example the UMAP/pseudotime and bacterial-order UMAP:

```bash
sbatch submit_steps.slurm 02 06
```

Each Slurm run writes into `results/slurm_<job-id>/` unless `MICROBIOME_OUTPUT` is set explicitly.

## Render figures only

Use the render-only workflow after an analysis has been accepted and only the layout, captions, internal titles, or export format need changing:

```bash
sbatch submit_render_only.slurm render_inputs/current
```

This workflow does **not** run MetaPhlAn, UMAP fitting, or Elastic Net fitting. It uses the saved abundance matrices, UMAP coordinates, pseudotime values, and held-out ENet predictions in `render_inputs/current/`. It writes 300-dpi PNG and vector PDF files to `results/render_<job-id>/`.

Figures intentionally have no internal titles to follow Frontiers formatting; their scientific descriptions belong in the numbered figure captions.

## GitHub guidance

The present package is approximately 31 MB and can be uploaded as a normal Git repository. Do not add raw sequencing reads, scratch directories, virtual environments, or Slurm log archives to Git. For larger future intermediate data, use an archival repository and cite its persistent accession in the manuscript.

## Reproducibility notes

- UMAP analyses use fixed random seeds.
- The saved outputs in `analysis_outputs/` and `render_inputs/current/` record the exact coordinates and held-out predictions behind the distributed figures.
- The final figures in `final_figures/` were produced by the render-only workflow from those saved values.
