# Render-only inputs

`render_figures_only.py` redraws figures using these saved outputs and does not
fit UMAP or Elastic Net models.  Populate `current/` with this structure:

```
01_figure1/figure1_order_abundance_qc3_duplicates_1pct.csv
01_figure1/qc3_sample_manifest.csv
02_figure2/dog_umap_coordinates.csv
02_figure2/pseudotime_main_component.csv
03_figure3/figure3_order_mean_composition_by_age_group_within_orders.csv
04_figure4/figure4_genus_relative_abundance_qc2.csv
05_figure5/figure5_cv_outcomes_no_1pct.csv
05_figure5/figure5_taxon_age_weight_spearman_qc4.csv
05_figure5/figure5_taxon_age_weight_spearman_labeled_qc4.csv
06_figure3_bacterial_umap/bacterial_order_umap_coordinates.csv
```

The clustermap is redrawn from its saved abundance matrix and sample manifest;
all other figures are drawn from saved coordinates or analytical outputs.
