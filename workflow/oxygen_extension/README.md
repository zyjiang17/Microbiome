# Species oxygen audit: first local analysis

Run date: 2026-10-20. Baseline: GitHub commit `8d42a6c`.

## Findings

The frozen GitHub scripts remain unchanged. The new audit is confined to this extension. No remote jobs, transfers, or GitHub pushes were performed.

The cleaned primary cohort contains 1,043 samples and 528 direct species features. The age >=1 sensitivity cohort contains 675 samples. All four frozen input cohorts have unique sample IDs and their hashes reproduce the archived hashes.

The proposed explanation needs refinement: Neisseria species do not share a uniform age trend.

| Species | Full cohort Spearman rho | BH q | Age >=1 rho | BH q |
|---|---:|---:|---:|---:|
| N. dumasiana | -0.052 | 0.303 | -0.338 | 3.17e-18 |
| N. shayeganii | 0.242 | 1.65e-14 | -0.064 | 0.274 |
| N. weaveri | 0.461 | 1.41e-54 | 0.177 | 2.20e-5 |

These are original-abundance correlations. BH adjustment covers all nonconstant direct-species tests separately in each cohort and denominator. Results include all observed species; rare detections require caution. Cohorts are overlapping sensitivity analyses, not independent replications.

The direct Neisseria genus and sum of its species give the same correlations: rho=-0.068 (p=0.0274) in the full cohort and rho=-0.333 (p=5.85e-19) at age >=1. These genus p-values are unadjusted checks. Thus a weak overall genus decline coexists with a strong increase in one aerobic species.

Conchiformibius steedae declines strongly (rho=-0.415, q=2.19e-43) in the full cohort. Bergeyella zoohelcum increases weakly (rho=0.095, q=0.0111). Porphyromonas cangingivalis increases (rho=0.421, q=1.03e-44). These findings do not yet establish the total community oxygen-group association.

Weight sensitivity uses the same 934 complete cases in QC3, and 618 in QC4, for adjusted/unadjusted comparisons. It regresses abundance ranks on age ranks and weight ranks using HC3 standard errors; this is an exploratory rank-linear adjustment, not a nonlinear final model. N. weaveri remains positively associated with age and N. dumasiana negatively associated after adjustment in both cohorts. N. shayeganii remains positive only in the full cohort. The adjustment q-values cover the nine Neisseria features within each cohort; they differ in testing family from the all-species table.

## Baseline reproduction

Executed the original 00, 03, 04 and 08 scripts in the existing isolated Python 3.11 environment. Composition age-group means, genus abundance/trend tables, Shannon sample tables/statistics, and all cohort sample manifests reproduce archived numerical values at rtol=1e-8, atol=1e-10. The cohort audit differs only in machine-specific input paths. The order-exclusion report matches after sorting rows by taxon; tied abundance rows changed order. UMAP, clustermap, and Elastic Net were not rerun in this audit. Therefore this is a partial baseline reproduction, not full workflow validation.

## Taxonomic aggregation issue

The existing composition script includes every feature containing an order token, summing order + family + genus + species and other descendant ranks. Parent abundance overlaps descendant abundance. This gives unequal weights to taxonomic branches and should not be treated as a direct-order composition.

The new comparison selects only columns ending at `o__` (73 columns versus 1,522 mixed-rank columns). In QC3 the mean per-sample total variation distance between the two normalized compositions is 0.0434; the maximum is 0.3810. The largest difference in age-group mean order abundance is 0.0447 (4.47 percentage points). Corresponding QC4 values are 0.0308, 0.2678, and 2.38 percentage points. These are diagnostic comparisons; existing figures have not been overwritten. The clustermap's aggregation helper also sums across descendants and needs a separate review before refitting.

## Denominator audit

Direct species sums average 0.9312 in QC3 (range begins at 0.0463) and 0.9468 in QC4 (minimum 0.6241). Their sums increase with age: rho=0.309 and 0.274 respectively. Direct-order sums average 0.99955 in QC3. Thus missing resolution to species is age-associated; it is not evidence of missing total bacterial abundance.

For primary oxygen-group estimates preserve original relative-abundance units and explicitly retain both unannotated species abundance and abundance unresolved to direct species. Also report sensitivity using the sum of all direct species as denominator. Never renormalize only over reviewed oxygen labels. The within-species normalization used here is a sensitivity analysis, not a fix for missing taxonomic resolution.

## Annotation status and evidence

`species_annotation_current.csv` contains all 528 species, prioritized by mean abundance. It is a review table, not a completed annotation. Only five species currently have non-unknown category seeds, covering on average 26.25% of direct species abundance in QC3. This limited and age-dependent coverage is insufficient for a community-wide claim. The sixth reviewed entry, C. steedae, remains unknown because aerobic culture conditions alone do not prove obligate oxygen dependence.

Evidence is recorded in `reviewed_evidence.csv`:

- N. weaveri: strictly aerobic phenotype in the primary species description (https://pubmed.ncbi.nlm.nih.gov/8240951/).
- N. dumasiana: observed facultative anaerobe entry, separate from model predictions (https://bacdive.dsmz.de/strain/158278; reference 65274). Primary description identified at https://doi.org/10.1099/ijsem.0.002148; its full text was not retrieved here.
- N. shayeganii: curated type-strain microaerophile description (https://bacdive.dsmz.de/strain/10496). Primary phenotypic confirmation remains desirable. Do not replace this with adjacent aerobic genome predictions.
- B. zoohelcum: aerobic description in primary case report (https://pubmed.ncbi.nlm.nih.gov/17660147/); this does not establish uniformity across all strains.
- P. cangingivalis: anaerobic isolates in the primary species description (https://pubmed.ncbi.nlm.nih.gov/7981097/). Preserve `anaerobic_unspecified` until obligate anaerobiosis is explicitly verified.
- C. steedae: aerobic incubation in the culture collection (https://catalogue-crbip.pasteur.fr/fiche_catalogue.xhtml?crbip=CIP+103435). Keep physiology unknown pending stronger evidence.

Periodontal associations must remain separate. A longitudinal dog study reports decreasing N. shayeganii with mild periodontitis (https://pubmed.ncbi.nlm.nih.gov/26507828/), whereas a habitat study reports a positive genus-level association with worse health scores in plaque (https://pubmed.ncbi.nlm.nih.gov/34001282/). These differ by design, sampling site and resolution; neither supplies an oxygen phenotype nor demonstrates that the current age association is mediated by disease.

## Remaining analysis

1. Complete species evidence review, beginning with dominant species and all nine Neisseria features. Resolve synonyms while preserving original feature IDs. Keep predicted traits separate. No genus-wide automatic oxygen assignment.
2. Verify C. steedae physiology and additional dominant taxa; unknown strain features must remain unknown unless evidence can be tied to them.
3. Fit oxygen-group models only after freezing annotations independently of age effects. Show annotation coverage, unresolved abundance, continuous age curves, nonlinear age sensitivity and matched-sample weight adjustment.
4. Run the no-Neisseria, no-Neisseriales and dominant-taxon leave-one-out analyses, using fixed denominators to distinguish taxon removal from renormalization.
5. Rebuild revised figures using the reviewed species labels and direct-rank aggregations. Complete the UMAP/ENet baseline checks before a full release.

## Reproduce

From the repository root with baseline dependencies installed:

```bash
python workflow/oxygen_extension/audit_species.py
python workflow/oxygen_extension/followup_checks.py
```

The first command regenerates an unreviewed queue; the second merges the preserved evidence file into the current annotation table. Outputs include input hashes, an environment snapshot, all-species tests, Neisseria plots, matched-case adjustment checks, denominator audits, annotation coverage and direct-order composition tables. No raw sequencing is rerun.

To compare the existing age-colored dog UMAP with the oxygen balance without
refitting the embedding, run:

```bash
python workflow/oxygen_extension/plot_umap_oxygen_balance.py
```

The resulting two-panel PNG and PDF use the frozen UMAP coordinates. The
oxygen panel colors matched dogs by the expanded oxygen-tolerant:anaerobe
log-ratio; unmatched samples are omitted from that panel.

First refit the bacterial-taxon UMAP on QC3, then render the revised combined
Figure 3, its standalone sample-level Figure 3B, and the separate
species-resolution coverage audit:

```bash
MICROBIOME_OUTPUT=workflow/oxygen_extension/analysis_results \
  python workflow/scripts/06_figure3_bacterial_umap.py
python workflow/oxygen_extension/plot_figure3_species_oxygen.py
```

Figure 3A is the bacterial-taxon UMAP colored by each taxon's Spearman
correlation with host age. Figure 3B shows one point per dog. Its y-axis is the oxygen-tolerant abundance
divided by the sum of reviewed oxygen-tolerant and anaerobic abundance. The
blue curve and band are a median and interquartile range calculated within a
fixed ±0.75-year age window, not a fitted inferential model. Abundance not resolved to species is reported separately in
the coverage figure.

Do not use the archived QC2 bacterial-taxon coordinates for Figure 3A. QC2 is
used only for the initial dog-level UMAP that identifies the excluded island
and DBSCAN noise; downstream analyses use QC3 (or QC4 where explicitly stated).

For sample-specific views, run:

```bash
python workflow/oxygen_extension/plot_sample_specific_oxygen.py
```

This creates an individual-dog oxygen-balance plot and a heatmap of the 12
most abundant reviewed oxygen-tolerant species across dogs ordered by age.
