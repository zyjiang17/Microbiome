# Exploratory HUMAnN pathway results: aging and the canine oral microbiome
**Analysis date:** October 2026 | **Status:** exploratory, hypothesis-generating; not for direct publication without further QC.

## Research framing and meeting alignment

The manuscript (*Age-associated decline of aerobic taxa in the canine oral microbiome*) emphasizes age-associated increases in *Porphyromonas*, *Tannerella* and other anaerobe-associated oral taxa, with declines in aerobic/oxygen-tolerant genera such as *Conchiformibius* and *Neisseria*. Historical Pellegrini lab meeting notes describe comparing continuous age and weight associations, nonlinear trends, log-transformed ordinations, species-level and higher-level taxonomic signals, and group-level shifts. HUMAnN addresses the complementary question: do encoded microbial **metabolic pathways** vary along this age gradient, and which taxa may contribute?

The October 2026 recorded meeting has **not** been transcribed/verified here. The priorities above are from the existing manuscript, earlier recorded project discussions and the lab's written meeting notes, **not claimed as quotations or decisions from the recent audio**.

## Inputs and essential cohort QC

- HUMAnN pathabundance: 695 unique samples × 458 community-level MetaCyc pathway columns, with 9,733 taxon-stratified pathway columns and age/weight metadata.
- HUMAnN pathcoverage: identical sample IDs; 458 matching MetaCyc pathways plus `UNMAPPED`, `UNINTEGRATED`, and species-stratified columns.
- Age recorded for all 695; weight recorded for 626. Age ranges from 0.019 to 17 (presumably years); units of weight are not explicitly documented in supplied TSVs.
- Community-level pathway CPM row sums are approximately 1 million as expected.
- **Unresolved: Madeline's October 15 email says `n=625`, whereas BOTH actual TSV files have `n=695`**. The two-million-read cutoff cannot be independently verified without read counts. The 695 samples should not automatically be equated with the cohort in the manuscript's earlier figures.
- No periodontal disease diagnosis, tooth exam, disease stage, sequencing depth or batch column is included.

## Per-pathway age statistics and sensitivity analysis

309 of 458 pathways appeared at nonzero CPM in ≥10% of samples. Age Spearman correlations were adjusted across 309 tests using the Benjamini–Hochberg FDR procedure.

| Analysis | FDR < .05 |
|---|---:|
| Raw community CPM vs age, 695 dogs | **257 / 309** |
| CLR(log(CPM+1)) vs age, 695 dogs | **242 / 309** |
| Raw CPM vs age, dogs ≥1 year, n=533 | **219 / 309** |
| CLR age associations, dogs ≥1 year, n=533 | **207 / 309** |
| Robust HC3 OLS log(CPM+1) ~ age + log(weight+1), n=626 | **222 / 309** |
| Same HC3 model using CLR outcome, n=626 | **225 / 309** |
| **FDR-significant same-direction in raw/CLR/all ages/adult-only/weight-adjusted raw age coefficient** | **136 / 309** |

Very widespread significance can reflect genuine biological gradients, but may also reflect highly correlated pathway definitions, compositional closure, sequencing-depth differences, sample/batch effects, and unmeasured variables. It does **not** imply that 257 independent metabolic systems change with age. The CLR +1 CPM pseudocount is a sensitivity check, not a replacement for a dedicated compositional approach and independent replication.

### Pathways matching the oral ecology hypothesis

| Pathway | Spearman ρ age (CPM) | FDR q | Pathway coverage ≥0.5 in dogs | Interpretation |
|---|---:|---:|---:|---|
| L-histidine degradation I (HISDEG-PWY) | +0.498 | 8.99e-44 | 50.2% | **Priority** amino-acid degradation signal; high attributed abundance from *Porphyromonas* spp. |
| L-histidine degradation III (PWY-5030) | +0.521 | 2.35e-48 | 12.9% | Strong association, lower coverage support; discuss as follow-up candidate |
| Anaerobic purine nucleobase degradation I (P164-PWY) | +0.515 | 2.99e-47 | approximately zero | Association only; insufficient coverage for pathway-completeness claim |
| L-lysine fermentation to acetate/butanoate (P163-PWY) | +0.400 | 1.86e-27 | zero | Low confidence for complete-pathway interpretation |
| Denitrification pathway (DENITRIFICATION-PWY) | +0.326 | 3.02e-18 | zero | Low confidence for complete-pathway interpretation |
| Fucose degradation (FUCCAT-PWY) | -0.530 | 2.73e-50 | check CSV | Younger-age-associated substrate utilization candidate; coverage relatively low |
| Mixed acid fermentation (FERMENTATION-PWY) | -0.493 | 1.07e-42 | check CSV | Counterexample: not all fermentation pathways rise with age |
| Gondoate biosynthesis (anaerobic) (PWY-7663) | -0.469 | 2.47e-38 | 95.0% | Strong *decrease* with age despite 'anaerobic' label; enzyme-level pathway classification is not a direct physiological oxygen phenotype |
| UMP biosynthesis III (PWY-7791) | +0.538 | 8.81e-52 | 73.4% | Well-covered candidate reflecting nucleic-acid precursor biosynthesis; not PD-specific |
| L-valine biosynthesis (VALSYN-PWY) | -0.524 | 6.41e-49 | 73.7% | Well-covered younger-age-associated pathway |

These are **cross-sectional age associations**, not evidence of periodontal disease or microbial pathway activity.

### Species attribution of higher-priority L-histidine degradation I

Across the full cohort, **67.9%** of cumulative community-level HISDEG-PWY abundance was assigned to named microbial species in the HUMAnN stratified output. Of that species-attributed portion, **69.8%** was assigned to *Porphyromonas* species collectively, and approximately **4.5%** to *Tannerella forsythia*. Top individual species fractions of the **attributed** portion were: *P. gulae* 29.5%, *P. cangingivalis* 24.7%, *Bergeyella zoohelcum* 17.3%, *P. gingivicanis* 12.2%, *T. forsythia* 4.5%. These cumulative attribution fractions do not imply direct enzyme expression, species causation, or disease status.

### Coverage-supported pathway prioritization

A provisional shortlist applies **two filters**: (a) concordant significant age association in raw CPM, CLR, age ≥1 subset, and weight-adjusted CPM models, and (b) coverage ≥0.5 in **at least 30% of dogs**. This yields **22 of 309** pathways (**8 age-positive, 14 age-negative**). This is an analysis heuristic, *not a validated universal coverage cutoff*. See `results/coverage_supported_priority_pathways.csv`.

The coverage table is especially important: **181/458 pathways have coverage identically zero across all 695 dogs**, even though pathway abundance may be nonzero. In all, 135/309 tested pathways never reach coverage ≥0.5, and 114 of these still showed raw age FDR <.05. Follow-up should ask why reconstruction of the full pathways is weak despite inferred abundance rather than treating `pathcoverage` as a normalized version of abundance.

## Pathway categories / MetaCyc

Provisional broad classes were generated using transparent, manually specified **keyword matching**, not MetaCyc hierarchy (see `results/pathway_category_keyword_audit.csv`). They are useful for organizing candidate pathways but not yet for an enrichment or formal superclass figure. A naive aggregate 'Fermentation and anaerobic energy' decreases with age, demonstrating how heterogeneous grouping can obscure individual mechanisms. Do **not** equate this class with anaerobic microbial abundance or an independently measured anaerobic function score.

`export_official_metacyc_mapping.R` is included to join the curated file2meco MetaCyc superclass annotation after package installation. A pathway may have multiple parents, and hierarchy rollups must avoid double-counting.

## Secondary weight associations

Among 626 dogs with weight, Spearman correlation detected 207/309 pathways associated with weight (raw CPM), versus 186/309 for CLR. HC3 log(CPM+1) ~ age + log(weight+1) retained 187/309 significant weight coefficients (181/309 for CLR). Age-weight Spearman ρ=+0.131 (p=0.001). These patterns require body-size/breed and read-depth adjustment before interpretation. They should not be compared naively with 695-dog age tests as though cohorts and confounding were identical.

## Functional ordination and pathway richness

- CLR (log(CPM+1), mean-centered) pathway PCA: PC1 explains approximately **14.9%** of profile variance; Spearman ρ(PC1,age) = **+0.574**.
- Number of detected community pathways per dog vs age: ρ = **+0.324**. This association could reflect read depth rather than biological functional diversity because per-sample raw read counts were unavailable.
- The sample distribution across age groups: <1 yr n=162; 1–<3 yr n=184; 3–<7 yr n=184; ≥7 yr n=165.

## Relation to previous CAZy analysis

**Recommendation:** Retain CAZy as carbohydrate-specific context and test HUMAnN as the broader functional analysis in the revised manuscript. Strong age associations and the *Porphyromonas*-linked histidine signal support this route. However, evidence is **insufficient** to state HUMAnN outperforms CAZy: no matched-sample CAZy comparison, independent outcome, or external validation is included in these files. Stronger p-values alone should never determine which analysis is preferred.

## Pending before manuscript integration

1. Ask Madeline for the **read-count/sample exclusion manifest** and reconcile 695 vs email 625.
2. Determine clinical periodontal severity availability; existing manuscript says no clinical disease status was collected. Do not run or claim disease association models without it.
3. Obtain MetaCyc hierarchy from `file2meco` and audit coverage and multi-parent pathway mappings; replace all provisional keyword categories.
4. Compare histidine-associated pathway profiles and contributor fractions on the exact same sample set as the MetaPhlAn taxonomic analysis, including oxygen-preference group definitions.
5. Validate abundance signals versus coverage, read depth, batch, and alternate transforms; consider examining UniRef/KO/eggNOG genes encoding pathway reactions when Madeline releases gene-family output.
6. Re-run relevant CAZy associations on the **same dogs**, age ranges, statistical thresholds and confounders to compare explanatory value fairly.
7. If pathway findings will be used to predict age, use nested held-out cross-validation and a properly matched taxonomic baseline; do not treat simple correlations as predictive accuracy.

## Draft Results paragraph (needs QC confirmation)

> Exploratory functional profiling of the canine oral microbiome identified 309 MetaCyc pathways detected in at least 10% of 695 samples. Of these, 257 exhibited age-associated abundance differences after false-discovery-rate adjustment, and 136 retained concordant significant associations across a series of analyses incorporating log-ratio transformation, exclusion of dogs younger than one year, and adjustment for body weight. L-histidine degradation I increased with age (Spearman ρ=0.498, FDR q=8.99×10⁻⁴⁴), with the majority of its taxonomically attributed abundance assigned to *Porphyromonas* species. However, several apparently age-associated anaerobic pathways showed little HUMAnN pathway coverage, underscoring the need for pathway-level validation. These results suggest age-associated functional restructuring that may complement previously identified taxonomic shifts; they do not establish a relationship to clinically diagnosed periodontal disease.

## Sources

- Data: attached HUMAnN abundance and coverage TSV tables supplied by Madeline Shono (Oct 15 2026); <https://github.com/zyjiang17/Microbiome/tree/main/humann_inputs>
- Existing draft: <https://docs.google.com/document/d/1L2C4WdQHHkFEG0KljJMz6dDVaQ9A29WpDvJpbdZA9-U/edit>
- Historical lab meeting notes: <https://docs.google.com/document/d/1rJNoHSYQr0bTUGrjPbMedOe8JytIlYjAE__s7MFmTPs/edit>
- file2meco: <https://github.com/ChiLiubio/file2meco>
- HUMAnN: <https://github.com/biobakery/humann>
