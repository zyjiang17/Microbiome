# Canine oral HUMAnN pathway analysis: methods, results, and cluster deployment

**Prepared:** October 2026
**Status:** Preliminary and exploratory. No clinical periodontal scores, raw read counts, or batch information are available in these inputs.
**Context:** The existing manuscript uses MetaPhlAn to describe age-associated shifts from oxygen-tolerant taxa (*Conchiformibius*, *Neisseria*) toward taxa often associated with canine periodontal disease (*Porphyromonas*, *Tannerella*, *Treponema*). HUMAnN adds *predicted functional potential*; it is not a measurement of pathway expression/activity or disease stage.

## A. Research questions

1. Is canine oral microbial pathway abundance associated with chronological age?
2. Do pathway profiles show a functional age gradient independently of individual taxonomic abundance trends?
3. Are effects sensitive to compositional transformation, exclusion of puppies, or adjustment for weight?
4. Are important age-associated pathways well supported by HUMAnN *coverage*?
5. Which species are assigned HUMAnN pathway abundance for mechanistically relevant pathways?
6. Once clinical metadata exist: which pathways associate with *measured* periodontal disease after controlling for age, weight/body size, breed, depth, and batch?

## B. Data provenance and reproducible preprocessing

Source: Madeline Shono HUMAnN3 outputs shared October 15, 2026. Two TSV files are included under `inputs/`. The tables are CPM-normalized community and species-stratified MetaCyc pathways, plus HUMAnN coverage scores. `UNMAPPED`/`UNINTEGRATED` were excluded from abundance but present in coverage. Keep these files frozen: do **not** rerun HUMAnN or replace original normalizations in this pipeline. A separate original taxonomic workflow exists at `/u/home/z/zyjiang/canine_oral_microbiome_repro/pipeline` and must not be overwritten.

- Match coverage rows to abundance sample IDs and verify metadata equality.
- Keep *unstratified* pathways (colon-containing MetaCyc IDs without `|`), excluding `UN...` records. Do not combine species-stratified with unstratified columns in totals or tests.
- Require observed CPM >0 in >=10% of samples. This retains **309 / 458** pathways.
- Age is recorded for **695 / 695** dogs; weight is present for **626 / 695**. Adults (age >=1 year) comprise **533** dogs.
- Age-group bins: <1 year **162**, 1-<3 **184**, 3-<7 **184**, >=7 **165**.
- Pathway CPM row totals are approximately one million; this closed composition motivates CLR sensitivity. Note: there are **695** unique sample IDs, whereas Madeline reported **625** after the >=2-million-read threshold. No raw read-count manifest is available for verification. The cohort discrepancy is a **publication blocker**.

## C. Per-pathway statistical workflow

**Primary exploratory screen:** Compute Spearman rho between each pathway's unstratified CPM and dog age. Apply Benjamini-Hochberg adjustment over all **309** tested pathways; significance threshold q<0.05. Use Spearman because age trends and abundance distributions are non-Gaussian and may be nonlinear. Zero CPMs represent lack of detected pathway evidence and produce tied ranks.

**Sensitivity 1 (composition):** Apply centered log-ratio to the **309-feature** filtered matrix after a +1 CPM pseudocount: z_ij = log(CPM_ij+1) - mean_k[log(CPM_ik+1)]. Compute age Spearman rho and BH q independently. This is a preliminary diagnostic, not a validated zero-aware compositional method. Pseudocount and feature filtering can alter the results.

**Sensitivity 2 (maturation):** Restrict to age>=1, n=533, and repeat raw CPM and CLR association tests. This distinguishes puppy-to-adult changes from associations across adulthood.

**Sensitivity 3 (weight):** For n=626 samples with weight, fit separate heteroskedasticity-robust HC3 OLS models, log1p(CPM) ~ age + log1p(weight), and CLR ~ age + log1p(weight), separately for each pathway. Apply BH correction to age coefficients. These models treat chronological age linearly and do not adjust for breed/body size, sequencing depth, batch, or periodontal status.

**Robustness heuristic:** A pathway is deemed 'supported_all_four' if its age effect is BH-significant and has the same direction in (i) raw CPM Spearman, (ii) CLR Spearman, (iii) adults-only CPM Spearman, and (iv) weight-adjusted log1p CPM regression. This is a *descriptive screening heuristic* and not statistical replication.

## D. Functional-community and pathway-level analyses

**Ordination:** PCA on CLR-transformed filtered pathways; inspect PC1/PC2 colored by age and calculate Spearman(PC1,age). PCA is unsupervised and PC-axis sign is arbitrary (the sign of PC1 might flip in another numerical environment).

**Pathway detection:** Count positive CPM features across all 458 unstratified pathways per dog and correlate with age. This is only a proxy for detected repertoire and may be confounded by library size.

**Coverage audit:** For each MetaCyc pathway compute mean HUMAnN coverage, fraction of dogs with coverage >=0.5 and >=0.9, and check discordance between nonzero abundance and zero coverage. HUMAnN coverage does not equal a read-depth-normalized abundance. We selected a provisional high-confidence shortlist by both concordant effects and coverage >=0.5 in >=30% of dogs; this is **not a validated threshold**.

**Species attribution:** For selected pathways, identify abundance columns of form `MetaCyc pathway|g__...s__...`. Sum each species' abundance across dogs. Report both (a) share among named species-attributed contributions and (b) the fraction of unstratified community total explained by those species. The first quantity must **not** be mistaken for the overall community fraction.

**Functional grouping:** Provisional categories based on keyword matching in pathway names are explicitly **not** official MetaCyc classes. A separate R script exports the `file2meco` curated hierarchy for auditing; superclass-level inference waits on this step.

**Weight:** As a secondary analysis, conduct raw/CLR Spearman tests for weight (n=626), and age-adjusted HC3 OLS. No interpretation of disease effects is permitted because clinical status is absent.

## E. Main quantitative results

| Analysis | Significant pathways (BH q < 0.05) |
|---|---:|
| Age vs raw CPM (n=695) | **257 / 309** |
| Age vs CLR abundance (n=695) | **242 / 309** |
| Age vs raw CPM (adults n=533) | **219 / 309** |
| Age vs CLR abundance (adults n=533) | **207 / 309** |
| Age coefficient adjusted for weight (n=626) | **222 / 309** |
| Age coefficient in CLR adjusted for weight (n=626) | **225 / 309** |
| Concordant/significant across four primary screens | **136 / 309** |
| Concordant + HUMAnN coverage >=0.5 in >=30% dogs | **22 pathways** (8 increasing, 14 decreasing) |

**Additional results:** PC1 of CLR pathway profiles explains **14.9%** of variance and correlates with age (Spearman rho **+0.574**, noting arbitrary PCA orientation). Pathway detection count correlates with age (**+0.324**). Age and weight are modestly correlated (rho **+0.131**, p=0.001). Weight associations are widespread too, but potentially confounded by breed, body size, age, depth and technical factors.

### Selected biologically relevant pathway associations

| MetaCyc pathway | Spearman rho(age, CPM) | BH q | Fraction dogs with coverage >=0.5 | Assessment |
|---|---:|---:|---:|---|
| HISDEG-PWY: L-histidine degradation I | **+0.498** | **8.99e-44** | **50.2%** | Priority mechanistic lead |
| PWY-5030: L-histidine degradation III | +0.521 | 2.35e-48 | 12.9% | Age signal, less coverage |
| P164-PWY: anaerobic purine nucleobase degradation I | +0.515 | 2.99e-47 | ~0% | Interpret cautiously |
| P163-PWY: lysine fermentation to acetate/butanoate | +0.400 | 1.86e-27 | 0% | Do not assert pathway completeness |
| DENITRIFICATION-PWY: denitrification | +0.326 | 3.02e-18 | 0% | Do not assert pathway completeness |
| FUCCAT-PWY: fucose degradation | -0.530 | 2.73e-50 | 6.0% | Lower-coverage signal |
| FERMENTATION-PWY: mixed acid fermentation | -0.493 | 1.07e-42 | Check QC table | Shows heterogeneous directionality |
| PWY-7663: gondoate biosynthesis (anaerobic) | -0.469 | 2.47e-38 | 95.0% | Strong counterexample to global anaerobic claim |
| PWY-7791: UMP biosynthesis III | +0.538 | 8.81e-52 | 73.4% | Well covered; not disease-specific |
| VALSYN-PWY: L-valine biosynthesis | -0.524 | 6.41e-49 | 73.7% | Well covered; younger-associated |

### Species-associated histidine I results

For HISDEG-PWY, the sum of stratified species contributions equals **67.9%** of summed unstratified pathway abundance in the cohort. Of that attributed fraction, *Porphyromonas* species account for **69.8%**; *Tannerella forsythia* for **4.5%**. Largest assigned species include *P. gulae* (29.5%), *P. cangingivalis* (24.7%), *Bergeyella zoohelcum* (17.3%), and *P. gingivicanis* (12.2%). These are **percentages of taxonomically attributed pathway CPM**, not percentages of all community pathway abundance. Species assignment does not demonstrate expression or disease causation.

### Important negative/limiting results

**181 of 458** pathways have coverage zero in every sample. **135 of 309** tested pathways never reach 0.5 coverage; **114** of these nonetheless have significant raw age associations. Thus very small p values alone are insufficient for confidently identifying complete pathways.

Mixed-acid fermentation and gondoate biosynthesis decrease with age even though the study's older-age taxonomic community is enriched for some anaerobe-associated taxa. It would be misleading to claim that all anaerobic metabolism increases with age.

## F. Figure atlas: 12 exploratory panels

The `reference_figures/` folder contains all 12 plots in PNG; `reports/FIGURE_CAPTIONS.md` links each plot with captions. First inspect **Figure 10** (well-supported pathway curves) and **Figure 11** (age-group pathway heatmap); then the age-volcano (Figure 02), raw/CLR robustness (Figure 03), and attribution CSV. Figure 05 uses unvalidated keyword classes and should **not** be treated as curated MetaCyc enrichment.

## G. Manuscript figure and study integration

Recommend testing a new **HUMAnN functional supplement or main figure**, combining: 1) pathways ranked by age correlation with coverage distinction, 2) well-supported representative abundance-age smooths, 3) overlap of functional and MetaPhlAn taxa using **exact sample-ID matching**, and 4) pathway-coverage diagnostics. Preserve the prior taxonomic plots. A CAZy replacement decision requires matched-sample pathway/CAZy comparison, not simply smaller p-values. No formal periodontal-disease models are possible with these data.

## H. Remaining QC and follow-up tasks

1. Madeline: resolve **695 samples in both TSV files vs 625 stated in her October 15 email** using >=2-million-read sample-level manifest. Rebuild analysis if any sample is ineligible.
2. Request library sizes/read counts, plates/batches, periodontal disease outcomes if they exist, and breed/sex metadata; reconcile against taxonomic QC1-QC4 sample sets in the existing pipeline.
3. Install/use `file2meco` to assign official MetaCyc hierarchy; report match rates and avoid double-counting multi-parent pathways.
4. Reconsider zero-aware compositional models and correlated pathway definitions; test age splines/GAM and adult age effects beyond the current Spearman screen.
5. Check pathway coverage and UniRef-to-KO/eggNOG gene evidence for HISDEG-PWY. Avoid mechanistic language about activity until independently supported.
6. Compare CAZy to HUMAnN on exactly matched dogs with the same age/weight covariates and a defined predictive or explanatory endpoint.
7. If prioritizing age prediction, benchmark MetaPhlAn-only vs HUMAnN-only vs combined features in nested CV, with no leakage from sample selection or preprocessing.

## I. Runtime and verification

- This bundle includes source TSVs in `inputs/`, versioned code in `code/`, the original reproducibility outputs in `reference_results/` and `reference_figures/`, and a Slurm job in `submit_humann.slurm`.
- User transfer is via the provided `deployment/upload_from_mac.sh` and local SSH authentication as `zyjiang` (no credentials stored).
- `sbatch submit_humann.slurm` reruns all 12 figures and outputs under `cluster_runs/slurm_JOBID` without modifying the existing taxonomic workflow. It also creates `preflight.json`, input SHA256 hashes, and a reference comparison file.
- See `README_DEPLOY.md` for transfer commands and credentials/security boundaries.
