# Oxygen analysis results

## Main finding

Across the 1,043-dog QC3 cohort, evidence-annotated oxygen-tolerant species decrease with age while anaerobic species increase. The pattern persists in the 675-dog age >=1 sensitivity cohort, after weight adjustment on complete cases, after removing Neisseria, and after removing all annotated Neisseriales.

| Annotation tier | Oxygen-tolerant rho | Anaerobic rho | Oxygen-tolerant/anaerobic log-ratio rho |
|---|---:|---:|---:|
| Strict species evidence | -0.289 | 0.576 | -0.572 |
| Expanded evidence | -0.534 | 0.620 | -0.616 |
| Broad evidence | -0.534 | 0.607 | -0.610 |

All correlations above use QC3 and have BH-adjusted p < 3e-21 within tier. In QC4, corresponding expanded-tier correlations are -0.504, 0.517, and -0.523. Weight-adjusted rank regressions retain the same directions in both cohorts. Exact results are in `oxygen_group_age_associations.csv`.

Expanded-tier mean original relative abundance changes across age groups:

| Group | Puppy (n=636) | Adult (n=272) | Mature (n=135) |
|---|---:|---:|---:|
| All oxygen-tolerant phenotypes | 67.5% | 44.4% | 29.1% |
| Anaerobe | 14.3% | 38.2% | 53.8% |
| Unassigned direct species | 9.8% | 12.4% | 14.0% |
| Unresolved below species | 8.4% | 5.0% | 3.6% |

These percentages use the original MetaPhlAn scale. They do not renormalize over annotated taxa.

## Robustness

The broad oxygen-tolerant association remains negative when important groups are removed without renormalization:

| Scenario | Spearman rho with age |
|---|---:|
| All broad oxygen-tolerant taxa | -0.534 |
| Without Neisseria | -0.556 |
| Without all Neisseriales | -0.303 |
| Without Conchiformibius steedae | -0.299 |
| Without Canicola haemoglobinophilus | -0.354 |

All listed leave-out associations have BH-adjusted p < 2e-23. The trend therefore extends beyond Neisseria and beyond any one dominant oxygen-tolerant species, although Neisseriales and C. steedae contribute materially to its magnitude.

## Important interpretation

The strict high-confidence categories do not support a uniform decline of organisms explicitly labeled aerobe or obligate aerobe. In QC3:

- strict `aerobe`: rho = 0.075;
- strict `obligate_aerobe` (currently N. weaveri): rho = 0.461;
- strict `microaerophile`: rho = -0.518;
- strict combined `oxygen-tolerant`: rho = -0.289.

The expanded `aerobe` category declines (rho = -0.434), but it includes moderate-confidence calls, especially C. steedae and unresolved Moraxella. Therefore the most defensible current claim is an **age-associated decline in oxygen-tolerant taxa accompanied by an increase in anaerobic taxa**. Retain “aerobic decline” only if the C. steedae and Moraxella assignments are upgraded and the manuscript defines aerobic broadly enough to include the intended taxa.

Neisseria is internally heterogeneous: N. weaveri increases, N. shayeganii increases in the full cohort but not the age >=1 cohort, N. dumasiana declines in the age >=1 cohort, and N. canis is approximately stable. A genus-level Neisseria decline should not be used as the mechanism for the community result.

## Scope

This is a cross-sectional relative-abundance analysis. It supports age-associated community composition differences, not within-dog change, absolute bacterial loss, oxygen measurements, or periodontal-disease mediation. Periodontal status was unavailable.
