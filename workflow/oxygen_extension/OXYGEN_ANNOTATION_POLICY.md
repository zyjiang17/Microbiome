# Oxygen annotation policy

## Recommended taxonomic resolution

The canonical annotation table is species-feature based. A label may be inherited from a genus only when the feature is unresolved below genus and the represented genus has a sufficiently consistent phenotype in observed records. Genus-inherited labels receive moderate or low confidence. Family- and order-level labels are not used as physiological annotations.

This is necessary because both between-species and within-species variation occur. The current data contain *Neisseria canis* (aerobe), *N. weaveri* (obligate aerobe), *N. dumasiana* (facultative anaerobe), and *N. shayeganii* (microaerophile). BacDive also reports aerobic, facultative-anaerobic, and microaerophilic records for different *N. animaloris* strains. A single Neisseria or Neisseriales oxygen label would erase observed biology.

The coarsest defensible primary resolution is therefore species. Genus is the coarsest permitted fallback for unresolved features, after checking that genus. Order is retained only for ecological visualization or a legacy sensitivity analysis and should be called `aerobe-associated`, `anaerobe-associated`, `mixed`, or `unreviewed`, rather than an organismal oxygen phenotype.

## Phenotype ontology

Keep these fields separate in the canonical table:

| Phenotype | Meaning |
|---|---|
| `obligate_aerobe` | Explicit evidence of required aerobic growth |
| `aerobe` | Reported aerobic phenotype without proof of obligacy |
| `facultative_anaerobe` | Growth with or without oxygen |
| `microaerophile` | Prefers reduced oxygen |
| `anaerobe` | Reported anaerobic phenotype; do not silently upgrade to obligate |
| `variable_*` | Multiple observed phenotypes across strains or sources |
| `unknown` | No adequate observed evidence |

For the broad two-group sensitivity analysis only, aerobe, obligate aerobe, facultative anaerobe, microaerophile, and their oxygen-tolerant variable combinations map to `oxygen_tolerant`; anaerobe maps to `anaerobic`. Unknown and unresolved variable records are excluded from the numerator but remain in the full-community denominator.

The manuscript's main result should retain the multicategory physiology. Calling all oxygen-tolerant organisms “aerobic” would conflate facultative anaerobes and microaerophiles with obligate aerobes. If the broad group drives the result, use “oxygen-tolerant taxa.” Use “aerobic taxa” only for an analysis restricted to aerobe and obligate-aerobe labels.

## Evidence hierarchy

1. Primary species description with an explicit oxygen phenotype.
2. Observed phenotype for the type strain or multiple curated strains in BacDive/culture collections.
3. Observed phenotype for another strain of the same species.
4. Reviewed species reference or consistent genus-level evidence for an unresolved feature.
5. Genome prediction, order/family convention, or culture atmosphere alone: hypothesis-generating only and excluded from the strict analysis.

Conflicting observed records are preserved as `variable`; they are not resolved by majority vote. Predictions are stored separately and never overwrite observed phenotypes.

## Analysis sets

- **Strict:** high-confidence, species-level observed evidence. Current coverage: 50.1% of mean direct-species abundance.
- **Expanded:** high- and moderate-confidence annotations, including vetted genus fallback. Low-confidence genus inheritance can be added as a separate sensitivity. Current broad base coverage including low-confidence calls: 92.9%.
- **Legacy:** original order labels, used only to show whether the prior conclusion changes under the old mapping.

Report unknown/unreviewed abundance and annotation coverage for every dog. Use original MetaPhlAn relative-abundance units for the main test and repeat with normalization over all direct species. Do not normalize over annotated taxa alone.

## Current base annotation

`base_oxygen_annotation_v1.csv` contains the initial evidence-backed calls for the dominant features. `outputs/base_oxygen_annotation_all_species.csv` preserves all 528 features, including unreviewed ones. The table is deliberately conservative and should be extended in descending abundance order before the final community model.
