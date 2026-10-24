#!/usr/bin/env Rscript
# Run where the R packages file2meco and microeco are installed.
# Converts file2meco's curated MetaCyc taxonomy into a CSV suitable for manual
# validation and subsequent joining by MetaCyc pathway ID.
suppressPackageStartupMessages(library(file2meco))
data('MetaCyc_pathway_map', package='file2meco')
map <- as.data.frame(MetaCyc_pathway_map)
map$metacyc_id <- rownames(map)
map <- map[,c('metacyc_id',setdiff(names(map),'metacyc_id')),drop=FALSE]
write.csv(map, file='results/file2meco_curated_MetaCyc_pathway_map.csv', row.names=FALSE)
cat('Exported',nrow(map),'MetaCyc mapped pathways.\n')

results <- read.csv('results/pathway_age_associations_sensitivity.csv', stringsAsFactors=FALSE)
results$metacyc_id <- sub(':.*$','',results$pathway)
matched <- sum(results$metacyc_id %in% map$metacyc_id)
cat(matched,'/',nrow(results),'tested pathway IDs match the hierarchy; audit the remainder.\n')
joined <- merge(results, map, by='metacyc_id', all.x=TRUE, sort=FALSE)
write.csv(joined, 'results/pathway_age_with_curated_MetaCyc_hierarchy.csv', row.names=FALSE)
# Inspect actual map column names before automatic superclass aggregation, and
# account for multi-parent pathways to prevent duplicative category totals.
print(names(map))
