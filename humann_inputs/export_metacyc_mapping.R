# Load the file2meco package and get the data (pathway map)
library(file2meco)
data("MetaCyc_pathway_map")

# What are the dimensions of the pathway_map file?
dim(MetaCyc_pathway_map)
# Column names of the mapping file
colnames(MetaCyc_pathway_map)
# Head of pathway file - unsure how the names of the pathways are formatted
head(MetaCyc_pathway_map)

# Make a column out of the row names (pathway_ids)
MetaCyc_pathway_map$pathway_id <- rownames(MetaCyc_pathway_map)

# Write the pathway map to a csv
write.csv(MetaCyc_pathway_map, "metacyc_pathway_map.csv", row.names = FALSE)