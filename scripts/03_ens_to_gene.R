# Reworking 4_clean_bulk to take input from 0_transcripts to genes
libs <- c("tidyverse", "biomaRt") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

# Load in bulk fractions
counts <- read.csv("data/processed/bulk/bulk_ensembl.csv", row.names = 1)

# Load the Ensembl dataset for mouse genes
ensembl <- useMart("ensembl")
ensembl <- useDataset("rnorvegicus_gene_ensembl",mart=ensembl)
annot <- getBM(attributes=c('external_gene_name','ensembl_gene_id'), mart = ensembl)
# Extract the Ensembl gene IDs from the row names
ensembl_ids <- rownames(counts)

# Join the gene symbols back to the counts data frame
# First, make sure the Ensembl IDs are a column in the counts data frame for merging
counts$ensembl_gene_id <- rownames(counts)

# Merge to add gene symbols
counts_with_symbols <- left_join(counts, annot, by = 'ensembl_gene_id')

# For each gene, keep the most abundant version
summarized_counts <- counts_with_symbols %>%
  #group_by(ensembl_gene_id_version) %>% 
  mutate(total = rowSums(dplyr::select(., where(is.numeric)), na.rm = TRUE)) |> 
  group_by(external_gene_name) |> 
  arrange(desc(total)) |> 
  filter(row_number() == 1) |> 
  ungroup() |> 
  subset(!is.na(external_gene_name)) |> 
  as.data.frame()

# If you wish to replace row names with gene symbols, where available
# Replace Ensembl IDs with gene symbols where available, keeping Ensembl ID where not
rownames(summarized_counts) <- summarized_counts$external_gene_name

summarized_counts <- summarized_counts |> dplyr::select(c(-total, -ensembl_gene_id, -external_gene_name))

write.csv(summarized_counts, "data/processed/bulk/bulk_gene.csv")

