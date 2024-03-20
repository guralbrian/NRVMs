# Test for differential gene expression with DESeq2 
# Include composition covariates to test if they mediate DE

# Load libs
libs <- c("tidyverse", "stats", "DESeq2") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

#### Loading and formatting of data ####
# Load expression and phenotype data
phenotypes <- read.csv("data/processed/bulk/phenotypes")
bulk <- read.csv("data/processed/bulk/bulk_gene.csv", row.names = 1, check.names = F)

# Subset to whole samples

# Prepare for DESeq2
bulk <- mutate_all(bulk, function(x) round(as.numeric(as.character(x)), digits = 0)) # round to integers

# Set reference factor levels for phenotypes
phenotypes <- phenotypes |> 
  mutate(sample.id = as.factor(paste0("BCJ_",sample.id)),
         NRVM = as.factor(NRVM),
         treatment = as.factor(treatment))

# Prepare sample information
sample_info <- data.frame(
  row.names = phenotypes$sample.id,
  NRVM = phenotypes$NRVM,
  treatment = phenotypes$treatment
)

sample_info$treatment <- relevel(sample_info$treatment, ref = "Ctl")
sample_info <- sample_info[colnames(bulk),]
#### Run DESeq2 ####

## DESeq without compositions
# Create a DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = bulk,
  colData = sample_info,
  design = ~ NRVM + treatment
)

# Run DESeq 
dds <- DESeq(dds)

# Save the results 
saveRDS(dds, "data/processed/deseq/deseq.RDS")
names <- resultsNames(dds)[-c(1:3)]

# make a list of each result
dds <- lapply(names, function(x){
  results(dds, name=x) |> 
    as.data.frame()
})

names_clean <- lapply(str_split(names, "treatment_"), "[[", 2) |> unlist()


for(i in 1:length(go.res)){
  write.csv(dds[[i]], paste0("results/04_deseq/", names_clean[[i]], "_DESeq.csv"))
}
