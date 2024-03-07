# This script is part of the initial processing pipeline for Jensen's NRVM RNAseq data
# In this script, we will compile the quant.sf files into a single count matrix
# Then, counts will be converted from transcript- to gene-leve anotations with tximport

# Load libs
#BiocManager::install("tximeta")
#BiocManager::install("BiocFileCache", version = "3.18")
libs <- c("tximport", "tximeta","BiocFileCache", "tidyverse", "SummarizedExperiment") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

### Use tximeta to make a linked transcriptome
# https://bioconductor.org/packages/release/bioc/vignettes/tximeta/inst/doc/tximeta.html

# Make linkedTxome for first run 
indexDir <- file.path("data/raw/anno/Rattus_norvegicus.mRatBN7.2.salmon")
fastaFTP <- "https://ftp.ensembl.org/pub/release-111/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz"
gtfPath  <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M34/gencode.vM34.annotation.gtf.gz"
suppressPackageStartupMessages(library(tximeta))
makeLinkedTxome(indexDir=indexDir,
                source="GENCODE",
                organism="Mus musculus",
                release="M34",
                genome="GRCm39",
                fasta=fastaFTP,
                gtf=gtfPath,
                write=FALSE)

# List sample quant.sf files
names <- list.files("data/raw/fastq")
names <- names[str_detect(names, "B6_")]
files <- file.path("data/raw/fastq", names, "quant.sf")

# Check that they exist
all(file.exists(files))

# Summarize gene-level expression
coldata <- data.frame(files = files, names = names)
se <- tximeta(coldata)
gse <- summarizeToGene(se)

# Save to processed data
save(gse, file="data/processed/bulk/rau_fractions_gse.RData")
write.csv(assay(gse),"data/processed/bulk/rau_fractions_ensembl.csv")