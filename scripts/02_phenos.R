library("tidyverse")
# Load and format the phenotype table from manually made csv from provided docx file
pheno.raw <- read.csv("data/raw/manually_made_pheno.csv", header = F)

pheno.split <- strsplit(pheno.raw$V2, " ")

pheno.clean <- data.frame(
  sample.id = pheno.raw$V1,
  NRVM = lapply(pheno.split, "[[", 1) |> unlist(),
  treatment = lapply(pheno.split, "[[", 2) |> unlist(),
  treatment.2 = lapply(pheno.split, "[[", 3) |> unlist()) |> 
  mutate(treatment = case_when(
    treatment == "Ctl" ~ treatment,
    .default = paste(treatment, treatment.2)
  )) |> 
  select(-treatment.2)


