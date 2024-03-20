# ClusterProfiler analysis
libs <- c("tidyverse", "clusterProfiler", "DESeq2", 
          "org.Rn.eg.db", "viridis", "stringr", "patchwork") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)


# Load the results and expression matrices
res <- readRDS("data/processed/deseq/deseq.RDS")  
names <- resultsNames(res)[-c(1:3)]
names_clean <- lapply(str_split(names, "treatment_"), "[[", 2) |> unlist()
# make a list of each result
res <- lapply(names, function(x){
  results(res, name=x) |> 
    as.data.frame()
})

# Function to perform GO enrichment with clusterProfiler
runGo <- function(data, onto){
# Get a list of significant genes
sig.genes <- data |> 
  filter(padj < 0.05 ) |> #& abs(log2FoldChange) >= 0.263) |> 
  row.names()
# Convert common gene names to ENSEMBLE IDs for clusterProfiler
gene.df <- bitr(sig.genes, fromType = "SYMBOL",
                toType = c("ENSEMBL"),
                OrgDb = org.Rn.eg.db)

# Check for enrichment within biological process gene clusters
ego <- enrichGO(gene          = gene.df$ENSEMBL,
                OrgDb         = org.Rn.eg.db,
                keyType = "ENSEMBL",
                ont           = onto,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE) |> 
        clusterProfiler::simplify(cutoff = 0.6)
}

# Apply the function
go.res   <- lapply(res, function(x){runGo(x, "BP")})

# Save GO Results

for(i in 1:length(go.res)){
  write.csv(go.res[[i]]@result, paste0("results/07_clusterProfiler/", names_clean[[i]], "_GO.csv"),
            row.names = F)
}

plotGO <- function(x, title){
df <- x |> 
  as.data.frame() |> 
  mutate(qscore = -log(p.adjust, base=10),
         desc.wrap = str_to_title(Description) |> 
                    str_wrap(width = 18) |> 
                     factor()) |> 
  arrange(desc(qscore)) |> 
  slice_head(n = 10) 
df$desc.wrap <- factor(df$desc.wrap, levels = rev(df$desc.wrap))

p.ego <- ggplot(df, aes(x = desc.wrap, y = qscore, fill = p.adjust)) +
    geom_bar(stat = "identity", color = "black") +
    coord_flip() +
    scale_fill_viridis(direction = -1) +
    theme(
      axis.text.y = element_text(hjust = 0.5, size = 12),
      axis.title.y = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_line(color = "darkgray"),
      legend.position = "bottom",
      legend.key.width=unit(1,"in"),
      title = element_text(size = 20)
    ) +
    labs(fill = "Adjusted\np-value",
         y = "Q Score",
         title = title)
}

# Make title lists
titles <- c("100nM A6", "1uM ISO", "1uM NE", "1uM PE", "50uM PE")

p.res <- lapply(1:length(go.res), function(n){plotGO(go.res[[n]], titles[[n]])})

# Save plot to results 
png(file = "results/07_clusterProfiler/go_clusters.png",
    width = 38, 
    height = 16,
    units = "in",
    res = 300)

wrap_plots(p.res)

dev.off()

