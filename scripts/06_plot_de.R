# Visualize DE 
libs <- c("tidyverse", "wesanderson", "ggrepel") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

# Load the results 
res <- readRDS("data/processed/deseq/deseq.RDS")  
names <- resultsNames(res)[-c(1:3)]

# make a list of each result
res <- lapply(names, function(x){
  results(res, name=x) |> 
    as.data.frame()
})

# Unadjusted results first 
plotGenes <- function(deseq.res, name){
  
#  deseq.res <- res[[1]]
# Find top DE genes
top.genes <- deseq.res |> 
  arrange(padj) |>
  slice_head(n = 10)
top.genes$gene <- row.names(top.genes)
#name <- names[[1]]

# Add conditional color formatting for significance
deseq.res <- deseq.res |> 
  mutate(significant = case_when(
    padj >  0.05 ~ FALSE,
    padj <= 0.05 ~ TRUE,
    .default = FALSE
  )) |> 
  subset(!is.na(padj))

plot <- ggplot(deseq.res, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 1, size = 6) + 
  scale_color_manual(values = c("#999999", "#ed9209")) +
  geom_text_repel(data = top.genes, aes(label = gene), size = 12, box.padding = unit(0.35, "lines"), 
                  force = 20, segment.linetype = 2, segment.size = 0.6, color = "black", force_pull = 0.01, min.segment.length = 0) + # label most sig genes 
  theme_minimal() +
  labs(x = "log2(Fold Change)", y = "-log10(adjusted p-value)", color = "p < 0.05", title = name) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "solid") + # add a line for p=0.05
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black", size = 25),
        axis.ticks = element_blank(),
        title = element_text(size = 25),
        legend.text = element_text(size = 25),
        legend.key.size = unit(0.7, 'in'),
        legend.title = element_text(size = 28, vjust = 0.7),
        axis.title = element_text(color = "black", size = 28),
        panel.background = element_rect(color="black"),
        plot.background = element_rect(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
}

de.plots <- lapply(1:length(res), function(n){plotGenes(res[[n]], names[[n]])})

  # Save 
  png(file = "results/06_plot_de/volcano_all.png",
      width = 40, 
      height = 27,
      units = "in",
      res = 300)
  
  wrap_plots(de.plots)
  
  dev.off()

bulk <- read.csv("data/processed/bulk/bulk_gene.csv", row.names = 1) 
phenotypes <- read.csv("data/processed/bulk/phenotypes")

phenotypes$sample.id <- paste0("BCJ_",phenotypes$sample.id)
# Prepare for DESeq2
bulk <- mutate_all(bulk, function(x) round(as.numeric(as.character(x)), digits = 0)) # round to integers

# Plot the PCA of samples
# Run PCA
pca <- prcomp(bulk)
  
PoV <- pca$sdev^2/sum(pca$sdev^2) * 100  
PoV <- round(PoV, digits = 1)
  # Merge with sample info, then plot
pca <- pca$rotation |> 
    as.data.frame() |> 
    dplyr::select(PC1, PC2) 
pca$sample.id <- row.names(pca)
  
  #my_palette <- c("#A6CEE3", "#1F78B4", "#FDBF6F", "#FF7F00")
  #legend.names <- c("Sham_1","Sham_2", "TAC_1", "TAC_2")
  
  pca <- pca |> left_join(phenotypes) 
  pca.plot <- pca |> 
    ggplot(aes(x = PC1, y = PC2, color = treatment, shape = NRVM)) +
    geom_point(size = 8, color = "black") +
    geom_point(size = 7) +
    #scale_color_manual(values = my_palette) +
    #scale_fill_manual(values = my_palette) +
    scale_x_continuous(expand = expansion(mult = 0.3), name = paste0("PC1", " (", PoV[1], " % of total variance)")) +
    scale_y_continuous(expand = expansion(mult = 0.3), name = paste0("PC2", " (", PoV[2], " % of total variance)")) +
    theme(axis.text.x = element_text(vjust = 0.5),
          axis.ticks = element_blank(),
          legend.position = "right",
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6),
          legend.text = element_text(size = 20),
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent', color=NA),
          panel.grid.major = element_line(color = "darkgrey"),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill='transparent'),
          legend.box.background = element_rect(fill='transparent'),
          plot.margin = unit(c(1,1,1,1), units = "cm"),
          text = element_text(size = 25))
  
  pca.plot
  # Save plot to results 
  png(file = "results/06_plot_de/pca.png",
      width = 12, 
      height = 9,
      units = "in",
      res = 300)
  
  pca.plot
  
  dev.off()
  
  
