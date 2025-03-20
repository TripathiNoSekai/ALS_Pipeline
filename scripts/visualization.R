library(ggplot2)
library(pheatmap)

# Load data
res <- read.csv(commandArgs(trailingOnly = TRUE)[1], row.names = 1)

# Volcano plot
png(paste0(commandArgs(trailingOnly = TRUE)[2], "/volcano.png"))
plot(res$log2FoldChange, -log10(res$pvalue), col = ifelse(res$padj < 0.05, 'red', 'black'), pch = 20)
dev.off()

# Heatmap
png(paste0(commandArgs(trailingOnly = TRUE)[2], "/heatmap.png"))
pheatmap(as.matrix(res), cluster_rows = TRUE, cluster_cols = TRUE)
dev.off()
