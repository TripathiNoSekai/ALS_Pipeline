library(DESeq2)

# Load data
data <- read.csv(commandArgs(trailingOnly = TRUE)[1], row.names = 1)
coldata <- data.frame(condition = rep(c('control', 'ALS'), each = 3), row.names = colnames(data))

# DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = data, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)

# Save results
write.csv(res, paste0(commandArgs(trailingOnly = TRUE)[2], "/DEGs.csv"))
