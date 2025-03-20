library(clusterProfiler)
library(org.Mm.eg.db)

# Load DEGs
res <- read.csv(commandArgs(trailingOnly = TRUE)[1], row.names = 1)

# Enrichment analysis
ego <- enrichGO(gene = rownames(res)[res$padj < 0.05],
                OrgDb = org.Mm.eg.db, 
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05)

# Save pathway results
write.csv(ego, paste0(commandArgs(trailingOnly = TRUE)[2], "/pathway.csv"))
