# ============================================================
# SCRIPT 02: RNA-SEQ QUALITY CONTROL AND DIFFERENTIAL ANALYSIS
# Load organized data from Script 01 if necessary ---------

# 1. Create DDS object

keytable <- bitr(rownames(counts),
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Mm.eg.db)
keytable <- na.omit(keytable)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ group)


dds <- DESeq(dds, fitType = "local")
result <- results(dds)

# Top table
tt <- data.frame(result)
tt$SYMBOL <- rownames(tt)

write.csv(
  tt,
  file = file.path(
    results_dir,
    "01_top_table.csv"
  ),
  quote = FALSE
)

tt2 <- merge(keytable, tt, by = "SYMBOL")

tt_filtered <- subset(tt, padj < 0.05 & abs(log2FoldChange) > 0.5)
tt_filtered_up <- subset(tt, padj < 0.05 & log2FoldChange > 0.5)
tt_filtered_down <- subset(tt, padj < 0.05 & log2FoldChange < - 0.5)

## Volcano plot
keyvals <- ifelse(tt$padj < 0.05 & tt$log2FoldChange >  0.5, "#CC0000",
                  ifelse(tt$padj < 0.05 & tt$log2FoldChange < -0.5, "#3399FF", "gray"))
keyvals[is.na(keyvals)] <- "gray"
names(keyvals) <- ifelse(keyvals == "#CC0000", "Up",
                         ifelse(keyvals == "#3399FF", "Down", "NS"))

v <- EnhancedVolcano(tt,
                     lab = tt$SYMBOL,
                     x = 'log2FoldChange',
                     y = 'padj',
                     xlab = bquote(~Log[2]~ 'fold change'),
                     ylab = bquote(~-Log[10]~ 'FDR'),
                     title = "Chemotherapy-induced gene expression change",
                     caption = paste("Total genes:", nrow(tt)),
                     pCutoff = 0.05,
                     FCcutoff = 0.5,
                     pointSize = 1.0,
                     labSize = 3,
                     xlim = c(-2, 4),
                     ylim = c(0, 8),
                     colCustom = keyvals,    
                     colAlpha = 0.5,
                     legendLabSize = 7,
                     legendIconSize = 7,
                     drawConnectors = TRUE,
                     widthConnectors = 0.3,
                     max.overlaps = 15)
v

ggsave("volcano.pdf", plot = v, width = 8, height = 7)
