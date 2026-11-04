# ============================================================
# SCRIPT 02: RNA-SEQ QUALITY CONTROL AND DIFFERENTIAL ANALYSIS
# Load organized data from Script 01 if necessary ---------

# 1. Create DDS object

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ group)


dds <- DESeq(dds, fitType = "local")
tt <- data.frame(results(dds))


plotPCA(vst(dds))

tt$SYMBOL <- rownames(tt)

write.csv(tt, "top_table.csv")

tt_entrez <- merge(keytable, tt, by = "SYMBOL")

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

## Heatmap
deg_counts <- counts
deg_counts <- deg_counts[rownames(counts) %in% tt_filtered$SYMBOL, ]
heat_legend <- colData
rownames(heat_legend) <- heat_legend$GEO_ID
heat_legend$GEO_ID <- NULL
heat_legend$condition <- NULL

heatmap <- pheatmap(deg_counts, 
                    annotation_col = heat_legend, 
                    show_rownames = T, 
                    cluster_cols = T, 
                    scale = "row",
                    filename = "heatmap.png")

## FDR histogram 
par(mar = c(4.1, 4.1, 3.1, 2.1))   # sem cex.axis = 0.3 e cex.lab = 2

histo <- ggplot(tt, aes(x = padj)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black",
                 boundary = 0, na.rm = TRUE) +
  labs(title = "FDR distribution", x = "padj", y = "Number of genes") +
  theme_classic()

ggsave("fdr_histogram.png", plot = histo,
       width = 8, height = 5, units = "in", dpi = 300)

## Box plot
df_long <- as.data.frame(counts) %>%
  rownames_to_column("SYMBOL") %>%
  pivot_longer(-SYMBOL, names_to = "sample", values_to = "expression")

boxplot <- ggplot(df_long, aes(x = sample, y = log2(expression + 1), fill = sample)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = c("#3399FF", "#3399FF", "#3399FF",
                               "#CC0000", "#CC0000", "#CC0000")) +
  labs(title = "Normalized expression",
       x = NULL, y = "log2(normalized counts + 1)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none")

ggsave("boxplot.png", plot = boxplot,
       width = 8, height = 5, units = "in", dpi = 300)

