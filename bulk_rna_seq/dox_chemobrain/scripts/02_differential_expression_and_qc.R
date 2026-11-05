# ============================================================
# SCRIPT 02: RNA-SEQ QUALITY CONTROL AND DIFFERENTIAL ANALYSIS
# Load organized data from Script 01 if necessary ---------
# ============================================================

### 1. Create DDS object and Top Table --------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ group)
dds <- DESeq(dds, fitType = "local")
resultsNames(dds)

keep <- rowSums(counts(dds) >= 10) >= 3   # >= 10 counts in at least 3 samples
dds <- dds[keep, ]

tt <- data.frame(results(dds)) # Top Table
tt$SYMBOL <- rownames(tt)

write.csv(tt, "data/top_table.csv") # Export TT

tt_entrez <- merge(keytable, tt, by = "SYMBOL")

tt_filtered <- subset(tt, padj < 0.05 & abs(log2FoldChange) > 0.5)
tt_filtered_up <- subset(tt, padj < 0.05 & log2FoldChange > 0.5)
tt_filtered_down <- subset(tt, padj < 0.05 & log2FoldChange < - 0.5)

### 1. Create DEG and QC plots ---------------------------------------

# PCA
vsd <- vst(dds, blind = TRUE)
pca_data   <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("#3399FF", "#CC0000", "#66A61E")) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"),
       color = "Group") +
  coord_fixed() +
  theme_bw(base_size = 12)

pca
ggsave("plots/pca_vst.png", plot = pca,
       width = 6, height = 5, dpi = 300)

# Volcano plot
keyvals <- ifelse(tt$padj < 0.05 & tt$log2FoldChange >  0.5, "#CC0000",
                  ifelse(tt$padj < 0.05 & tt$log2FoldChange < -0.5, "#3399FF", "gray"))
keyvals[is.na(keyvals)] <- "gray"
names(keyvals) <- ifelse(keyvals == "#CC0000", "Up",
                         ifelse(keyvals == "#3399FF", "Down", "NS"))

vp <- EnhancedVolcano(tt,
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
                     xlim = c(-5, 5),
                     #ylim = c(0, 10),
                     colCustom = keyvals,    
                     colAlpha = 0.5,
                     legendLabSize = 7,
                     legendIconSize = 7,
                     drawConnectors = TRUE,
                     widthConnectors = 0.3,
                     max.overlaps = 15)

vp
ggsave("plots/volcano.pdf", plot = vp, width = 8, height = 7,
       units = "in", dpi = 300)

# Heatmap
deg_counts <- counts
deg_counts <- deg_counts[rownames(counts) %in% tt_filtered$SYMBOL, ]
heat_legend <- colData
rownames(heat_legend) <- heat_legend$GEO_ID
heat_legend$GEO_ID <- NULL
heat_legend$condition <- NULL

hm <- pheatmap(deg_counts, 
                    annotation_col = heat_legend, 
                    show_rownames = T, 
                    cluster_cols = T, 
                    scale = "row",
                    fontsize = 8,
                    cellwidth = 10,
                    cellheight = 12)

hm
ggsave("plots/heatmap.pdf", plot = hm, 
       width = 4, height = 12, units = "in", dpi = 300)

# FDR histogram 
histo <- ggplot(tt, aes(x = padj)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black",
                 boundary = 0, na.rm = TRUE) +
  labs(title = "FDR distribution", x = "padj", y = "Number of genes") +
  theme_classic()

histo
ggsave("plots/fdr_histogram.png", plot = histo,
       width = 8, height = 5, units = "in", dpi = 300)

# Box plot
df_long <- as.data.frame(counts) %>%
  rownames_to_column("SYMBOL") %>%
  pivot_longer(-SYMBOL, names_to = "sample", values_to = "expression")

bp <- ggplot(df_long, aes(x = sample, y = log2(expression + 1), fill = sample)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = c("#3399FF", "#3399FF", "#3399FF",
                               "#CC0000", "#CC0000", "#CC0000")) +
  labs(title = "Normalized expression",
       x = NULL, y = "log2(normalized counts + 1)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none")

bp
ggsave("plots/boxplot.png", plot = bp,
       width = 8, height = 5, units = "in", dpi = 300)

# QQ-Plot
png("plots/QQPlot.png", 
    width = 6, height = 6, units = "in", res = 300)
qqman::qq(tt$pvalue)
dev.off()

# MA plot
resultsNames(dds)
res_shr <- lfcShrink(dds, coef = "group_CHEM_vs_VEH", type = "apeglm")

ma_df <- as.data.frame(res_shr)
ma_df <- ma_df[ma_df$baseMean > 0, ]

ma_df$sig <- ifelse(is.na(ma_df$padj) | ma_df$padj >= 0.05, "ns",
                    ifelse(ma_df$log2FoldChange > 0, "Up", "Down"))
ma_df$sig <- factor(ma_df$sig, levels = c("Up", "Down", "ns"))

p_ma <- ggplot(ma_df, aes(x = baseMean, y = log2FoldChange, color = sig)) +
  geom_point(size = 0.8, alpha = 0.6, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_log10() +
  scale_color_manual(values = c("Up" = "#CC0000", "Down" = "#0066CC", "ns" = "grey70")) +
  labs(x = "Mean of normalized counts (log10)",
       y = "log2 fold change",
       color = NULL,
       title = "MA plot") +
  coord_cartesian(ylim = c(-3, 3)) +
  theme_classic(base_size = 12)

p_ma
ggsave("plots/MAplot.png", plot = p_ma,
       width = 6, height = 5, dpi = 300)

# Dispersion plot
png("plots/dispersion_plot.png",
    width = 6, height = 5, units = "in", res = 300)
plotDispEsts(dds, main = "Dispersion estimates")
dev.off()

# Normalized counts density
norm_counts <- counts(dds, normalized = TRUE)

df_dens <- as.data.frame(log2(norm_counts + 1)) %>%
  pivot_longer(everything(), names_to = "amostra", values_to = "log2_expr")

norm_counts <- counts(dds, normalized = TRUE)

df_dens <- as.data.frame(log2(norm_counts + 1)) %>%
  pivot_longer(everything(), names_to = "amostra", values_to = "log2_expr")
meta <- data.frame(amostra = colnames(dds), group = colData(dds)$group)
df_dens <- left_join(df_dens, meta, by = "amostra")

p_dens <- ggplot(df_dens, aes(x = log2_expr, color = amostra, linetype = group)) +
  geom_density(linewidth = 0.7) +
  labs(x = "log2(normalized counts + 1)", y = "Density",
       title = "Density of normalized counts") +
  theme_classic(base_size = 12)

p_dens

ggsave("plots/density_counts.png", plot = p_dens,
       width = 7, height = 5, dpi = 300)

# Cooks Distance
png("plots/cooks_distance.png",
    width = 7, height = 5, units = "in", res = 300)
par(mar = c(8, 4.5, 3, 1))
boxplot(log10(assays(dds)[["cooks"]]), range = 0, las = 2,
        ylab = "log10(Cook's distance)", main = "Cook's distance")
dev.off()
