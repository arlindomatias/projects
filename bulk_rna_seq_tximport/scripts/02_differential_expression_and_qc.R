# ============================================================
# SCRIPT 02: RNA-SEQ QUALITY CONTROL AND DIFFERENTIAL ANALYSIS
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))
suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(EnhancedVolcano)
  library(qqman)
})
out <- snakemake@output

dds <- readRDS(snakemake@input$dds)
padj_cut <- cfg$de$padj
lfc_cut  <- cfg$de$lfc

### 1. Filtro e ajuste -------------------------------------------------------
# tximport traz contagens estimadas (fracionárias); o filtro usa as contagens arredondadas do dds.
min_samples <- cfg$filter$min_samples %||% min(table(dds$group))
keep <- rowSums(counts(dds) >= cfg$filter$min_count) >= min_samples
message("Filter: >= ", cfg$filter$min_count, " counts in >= ", min_samples,
        " samples -> ", sum(keep), " of ", nrow(dds), " genes")
dds <- dds[keep, ]
dds <- DESeq(dds, fitType = cfg$de$fit_type)
saveRDS(dds, out$dds)

coef_name <- paste0("group_", trt_level, "_vs_", ref_level)
message("Coefficients: ", paste(resultsNames(dds), collapse = ", "))
stopifnot(coef_name %in% resultsNames(dds))

res     <- results(dds, contrast = c("group", trt_level, ref_level), alpha = padj_cut)
res_shr <- lfcShrink(dds, coef = coef_name, res = res, type = "apeglm")
summary(res)

ann <- as.data.frame(mcols(dds)[, c("gene_id", "SYMBOL")])
tt <- data.frame(ann, as.data.frame(res), log2FoldChange_shrunk = res_shr$log2FoldChange,
                 row.names = NULL, check.names = FALSE)
tt$label <- gene_label(tt$SYMBOL, tt$gene_id)
tt <- tt[order(tt$padj, tt$pvalue), ]
tt_filtered <- subset(tt, padj < padj_cut & abs(log2FoldChange) > lfc_cut)

write_csv(tt, out$top_table)

de_summary <- tibble(
  contrast   = paste(trt_level, "vs", ref_level),
  design     = deparse(design(dds)),
  samples    = ncol(dds),
  genes_tested = nrow(tt),
  padj_cutoff  = padj_cut,
  lfc_cutoff   = lfc_cut,
  deg_up   = sum(tt_filtered$log2FoldChange > 0),
  deg_down = sum(tt_filtered$log2FoldChange < 0),
  padj_only = sum(tt$padj < padj_cut, na.rm = TRUE)
)
write_tsv(de_summary, out$summary)
print(as.data.frame(de_summary))

norm_counts <- counts(dds, normalized = TRUE)
write_csv(data.frame(gene_id = rownames(norm_counts), SYMBOL = mcols(dds)$SYMBOL,
                     round(norm_counts, 3), check.names = FALSE),
          out$norm_counts)

### 2. QC plots ---------------------------------------------------------------
message("Creating QC and DEG plots")
sample_labels <- setNames(dds$label, colnames(dds))
vsd <- vst(dds, blind = TRUE)

# PCA
pca_data   <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
pca_data$label <- sample_labels[pca_data$name]

pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group, label = label)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
  scale_color_manual(values = group_colors) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"),
       color = "Group") +
  coord_fixed() +
  theme_bw(base_size = 12)
ggsave(out$pca, plot = pca, width = 6, height = 5, dpi = 300)

# Sample-to-sample distances (VST)
annot_col  <- data.frame(Group = dds$group, row.names = colnames(dds))
annot_cols <- list(Group = group_colors)
dist_mat <- as.matrix(dist(t(assay(vsd))))
dimnames(dist_mat) <- list(sample_labels[rownames(dist_mat)], sample_labels[colnames(dist_mat)])
rownames(annot_col) <- sample_labels[rownames(annot_col)]
pheatmap(dist_mat,
         clustering_distance_rows = as.dist(dist_mat),
         clustering_distance_cols = as.dist(dist_mat),
         annotation_col = annot_col, annotation_colors = annot_cols,
         color = colorRampPalette(c("#08306B", "#F7FBFF"))(100),
         main = "Euclidean distance between samples (VST)",
         filename = out$sample_dist, width = 7, height = 6)

# Marker genes (e.g. sex and cell type markers) per sample
markers <- unlist(cfg$plots$marker_genes)
mk_idx  <- which(mcols(dds)$SYMBOL %in% markers)
mk_df <- as.data.frame(norm_counts[mk_idx, , drop = FALSE])
mk_df$SYMBOL <- mcols(dds)$SYMBOL[mk_idx]
mk_df <- mk_df |>
  pivot_longer(-SYMBOL, names_to = "sample", values_to = "count") |>
  mutate(label = factor(sample_labels[sample], levels = sample_labels),
         group = dds$group[match(sample, colnames(dds))])
not_found <- setdiff(markers, mk_df$SYMBOL)
if (length(not_found)) message("Marker genes absent after filtering: ",
                               paste(not_found, collapse = ", "))
p_mk <- ggplot(mk_df, aes(x = label, y = log2(count + 1), fill = group)) +
  geom_col() +
  facet_wrap(~ SYMBOL, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  labs(x = NULL, y = "log2(normalized counts + 1)", fill = "Group",
       title = "Marker genes",
       caption = if (length(not_found)) paste("Filtered out / absent:",
                                              paste(not_found, collapse = ", "))) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(out$markers, plot = p_mk, width = 9, height = 6, dpi = 300)

# Volcano plot
keyvals <- ifelse(tt$padj < padj_cut & tt$log2FoldChange >  lfc_cut, "#CC0000",
                  ifelse(tt$padj < padj_cut & tt$log2FoldChange < -lfc_cut, "#3399FF", "gray"))
keyvals[is.na(keyvals)] <- "gray"
names(keyvals) <- ifelse(keyvals == "#CC0000", "Up",
                         ifelse(keyvals == "#3399FF", "Down", "NS"))

vp <- EnhancedVolcano(tt,
                      lab = tt$label,
                      x = "log2FoldChange",
                      y = "padj",
                      xlab = bquote(~Log[2]~ "fold change"),
                      ylab = bquote(~-Log[10]~ "FDR"),
                      title = cfg$title,
                      subtitle = paste(trt_level, "vs", ref_level),
                      caption = paste("Total genes:", nrow(tt)),
                      pCutoff = padj_cut,
                      FCcutoff = lfc_cut,
                      pointSize = 1.0,
                      labSize = 3,
                      colCustom = keyvals,
                      colAlpha = 0.5,
                      legendLabSize = 7,
                      legendIconSize = 7,
                      drawConnectors = TRUE,
                      widthConnectors = 0.3,
                      max.overlaps = 15)
ggsave(out$volcano, plot = vp, width = 8, height = 7, units = "in", dpi = 300)

# Heatmap: top genes by FDR (DEGs first; capped so the plot stays legible)
n_max <- cfg$plots$heatmap_max_genes
hm_genes <- head(tt_filtered$gene_id, n_max)
hm_title <- paste0("Top ", length(hm_genes), " of ", nrow(tt_filtered), " DEGs")
if (length(hm_genes) < 2) {
  hm_genes <- head(tt$gene_id[!is.na(tt$padj)], n_max)
  hm_title <- paste0("No DEGs: top ", length(hm_genes), " genes by FDR")
}
hm_mat <- assay(vsd)[hm_genes, , drop = FALSE]
rownames(hm_mat) <- tt$label[match(hm_genes, tt$gene_id)]
colnames(hm_mat) <- sample_labels[colnames(hm_mat)]
pheatmap(hm_mat,
         annotation_col    = annot_col,
         annotation_colors = annot_cols,
         show_rownames = TRUE,
         cluster_rows  = TRUE,
         cluster_cols  = TRUE,
         scale         = "row",
         fontsize      = 8,
         main          = hm_title,
         filename = out$heatmap, width = 6,
         height = min(4 + 0.13 * length(hm_genes), 14))

# FDR histogram
histo <- ggplot(tt, aes(x = padj)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black",
                 boundary = 0, na.rm = TRUE) +
  labs(title = "FDR distribution", x = "padj", y = "Number of genes") +
  theme_classic()
ggsave(out$fdr, plot = histo, width = 8, height = 5, units = "in", dpi = 300)

# Box plot of normalized counts, colored by group
meta <- data.frame(sample = colnames(dds), group = dds$group,
                   label = factor(sample_labels, levels = sample_labels))
df_long <- as.data.frame(norm_counts) |>
  rownames_to_column("gene_id") |>
  pivot_longer(-gene_id, names_to = "sample", values_to = "expression") |>
  left_join(meta, by = "sample")

bp <- ggplot(df_long, aes(x = label, y = log2(expression + 1), fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = group_colors) +
  labs(title = "Normalized expression",
       x = NULL, y = "log2(normalized counts + 1)", fill = "Group") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(out$boxplot, plot = bp, width = 8, height = 5, units = "in", dpi = 300)

# QQ-Plot
png(out$qqplot, width = 6, height = 6, units = "in", res = 300)
qqman::qq(tt$pvalue[!is.na(tt$pvalue)])
dev.off()

# MA plot (shrunken log2FC)
ma_df <- as.data.frame(res_shr)
ma_df <- ma_df[ma_df$baseMean > 0, ]
ma_df$sig <- ifelse(is.na(ma_df$padj) | ma_df$padj >= padj_cut, "ns",
                    ifelse(ma_df$log2FoldChange > 0, "Up", "Down"))
ma_df$sig <- factor(ma_df$sig, levels = c("Up", "Down", "ns"))

p_ma <- ggplot(ma_df, aes(x = baseMean, y = log2FoldChange, color = sig)) +
  geom_point(size = 0.8, alpha = 0.6, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_log10() +
  scale_color_manual(values = c("Up" = "#CC0000", "Down" = "#0066CC", "ns" = "grey70")) +
  labs(x = "Mean of normalized counts (log10)",
       y = "log2 fold change (apeglm)",
       color = NULL,
       title = "MA plot") +
  coord_cartesian(ylim = c(-3, 3)) +
  theme_classic(base_size = 12)
ggsave(out$ma, plot = p_ma, width = 6, height = 5, dpi = 300)

# Dispersion plot
png(out$dispersion, width = 6, height = 5, units = "in", res = 300)
plotDispEsts(dds, main = "Dispersion estimates")
dev.off()

# Normalized counts density
df_dens <- df_long |> mutate(log2_expr = log2(expression + 1))
p_dens <- ggplot(df_dens, aes(x = log2_expr, color = label, linetype = group)) +
  geom_density(linewidth = 0.7) +
  labs(x = "log2(normalized counts + 1)", y = "Density", color = "Sample",
       title = "Density of normalized counts") +
  theme_classic(base_size = 12)
ggsave(out$density, plot = p_dens, width = 7, height = 5, dpi = 300)

# Cook's distance
cooks <- assays(dds)[["cooks"]]
colnames(cooks) <- sample_labels[colnames(cooks)]
png(out$cooks, width = 7, height = 5, units = "in", res = 300)
par(mar = c(8, 4.5, 3, 1))
boxplot(log10(cooks), range = 0, las = 2,
        ylab = "log10(Cook's distance)", main = "Cook's distance")
dev.off()
