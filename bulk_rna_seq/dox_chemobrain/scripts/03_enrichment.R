# ============================================================
# SCRIPT 03: ENRICHMENT ANALYSIS
# Load organized data from Script 02 if necessary ---------
# ============================================================

### 1. Prepare count matrix with SYMBOLS --------------------
counts_raw <- counts
counts_raw$"SYMBOL" <- rownames(counts)
counts_raw <- merge(tt, counts_raw, by = "SYMBOL")

counts_sig <- counts_raw[counts_raw$SYMBOL %in% tt_filtered$SYMBOL, ]

# Prepare ranking list for GSEA
gene_list <- counts_raw$log2FoldChange
names(gene_list) <- counts_raw$SYMBOL

gene_list <- gene_list[!is.na(gene_list) & !is.na(names(gene_list))]
gene_list <- gene_list[!duplicated(names(gene_list))]
gene_list <- sort(gene_list, decreasing = TRUE)

### 2. GSEA Analysis ----------------------------------------
# GSEA with GO
gse_GO_all <- gseGO(geneList= gene_list, 
                    ont ="ALL", # BP, CC, MF
                    keyType = "SYMBOL", 
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    pvalueCutoff = 0.05, 
                    verbose = TRUE, 
                    OrgDb = "org.Mm.eg.db",
                    pAdjustMethod = "fdr")
gse_GO_all <- pairwise_termsim(gse_GO_all) 

# Dotplot
dot <- dotplot(gse_GO_all,
                        showCategory = 10,
                        split = "ONTOLOGY",
                        font.size = 12) +
  facet_grid(ONTOLOGY ~ .sign, scales = "free_y", space = "free_y")

dot
ggsave("plots/dotplot_GO.png", plot = dot, width = 7, height = 12, dpi = 300)

# Ridgeplot
make_ridge <- function(obj, ont, n = 10) {
  sub <- obj
  sub@result <- obj@result[obj@result$ONTOLOGY == ont, ]
  ridgeplot(sub,
            showCategory = n,
            fill = "p.adjust",
            core_enrichment = TRUE,
            label_format = 50,
            decreasing = TRUE) +
    ggtitle(ont)
}

onts <- c("BP", "CC", "MF")
onts <- onts[onts %in% gse_GO_all@result$ONTOLOGY]
ridge_list <- lapply(onts, make_ridge, obj = gse_GO_all)
names(ridge_list) <- onts

ridge_all <- wrap_plots(ridge_list, ncol = 1)
ridge_all
ggsave("plots/ridge.png", plot = ridge_all, width = 7, height = 12, dpi = 300)

# Cnetplot
make_cnet <- function(obj, ont, n = 20) {
  sub <- obj
  sub@result <- obj@result[obj@result$ONTOLOGY == ont, ]
  cnetplot(sub,
           layout = igraph::layout_nicely,
           showCategory = n,
           color_category = "#E5C494",
           size_category = 1,
           color_item = "#B3B3B3",
           size_item = 1,
           color_edge = "grey",
           size_edge = 0.5,
           node_label = "category",
           foldChange = gene_list,
           hilight = "none",
           hilight_alpha = 0.3) +
    ggtitle(ont)
}
cnet_list <- lapply(onts, make_cnet, obj = gse_GO_all)
names(cnet_list) <- onts
cnet_all <- wrap_plots(cnet_list, ncol = 1)

cnet_all
ggsave("plots/cnet.png", plot = cnet_all, width = 10, height = 10, dpi = 300)

# Emaplot
make_ema <- function(obj, ont, n = 20) {
  sub <- obj
  sub@result <- obj@result[obj@result$ONTOLOGY == ont, ]
  emapplot(sub,
           layout = igraph::layout_with_kk,
           showCategory = n,
           color = "p.adjust",
           size_category = 1,
           min_edge = 0.2,
           color_edge = "grey",
           size_edge = 0.5,
           node_label = "category",
           group = TRUE,
           label_format = 30,
           clusterFunction = stats::kmeans,
           nWords = 4,
           nCluster = NULL) +
    ggtitle(ont)
}
ema_list <- lapply(onts, make_ema, obj = gse_GO_all)
names(ema_list) <- onts
ema_all <- wrap_plots(ema_list, ncol = 1)

ema_all
ggsave("plots/emaplot.png", plot = ema_all, width = 9, height = 40, dpi = 300)

# Single pathway GSEA plot
go_sets <- data.frame(gse_GO_all@result[["Description"]])
id_go = 4 # Choose pathway of interest

gsea <- gseaplot2(gse_GO_all, 
          geneSetID = id_go,
          pvalue_table = TRUE,
          title = gse_GO_all$Description[id_go],
          color = c("#E495A5", "#86B875", "#7DB0DD"), 
          ES_geom = "line")
gsea
ggsave("plots/gsea.png", plot = gsea, width = 9, height = 7, dpi = 300)

### 3. Convert gene_list to Entrez ID for further analyzes --------
map <- bitr(names(gene_list),
            fromType = "SYMBOL",
            toType   = "ENTREZID",
            OrgDb    = org.Mm.eg.db)

gene_list_entrez <- gene_list[map$SYMBOL]
names(gene_list_entrez) <- map$ENTREZID

gene_list_entrez <- gene_list_entrez[!is.na(names(gene_list_entrez))]
gene_list_entrez <- gene_list_entrez[!duplicated(names(gene_list_entrez))]
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

gse_kegg <- gseKEGG(geneList = gene_list_entrez,
                    keyType = "kegg",
                    nPerm = 10000,
                    organism = "mmu",
                    minGSSize = 3,
                    maxGSSize = 800,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "none")
gse_kegg <- pairwise_termsim(gse_kegg)

# Dotplot KEGG
dot_k <- dotplot(gse_kegg,
                 showCategory = 10,
                 split = ".sign",
                 font.size = 12) +
  facet_grid(. ~ .sign)

dot_k
ggsave("plots/dotplot_kegg.png", plot = dot, width = 7, height = 12, dpi = 300)

# Emaplot KEGG
emakegg <- emapplot(gse_kegg)

emakegg
ggsave("plots/ema_kegg.png", plot = emakegg, width = 8, height = 9, 
       dpi = 300, bg = "white")

# Cnetplot KEGG
cnetkegg <- cnetplot(gse_kegg,
                     layout = igraph::layout_nicely,
                     showCategory = 10,
                     color_category = "#E5C494",
                     size_category = 1,
                     size_item = 1,
                     color_edge = "grey",
                     size_edge = 0.5,
                     node_label = "category",
                     foldChange = gene_list_entrez,
                     hilight = "none",
                     hilight_alpha = 0.3)
cnetkegg
ggsave("plots/cnet_kegg.png", plot = cnetkegg, width = 7, height = 5, 
       dpi = 300, bg = "white")

# Ridge KEGG
ridgekegg <- ridgeplot(gse_kegg,
          showCategory = 10,
          fill = "p.adjust",
          core_enrichment = TRUE,
          label_format = 30
) + labs(x = "enrichment distribution")

ridgekegg
ggsave("plots/ridge_kegg.png", plot = ridgekegg, width = 7, height = 5, dpi = 300)

# Pathview
kegg_sets <- data.frame(Description = gse_kegg@result[["Description"]],
                        ID = gse_kegg@result[["ID"]])

run_pathview <- function(pathway_id, out_dir = file.path(project_dir, "plots", "pathview"), ...) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  old_wd <- setwd(out_dir)
  on.exit(setwd(old_wd))  
  pathview(gene.data = gene_list_entrez, pathway.id = pathway_id,
           species = "mmu", gene.idtype = "ENTREZ", kegg.dir = out_dir, ...)
}
run_pathview("mmu04668", limit = list(gene = 1, cpd = 1),
             out.suffix = "log2FC")

id_kegg = 8
gsea_kegg <- gseaplot2(gse_kegg, 
          geneSetID = id_kegg,
          pvalue_table = TRUE,
          title = gse_kegg$Description[id_kegg],
          color = c("#E495A5", "#86B875", "#7DB0DD"), 
          ES_geom = "line")

gsea_kegg
ggsave("plots/gsea_kegg.png", plot = gsea_kegg, width = 9, height = 7, dpi = 300)
