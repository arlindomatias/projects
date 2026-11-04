# ============================================================
# SCRIPT 03: ENRICHMENT ANALYSIS
# Load organized data from Script 02 if necessary ---------


## Prepare count matrix for entrez IDs

entrez_count <- counts
entrez_count$"SYMBOL" <- rownames(entrez_count)
entrez_count <- merge(tt_entrez, entrez_count, by = "SYMBOL")

entrez_sig <- entrez_count[entrez_count$SYMBOL %in% tt_filtered$SYMBOL, ]

entrez_single <- entrez_count %>%
  group_by(ENTREZID) %>%
  summarise(
    SYMBOL = dplyr::first(SYMBOL),
    across(where(is.numeric), mean)
  )
entrez_single <- entrez_single[!is.na(entrez_single$log2FoldChange), ]

## Prepare ranking list for GSEA
gene_list <- entrez_single$log2FoldChange
names(gene_list) <- entrez_single$ENTREZID

gene_list <- gene_list[!is.na(gene_list) & !is.na(names(gene_list))]
gene_list <- gene_list[!duplicated(names(gene_list))]
gene_list <- sort(gene_list, decreasing = TRUE)

## GSEA with GO
gse_GO_all <- gseGO(geneList= gene_list, 
                    ont ="ALL", # BP, CC, MF
                    keyType = "ENTREZID", 
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    pvalueCutoff = 0.05, 
                    verbose = TRUE, 
                    OrgDb = "org.Mm.eg.db",
                    pAdjustMethod = "none")
gse_GO_all <- pairwise_termsim(gse_GO_all) 

dot <- dotplot(gse_GO_all,
        font.size = 12,
        showCategory=5, 
        split=".sign") + facet_grid(.~.sign)
dot
ggsave("dotplot_GO.png", plot = dot, width = 7, height = 5, dpi = 300)

ridge <- ridgeplot(gse_GO_all,
          showCategory = 10,
          fill = "p.adjust",
          core_enrichment = TRUE,
          label_format = 50,
          decreasing = TRUE)
ridge
ggsave("ridge.png", plot = ridge, width = 7, height = 5, dpi = 300)

cnet <- cnetplot(gse_GO_all,
         layout = igraph::layout_nicely,
         showCategory = 20,
         color_category = "#E5C494",
         size_category = 1,
         color_item = "#B3B3B3",
         size_item = 1,
         color_edge = "grey",
         size_edge = 0.5,
         node_label = "category",
         foldChange = gene_list,
         hilight = "none",
         hilight_alpha = 0.3)
cnet
ggsave("cnet.png", plot = cnet, width = 7, height = 5, dpi = 300)


ema <- emapplot(gse_GO_all,
                layout = igraph::layout_with_kk,
                showCategory = 20,
                color = "p.adjust",
                size_category = 1,
                min_edge = 0.2,
                color_edge = "grey",
                size_edge = 0.5,
                node_label = "category",
                pie = "equal",
                group = TRUE,
                label_format = 30,
                clusterFunction = stats::kmeans,
                nWords = 4,
                nCluster = NULL)
ema
ggsave("emapplot_GSEA_GO.png", plot = ema, width = 9, height = 7, dpi = 300)



go_sets <- data.frame(gse_GO_all@result[["Description"]])
id_go = 5

gsea <- gseaplot2(gse_GO_all, 
          geneSetID = id_go,
          pvalue_table = TRUE,
          title = gse_GO_all$Description[id_go],
          color = c("#E495A5", "#86B875", "#7DB0DD"), 
          ES_geom = "line")
gsea
ggsave("gsea.png", plot = gsea, width = 9, height = 7, dpi = 300)


gse_kegg <- gseKEGG(geneList = gene_list,
                    keyType = "kegg",
                    nPerm = 10000,
                    organism = "mmu",
                    minGSSize = 3,
                    maxGSSize = 800,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "none")
gse_kegg <- pairwise_termsim(gse_kegg)


kegg_sets <- data.frame(gse_kegg@result[["Description"]])


dotkegg <- dotplot(gse_kegg, 
        showCategory = 10, 
        title = "Enriched Pathways" , 
        split=".sign") + facet_grid(.~.sign)
dotkegg
ggsave("dotplot_kegg.png", plot = dotkegg, width = 7, height = 5, dpi = 300)

emakegg <- emapplot(gse_kegg)
emakegg
ggsave("ema_kegg.png", plot = emakegg, width = 7, height = 5, dpi = 300)


cnetkegg <- cnetplot(gse_kegg,
         layout = igraph::layout_nicely,
         showCategory = 10,
         color_category = "#E5C494",
         size_category = 1,
         color_item = "#B3B3B3",
         size_item = 1,
         color_edge = "grey",
         size_edge = 0.5,
         node_label = "category",
         foldChange = gene_list,
         hilight = "none",
         hilight_alpha = 0.3)
cnetkegg
ggsave("cnet_kegg.png", plot = cnetkegg, width = 7, height = 5, dpi = 300)

ridgekegg <- ridgeplot(gse_kegg,
          #showCategory = 10,
          #fill = "p.adjust",
          #core_enrichment = TRUE,
          #label_format = 30
) + labs(x = "enrichment distribution")
ridgekegg
ggsave("ridge_kegg.png", plot = ridgekegg, width = 7, height = 5, dpi = 300)

path <- pathview(gene.data = gene_list,
                 pathway.id="mmu00071", 
                 species = "mmu")






