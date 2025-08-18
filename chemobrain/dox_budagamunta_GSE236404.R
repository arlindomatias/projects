######################## Packages ########################
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.21")
.libPaths("~/R/library")
setwd("~/R/wd")

library("tidyverse")
library("GEOquery")
library("limma")
library("dplyr")
library("ggplot2")
library("clusterProfiler")
library("org.Mm.eg.db")
library("org.Rn.eg.db")
library("enrichplot")
library("EnhancedVolcano")
library("pheatmap")
library("enrichplot")
library("DOSE")
library("escape")
library("readxl")
library("Cairo")
library("DESeq2")
library("apeglm")
library("pathview")

######################## Adjusts ########################
par(mar=c(4.1, 4.1, 3.1, 2.1), cex.axis = 2, cex.lab = 2) 
palette(c("#3399FF", "#CC0000", "#66A61E"))
dev.new(width= 15, height=15, noRStudioGD = TRUE) 
options(scipen=0) 

######################## Functions ######################
# Save result tables
save_table <- function(table, file_name){
  if (!grepl("\\.csv$", file_name, ignore.case = TRUE)) {
    file_name <- paste0(file_name, ".csv")
  }
  write.csv(tabela, file = file_name, row.names = TRUE)
}

# Fast save of the last plot
fast_save <- ggsave(paste0("graph_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
             plot = last_plot(),
             dpi = 300,
             width = 8, height = 8, 
             units = "in",    
             device = "png",    
             bg = "white")   

## Função aprimorada para salvar gráficos
picture_save <- function(file_name = paste0("graph_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"), 
               plot_obj = last_plot(), largura = 2000, altura = 2000, resolucao = 300, formato = "png") {
  if (formato == "png") {
    CairoPNG(filename = file_name, width = largura, height = altura, res = resolucao)
  } else if (formato == "eps") {
    CairoPS(file = file_name, width = largura / resolucao, height = altura / resolucao, onefile = FALSE)
  } else if (formato == "pdf") {
    CairoPDF(file = file_name, width = largura / resolucao, height = altura / resolucao)
  } else {
    stop("Unsurported format. Use 'png', 'eps' ou 'pdf'.")
  }
  on.exit(dev.off())  
  print(plot_obj)     
}

######################## Import Data ###################
gse <- getGEO("GSE236404", GSEMatrix = TRUE)[[1]]
fvarLabels(gse) <- make.names(fvarLabels(gse))
info <- pData(gse)
fdata <- fData(gse)
info <- info[info$`treatment:ch1` != "Dox + ABT-263",]

counts <- read_csv("GSE236404_DOX_Normalized_Count-6-23-23.csv.gz")
counts <- counts[ , c(1:2, order(colnames(counts)[3:ncol(counts)]) + 2)]
counts <- counts[ , -c(11:18)]

# Create groups
group <- make.names(info$`treatment:ch1`) %>% factor()
levels(group)[levels(group) == "Dox"] <- "CHEM"
levels(group)[levels(group) == "Vehicle"] <- "VEH"

# Criar colData e garantir os nomes das amostras
colData <- data.frame(condition = group)
rownames(colData) <- rownames(info)

# Ajustar colnames de counts
counts_raw <- counts %>%
  group_by(gene_id) %>%
  summarise(
    gene_id = dplyr::first(gene_id),
    across(where(is.numeric), mean)
  ) %>% data.frame()

rownames(counts_raw) <- counts_raw$gene_id
counts_raw$gene_id <- NULL
colnames(counts_raw) <- rownames(info)

# Map ENSEMBL -> ENTREZ
keytable <- bitr(counts$gene_id,
                          fromType = "ENSEMBL",
                          toType = "ENTREZID",
                          OrgDb = org.Rn.eg.db)
keytable <- na.omit(keytable)


# Analysis
limma_table <- normalizeBetweenArrays(counts_raw) %>% data.frame()
design <- model.matrix(~colData$condition + 0, limma_table)
colnames(design) <- levels(group)

vooma_graph <- vooma(limma_table, design, plot = T)
vooma_graph$genes <- fdata


# Filtro: genes com FDR < 0.02 e |logFC| > 1
tt_filtered <- subset(tt2, pvalue < 0.01 & abs(log2FoldChange) > 0.5)
tt_filtered_up <- subset(tt2, pvalue < 0.01 & log2FoldChange > 0.5)
tt_filtered_down <- subset(tt2, pvalue < 0.01 & log2FoldChange < - 0.5)

plotPCA(vst(dds))

write.table(tt, "C:/Users/AsRock/Documents/cavalier.txt")

## Volcano plot
# Adicionar coluna para categorização
tt2$threshold <- "Not significant"
tt2$threshold[tt2$pvalue < 0.01 & tt2$log2FoldChange > 0.2] <- "Upregulated"
tt2$threshold[tt2$pvalue < 0.01 & tt2$log2FoldChange < - 0.2] <- "Downregulated"

# Renderizar
v <- EnhancedVolcano(tt2, 
                      lab = tt2$SYMBOL,   # Labels dos pontos (nomes dos genes)
                      x = 'log2FoldChange',    # Coluna com os dados de log2 fold change
                      y = 'pvalue',             # Coluna com os dados de p-valor ajustado (FDR)
                      xlab = bquote(~Log[2]~ 'fold change'),
                      title = "Chemotherapy-induced gene expression change", # Título do gráfico
                      caption = paste("Total genes:", nrow(tt)),  # Legenda com número de genes
                      pCutoff = 0.01,      # Limite de p-valor (ou FDR) para destacar genes
                      FCcutoff = 1,       # Limite mínimo de log2FC para considerar relevante
                      pointSize = 3.0,      # Tamanho dos pontos (bolinhas) no gráfico
                      labSize = 3,
                      xlim = c(-6, 6),    # Limites do eixo X
                      ylim = c(0, 6),    # Limites do eixo Y
                      col = c("gray", "#3399FF", "#009933", "#CC0000"),  # Cores: NS, Down, Up, ambos
                      colAlpha = 0.5,              # Transparência dos pontos
                      legendLabels = c("NS", "Log2FC", "FDR", "FDR & Log2FC"), # Legenda das cores
                      legendLabSize = 7,          # Tamanho do texto da legenda
                      legendIconSize = 7,        # Tamanho dos ícones da legenda
                      drawConnectors = TRUE,       # Desenhar linhas ligando labels aos pontos
                      widthConnectors = 0.3,      # Espessura das linhas que ligam os labels
                      max.overlaps = 30
)

#heatmap
deg_counts <- counts_raw[1:6]
deg_counts <- deg_counts[rownames(counts_raw) %in% tt_filtered$SYMBOL, ]
heatmap <- pheatmap(deg_counts, 
                    annotation_col = colData, 
                    show_rownames = T, 
                    cluster_cols = F, 
                    scale = "row")


boxplot <- boxplot(deg_counts, 
                   boxwex=0.6, 
                   notch=T, 
                   main= "Normalized expression", 
                   outline=FALSE, 
                   las=2, 
                   col= deg_counts) 

## Histograma de valores p e fdr 
hist(tt$padj, col = "grey", border = "white", xlab = "p value",
     ylab = "Number of genes", main = "FDR distribution")

ego_total <- enrichGO(gene = tt_filtered$SYMBOL,
                #universe = counts$SYMBOL,
                OrgDb         = org.Mm.eg.db,   
                keyType       = "SYMBOL",
                ont           = "CC", #BP, MF, CC
                pAdjustMethod = "none",
                pvalueCutoff  = 0.05)

ego_up <- enrichGO(gene = tt_filtered_up$SYMBOL,
                      #universe = tt$SYMBOL,
                      OrgDb         = org.Mm.eg.db,   
                      keyType       = "SYMBOL",
                      ont           = "CC",
                      pAdjustMethod = "none",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

ego_down <- enrichGO(gene = tt_filtered_down$SYMBOL,
                   #universe = tt$SYMBOL,
                   OrgDb         = org.Mm.eg.db,   
                   keyType       = "SYMBOL",
                   ont           = "CC",
                   pAdjustMethod = "none",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
# Gráficos
dotplot(ego_down,
        x = "geneRatio",
        font.size = 12,
        showCategory = 10) + scale_color_gradientn(colors = c("#0099ff", "#6600ff", "#ff0000"))


barplot(ego_total, 
        showCategory=20,
        font.size = 8,
        title = "Gene Enrichment Analysis",
        label_format = 60,) + scale_color_gradientn(colors = c("#0099ff", "#6600ff", "#ff0000"))


# Preparar matriz de contagem com ENTREZ
entrez_count <- counts_raw[1:6]
entrez_count$"SYMBOL" <- rownames(entrez_count)
entrez_count <- merge(tt2, entrez_count, by = "SYMBOL")

# Filtrar apenas os genes diferencialmente expressos
entrez_sig <- entrez_count[entrez_count$SYMBOL %in% tt_filtered$SYMBOL, ]

# Remover duplicatas por ENTREZID
entrez_single <- entrez_count %>%
  group_by(ENTREZID) %>%
  summarise(
    SYMBOL = dplyr::first(SYMBOL),
    across(where(is.numeric), mean)
  )
entrez_single <- entrez_single[!is.na(entrez_single$log2FoldChange), ]

# Criar lista
gene_list <- sort(entrez_single$log2FoldChange, decreasing = TRUE)
names(gene_list) <- entrez_single$ENTREZID

## GSEA usando GO
gse_GO_all <- gseGO(geneList= gene_list, 
                ont ="ALL", 
                keyType = "ENTREZID", 
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE, 
                OrgDb = "org.Mm.eg.db",
                pAdjustMethod = "none")

gse_GO_bp <- gseGO(geneList= gene_list, 
                    ont ="BP", 
                    keyType = "ENTREZID", 
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    pvalueCutoff = 0.05, 
                    verbose = TRUE, 
                    OrgDb = "org.Mm.eg.db",
                    pAdjustMethod = "none")

gse_GO_mf <- gseGO(geneList= gene_list, 
                   ont ="MF", 
                   keyType = "ENTREZID", 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 0.05, 
                   verbose = TRUE, 
                   OrgDb = "org.Mm.eg.db",
                   pAdjustMethod = "none")

gse_GO_cc <- gseGO(geneList= gene_list, 
                   ont ="CC", 
                   keyType = "ENTREZID", 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 0.05, 
                   verbose = TRUE, 
                   OrgDb = "org.Mm.eg.db",
                   pAdjustMethod = "none")

gse_GO_all <- pairwise_termsim(gse_GO_all) # Adicionar matriz de similaridade
gse_GO_bp <- pairwise_termsim(gse_GO_bp) # Adicionar matriz de similaridade
gse_GO_cc <- pairwise_termsim(gse_GO_cc) # Adicionar matriz de similaridade
gse_GO_mf <- pairwise_termsim(gse_GO_mf) # Adicionar matriz de similaridade


# Gráficos
dotplot(gse_GO_mf,
        font.size = 12,
        showCategory=5, 
        split=".sign") + facet_grid(.~.sign)

ridgeplot(gse_GO_bp,
  showCategory = 10,
  fill = "p.adjust",
  core_enrichment = TRUE,
  label_format = 50,
  decreasing = TRUE
)

cnetplot(gse_GO_mf,
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

emapplot(gse_GO_mf,
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
         group_style = "ggforce",
         label_group_style = "shawdowtext",
         label_format = 30,
         clusterFunction = stats::kmeans,
         nWords = 4,
         nCluster = NULL
)

go_sets <- data.frame(gse_GO_all@result[["Description"]])

id_go = 704
gseaplot2(gse_GO_all, 
          geneSetID = id_go,
          pvalue_table = TRUE,
          title = gse_GO_all$Description[id_go],
          color = c("#E495A5", "#86B875", "#7DB0DD"), 
          ES_geom = "line")


gsearank(gse_GO_all, 711, title = gse_GO_all[711, "Description"])


## GSEA KEGG (Precisa ser ENTREZ ID) INCOMPLETO
gse_kegg <- gseKEGG(geneList = gene_list,
                    keyType = "kegg",
                    nPerm = 10000,
                    organism = "mmu",
                    minGSSize = 3,
                    maxGSSize = 800,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "none",
                    )

kegg_sets <- data.frame(gse_kegg@result[["Description"]])

dotplot(gse_kegg, 
        showCategory = 10, 
        title = "Enriched Pathways" , 
        split=".sign") + facet_grid(.~.sign)

gse_kegg <- pairwise_termsim(gse_kegg)

emapplot(gse_kegg)

cnetplot(gse_kegg,
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

ridgeplot(gse_kegg,
          #showCategory = 10,
          #fill = "p.adjust",
          #core_enrichment = TRUE,
          #label_format = 30
          ) + labs(x = "enrichment distribution")


path <- pathview(gene.data = gene_list,
         pathway.id="mmu00071", 
         species = "mmu")

id = 67
gseaplot2(gse_kegg, 
          geneSetID = id,
          pvalue_table = TRUE,
          title = gse_kegg$Description[id],
          color = c("#E495A5", "#86B875", "#7DB0DD"), 
          ES_geom = "line")


