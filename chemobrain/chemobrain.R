######################## Pacotes ########################
library(tidyverse)
library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(EnhancedVolcano)
library(pheatmap)
library(enrichplot)
library(DOSE)
library(escape)
library(readxl)
library(Cairo)
library(DESeq2)
library(apeglm)
library(pathview)

# Função para salvar tabela com verificação de extensão
save <- function(tabela, nome_arquivo){
  if (!grepl("\\.csv$", nome_arquivo, ignore.case = TRUE)) {
    nome_arquivo <- paste0(nome_arquivo, ".csv")
  }
  write.csv(tabela, file = nome_arquivo, row.names = TRUE)
}

par(mar=c(4.1, 4.1, 3.1, 2.1), cex.axis = 0.3, cex.lab = 2) 
palette(c("#3399FF", "#CC0000", "#66A61E"))
dev.new(width= 100, height=100, noRStudioGD = TRUE)
options(scipen=0)

# Baixa metadados e arquivos associados
gse <- getGEO("GSE295281", GSEMatrix = TRUE)[[1]]
fvarLabels(gse) <- make.names(fvarLabels(gse))
info <- pData(gse) # Metadados das amostras
counts <- read.csv("GSE295281_ACT_combined_count_matrix_v02.csv.gz", 
                   row.names = 2)
colnames(counts)[colnames(counts) == "Gene.Name"] <- "SYMBOL"

# Extrair grupos dos metadados
group <- make.names(gse$`characteristics_ch1.4`)
group <- factor(group)
levels(group)[levels(group) == "treatment..2.mg.kg.doxorubicin..A...50.mg.kg.cyclophosphamide..C...5.mg.kg.paclitaxel..T."] <- "CHEM"
levels(group)[levels(group) == "treatment..saline..0.9..sodium.chloride."] <- "VEH"

# Criar colData e garantir os nomes das amostras
colData <- data.frame(condition = group)
rownames(colData) <- rownames(info)

# Ajustar colnames de counts
counts_raw <- counts %>%
  group_by(SYMBOL) %>%
  summarise(
    SYMBOL = dplyr::first(SYMBOL),
    across(where(is.numeric), mean)
  ) %>% data.frame()

rownames(counts_raw) <- counts_raw$SYMBOL
counts_raw$SYMBOL <- NULL
colnames(counts_raw) <- info$geo_accession


# Obter IDs ENSEMBL
keytable <- data.frame("ENSEMBL" = rownames(counts), "SYMBOL" = counts$SYMBOL)


# Mapear ENSEMBL -> ENTREZ
ensembl_to_entrez <- bitr(keytable$ENSEMBL,
                          fromType = "ENSEMBL",
                          toType = "ENTREZID",
                          OrgDb = org.Mm.eg.db)
ensembl_to_entrez <- na.omit(ensembl_to_entrez)

keytable <- merge(ensembl_to_entrez, keytable, by = "ENSEMBL")

# Garantir que são inteiros
id_count <- round(as.matrix(counts_raw))
storage.mode(id_count) <- "integer"
id_count <- data.frame(id_count)

# Criar objeto DESeq2
dds <- DESeqDataSetFromMatrix(countData = id_count,
                              colData = colData,
                              design = ~ condition)

dds <- DESeq(dds, fitType = "local")


result <- results(dds)
tt <- data.frame(result)
tt$SYMBOL <- rownames(tt)
tt2 <- merge(keytable, tt, by = "SYMBOL")

# Filtro: genes com FDR < 0.05 e |logFC| > 1
tt_filtered <- subset(tt2, padj < 0.05 & abs(log2FoldChange) > 0.5)
tt_filtered_up <- subset(tt2, padj < 0.05 & log2FoldChange > 0.5)
tt_filtered_down <- subset(tt2, padj < 0.05 & log2FoldChange < -0.5)
tt_filtered2 <- subset(tt2, padj < 0.01 & abs(log2FoldChange) > 2)

plotPCA(vst(dds))

write.table(tt, "C:/Users/AsRock/Documents/tt.txt")

## Volcano plot
# Adicionar coluna para categorização
tt2$threshold <- "Not significant"
tt2$threshold[tt2$padj < 0.05 & tt2$log2FoldChange > 1] <- "Upregulated"
tt2$threshold[tt2$padj < 0.05 & tt2$log2FoldChange < -1] <- "Downregulated"

# Renderizar
v           <- EnhancedVolcano(tt2, 
                               lab = tt2$SYMBOL,         # Labels dos pontos (nomes dos genes)
                               x = 'log2FoldChange',                 # Coluna com os dados de log2 fold change
                               y = 'padj',             # Coluna com os dados de p-valor ajustado (FDR)
                               xlab = bquote(~Log[2]~ 'fold change'),
                               title = "Chemotherapy-induced gene expression change", # Título do gráfico
                               caption = paste("Total genes:", nrow(tt)),  # Legenda com número de genes
                               pCutoff = 0.01,              # Limite de p-valor (ou FDR) para destacar genes
                               FCcutoff = 2,                # Limite mínimo de log2FC para considerar relevante
                               pointSize = 3.0,             # Tamanho dos pontos (bolinhas) no gráfico
                               labSize = 3,
                               xlim = c(-25, 25),             # Limites do eixo X
                               ylim = c(0, 40),              # Limites do eixo Y
                               col = c("gray", "#3399FF", "#009933", "#CC0000"),  # Cores: NS, Down, Up, ambos
                               colAlpha = 0.5,              # Transparência dos pontos
                               legendLabels = c("NS", "Log2FC", "FDR", "FDR & Log2FC"), # Legenda das cores
                               legendLabSize = 10,          # Tamanho do texto da legenda
                               legendIconSize = 10,        # Tamanho dos ícones da legenda
                               drawConnectors = TRUE,       # Desenhar linhas ligando labels aos pontos
                               widthConnectors = 0.3,        # Espessura das linhas que ligam os labels
                               max.overlaps = 20
)

#heatmap
deg_counts <- counts_raw[rownames(counts_raw) %in% tt_filtered2$SYMBOL, ]
heatmap <- pheatmap(deg_counts, 
                    annotation_col = colData, 
                    show_rownames = T, 
                    cluster_cols = F, 
                    scale = "row")




boxplot <- boxplot(tt_counts, 
                   boxwex=0.6, 
                   notch=T, 
                   main= "Normalized expression", 
                   outline=FALSE, 
                   las=2, 
                   col= tt_counts) 

## Histograma de valores p e fdr 
hist(tt$padj, col = "grey", border = "white", xlab = "FDR",
     ylab = "Number of genes", main = "FDR distribution")

tt_filtered <- merge(keytable, tt_filtered, by = "ENSEMBL")

ego_total <- enrichGO(gene = tt_filtered$SYMBOL,
                #universe = counts$SYMBOL,
                OrgDb         = org.Mm.eg.db,   
                keyType       = "SYMBOL",
                ont           = "CC", #BP, MF, CC
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05)

ego_up <- enrichGO(gene = tt_filtered_up$SYMBOL,
                      #universe = tt$SYMBOL,
                      OrgDb         = org.Mm.eg.db,   
                      keyType       = "SYMBOL",
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

ego_down <- enrichGO(gene = tt_filtered_down$SYMBOL,
                   #universe = tt$SYMBOL,
                   OrgDb         = org.Mm.eg.db,   
                   keyType       = "SYMBOL",
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
# Gráficos
dotplot(ego_down,
        x = "geneRatio",
        font.size = 12,
        showCategory = 10) + scale_color_gradientn(colors = c("#0099ff", "#6600ff", "#ff0000"))


barplot(ego_total, 
        showCategory=50,
        font.size = 8,
        title = "Gene Enrichment Analysis",
        label_format = 60,) + scale_color_gradientn(colors = c("#0099ff", "#6600ff", "#ff0000"))




# Preparar matriz de contagem com ENTREZ
entrez_count <- counts
entrez_count$"ENSEMBL" <- rownames(entrez_count)
entrez_count <- merge(ensembl_to_entrez, entrez_count, by = "ENSEMBL")

# Remover duplicatas por ENTREZID
entrez_single <- entrez_count %>%
  group_by(ENTREZID) %>%
  summarise(
    SYMBOL = dplyr::first(SYMBOL),
    across(where(is.numeric), mean)
  )
keytable <- entrez_count[!duplicated(entrez_count$SYMBOL), 1:3]


# Filtrar apenas os genes diferencialmente expressos
entrez_sig <- entrez_single[entrez_single$SYMBOL %in% tt_filtered$SYMBOL, ]


## Criar lista de ranking dos genes
tt_entrez_single <- tt %>%
  group_by(ENTREZID) %>%
  summarise(
    ENTREZID = dplyr::first(ENTREZID),
    across(where(is.numeric), mean)
  )

tt_entrez_single  <- tt_entrez_single[complete.cases(tt_entrez_single), ]
gene_list <- sort(tt_entrez_single$log2FoldChange, decreasing = TRUE)
names(gene_list) <- tt_entrez_single$ENTREZID


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
  label_format = 100,
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

sets <- data.frame(gse_GO_all@result[["Description"]])

gseaplot(gse_GO_all,
         geneSetID = 21,
         by = "runningScore", 
         title = gse_GO_all$Description[21],
         color = "black",
         color.line = "green",
         color.vline = "#FA5860",
         )

gseaplot2(gse_GO_all, 
          geneSetID = 648,
          pvalue_table = TRUE,
          title = gse_GO_all$Description[648],
          color = c("#E495A5", "#86B875", "#7DB0DD"), 
          ES_geom = "line")


gsearank(gse_GO_all, 20, title = gse_GO_all[20, "Description"])


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
         showCategory = 30,
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

ridgeplot(gse_kegg) + labs(x = "enrichment distribution")


path <- pathview(gene.data = gene_list,
         pathway.id="mmu04979", 
         species = "mmu")

path2 <- pathview(gene.data= gene_list, 
                  pathway.id="mmu04979", 
                  species = "mmu", 
                  kegg.native = F)

knitr::include_graphics("mmu04979.pathview.png")


gseaplot2(gse_kegg, 
          geneSetID = 23,
          pvalue_table = TRUE,
          title = gse_kegg$Description[23],
          color = c("#E495A5", "#86B875", "#7DB0DD"), 
          ES_geom = "line")


