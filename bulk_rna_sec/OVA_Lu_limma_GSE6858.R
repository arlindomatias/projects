######################## Packages ########################
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.21")
.libPaths("~/R/library")
setwd("~/R/wd")
BiocManager::install("mouse4302.db")
library(mouse4302.db)
library(AnnotationDbi)
library("tidyverse")
library("GEOquery")
library("limma")
library("dplyr")
library("ggplot2")
library("clusterProfiler")
library("org.Mm.eg.db")
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
library("grid")
library("tibble")

######################## Adjusts ########################
par(mar=c(4.1, 4.1, 3.1, 2.1), cex.axis = 2, cex.lab = 2) 
palette(c("#3399FF", "#CC0000", "#66A61E"))
x11(width= 15, height=15) 
options(scipen=0) 

######################## Import Data ###################
# Obter informações dos dados pelo Geoquery
gset <- getGEO("GSE6858")[[1]]
fvarLabels(gset) <- make.names(fvarLabels(gset))
counts <- exprs(gset) # Matriz de expressão / Nem sempre disponível, baixar manualmente
metadados <- pData(gset) # Metadados das amostras
features <- fData(gset) # Anotação das features (genes, probes, etc.)
probes <- rownames(counts)

######
# 1. Garantir que counts tenha Probeset como rownames
counts_df <- counts %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ProbesetID")

# 2. Obter os probesets (IDs da plataforma)
probes <- counts_df$ProbesetID

# 3. Mapear probes para símbolos de genes
gene_symbols <- mapIds(mouse4302.db,
                           keys = probes,
                           column = "SYMBOL",
                           keytype = "PROBEID",
                           multiVals = "list")

# 4. Adicionar coluna de símbolo ao counts_df
counts_df$gene_name <- gene_symbols[counts_df$ProbesetID]

# 5. Remover probes sem símbolo (opcional)
counts_df <- counts_df %>% filter(!is.na(gene_name))

# 6. Resolver duplicatas por gene 
counts_agg <- counts_df %>%
  group_by(gene_name) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  as.data.frame()

counts_agg <- counts_agg %>%
  dplyr::select(gene_name, dplyr::everything(), -ProbesetID)

metadados_sub <- metadados %>%
  filter(grepl("BALB/c_OVA|BALB/c_PBS", title))

metadados_sub$grupo <- ifelse(grepl("BALB/c_OVA", metadados_sub$title),
                              "OVA", "Control")

metadados_sub$grupo <- factor(metadados_sub$grupo, levels = c("Control", "OVA"))

################### Cálculo de expressão diferencial ###################
analise <- function(data, sample_metadata, group_column = "grupo") {
  # Cria um fator com os grupos "Control" e "OVA"
  group <- factor(sample_metadata[[group_column]], levels = c("Control", "OVA"))
  
  # Design matrix
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  
  # Ajuste do modelo linear
  fit <- lmFit(data, design)
  
  # Contraste OVA vs Control
  contrast_matrix <- makeContrasts(OVA_vs_Control = OVA - Control, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contrast_matrix))
  
  # Top genes
  top_genes <- topTable(fit2, coef = "OVA_vs_Control", number = Inf)
  
  # Os rownames de top_genes já vêm do objeto 'data'
  list(
    top_genes = top_genes,
    design = design,
    groups = colnames(design),
    fit2 = fit2
  )
}

# Chamada
results <- analise(counts, metadados_sub, group_column = "grupo")
tt <- results$top_genes

# Filter
tt_filtered <- subset(tt, adj.P.Val < 0.05 & abs(logFC) > 0.5)

# Save result table
save_table <- function(table, file_name){
  if (!grepl("\\.csv$", file_name, ignore.case = TRUE)) {
    file_name <- paste0(file_name, ".csv")
  }
  write.csv(table, file = file_name, row.names = TRUE)
}

######################## Visualização ########################
## Volcano plot
# Adicionar coluna para categorização
tt$threshold <- "Not significant"
tt$threshold[tt$P.Value < 0.05 & tt$logFC > 1] <- "Upregulated"
tt$threshold[tt$P.Value < 0.05 & tt$logFC < -1] <- "Downregulated"
CairoPNG("volcano.png", width = 2000, height = 2000, res = 300)

EnhancedVolcano(
  tt, 
  lab = rownames(tt),         # Labels dos pontos (nomes dos genes)
  x = 'logFC',                 # Coluna com os dados de log2 fold change
  y = 'adj.P.Val',             # Coluna com os dados de p-valor ajustado (FDR)
  title = "OVA-induced gene expression change", # Título do gráfico
  caption = paste("Total genes:", nrow(tt)),  # Legenda com número de genes
  pCutoff = 0.05,              # Limite de p-valor (ou FDR) para destacar genes
  FCcutoff = 2,                # Limite mínimo de log2FC para considerar relevante
  pointSize = 3.0,             # Tamanho dos pontos (bolinhas) no gráfico
  xlim = c(-3,3),             # Limites do eixo X
  ylim = c(0, 6),              # Limites do eixo Y
  col = c("gray", "#3399FF", "#009933", "#CC0000"),  # Cores: NS, Down, Up, ambos
  colAlpha = 0.5,              # Transparência dos pontos
  legendLabels = c("NS", "Log2FC", "FDR", "FDR & Log2FC"), # Legenda das cores
  legendLabSize = 10,          # Tamanho do texto da legenda
  legendIconSize = 3.5,        # Tamanho dos ícones da legenda
  drawConnectors = TRUE,       # Desenhar linhas ligando labels aos pontos
  widthConnectors = 0.5,        # Espessura das linhas que ligam os labels
  max.overlaps = 30
)


grid::grid.newpage()
grid::grid.draw(volcano)  # desenha o gráfico do objeto
dev.off()

## Heatmap
deg_counts <- counts[rownames(counts) %in% rownames(tt_filtered), ]
CairoPNG("heatmap.png", width = 2000, height = 2000, res = 300)

ann <- data.frame(
  Group = factor(metadados_sub$grupo)
)
rownames(ann) <- colnames(deg_counts)

pheatmap(
  deg_counts,
  annotation_col = ann,
  show_rownames = FALSE,
  cluster_cols = TRUE,
  scale = "row"
)

dev.off()

## Boxplot
CairoPNG("boxplot.png", width = 2000, height = 2000, res = 300)

boxplot(counts, 
  boxwex=0.6, 
  notch=T, 
  main= "Normalized expression", 
  outline=FALSE, 
  las=2, 
  col= deg_counts) 

dev.off()

## Histograma de valores p e fdr 
CairoPNG("histogram.png", width = 2000, height = 2000, res = 300)

hist(
  tt$P.Value, 
  col = "grey", 
  border = "white", 
  xlab = "P-Value",
  ylab = "Number of genes", main = "P-value distribution")

dev.off()

## Densidade de expressão
CairoPNG("density.png", width = 2000, height = 2000, res = 300)

plotDensities(counts,
  group= results$groups, 
  main= "Exploression density",
  legend ="topright")

dev.off()

## QQ Plot de estatística t
t.good <- which(!is.na(results$fit2$F)) # Filtrar sondas ruins
CairoPNG("qqplot.png", width = 2000, height = 2000, res = 300)

qqt(
  results$fit2$t[t.good], 
  results$fit2$df.total[t.good], 
  main="Q-Q Plot")

dev.off()

## MD plot (log fold change vs mean log expression)
dT = decideTests(results$fit2, 
                 adjust.method="none", 
                 p.value=0.01, 
                 lfc=0)
ct <- 1 # contraste
CairoPNG("mdplot.png", width = 2000, height = 2000, res = 300)

plotMD(results$fit2, 
       column=ct, 
       status=dT[,ct], 
       legend=F, 
       pch=20, 
       cex=1)

dev.off()

######################## Enrichment ######################## parei aqui
ego_total <- enrichGO(gene = rownames(tt_filtered),
                #universe = counts$SYMBOL,
                OrgDb         = org.Rn.eg.db,   
                keyType       = "SYMBOL",
                ont           = "BP", #BP, MF, CC
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
dotplot(ego_total,
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


