######################## Pacotes ########################
library(tidyverse)
library(GEOquery)
library(limma)
library(dplyr)
library(readr)
library(ggplot2)
library(clusterProfiler)
library(org.Rn.eg.db)
library(enrichplot)
library(EnhancedVolcano)
library(pheatmap)
library(enrichplot)
library(DOSE)
library(escape)
library(readxl)
library(Cairo)
library(apeglm)
library(pathview)
library(janitor)
library(umap)
#### Lembrar que logFC positivo indica maior expressão no primeiro grupo

# Save result table
save_table <- function(table, file_name){
  if (!grepl("\\.csv$", file_name, ignore.case = TRUE)) {
    file_name <- paste0(file_name, ".csv")
  }
  write.csv(table, file = file_name, row.names = TRUE)
}

######################## Adjusts ########################
par(mar=c(4.1, 4.1, 3.1, 2.1), cex.axis = 2, cex.lab = 2) 
palette(c("#3399FF", "#CC0000", "#66A61E"))
x11(width= 15, height=15) 
options(scipen=0) 

# load series and platform data from GEO

gset <- getGEO("GSE60302", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6101", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]


# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
fxm_stress <- "XXXXXXXXXX00000XXXXXXXXXX11111"
ctrlxstress_f <- "00000XXXXX11111XXXXXXXXXXXXXXX"
ctrlxstress_m <- "XXXXXXXXXXXXXXX00000XXXXX11111"
compare <- fxm_stress
sml <- strsplit(compare, split="")[[1]]


# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("Female","Male"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

# calculate precision weights and show plot of mean-variance trend
v <- vooma(gset, design, plot=F)
# OR weights by group
#v <- voomaByGroup(gset, group=groups, design, plot=T, cex=0.1, pch=".", col=1:nlevels(gs))
v$genes <- fData(gset) # attach gene annotations

# fit linear model
fit  <- lmFit(v)

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tt <- topTable(fit2, adjust="fdr", sort.by="B", number = Inf)
tt <- subset(tt, select=c("Gene.symbol","Gene.title","Gene.ID","adj.P.Val","P.Value","t","B","logFC"))
tt <- tt[!is.na(tt$Gene.symbol) & tt$Gene.symbol != "", ]
tt_filtered <- subset(tt, adj.P.Val < 0.05 & abs(logFC) > 0.5)

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
png("hist.png", width = 2000, height = 2000, res = 300)
hist(tt$adj.P.Val, 
     col = "grey", border = "white",
     xlab = "P-adj", ylab = "Number of genes",
     main = "P-adj value distribution")
dev.off()

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=1)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
ct = 1
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

# box-and-whisker plot
png("box.png", width = 2000, height = 2000, res = 300)
ord <- order(gs)  # order samples by group
par(mar=c(7,4,2,1))
title <- paste ("GSE60302", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")
dev.off()

# expression value distribution
png("dist.png", width = 2000, height = 2000, res = 300)
par(mar=c(4,4,2,1))
title <- paste ("GSE60302", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")
dev.off()


# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 5, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=5", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)

######################## Visualização ########################
## Volcano plot
# Adicionar coluna para categorização
tt$threshold <- "Not significant"
tt$threshold[tt$adj.P.Val < 0.05 & tt$logFC > 1] <- "Upregulated"
tt$threshold[tt$adj.P.Val < 0.05 & tt$logFC < -1] <- "Downregulated"
CairoPNG("volcano.png", width = 3000, height = 2000, res = 300)

v <- EnhancedVolcano(
  tt, 
  lab = tt$Gene.symbol,         # Labels dos pontos (nomes dos genes)
  x = 'logFC',                 # Coluna com os dados de log2 fold change
  y = 'adj.P.Val',             # Coluna com os dados de p-valor ajustado (FDR)
  title = "Female vs Male", # Título do gráfico
  caption = paste("Total genes:", nrow(tt)),  # Legenda com número de genes
  pCutoff = 0.05,              # Limite de p-valor (ou FDR) para destacar genes
  FCcutoff = 1,                # Limite mínimo de log2FC para considerar relevante
  pointSize = 1.5,             # Tamanho dos pontos (bolinhas) no gráfico
  xlim = c(-3,3),             # Limites do eixo X
  ylim = c(0, 4),              # Limites do eixo Y
  col = c("gray", "#3399FF", "#009933", "#CC0000"),  # Cores: NS, Down, Up, ambos
  colAlpha = 0.5,              # Transparência dos pontos
  legendLabels = c("NS", "Log2FC", "FDR", "FDR & Log2FC"), # Legenda das cores
  legendLabSize = 9,          # Tamanho do texto da legenda
  legendIconSize = 3,        # Tamanho dos ícones da legenda
  labSize = 4,
  drawConnectors = TRUE,       # Desenhar linhas ligando labels aos pontos
  widthConnectors = 0.5,        # Espessura das linhas que ligam os labels
  max.overlaps = 20
)

grid::grid.newpage()
grid::grid.draw(v)  # desenha o gráfico do objeto
dev.off()

## Heatmap
top <- tt_filtered[order(tt_filtered$B, decreasing = TRUE)[1:50], ]
rowData <- subset(top, select=c("Gene.symbol"))
rowData[] <- lapply(rowData, function(x) {
  if(is.numeric(x)) factor(x) else x
})

deg_counts <- ex[rownames(ex) %in% rownames(top), ]
deg_counts <- as.data.frame(deg_counts)

deg_counts <- merge(
  x = deg_counts,
  y = rowData,
  by = "row.names",  # transforma os rownames em coluna temporária
  all = TRUE          # TRUE = mantém todos (equivalente a full join)
)
rownames(deg_counts) <- deg_counts$Row.names
deg_counts$Row.names <- NULL
deg_counts <- deg_counts %>%
  group_by(Gene.symbol) %>%
  summarise(
    Gene.symbol = dplyr::first(Gene.symbol),
    across(where(is.numeric), mean)
  ) %>% as.data.frame()

rownames(deg_counts) <- deg_counts$Gene.symbol
deg_counts <- deg_counts[,2:11]
metadados <- pData(gset) # Metadados das amostras
colData <- subset(metadados, select=c("group"))


png("heatmap.png", width = 2000, height = 2000, res = 300)
heatmap <- pheatmap(deg_counts, 
                    annotation_col = colData, 
                    show_rownames = T, 
                    cluster_cols = F, 
                    scale = "row")
dev.off()


######################## Enrichment ######################## 
up_genes <- tt_filtered[tt_filtered$logFC > 0, ]
down_genes <- tt_filtered[tt_filtered$logFC < 0, ]

ego_total <- enrichGO(gene = tt_filtered$Gene.symbol,
                      #universe = counts$SYMBOL,
                      OrgDb         = org.Rn.eg.db,   
                      keyType       = "SYMBOL",
                      ont           = "MF", #BP, MF, CC
                      pAdjustMethod = "none",
                      pvalueCutoff  = 0.05)


# Gráficos
png("egomf.png", width = 2000, height = 2000, res = 300)
dotplot(ego_total,
          x = "GeneRatio",   # confirme a capitalização correta
          font.size = 10,
          showCategory = 20) +
    scale_color_gradientn(colors = c("#0099ff", "#6600ff", "#ff0000"))
dev.off()



# Preparar matriz de contagem com ENTREZ
genes <- topTable(fit2, adjust="fdr", sort.by="B", number = Inf)
genes <- subset(genes, select=c("Gene.symbol","Gene.ID"))
genes <- genes[!is.na(genes$Gene.ID) & genes$Gene.ID != "", ]
entrez_counts <- merge(
  x = ex,
  y = genes,
  by = "row.names",  # transforma os rownames em coluna temporária
  all = TRUE          # TRUE = mantém todos (equivalente a full join)
)
entrez_counts$Row.names <- NULL

# Remover duplicatas por ENTREZID
entrez_single <- merge(entrez_counts, tt, by = c("Gene.ID", "Gene.symbol"))
entrez_single <- entrez_single[!is.na(entrez_single$logFC), ]

entrez_single <- entrez_single %>%
  filter(!is.na(Gene.ID)) %>%  # remove NA
  group_by(Gene.ID) %>%
  summarise(
    Gene.symbol = dplyr::first(Gene.symbol),
    across(where(is.numeric), mean)
  )
entrez_single <- as.data.frame(entrez_single)
rownames(entrez_single) <- entrez_single$Gene.ID

# Criar lista
gene_list <- sort(entrez_single$logFC, decreasing = TRUE)
names(gene_list) <- entrez_single$Gene.ID

## GSEA usando GO
gse_GO_all <- gseGO(geneList= gene_list, 
                    ont ="ALL", 
                    keyType = "ENTREZID", 
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    pvalueCutoff = 0.05, 
                    verbose = TRUE, 
                    OrgDb = "org.Rn.eg.db",
                    pAdjustMethod = "BH")
gse_GO_all <- pairwise_termsim(gse_GO_all) # Adicionar matriz de similaridade

# Filtrar resultados da ontologia BP (Biological Process)
gse_BP <- gse_GO_all
gse_BP@result <- gse_BP@result[gse_BP@result$ONTOLOGY == "BP", ]

# Para CC (Cellular Component)
gse_CC <- gse_GO_all
gse_CC@result <- gse_CC@result[gse_CC@result$ONTOLOGY == "CC", ]

# Para MF (Molecular Function)
gse_MF <- gse_GO_all
gse_MF@result <- gse_MF@result[gse_MF@result$ONTOLOGY == "MF", ]

# Gráficos
png("gsecc.png", width = 2000, height = 2000, res = 300)
dotplot(gse_CC,
        font.size = 12,
        showCategory=10, 
        split=".sign") + facet_grid(.~.sign)
dev.off()

png("gsebp.png", width = 2000, height = 2000, res = 300)
dotplot(gse_BP,
        font.size = 12,
        showCategory=10, 
        split=".sign") + facet_grid(.~.sign)
dev.off()

png("gsemf.png", width = 2000, height = 2000, res = 300)
dotplot(gse_MF,
        font.size = 12,
        showCategory=10, 
        split=".sign") + facet_grid(.~.sign)
dev.off()

png("cnet.png", width = 2000, height = 2000, res = 300)
cnetplot(gse_GO_all,
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
dev.off()


png("ema.png", width = 2000, height = 2000, res = 300)
emapplot(gse_GO_all,
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
dev.off()

go_sets <- data.frame(gse_GO_all@result[["Description"]])


id_go = 52
png("gsea.png", width = 2500, height = 1000, res = 300)
gseaplot2(gse_GO_all, 
          geneSetID = id_go,
          pvalue_table = TRUE,
          title = gse_GO_all$Description[id_go],
          color = c("#E495A5", "#86B875", "#7DB0DD"), 
          ES_geom = "line")
dev.off()

png("gsearank52.png", width = 2500, height = 1000, res = 300)
gsearank(gse_GO_all, id_go, title = gse_GO_all[id_go, "Description"])
dev.off()

## GSEA KEGG (Precisa ser ENTREZ ID) INCOMPLETO
gse_kegg <- gseKEGG(geneList = gene_list,
                    keyType = "kegg",
                    nPerm = 10000,
                    organism = "rno",
                    minGSSize = 3,
                    maxGSSize = 800,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
)

kegg_sets <- data.frame(gse_kegg@result[["Description"]])

png("gsekegg.png", width = 2000, height = 2000, res = 300)
dotplot(gse_kegg, 
        showCategory = 20, 
        title = "Enriched Pathways" , 
        split=".sign") + facet_grid(.~.sign)
dev.off()

gse_kegg <- pairwise_termsim(gse_kegg)

png("emakegg.png", width = 2000, height = 2000, res = 300)
emapplot(gse_kegg)
dev.off()

png("cnetkegg.png", width = 2000, height = 2000, res = 300)
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
dev.off()

png("ridgekegg.png", width = 2000, height = 2000, res = 300)
ridgeplot(gse_kegg,
          #showCategory = 10,
          #fill = "p.adjust",
          #core_enrichment = TRUE,
          #label_format = 30
) + labs(x = "enrichment distribution")
dev.off()

png("pathview.png", width = 2000, height = 2000, res = 300)
path <- pathview(gene.data = gene_list,
                 pathway.id="rno01212 ", 
                 species = "rno")
dev.off()

id = 37
png("gseakegg37.png", width = 2500, height = 1000, res = 300)
gseaplot2(gse_kegg, 
          geneSetID = id,
          pvalue_table = TRUE,
          title = gse_kegg$Description[id],
          color = c("#E495A5", "#86B875", "#7DB0DD"), 
          ES_geom = "line")
dev.off()

save_table(tt, 'tt')
