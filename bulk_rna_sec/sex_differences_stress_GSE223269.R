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
gse <- getGEO("GSE223269", GSEMatrix = TRUE)[[1]]
fvarLabels(gse) <- make.names(fvarLabels(gse))

counts <- read_delim("GSE223269_Raw_gene_counts_matrix.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)
colnames(counts)[1] <- "SYMBOL"

# lê sua tabela de metadados baixada do SRA (provavelmente um .tsv)
meta <- read_csv("SraRunTable.csv")

# cria uma versão simplificada só com o que interessa
info <- meta %>%
  select(Run, `Sample Name`, phenotype, tissue) %>%
  mutate(Sample_ID = paste0("sample_", row_number()))
