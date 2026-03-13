# Pacotes
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Bioconductor")
library(BiocManager)
library(tximport)
library(readr)
library(DESeq2)
library(GenomicFeatures)
library(readr)
library(dplyr)

setwd("/home/arlindo/")

# Identificar amostras
Metadados <- read_delim("Metadados.tsv", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)
Metadados <- as.data.frame(Metadados)
Metadados <- arrange(Metadados, Metadados$run_accession)
rownames(Metadados) <- Metadados$run_accession
Metadados$run_accession <- NULL

delineamento <- Metadados[ncol(Metadados)]


# Listar os quant.sf
files <- file.path("salmon_quant", rownames(delineamento), "quant.sf")
names(files) <- rownames(delineamento)
all(file.exists(files))

txdb <- makeTxDbFromGFF("ref/gencode.vM35.annotation.gtf", format = "gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)


dds <- DESeqDataSetFromTximport(txi, colData = sample_table, design = ~ condition)

dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "KO", "WT"))
summary(res)

