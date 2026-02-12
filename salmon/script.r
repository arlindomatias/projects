# Pacotes
library(tximport)
library(readr)
library(DESeq2)
library(GenomicFeatures)
library(txdbmaker)

# Identificar amostras
samples <- c(
  "SRR33247676","SRR33247677"
)

setwd("/home/arlindo/github/teste")

# Listar os quant.sf
files <- file.path("salmon_quant", samples, "quant.sf")
names(files) <- samples
files
all(file.exists(files))

# Mapeamento transcrito
txdb <- txdbmaker::makeTxDbFromGFF("gencode.vM35.annotation.gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")
head(tx2gene)

# Importar os dados
txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM"
)

str(txi)

# Criar DF
colData <- data.frame(
  row.names = samples,
  condition = as.factor(c("control", "treated"))
)


# Análise DE
dds <- DESeqDataSetFromTximport(
  txi,
  colData = colData,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "treated", "control"))

res <- res[order(res$padj), ]
head(res)

summary(res)
plotMA(res)




