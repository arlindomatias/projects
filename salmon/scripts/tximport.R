# Agrega quant.sf do Salmon para nível de gene com tximport e grava as tabelas.
# Uso: Rscript tximport.R samples.tsv <dir_salmon> tx2gene.tsv \
#        counts.tsv counts_lengthScaledTPM.tsv tpm.tsv length.tsv txi.rds
suppressPackageStartupMessages({
  library(tximport)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 8)
samples_tsv <- args[1]; salmon_dir <- args[2]; tx2gene_tsv <- args[3]
out_counts <- args[4]; out_scaled <- args[5]; out_tpm <- args[6]; out_length <- args[7]
out_rds <- args[8]

samples <- read_tsv(samples_tsv, show_col_types = FALSE)
files <- file.path(salmon_dir, samples$run, "quant.sf")
names(files) <- samples$run
stopifnot(all(file.exists(files)))

tx2gene <- read_tsv(tx2gene_tsv, show_col_types = FALSE)
genes <- unique(tx2gene[, c("gene_id", "gene_name")])

txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, c("transcript_id", "gene_id")])
txi_scaled <- tximport(files, type = "salmon", tx2gene = tx2gene[, c("transcript_id", "gene_id")],
                       countsFromAbundance = "lengthScaledTPM")

write_matrix <- function(m, path, digits) {
  df <- data.frame(gene_id = rownames(m), round(m, digits), check.names = FALSE)
  df <- merge(genes, df, by = "gene_id", sort = FALSE)
  df <- df[order(df$gene_id), ]
  write_tsv(df, path)
}
write_matrix(txi$counts, out_counts, 3)
write_matrix(txi_scaled$counts, out_scaled, 3)
write_matrix(txi$abundance, out_tpm, 4)
write_matrix(txi$length, out_length, 3)

# Objeto completo para DESeqDataSetFromTximport(txi, colData, design).
saveRDS(txi, out_rds)

cat(sprintf("%d genes x %d amostras\n", nrow(txi$counts), ncol(txi$counts)))
cat("Reads atribuídos por amostra:\n")
print(round(colSums(txi$counts)))
