# Naive Bayesian classifier against SILVA 138.2 (to genus), then exact-match species.
log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages(library(dada2))

set.seed(snakemake@params[["seed"]])
seqtab <- readRDS(snakemake@input[["seqtab"]])
seqs <- colnames(seqtab)

taxa <- assignTaxonomy(seqs, snakemake@input[["trainset"]],
                       multithread = snakemake@threads, tryRC = TRUE)
taxa <- addSpecies(taxa, snakemake@input[["species"]], tryRC = TRUE)

out <- data.frame(asv = paste0("ASV", seq_along(seqs)), taxa, sequence = seqs,
                  row.names = NULL, check.names = FALSE)
cat("Assigned fraction per rank:\n")
print(colMeans(!is.na(taxa)))
write.table(out, snakemake@output[["taxa"]], sep = "\t", quote = FALSE, row.names = FALSE)
