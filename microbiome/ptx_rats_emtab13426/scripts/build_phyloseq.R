# Assemble the phyloseq object (ASV counts, taxonomy, midpoint-rooted tree, sample metadata).
log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({
  library(phyloseq)
  library(Biostrings)
  library(ape)
  library(phangorn)
})

seqtab <- readRDS(snakemake@input[["seqtab"]])
taxa <- read.delim(snakemake@input[["taxa"]], check.names = FALSE)
meta <- read.delim(snakemake@input[["metadata"]], colClasses = "character")

asv_ids <- paste0("ASV", seq_len(ncol(seqtab)))
stopifnot(identical(taxa$sequence, colnames(seqtab)), identical(taxa$asv, asv_ids))
seqs <- DNAStringSet(colnames(seqtab)); names(seqs) <- asv_ids
colnames(seqtab) <- asv_ids

rownames(meta) <- meta$sample
meta$group <- factor(meta$group, levels = c("Control", "PTX_WK2", "PTX_WK4"))
meta$week <- suppressWarnings(as.integer(meta$week))
stopifnot(setequal(rownames(seqtab), rownames(meta)))

tax_mat <- as.matrix(taxa[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")])
rownames(tax_mat) <- taxa$asv

tree <- midpoint(read.tree(snakemake@input[["tree"]]))

ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE), tax_table(tax_mat),
               sample_data(meta[rownames(seqtab), ]), phy_tree(tree), refseq(seqs))
print(ps)
saveRDS(ps, snakemake@output[["ps"]])

counts <- data.frame(asv = asv_ids, t(seqtab), taxa[, -1], check.names = FALSE)
write.table(counts, snakemake@output[["counts"]], sep = "\t", quote = FALSE, row.names = FALSE)
