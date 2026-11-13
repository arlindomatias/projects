# DADA2: error models, sample inference (pooled = FALSE), merging, chimera removal.
# Writes the chimera-free sequence table, ASV FASTA, and a per-sample read tracking table.
log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({
  library(dada2)
  library(ggplot2)
  library(jsonlite)
})

set.seed(snakemake@params[["seed"]])
threads <- snakemake@threads
samples <- unlist(snakemake@params[["samples"]])
filt_f <- setNames(unlist(snakemake@input[["fwd"]]), samples)
filt_r <- setNames(unlist(snakemake@input[["rev"]]), samples)

# Error profiles differ between sequencing runs (2x300 vs 2x250) and between submitted R1 and
# R2. Mate orientation differs between samples, so the "fwd" file comes from R1 in some and
# from R2 in others. Error models are therefore learned, and samples denoised, separately per
# batch = sequencing run x orientation. Merged ASVs span the same amplicon in every batch.
cut_json <- setNames(lapply(unlist(snakemake@input[["cutadapt"]]), fromJSON), samples)
orientation <- sapply(cut_json, `[[`, "orientation")
seq_run <- setNames(unlist(snakemake@params[["seq_runs"]]), samples)
batch <- paste0("run ", seq_run, ", ", orientation)
cat("Samples per batch:\n"); print(table(batch))

# dada() returns a bare object instead of a list when given a single sample.
as_list <- function(x, s) if (inherits(x, "dada")) setNames(list(x), s) else x
dada_f <- dada_r <- list()
err_plots <- list()
for (b in sort(unique(batch))) {
  s <- samples[batch == b]
  err_f <- learnErrors(filt_f[s], multithread = threads, randomize = TRUE)
  err_r <- learnErrors(filt_r[s], multithread = threads, randomize = TRUE)
  err_plots[[paste(b, "fwd")]] <- plotErrors(err_f, nominalQ = TRUE) +
    ggtitle(sprintf("Forward reads, %s batch (n = %d)", b, length(s)))
  err_plots[[paste(b, "rev")]] <- plotErrors(err_r, nominalQ = TRUE) +
    ggtitle(sprintf("Reverse reads, %s batch (n = %d)", b, length(s)))
  dada_f[s] <- as_list(dada(filt_f[s], err = err_f, multithread = threads), s)
  dada_r[s] <- as_list(dada(filt_r[s], err = err_r, multithread = threads), s)
}
dada_f <- dada_f[samples]
dada_r <- dada_r[samples]

png(snakemake@output[["errors"]], width = 12, height = 6 * ceiling(length(err_plots) / 2),
    units = "in", res = 100)
gridExtra::grid.arrange(grobs = err_plots, ncol = 2)
dev.off()

mergers <- mergePairs(dada_f, filt_f, dada_r, filt_r, verbose = TRUE)

seqtab <- makeSequenceTable(mergers)
cat("Merged sequence length distribution:\n")
print(table(nchar(getSequences(seqtab))))

seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = threads, verbose = TRUE)
cat(sprintf("Non-chimeric fraction of reads: %.4f\n", sum(seqtab_nochim) / sum(seqtab)))
cat(sprintf("ASVs: %d (before chimera removal: %d)\n", ncol(seqtab_nochim), ncol(seqtab)))

# Name ASVs by decreasing total abundance; keep sequences in the FASTA.
seqtab_nochim <- seqtab_nochim[, order(-colSums(seqtab_nochim)), drop = FALSE]
saveRDS(seqtab_nochim, snakemake@output[["seqtab"]])
seqs <- colnames(seqtab_nochim)
writeLines(paste0(">ASV", seq_along(seqs), "\n", seqs), snakemake@output[["fasta"]])

# Read tracking: raw pairs -> primers found -> filtered -> denoised -> merged -> non-chimeric.
cut_stats <- t(sapply(cut_json, function(j)
  c(raw = j$read_counts$input, primers = j$read_counts$output)))
filt_stats <- do.call(rbind, lapply(unlist(snakemake@input[["stats"]]), read.delim))
get_n <- function(x) sum(getUniques(x))
track <- data.frame(
  sample = samples,
  seq_run = seq_run,
  orientation = orientation,
  raw = cut_stats[, "raw"],
  primers_found = cut_stats[, "primers"],
  filtered = filt_stats$reads.out[match(samples, filt_stats$sample)],
  denoised_fwd = sapply(dada_f, get_n),
  denoised_rev = sapply(dada_r, get_n),
  merged = sapply(mergers, get_n),
  nonchim = rowSums(seqtab_nochim)[samples],
  row.names = NULL
)
track$retained_pct <- round(100 * track$nonchim / track$raw, 1)
print(track)
write.table(track, snakemake@output[["track"]], sep = "\t", quote = FALSE, row.names = FALSE)
