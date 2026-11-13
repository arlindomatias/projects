# Aggregate quality profiles of primer-trimmed reads (fwd = 338F side, rev = 806R side).
log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({
  library(dada2)
  library(ggplot2)
})

p_fwd <- plotQualityProfile(snakemake@input[["fwd"]], aggregate = TRUE) +
  ggtitle("Forward reads (after 338F)")
p_rev <- plotQualityProfile(snakemake@input[["rev"]], aggregate = TRUE) +
  ggtitle("Reverse reads (after 806R)")

png(snakemake@output[[1]], width = 12, height = 5, units = "in", res = 150)
gridExtra::grid.arrange(p_fwd, p_rev, ncol = 2)
dev.off()
