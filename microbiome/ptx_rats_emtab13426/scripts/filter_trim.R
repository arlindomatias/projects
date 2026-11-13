# DADA2 quality filtering and truncation for one sample.
log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages(library(dada2))

out <- filterAndTrim(
  fwd = snakemake@input[["fwd"]], filt = snakemake@output[["fwd"]],
  rev = snakemake@input[["rev"]], filt.rev = snakemake@output[["rev"]],
  truncLen = unlist(snakemake@params[["trunc_len"]]),
  maxEE = unlist(snakemake@params[["max_ee"]]),
  maxN = 0, truncQ = 2, rm.phix = TRUE,
  compress = TRUE, multithread = FALSE
)
print(out)
write.table(
  data.frame(sample = snakemake@wildcards[["sample"]],
             reads.in = out[1, "reads.in"], reads.out = out[1, "reads.out"]),
  snakemake@output[["stats"]], sep = "\t", quote = FALSE, row.names = FALSE
)
