# ============================================================
# SCRIPT 00: configuração comum (carregado por todos os passos)
# Os scripts rodam pela diretiva `script:` do Snakemake, que cria o objeto `snakemake`.
# ============================================================
if (!exists("snakemake")) {
  stop("Run this script through Snakemake (`snakemake --use-conda --cores 1`).")
}

# Mensagens e saída do R vão para o log da regra.
log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

cfg <- snakemake@config

for (d in c("results", "plots")) dir.create(d, showWarnings = FALSE)

options(timeout = 600, lifecycle_verbosity = "quiet")
set.seed(cfg$seed %||% 123)
if (!interactive()) pdf(NULL)

ref_level <- cfg$groups$reference
trt_level <- cfg$groups$treatment
group_colors <- unlist(cfg$groups$colors)

# Rótulo legível do gene: símbolo, ou o Ensembl ID quando não houver símbolo.
gene_label <- function(symbol, gene_id) {
  ifelse(is.na(symbol) | symbol == "", gene_id, symbol)
}
