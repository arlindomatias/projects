# Install BiocManager if necessary
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# CRAN packages
cran_packages <- c(
  "readr", "dplyr", "tidyr", "ggplot2",
  "pheatmap", "ggrepel", "readxl"
)

# Bioconductor packages
bioc_packages <- c(
  "GEOquery", "DESeq2", "clusterProfiler", "org.Mm.eg.db",
  "EnhancedVolcano"
)

# Load packages (each one only once)
invisible(lapply(
  c(cran_packages, bioc_packages),
  library,
  character.only = TRUE
))

# Project configuration
# 1. Project configuration -----------------------------------
project_dir <- path.expand(
  "~/github/private/portfolio/dox_chemobrain"
)

setwd(project_dir)

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

      