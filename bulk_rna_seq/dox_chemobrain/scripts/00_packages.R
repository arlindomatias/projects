# Install BiocManager if necessary
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# CRAN packages
cran_packages <- c(
  "readr", "dplyr", "tidyr", "ggplot2", "gridExtra",
  "pheatmap", "ggrepel", "readxl", "tibble", "ggforce", "shadowtext"
)

# Bioconductor packages
bioc_packages <- c(
  "GEOquery", "DESeq2", "clusterProfiler", "org.Mm.eg.db",
  "EnhancedVolcano", "enrichplot", "org.Hs.eg.db", "pathview"
)

if (!require(bioc_packages, quietly = TRUE))
  BiocManager::install(bioc_packages)

if (!require(cran_packages, quietly = TRUE))
  install.packages(cran_packages)

# Load packages (each one only once)
invisible(lapply(
  c(cran_packages, bioc_packages),
  library,
  character.only = TRUE
))

# Project configuration
# 1. Project configuration -----------------------------------
project_dir <- path.expand("~/github/projects/bulk_rna_seq/dox_chemobrain")

setwd(project_dir)

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

par(mar=c(4.1, 4.1, 3.1, 2.1), cex.axis = 0.3, cex.lab = 2) 
palette(c("#3399FF", "#CC0000", "#66A61E"))
dev.new(width= 10, height=10, noRStudioGD = TRUE)
options(scipen=0)
