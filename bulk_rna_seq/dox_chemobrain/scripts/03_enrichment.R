# ============================================================
# SCRIPT 03: FUNCTIONAL ENRICHMENT AND GSEA
# Dataset: GSE295281
# Groups: VEH (Saline) vs CHEM (AC-T chemotherapy)
#
# Purpose:
# 1. Load differential expression results from Script 02.
# 2. Map gene symbols to Entrez Gene IDs.
# 3. Perform GO over-representation analysis.
# 4. Perform GO Gene Set Enrichment Analysis (GSEA).
# 5. Perform KEGG Gene Set Enrichment Analysis.
# 6. Generate enrichment plots and export results.
#
# IMPORTANT:
# The differential expression analysis in Script 02 is
# exploratory because the original expression matrix scale
# has not yet been verified.
# ============================================================


# 1. Package installation and loading -------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c(
  "ggplot2",
  "dplyr",
  "readr"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

bioc_packages <- c(
  "clusterProfiler",
  "org.Mm.eg.db",
  "enrichplot",
  "AnnotationDbi",
  "pathview"
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(AnnotationDbi)
library(ggplot2)
library(dplyr)
library(readr)


# 2. Project configuration -----------------------------------

project_dir <- path.expand(
  "~/github/private/portfolio/rna_seq"
)

results_dir <- file.path(
  project_dir,
  "results"
)

input_dir <- file.path(
  results_dir,
  "02_quality_control"
)

enrichment_dir <- file.path(
  results_dir,
  "03_enrichment"
)

plot_dir <- file.path(
  enrichment_dir,
  "plots"
)

dir.create(
  enrichment_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 3. Load differential expression results --------------------

message("Loading Script 02 results...")

all_genes <- read.csv(
  file.path(input_dir, "02_limma_all_genes.csv"),
  stringsAsFactors = FALSE
)

significant_genes <- read.csv(
  file.path(
    input_dir,
    "02_limma_significant_genes_FDR05.csv"
  ),
  stringsAsFactors = FALSE
)

# Validate required columns.
required_columns <- c(
  "ENSEMBL",
  "SYMBOL",
  "logFC",
  "P.Value",
  "adj.P.Val"
)

if (!all(required_columns %in% colnames(all_genes))) {
  stop("Required columns are missing from the differential expression table.")
}

if (!all(required_columns %in% colnames(significant_genes))) {
  stop("Required columns are missing from the significant gene table.")
}

# Keep rows with valid gene symbols.
all_genes <- all_genes[
  !is.na(all_genes$SYMBOL) &
    all_genes$SYMBOL != "" &
    !is.na(all_genes$logFC),
]

significant_genes <- significant_genes[
  !is.na(significant_genes$SYMBOL) &
    significant_genes$SYMBOL != "" &
    !is.na(significant_genes$logFC),
]

# Separate upregulated and downregulated genes.
genes_up <- significant_genes[
  significant_genes$logFC > 0,
  ,
  drop = FALSE
]

genes_down <- significant_genes[
  significant_genes$logFC < 0,
  ,
  drop = FALSE
]

cat("\n========== INPUT SUMMARY ==========\n")
cat("Genes in differential expression table:", nrow(all_genes), "\n")
cat("Significant genes:", nrow(significant_genes), "\n")
cat("Upregulated genes:", nrow(genes_up), "\n")
cat("Downregulated genes:", nrow(genes_down), "\n")


# 4. Map gene symbols to Entrez IDs ---------------------------

message("Mapping mouse gene symbols to Entrez IDs...")

# Map symbols using the mouse annotation database.
gene_mapping <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(all_genes$SYMBOL),
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENTREZID")
)

gene_mapping <- gene_mapping[
  !is.na(gene_mapping$ENTREZID),
  ,
  drop = FALSE
]

# Export the mapping table.
write.csv(
  gene_mapping,
  file = file.path(
    enrichment_dir,
    "03_gene_symbol_to_entrez_mapping.csv"
  ),
  row.names = FALSE
)

# Keep one Entrez ID per symbol for over-representation analysis.
# Ambiguous one-to-many mappings are excluded from these lists.
gene_mapping_unique <- gene_mapping[
  !duplicated(gene_mapping$SYMBOL) &
    !duplicated(gene_mapping$ENTREZID),
  ,
  drop = FALSE
]

# Helper function for mapping a list of gene symbols.
map_symbols <- function(data, mapping) {
  
  result <- merge(
    data,
    mapping,
    by = "SYMBOL",
    all = FALSE,
    sort = FALSE
  )
  
  result <- result[
    !is.na(result$ENTREZID),
    ,
    drop = FALSE
  ]
  
  result
}

all_genes_mapped <- map_symbols(
  all_genes,
  gene_mapping_unique
)

significant_mapped <- map_symbols(
  significant_genes,
  gene_mapping_unique
)

up_mapped <- map_symbols(
  genes_up,
  gene_mapping_unique
)

down_mapped <- map_symbols(
  genes_down,
  gene_mapping_unique
)

cat("\n========== GENE MAPPING ==========\n")
cat("Mapped genes in full table:", nrow(all_genes_mapped), "\n")
cat("Mapped significant genes:", nrow(significant_mapped), "\n")
cat("Mapped upregulated genes:", nrow(up_mapped), "\n")
cat("Mapped downregulated genes:", nrow(down_mapped), "\n")


# 5. Prepare gene lists for GO enrichment ---------------------

message("Preparing gene lists for GO enrichment...")

# Remove duplicate identifiers.
gene_symbols_total <- unique(significant_mapped$SYMBOL)
gene_symbols_up <- unique(up_mapped$SYMBOL)
gene_symbols_down <- unique(down_mapped$SYMBOL)

# Use all tested, mapped genes as the enrichment universe.
gene_universe <- unique(all_genes_mapped$SYMBOL)

write.csv(
  data.frame(SYMBOL = gene_universe),
  file = file.path(
    enrichment_dir,
    "03_GO_gene_universe.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(SYMBOL = gene_symbols_total),
  file = file.path(
    enrichment_dir,
    "03_significant_gene_symbols.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(SYMBOL = gene_symbols_up),
  file = file.path(
    enrichment_dir,
    "03_upregulated_gene_symbols.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(SYMBOL = gene_symbols_down),
  file = file.path(
    enrichment_dir,
    "03_downregulated_gene_symbols.csv"
  ),
  row.names = FALSE
)


# 6. GO over-representation analysis --------------------------

message("Running GO Cellular Component enrichment...")

# This helper safely handles empty gene lists or analyses
# that return no enriched terms.
run_enrichGO <- function(genes, universe, ontology = "CC") {
  
  if (length(genes) == 0) {
    message("Skipping GO ", ontology, ": empty gene list.")
    return(NULL)
  }
  
  tryCatch(
    enrichGO(
      gene = genes,
      universe = universe,
      OrgDb = org.Mm.eg.db,
      keyType = "SYMBOL",
      ont = ontology,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.20,
      readable = TRUE
    ),
    error = function(e) {
      message("GO enrichment failed: ", conditionMessage(e))
      return(NULL)
    }
  )
}

ego_total <- run_enrichGO(
  gene_symbols_total,
  gene_universe,
  ontology = "CC"
)

ego_up <- run_enrichGO(
  gene_symbols_up,
  gene_universe,
  ontology = "CC"
)

ego_down <- run_enrichGO(
  gene_symbols_down,
  gene_universe,
  ontology = "CC"
)


# 7. Export GO enrichment results -----------------------------

export_enrichment <- function(object, filename) {
  
  if (is.null(object)) {
    message("No result object to export: ", filename)
    return(invisible(NULL))
  }
  
  result_table <- as.data.frame(object)
  
  write.csv(
    result_table,
    file = file.path(enrichment_dir, filename),
    row.names = FALSE
  )
  
  invisible(result_table)
}

ego_total_df <- export_enrichment(
  ego_total,
  "03_GO_CC_total.csv"
)

ego_up_df <- export_enrichment(
  ego_up,
  "03_GO_CC_upregulated.csv"
)

ego_down_df <- export_enrichment(
  ego_down,
  "03_GO_CC_downregulated.csv"
)


# 8. GO enrichment plots --------------------------------------

message("Generating GO enrichment plots...")

if (!is.null(ego_total) && nrow(as.data.frame(ego_total)) > 0) {
  
  p_go_total <- barplot(
    ego_total,
    showCategory = 20,
    font.size = 8,
    title = "GO Cellular Component: significant genes",
    label_format = 60
  )
  
  ggsave(
    file.path(plot_dir, "03_GO_CC_total_barplot.png"),
    p_go_total,
    width = 11,
    height = 8,
    dpi = 300
  )
}

if (!is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0) {
  
  p_go_down <- dotplot(
    ego_down,
    x = "GeneRatio",
    font.size = 12,
    showCategory = 10
  )
  
  ggsave(
    file.path(plot_dir, "03_GO_CC_downregulated_dotplot.png"),
    p_go_down,
    width = 10,
    height = 7,
    dpi = 300
  )
}

if (!is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) {
  
  p_go_up <- dotplot(
    ego_up,
    x = "GeneRatio",
    font.size = 12,
    showCategory = 10
  )
  
  ggsave(
    file.path(plot_dir, "03_GO_CC_upregulated_dotplot.png"),
    p_go_up,
    width = 10,
    height = 7,
    dpi = 300
  )
}


# 9. Prepare ranked gene list for GSEA ------------------------

message("Preparing ranked gene list for GSEA...")

# GSEA uses all tested genes, not only significant genes.
# Collapse duplicated Entrez IDs by retaining the gene with
# the largest absolute log fold change.
ranked_genes <- all_genes_mapped[
  order(abs(all_genes_mapped$logFC), decreasing = TRUE),
  ,
  drop = FALSE
]

ranked_genes <- ranked_genes[
  !duplicated(ranked_genes$ENTREZID),
  ,
  drop = FALSE
]

gene_list <- ranked_genes$logFC
names(gene_list) <- ranked_genes$ENTREZID

gene_list <- sort(
  gene_list,
  decreasing = TRUE
)

# Remove invalid values and ensure unique Entrez IDs.
gene_list <- gene_list[
  is.finite(gene_list) &
    !is.na(names(gene_list)) &
    names(gene_list) != ""
]

stopifnot(!anyDuplicated(names(gene_list)))

write.csv(
  data.frame(
    ENTREZID = names(gene_list),
    logFC = as.numeric(gene_list)
  ),
  file = file.path(
    enrichment_dir,
    "03_GSEA_ranked_gene_list.csv"
  ),
  row.names = FALSE
)

cat("Genes in ranked GSEA list:", length(gene_list), "\n")


# 10. GO GSEA: ALL, BP, MF and CC -----------------------------

message("Running GO GSEA...")

run_gseGO <- function(ranked_list, ontology) {
  
  if (length(ranked_list) == 0) {
    return(NULL)
  }
  
  tryCatch(
    gseGO(
      geneList = ranked_list,
      OrgDb = org.Mm.eg.db,
      keyType = "ENTREZID",
      ont = ontology,
      minGSSize = 3,
      maxGSSize = 800,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = TRUE
    ),
    error = function(e) {
      message(
        "GSEA GO failed for ", ontology, ": ",
        conditionMessage(e)
      )
      return(NULL)
    }
  )
}

gse_GO_all <- run_gseGO(gene_list, "ALL")
gse_GO_bp <- run_gseGO(gene_list, "BP")
gse_GO_mf <- run_gseGO(gene_list, "MF")
gse_GO_cc <- run_gseGO(gene_list, "CC")


# 11. Export GO GSEA results ----------------------------------

gse_go_objects <- list(
  ALL = gse_GO_all,
  BP = gse_GO_bp,
  MF = gse_GO_mf,
  CC = gse_GO_cc
)

for (ontology_name in names(gse_go_objects)) {
  
  object <- gse_go_objects[[ontology_name]]
  
  if (!is.null(object)) {
    
    write.csv(
      as.data.frame(object),
      file = file.path(
        enrichment_dir,
        paste0("03_GSEA_GO_", ontology_name, ".csv")
      ),
      row.names = FALSE
    )
  }
}


# 12. GO GSEA plots -------------------------------------------

message("Generating GO GSEA plots...")

if (!is.null(gse_GO_mf) && nrow(as.data.frame(gse_GO_mf)) > 0) {
  
  p_gsea_mf <- dotplot(
    gse_GO_mf,
    showCategory = 5,
    split = ".sign"
  ) +
    facet_grid(. ~ .sign)
  
  ggsave(
    file.path(plot_dir, "03_GSEA_GO_MF_dotplot.png"),
    p_gsea_mf,
    width = 11,
    height = 7,
    dpi = 300
  )
}

if (!is.null(gse_GO_bp) && nrow(as.data.frame(gse_GO_bp)) > 0) {
  
  p_gsea_bp <- ridgeplot(
    gse_GO_bp,
    showCategory = 10,
    fill = "p.adjust"
  )
  
  ggsave(
    file.path(plot_dir, "03_GSEA_GO_BP_ridgeplot.png"),
    p_gsea_bp,
    width = 11,
    height = 8,
    dpi = 300
  )
}

if (!is.null(gse_GO_mf) && nrow(as.data.frame(gse_GO_mf)) > 0) {
  
  gse_GO_mf_sim <- pairwise_termsim(gse_GO_mf)
  
  saveRDS(
    gse_GO_mf_sim,
    file = file.path(
      enrichment_dir,
      "03_GSEA_GO_MF_similarity.rds"
    )
  )
  
  p_emap_mf <- emapplot(
    gse_GO_mf_sim,
    showCategory = 20
  )
  
  ggsave(
    file.path(plot_dir, "03_GSEA_GO_MF_emapplot.png"),
    p_emap_mf,
    width = 12,
    height = 10,
    dpi = 300
  )
}

if (!is.null(gse_GO_mf) && nrow(as.data.frame(gse_GO_mf)) > 0) {
  
  p_cnet_mf <- cnetplot(
    gse_GO_mf,
    showCategory = 10,
    foldChange = gene_list
  )
  
  ggsave(
    file.path(plot_dir, "03_GSEA_GO_MF_cnetplot.png"),
    p_cnet_mf,
    width = 12,
    height = 10,
    dpi = 300
  )
}

# Export a representative GSEA running enrichment plot.
if (!is.null(gse_GO_all) && nrow(as.data.frame(gse_GO_all)) > 0) {
  
  gsea_all_df <- as.data.frame(gse_GO_all)
  
  top_go_id <- gsea_all_df$ID[1]
  
  p_gsea_curve <- gseaplot2(
    gse_GO_all,
    geneSetID = top_go_id,
    pvalue_table = TRUE,
    title = gsea_all_df$Description[1]
  )
  
  ggsave(
    file.path(plot_dir, "03_GSEA_GO_top_term_curve.png"),
    p_gsea_curve,
    width = 10,
    height = 7,
    dpi = 300
  )
}


# 13. GSEA KEGG -----------------------------------------------

message("Running KEGG GSEA...")

gse_kegg <- tryCatch(
  gseKEGG(
    geneList = gene_list,
    organism = "mmu",
    keyType = "ncbi-geneid",
    minGSSize = 3,
    maxGSSize = 800,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = TRUE
  ),
  error = function(e) {
    message("KEGG GSEA failed: ", conditionMessage(e))
    return(NULL)
  }
)

if (!is.null(gse_kegg)) {
  
  write.csv(
    as.data.frame(gse_kegg),
    file = file.path(
      enrichment_dir,
      "03_GSEA_KEGG.csv"
    ),
    row.names = FALSE
  )
}


# 14. KEGG GSEA plots -----------------------------------------

message("Generating KEGG plots...")

if (!is.null(gse_kegg) && nrow(as.data.frame(gse_kegg)) > 0) {
  
  p_kegg_dot <- dotplot(
    gse_kegg,
    showCategory = 10,
    title = "KEGG GSEA",
    split = ".sign"
  ) +
    facet_grid(. ~ .sign)
  
  ggsave(
    file.path(plot_dir, "03_GSEA_KEGG_dotplot.png"),
    p_kegg_dot,
    width = 11,
    height = 7,
    dpi = 300
  )
  
  gse_kegg_sim <- pairwise_termsim(gse_kegg)
  
  saveRDS(
    gse_kegg_sim,
    file = file.path(
      enrichment_dir,
      "03_GSEA_KEGG_similarity.rds"
    )
  )
  
  p_kegg_emap <- emapplot(
    gse_kegg_sim,
    showCategory = 20
  )
  
  ggsave(
    file.path(plot_dir, "03_GSEA_KEGG_emapplot.png"),
    p_kegg_emap,
    width = 12,
    height = 10,
    dpi = 300
  )
  
  p_kegg_ridge <- ridgeplot(
    gse_kegg,
    showCategory = 10
  ) +
    labs(x = "Enrichment distribution")
  
  ggsave(
    file.path(plot_dir, "03_GSEA_KEGG_ridgeplot.png"),
    p_kegg_ridge,
    width = 11,
    height = 8,
    dpi = 300
  )
}


# 15. Export KEGG pathway visualization -----------------------

# This section is optional and runs only if KEGG has results.
# The pathway ID below is an example, not a predicted result.

if (!is.null(gse_kegg) && nrow(as.data.frame(gse_kegg)) > 0) {
  
  example_pathway <- "mmu00071"
  
  pathview_result <- tryCatch(
    pathview(
      gene.data = gene_list,
      pathway.id = sub("^mmu", "", example_pathway),
      species = "mmu",
      gene.idtype = "entrez"
    ),
    error = function(e) {
      message("Pathview failed: ", conditionMessage(e))
      return(NULL)
    }
  )
}


# 16. Save analysis objects -----------------------------------

saveRDS(
  list(
    ego_total = ego_total,
    ego_up = ego_up,
    ego_down = ego_down,
    gse_GO_all = gse_GO_all,
    gse_GO_bp = gse_GO_bp,
    gse_GO_mf = gse_GO_mf,
    gse_GO_cc = gse_GO_cc,
    gse_kegg = gse_kegg,
    gene_list = gene_list,
    gene_mapping = gene_mapping,
    gene_universe = gene_universe
  ),
  file = file.path(
    enrichment_dir,
    "03_enrichment_analysis_objects.rds"
  )
)


# 17. Record session information -----------------------------

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    enrichment_dir,
    "03_sessionInfo.txt"
  )
)


# 18. Final report -------------------------------------------

cat("\n")
cat("====================================================\n")
cat("SCRIPT 03: FUNCTIONAL ENRICHMENT AND GSEA\n")
cat("====================================================\n")

cat("Dataset: GSE295281\n")
cat("Organism: Mus musculus\n")
cat("GO ontology for over-representation: CC\n")
cat("Significant genes:", nrow(significant_genes), "\n")
cat("Upregulated genes:", nrow(genes_up), "\n")
cat("Downregulated genes:", nrow(genes_down), "\n")
cat("Genes in ranked GSEA list:", length(gene_list), "\n")
cat("Output directory:", enrichment_dir, "\n")
cat("Status: EXPLORATORY; differential expression scale unverified\n")
cat("====================================================\n")

