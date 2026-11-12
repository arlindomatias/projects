# ============================================================
# SCRIPT 01: tximport + metadados -> DESeqDataSet
# Substitui o 01_download_and_prepare.R do bulk-rna-seq (que lia a planilha do GEO).
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))
suppressPackageStartupMessages(library(DESeq2))

### 1. Importar ------------------------------------------------------------
txi <- readRDS(snakemake@input$txi)
stopifnot(all(c("counts", "abundance", "length") %in% names(txi)))
if (!identical(txi$countsFromAbundance, "no")) {
  stop("txi must be imported with countsFromAbundance = 'no' for DESeqDataSetFromTximport ",
       "(found: ", txi$countsFromAbundance, ").")
}

samples <- as.data.frame(read_tsv(snakemake@input$samples, show_col_types = FALSE))
tx2gene <- read_tsv(snakemake@input$tx2gene, show_col_types = FALSE)

id_col    <- cfg$samples$id_column
cond_col  <- cfg$samples$condition_column
label_col <- cfg$samples$label_column %||% id_col
labels    <- unlist(cfg$groups$labels)

### 2. Metadados das amostras ---------------------------------------------
missing_in_txi <- setdiff(samples[[id_col]], colnames(txi$counts))
missing_in_tsv <- setdiff(colnames(txi$counts), samples[[id_col]])
if (length(missing_in_txi) || length(missing_in_tsv)) {
  stop("Samples differ between samples table and txi. Only in table: ",
       paste(missing_in_txi, collapse = ", "), " | only in txi: ",
       paste(missing_in_tsv, collapse = ", "))
}

unknown <- setdiff(unique(samples[[cond_col]]), names(labels))
if (length(unknown)) stop("Conditions without a label in config groups$labels: ",
                          paste(unknown, collapse = ", "))
bad_labels <- labels[make.names(labels) != labels]
if (length(bad_labels)) stop("Group labels must be syntactic R names: ",
                             paste(bad_labels, collapse = ", "))
stopifnot(ref_level %in% labels, trt_level %in% labels)

exclude <- unlist(cfg$exclude_samples)
if (length(exclude)) {
  stopifnot(all(exclude %in% samples[[id_col]]))
  message("Excluding samples: ", paste(exclude, collapse = ", "))
}

coldata <- samples |>
  filter(!.data[[id_col]] %in% exclude) |>
  dplyr::select(!matches("fastq|md5|ftp|url", ignore.case = TRUE)) |>
  mutate(
    sample_id = .data[[id_col]],
    label     = as.character(.data[[label_col]]),
    group     = factor(labels[.data[[cond_col]]], levels = c(ref_level, trt_level))
  ) |>
  filter(!is.na(group)) |>          # condições fora do contraste
  arrange(group, sample_id) |>
  as.data.frame()
rownames(coldata) <- coldata$sample_id

covariates <- unlist(cfg$covariates)
stopifnot(all(covariates %in% names(coldata)))
for (cv in covariates) if (is.character(coldata[[cv]])) coldata[[cv]] <- factor(coldata[[cv]])

message("Samples per group:")
print(table(coldata$group))

### 3. Recortar o txi nas amostras usadas, na ordem do colData -------------
for (m in c("counts", "abundance", "length")) {
  txi[[m]] <- txi[[m]][, coldata$sample_id, drop = FALSE]
}

### 4. Anotação dos genes --------------------------------------------------
genes <- unique(tx2gene[, c("gene_id", "gene_name")])
stopifnot(!anyDuplicated(genes$gene_id))

annot <- data.frame(gene_id_version = rownames(txi$counts))
annot$gene_id <- sub("\\.[0-9]+(_PAR_Y)?$", "\\1", annot$gene_id_version)
annot$SYMBOL  <- genes$gene_name[match(annot$gene_id_version, genes$gene_id)]
stopifnot(!anyDuplicated(annot$gene_id))
message("Genes: ", nrow(annot), " | without symbol: ", sum(is.na(annot$SYMBOL)),
        " | symbols shared by >1 gene: ",
        sum(duplicated(annot$SYMBOL) | duplicated(annot$SYMBOL, fromLast = TRUE)))

### 5. DESeqDataSet --------------------------------------------------------
design <- as.formula(paste("~", paste(c(covariates, "group"), collapse = " + ")))
message("Design: ", deparse(design))

dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = design)
rownames(dds) <- annot$gene_id
mcols(dds) <- cbind(mcols(dds), annot)

saveRDS(dds, snakemake@output$dds)
write_csv(coldata, snakemake@output$coldata)
message("DESeqDataSet: ", nrow(dds), " genes x ", ncol(dds), " samples")
