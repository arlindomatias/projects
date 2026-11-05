# ============================================================
# SCRIPT 01: DATA ACQUISITION AND ORGANIZATION
# Dataset: GSE178812
# ============================================================

### 1. Configure used files -----------------------------------
matrix_file <- file.path(
  data,
  "GSE178812_Doxo-counts.xlsx"
)

matrix_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178812/suppl/GSE178812%5FDoxo%2Dcounts%2Exlsx"

if (!file.exists(matrix_file)) {
  
  message("Downloading GEO supplementary expression matrix...")
  
  download.file(
    url = matrix_url,
    destfile = matrix_file,
    mode = "wb"
  )
  
} else {
  
  message("Expression matrix already exists. Skipping download.")
  
}

if (!file.exists(matrix_file)) {
  stop("Expression matrix file was not downloaded successfully.")
}

counts_raw <- read_excel( # Import count matrix
  file.path(data, "GSE178812_Doxo-counts.xlsx"))
colnames(counts_raw)[1] <- c("SYMBOL")

### 2. Retrieve GEO metadata -------------------------------------

local_file <- file.path(project_dir, paste0(geo_accession, "_series_matrix.txt.gz"))

# Use the local file if it exists; otherwise download from GEO
if (file.exists(local_file)) {
  message("Local file found: ", local_file)
  gse <- getGEO(filename = local_file, GSEMatrix = TRUE)
} else {
  message("Local file not found. Downloading from GEO: ", geo_accession)
  gse <- getGEO(geo_accession, GSEMatrix = TRUE, destdir = project_dir)
  if (is.list(gse)) gse <- gse[[1]]   # online retrieval returns a list
}

# Extract sample metadata according to the object class
if (inherits(gse, "SummarizedExperiment")) {
  colnames(rowData(gse)) <- make.names(colnames(rowData(gse)))
  info <- as.data.frame(colData(gse))
} else {
  fvarLabels(gse) <- make.names(fvarLabels(gse))
  info <- pData(gse)
}

# Remove non-used group
info <- info[!grepl("DMQ", info$title), ]

colData <- data.frame(
  GEO_ID = rownames(info),
  condition = sub(":.*", "", info$title))


### 3. Define experimental groups -------------------------------
group <- ifelse(
  grepl("Doxo-treated", info$title),
  "CHEM",
  ifelse(
    grepl("Sham", info$title),
    "VEH",
    NA_character_
  )) 

colData$group <- factor(group) 
colData$group <- factor(colData$group, levels = c("VEH", "CHEM"))
counts_raw <- counts_raw %>% dplyr::select(SYMBOL, all_of(colData$condition))

### 4. Inspect gene identifiers ---------------------------------
# Missing identifiers.
missing_symbols <- sum(is.na(counts_raw$SYMBOL) | counts_raw$SYMBOL == "")

# Duplicated gene symbols.
duplicated_symbol_values <- length(unique(
  counts_raw$SYMBOL[
    duplicated(counts_raw$SYMBOL)
  ]
))

cat("Missing symbols:", missing_symbols, "|" ,"Duplicated symbols:", duplicated_symbol_values)

counts <- counts_raw[
  !duplicated(counts_raw$SYMBOL),
]

counts <- as.data.frame(counts)

rownames(counts) <- counts$SYMBOL

counts$SYMBOL <- NULL

colnames(counts) <- colData$GEO_ID

all(counts == round(counts))  # Should return True
any(counts < 0)               # Should return False
anyNA(counts)                 # Should return False

### 5. Create keytable to convert gene IDs ----------------------
keytable <- bitr(rownames(counts),
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Mm.eg.db)
keytable <- na.omit(keytable)

# 6. Export organized data ------------------------------------
write.csv(counts, "data/counts.csv", row.names = F)
write.csv(colData, "data/colData.csv", row.names = F)
write.csv(keytable, "data/keytable.csv", row.names = F)

gc()
