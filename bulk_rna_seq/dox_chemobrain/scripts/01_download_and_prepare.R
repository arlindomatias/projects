# ============================================================
# SCRIPT 01: DATA ACQUISITION AND ORGANIZATION
# Dataset: GSE178812

# 1. Project configuration -----------------------------------
geo_accession <- "GSE178812"
options(timeout = 600)

matrix_file <- file.path(
  project_dir,
  "GSE178812_Doxo-counts.xlsx"
)

matrix_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178812/suppl/GSE178812%5FDoxo%2Dcounts%2Exlsx"

# 2. Download expression matrix -------------------------------

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


# 3. Import expression matrix ----------------------------------
counts_raw <- read_excel("GSE178812_Doxo-counts.xlsx")
colnames(counts_raw)[1] <- c("SYMBOL")

# 4. Retrieve GEO metadata -------------------------------------
gse <- getGEO(geo_accession, GSEMatrix = TRUE, returnType = "ExpressionSet")[[1]]
fvarLabels(gse) <- make.names(fvarLabels(gse))
info <- pData(gse)

info <- info[!grepl("DMQ", info$title), ]


colData <- data.frame(
  GEO_ID = rownames(info),
  condition = sub(":.*", "", info$title))


# 5. Define experimental groups -------------------------------
group <- ifelse(
  grepl("Doxo-treated", info$title),
  "CHEM",
  ifelse(
    grepl("Sham", info$title),
    "VEH",
    NA_character_
  )) 

colData$group <- factor(group) 
counts_raw <- counts_raw %>% dplyr::select(SYMBOL, all_of(colData$condition))

# 6. Inspect gene identifiers ---------------------------------
# Check missing identifiers.
missing_symbols <- sum(is.na(counts_raw$SYMBOL) | counts_raw$SYMBOL == "")

# Identify duplicated gene symbols.
duplicated_symbol_values <- length(unique(
  counts_raw$SYMBOL[
    duplicated(counts_raw$SYMBOL)
  ]
))

cat("Missing symbols:", missing_symbols, "||" ,"Duplicated symbols:", duplicated_symbol_values)

counts <- counts_raw[
  !duplicated(counts_raw$SYMBOL),
]

counts <- as.data.frame(counts)

rownames(counts) <- counts$SYMBOL

counts$SYMBOL <- NULL

colnames(counts) <- colData$GEO_ID


keytable <- bitr(rownames(counts),
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Mm.eg.db)
keytable <- na.omit(keytable)

# 7. Export organized data ------------------------------------

write.csv(counts, "counts.csv")

write.csv(colData, "colData.csv")

write.csv(keytable, "keytable.csv")
