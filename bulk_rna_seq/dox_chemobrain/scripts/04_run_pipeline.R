# Main RNA-seq pipeline

scripts <- c(
  "scripts/01_download_and_prepare.R",
  "scripts/02_differential_expression_and_qc.R",
  "scripts/03_enrichment.R"
)

for (script in scripts) {
  message("\nRunning: ", script)
  
  source(script)
  
  message("Completed: ", script)
}

message("\nRNA-seq pipeline completed successfully.")