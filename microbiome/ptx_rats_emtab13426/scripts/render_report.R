# Render report.Rmd. The Rmd reads inputs, params and outputs from the `snakemake` object,
# which stays visible because rendering happens in a child of the global environment.
log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

html <- normalizePath(snakemake@output[["html"]], mustWork = FALSE)
rmarkdown::render(
  snakemake@input[["rmd"]],
  output_file = basename(html),
  output_dir = dirname(html),
  intermediates_dir = tempdir(),
  knit_root_dir = getwd(),
  envir = new.env(parent = globalenv())
)
