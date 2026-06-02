args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 09_summary_report.R <results_dir> <output_md>", call. = FALSE)
}

results_dir <- args[[1]]
output_md <- args[[2]]

if (!dir.exists(results_dir)) {
  stop("Results directory not found: ", results_dir)
}

dataset_dirs <- list.dirs(results_dir, full.names = TRUE, recursive = FALSE)
if (length(dataset_dirs) == 0) {
  stop("No dataset directories found in ", results_dir)
}

lines <- c(
  "# MuSiC deconvolution summary",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
  ""
)

for (dataset_dir in dataset_dirs) {
  dataset <- basename(dataset_dir)
  qc_dir <- file.path(dataset_dir, "qc")
  cell_summary_path <- file.path(qc_dir, paste0(dataset, "_celltype_summary.tsv"))
  sample_summary_path <- file.path(qc_dir, paste0(dataset, "_sample_summary.tsv"))

  if (!file.exists(cell_summary_path) || !file.exists(sample_summary_path)) {
    next
  }

  cell_summary <- read.delim(cell_summary_path, check.names = FALSE)
  sample_summary <- read.delim(sample_summary_path, check.names = FALSE)

  n_samples <- nrow(sample_summary)
  n_celltypes <- nrow(cell_summary)
  mean_total <- mean(sample_summary$total_fraction, na.rm = TRUE)
  top_celltypes <- head(cell_summary$celltype, 5)

  lines <- c(
    lines,
    paste0("## ", dataset),
    "",
    paste0("- Samples: ", n_samples),
    paste0("- Cell types: ", n_celltypes),
    paste0("- Mean total fraction: ", sprintf("%.3f", mean_total)),
    paste0("- Top cell types: ", paste(top_celltypes, collapse = ", ")),
    ""
  )
}

if (length(lines) <= 4) {
  stop("No QC summaries found in ", results_dir)
}

output_dir <- dirname(output_md)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

writeLines(lines, output_md)
cat("Summary report written to", output_md, "\n")
