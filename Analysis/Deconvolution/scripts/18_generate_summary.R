#!/usr/bin/env Rscript
# Generate Final Summary Report for Deconvolution Pipeline
# Collates QC metrics and Deconvolution Composition

suppressPackageStartupMessages({
  library(knitr)
  library(rmarkdown)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
RESULTS_BASE <- ifelse(length(args) > 0, args[1], "RNA-seq/deconvolution/results")
OUTPUT_FILE <- ifelse(length(args) > 1, args[2], file.path(RESULTS_BASE, "summary_report.md"))

# Generate Markdown Report
report_content <- c(
  "# Deconvolution Pipeline Summary Report",
  paste("Date:", Sys.Date()),
  "",
  "## 1. Overview",
  "This report summarizes the deconvolution results for Cas13 bulk RNA-seq samples using:",
  "- **MuSiC**: Multi-subject Single Cell deconvolution",
  "- **InstaPrism**: Fast implementation of BayesPrism",
  "",
  "## 2. Dataset Summary",
  "| Dataset | Species | Samples | Status |",
  "|:--------|:--------|:--------|:-------|"
)

# Auto-discover datasets from results directory; fall back to known list
discovered <- list.dirs(RESULTS_BASE, full.names = FALSE, recursive = FALSE)
discovered <- discovered[!discovered %in% c("summary_plots", "cibersortx_sigmatrix", "")]
known_order <- c("inhouse_MCD", "GSE156918", "GSE205974",
                 "GSE130970", "GSE135251", "GSE126848", "PRJNA512027", "GSE167523")
datasets <- c(known_order[known_order %in% discovered],
              setdiff(discovered, known_order))

mouse_datasets <- c("inhouse_MCD", "GSE156918", "GSE205974")

for (ds in datasets) {
  music_file  <- file.path(RESULTS_BASE, ds, paste0(ds, "_music_prop_weighted.tsv"))
  bp_file     <- file.path(RESULTS_BASE, ds, paste0(ds, "_bayesprism_proportions.tsv"))
  music_done  <- file.exists(music_file)
  bp_done     <- file.exists(bp_file)
  if (music_done) {
    props <- read.delim(music_file, row.names = 1)
    n_samples <- nrow(props)
  } else {
    n_samples <- 0
  }
  status <- dplyr::case_when(
    music_done & bp_done  ~ "✅ MuSiC + BayesPrism",
    music_done & !bp_done ~ "⚠️ MuSiC only",
    !music_done & bp_done ~ "⚠️ BayesPrism only",
    TRUE                  ~ "❌ Missing"
  )
  species <- ifelse(ds %in% mouse_datasets, "Mouse", "Human")
  report_content <- c(report_content, paste("|", ds, "|", species, "|", n_samples, "|", status, "|"))
}

report_content <- c(report_content, 
  "",
  "## 3. Method Comparison",
  "Detailed comparison plots are available in `Analysis/Deconvolution/comparison_plots/`.",
  "",
  "### Correlation Analysis",
  "Correlation between MuSiC and InstaPrism estimated proportions for each dataset.",
  "![Correlation Plots](../comparison_plots/method_correlation.pdf)",
  "",
  "### Composition Analysis",
  "Stacked bar plots showing cell type composition across all samples.",
  "![Composition Comparison](../comparison_plots/method_comparison_barplots.pdf)",
  "",
  "## 4. Key Findings",
  "- **Consistency**: Check correlation plots. High correlation (>0.8) indicates robust deconvolution.",
  "- **Discrepancies**: Significant differences may arise for rare cell types or those with similar expression profiles.",
  "",
  "## 5. Next Steps",
  "- Use CIBERSORTx web results to further validate these findings.",
  "- Proceed with differential expression or other downstream analyses using the most robust method (or consensus)."
)

writeLines(report_content, OUTPUT_FILE)
cat("Summary report generated at:", OUTPUT_FILE, "\n")
