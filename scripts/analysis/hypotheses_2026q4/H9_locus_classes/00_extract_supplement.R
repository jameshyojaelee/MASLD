#!/usr/bin/env Rscript
# H9 step 0: convert every sheet of the downloaded Jamialahmadi 2024 supplement to plain TSV.
# The xlsx is untrusted data: it is only read as text cells (col_types = "text"), never executed.
suppressPackageStartupMessages({library(readxl); library(data.table)})
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 2L) stop("Usage: 00_extract_supplement.R xlsx out_dir")
f <- argv[[1]]; out <- argv[[2]]
if (dir.exists(out)) stop("Refusing to overwrite ", out)
dir.create(out, recursive = TRUE)
set.seed(20260923)
for (s in excel_sheets(f)) {
  x <- as.data.table(suppressMessages(read_excel(f, sheet = s, col_names = FALSE, col_types = "text")))
  fwrite(x, file.path(out, paste0(s, ".tsv")), sep = "\t", col.names = FALSE, na = "NA", quote = FALSE)
  cat(s, nrow(x), "x", ncol(x), "\n")
}
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
