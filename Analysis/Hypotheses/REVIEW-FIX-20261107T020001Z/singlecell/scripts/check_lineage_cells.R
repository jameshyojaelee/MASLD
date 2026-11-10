suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
root <- Sys.getenv("MASLD_PROJECT_ROOT")
source(file.path(root, "Analysis/SingleCell/scripts/lib_source_diagnosis.R"))
libraries <- data.table(sample = c("run1", "run2", "excluded_sort", "absent_run"),
  donor = c("donorA", "donorA", "donorA", "donorB"),
  n_Hepatocytes = c(30L, 25L, 500L, 100L),
  `n_Circulating_NK/NKT` = c(20L, 35L, 500L, 100L))
mapping <- setNames(libraries$donor, libraries$sample)
observed <- source_lineage_cell_counts(libraries, "Hepatocytes", c("run1", "run2"), mapping)
stopifnot(nrow(observed) == 1L, observed$sample == "donorA", observed$n_cells == 55L,
          observed$n_libraries == 2L)
stopifnot(source_lineage_cell_counts(libraries, "Circulating_NK_NKT", c("run1", "run2"), mapping)$n_cells == 55L)
stopifnot(source_lineage_cell_counts(libraries, "Hepatocytes", "run1", mapping)$n_cells < 50L)
mapping["run1"] <- "wrong_donor"
stopifnot(inherits(try(source_lineage_cell_counts(libraries, "Hepatocytes", "run1", mapping), silent = TRUE), "try-error"))
message("PASS: lineage counts include only present eligible columns, sum split runs, preserve lineage names and reject mismatched donors")
for (script in c("pseudobulk_de.R", "run_pseudobulk_de_coarse_stages.R", "lib_source_diagnosis.R")) {
  parse(file.path(root, "Analysis/SingleCell/scripts", script))
}
message("PASS: three modified R sources parse")
