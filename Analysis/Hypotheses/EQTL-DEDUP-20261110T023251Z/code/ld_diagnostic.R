#!/usr/bin/env Rscript
# Minimum eigenvalue of eQTL LD (PolyFun EUR) for affected and control genes.
suppressPackageStartupMessages(library(data.table))
run <- Sys.getenv("RUN_DIR")
fm <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
setwd(fm)
source(file.path(fm, "src/finemapping_functions.R"))
orig_match <- get_ld_per_locus
source(file.path(run, "code/unique_ld_match.R"))
Sys.unsetenv("IDENTITY_LOSS_DIR")
aff <- fread(file.path(run, "inputs/affected_primary_genes.tsv"))[in_refit == TRUE, .(ensembl, chr, set = "affected_unique_match")]
ctl <- fread(file.path(run, "inputs/control_unaffected_genes.tsv"))[, .(ensembl, chr, set = "control_production_match")]
todo <- rbind(aff, ctl)
res <- list()
for (ch in sort(unique(todo$chr))) {
  e <- fread(sprintf("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/broadaway_eqtl/chr%d_marginal_summary_results.tsv", ch))
  e <- e[!is.na(POS) & !is.na(Beta) & !is.na(SE) & SE > 0 & ENSG != ""]
  for (i in which(todo$chr == ch)) {
    g <- todo[i]; gd <- e[ENSG == g$ensembl]
    START <- max(min(gd$POS) - 500000L, 1L); END <- max(gd$POS) + 500000L
    ss <- data.frame(chromosome = gd$CHR, position = gd$POS, allele1 = gd$NEA, allele2 = gd$EA,
                     beta = gd$Beta, standard_error = gd$SE)
    f <- if (g$set == "affected_unique_match") get_ld_per_locus_unique else orig_match
    r <- tryCatch(f(ss, LOCUS = g$ensembl, CHR = ch, START = START, END = END, ancestry = "EUR"),
                  error = function(err) NULL)
    if (is.null(r)) { res[[length(res) + 1]] <- data.table(g, n = NA_integer_); next }
    R <- as.matrix(r[[2]])
    ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
    res[[length(res) + 1]] <- data.table(g, n = nrow(R), min_eigen = min(ev),
      n_eigen_below_minus_0.001 = sum(ev < -0.001),
      chol_after_0.001 = tryCatch({chol(R + diag(1e-3, nrow(R))); TRUE}, error = function(e) FALSE),
      ld_value_decimals = max(nchar(sub("^-?[0-9]*\\.?", "", format(R[1, 1:min(50, ncol(R))], scientific = FALSE, digits = 15)))))
    cat(g$set, g$ensembl, nrow(R), min(ev), "\n")
    rm(R, r); gc(FALSE)
  }
}
out <- rbindlist(res, fill = TRUE)
dir.create(file.path(run, "results"), showWarnings = FALSE)
fwrite(out, file.path(run, "results/eqtl_ld_matrix_check.tsv"), sep = "\t")
print(out[, .(genes = .N, chol_pass = sum(chol_after_0.001, na.rm = TRUE), median_min_eigen = median(min_eigen, na.rm = TRUE),
              max_min_eigen = max(min_eigen, na.rm = TRUE)), by = set])
