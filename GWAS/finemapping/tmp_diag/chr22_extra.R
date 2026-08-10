suppressPackageStartupMessages(library(data.table))
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
nw <- basename(list.files(file.path(FM, "results/eqtl_susie_polyfun/chr22"), pattern = "[.]rds$"))
od <- basename(list.files(file.path(FM, "results/eqtl_susie/chr22"),        pattern = "[.]rds$"))
extra <- setdiff(nw, od)
cat(sprintf("chr22 NEW-only: %d | OLD-only: %d\n", length(extra), length(setdiff(od, nw))))
grab <- function(dir, f) {
  s <- tryCatch(readRDS(file.path(FM, dir, f)), error = function(e) NULL)
  if (is.null(s)) return(NULL)
  data.table(gene = sub("_susie.rds", "", f),
             converged = isTRUE(s$converged),
             n_cs = if (is.null(s$sets$cs)) 0L else length(s$sets$cs))
}
e <- rbindlist(Filter(Negate(is.null), lapply(head(extra, 25), function(f) grab("results/eqtl_susie_polyfun/chr22", f))))
cat("\n=== the 57 NEW-only fits (first 25) ===\n")
cat(sprintf("  converged: %d/%d | with >=1 credible set: %d\n",
            sum(e$converged), nrow(e), sum(e$n_cs > 0)))
# baseline: shared genes, same directory
sh <- head(intersect(nw, od), 25)
s2 <- rbindlist(Filter(Negate(is.null), lapply(sh, function(f) grab("results/eqtl_susie_polyfun/chr22", f))))
cat(sprintf("  BASELINE shared genes: converged %d/%d | with >=1 CS: %d\n",
            sum(s2$converged), nrow(s2), sum(s2$n_cs > 0)))
cat("\nOnly fits with a credible set can enter coloc.susie(), so that is the\n")
cat("number that matters for the Fig2 SuSiE arm.\n")
