# Audit H80: read converged flag, CS count and LD provenance from the production eQTL rds of the
# 111 positive/tested genes whose last summary row says not converged, plus ENSG00000158863.
suppressPackageStartupMessages(library(data.table))
`%||%` <- function(a, b) if (is.null(a)) b else a
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/eqtl_susie_polyfun"
D <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/audit/H80"
g <- fread(file.path(D, "notconverged_summary_genes.tsv"))
g <- rbind(g, data.table(ensembl = "ENSG00000158863", susie_state_t12 = "pilot_gene"))
paths <- Sys.glob(file.path(FM, "chr*", paste0(g$ensembl, "_susie.rds")))
names(paths) <- sub("_susie.rds", "", basename(paths))
out <- rbindlist(lapply(g$ensembl, function(id) {
  p <- paths[id]; if (is.na(p)) return(data.table(ensembl = id, rds = FALSE))
  s <- readRDS(p)
  data.table(ensembl = id, rds = TRUE, converged = isTRUE(s$converged), n_cs = length(s$sets$cs),
             niter = s$niter, ld_panel = s$ld_provenance$ld_panel %||% NA_character_,
             ld_base_dir = s$ld_provenance$ld_base_dir %||% NA_character_,
             generated = s$ld_provenance$generated_at_utc %||% NA_character_,
             mtime = format(file.mtime(p), "%Y-%m-%d"))
}), fill = TRUE)
`%||%` <- function(a, b) if (is.null(a)) b else a
out <- merge(g, out, by = "ensembl")
fwrite(out, file.path(D, "rds_check.tsv"), sep = "\t")
print(out[, .N, by = .(susie_state_t12, rds, converged)])
print(out[, .N, by = .(ld_panel, substr(generated, 1, 10))])
print(out[ensembl == "ENSG00000158863"])
