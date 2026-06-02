#!/usr/bin/env Rscript
# ============================================================================
# 55e_combine_regulatory_atlas.R
# Combiner for the Tier-2 regulatory variant-to-gene CROSS-VALIDATION layers.
# Left-joins the three per-layer gene-level tables on human_symbol into ONE
# RNA-seq/results/multi_evidence/regulatory_atlas_columns.tsv.
#
# Does NOT edit 27a or any shared pipeline file. 27a's post-assembly merge block
# can pick this file up by the same drop-then-merge-on-human_symbol convention
# used for sqtl_atlas_columns.tsv etc.
#
# Inputs (each produced by its own script; missing = treated as empty):
#   abc_atlas_columns.tsv   (55b)
#   caqtl_atlas_columns.tsv (55c)
#   ccre_atlas_columns.tsv  (55d)
# ============================================================================

suppressPackageStartupMessages({ library(data.table) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ME  <- file.path(BASE, "RNA-seq/results/multi_evidence")
OUT <- file.path(ME, "regulatory_atlas_columns.tsv")

read_layer <- function(f) {
  p <- file.path(ME, f)
  if (!file.exists(p)) { cat("  MISSING (skip): ", f, "\n", sep = ""); return(NULL) }
  dt <- fread(p)
  if (!"human_symbol" %in% names(dt) || nrow(dt) == 0) {
    cat("  EMPTY (skip): ", f, "\n", sep = ""); return(NULL) }
  dt <- dt[!is.na(human_symbol) & human_symbol != ""]
  dt <- dt[!duplicated(human_symbol)]
  cat(sprintf("  %-26s %5d genes, %2d cols\n", f, nrow(dt), ncol(dt)))
  dt
}

cat("== 55e_combine_regulatory_atlas ==\n")
abc   <- read_layer("abc_atlas_columns.tsv")
caqtl <- read_layer("caqtl_atlas_columns.tsv")
ccre  <- read_layer("ccre_atlas_columns.tsv")

layers <- Filter(Negate(is.null), list(abc, caqtl, ccre))
if (length(layers) == 0) stop("No regulatory layer outputs found -- run 55b/55c/55d first.")

# Full outer join on human_symbol (Reduce + merge all=TRUE).
combined <- Reduce(function(x, y) merge(x, y, by = "human_symbol", all = TRUE), layers)

# Derived: how many of the 3 layers cross-validate this gene.
# Each hit flag is coerced to a non-NA logical; layers absent from the join
# (their script didn't run / found nothing) are treated as all-FALSE.
combined[, abc_v2g_hit       := if ("abc_v2g_hit"       %in% names(combined)) (!is.na(abc_v2g_hit)       & abc_v2g_hit       == TRUE) else FALSE]
combined[, caqtl_credset_hit := if ("caqtl_credset_hit" %in% names(combined)) (!is.na(caqtl_credset_hit) & caqtl_credset_hit == TRUE) else FALSE]
combined[, ccre_credset_hit  := if ("ccre_credset_hit"  %in% names(combined)) (!is.na(ccre_credset_hit)  & ccre_credset_hit  == TRUE) else FALSE]
combined[, n_regulatory_layers := as.integer(abc_v2g_hit) + as.integer(caqtl_credset_hit) + as.integer(ccre_credset_hit)]

# Any MASLD-disease-driven support across layers (vs liver-enzyme-proxy only).
md_cols <- intersect(c("abc_masld_gwas_driven", "caqtl_masld_gwas_driven", "ccre_masld_gwas_driven"),
                     names(combined))
if (length(md_cols) > 0) {
  combined[, regulatory_masld_gwas_driven :=
    Reduce(`|`, lapply(md_cols, function(c) !is.na(get(c)) & get(c) == TRUE))]
}

setcolorder(combined, c("human_symbol", "n_regulatory_layers",
                        intersect(c("abc_v2g_hit","caqtl_credset_hit","ccre_credset_hit"), names(combined))))
combined <- combined[order(-n_regulatory_layers, human_symbol)]

fwrite(combined, OUT, sep = "\t")
cat(sprintf("\nWROTE %s : %d genes x %d cols\n", OUT, nrow(combined), ncol(combined)))
cat("  genes per layer:\n")
cat(sprintf("    ABC v2g           : %d\n", sum(combined$abc_v2g_hit)))
cat(sprintf("    caQTL credset     : %d\n", sum(combined$caqtl_credset_hit)))
cat(sprintf("    cCRE credset      : %d\n", sum(combined$ccre_credset_hit)))
cat("  n_regulatory_layers distribution:\n")
print(combined[, .N, by = n_regulatory_layers][order(-n_regulatory_layers)])
cat(sprintf("  ACTG1: %s\n", ifelse("ACTG1" %in% combined$human_symbol, "PRESENT", "absent")))
if ("ACTG1" %in% combined$human_symbol) {
  show <- intersect(c("human_symbol","n_regulatory_layers","abc_v2g_hit","abc_max_score",
                      "caqtl_credset_hit","ccre_credset_hit","ccre_class"), names(combined))
  print(combined[human_symbol == "ACTG1", ..show])
}
cat("\nGenes cross-validated by >=2 layers (top 20):\n")
print(head(combined[n_regulatory_layers >= 2], 20))
