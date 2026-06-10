#!/usr/bin/env Rscript
# lvqw_divergence_report.R
# ---------------------------------------------------------------------------
# DIVERGENCE REPORT: LVQW-fixed staging sex calls vs canonical dream sex calls.
# READ-ONLY. Compares:
#   STAGING  : results/integration/sex_v3/lvqw_staging/{sex_deg_classification_v3,
#              sex_interaction_dream_v3}.csv
#   CANONICAL: results/integration/sex_v3/{sex_deg_classification_v3,
#              sex_interaction_dream_v3}.csv
# Emits a console report + writes lvqw_staging/divergence_report.csv summaries.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table) })

BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RDIR  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                   "results/integration")
CANON <- file.path(RDIR, "sex_v3")
STAGE <- file.path(CANON, "lvqw_staging")

f_canon_cls <- file.path(CANON, "sex_deg_classification_v3.csv")
f_stage_cls <- file.path(STAGE, "sex_deg_classification_v3.csv")
f_canon_int <- file.path(CANON, "sex_interaction_dream_v3.csv")
f_stage_int <- file.path(STAGE, "sex_interaction_dream_v3.csv")
for (f in c(f_canon_cls, f_stage_cls, f_canon_int, f_stage_int))
  if (!file.exists(f)) stop("missing input: ", f)

cc <- fread(f_canon_cls); sc <- fread(f_stage_cls)
ci <- fread(f_canon_int); si <- fread(f_stage_int)

cat("=====================================================================\n")
cat("LVQW (dataset fixed) vs DREAM (random slope) — sex divergence report\n")
cat("=====================================================================\n")
cat("canonical classification rows:", nrow(cc), "  staging rows:", nrow(sc), "\n\n")

## (a) per-class counts ------------------------------------------------------
cat("---- (a) sex_class distribution ----\n")
lvls <- c("Concordant","Female_biased","Male_biased","Divergent")
tab_c <- table(factor(cc$sex_class, levels = lvls))
tab_s <- table(factor(sc$sex_class, levels = lvls))
cmp <- data.table(sex_class = lvls,
                  dream_canonical = as.integer(tab_c),
                  lvqw_staging    = as.integer(tab_s))
cmp[, delta := lvqw_staging - dream_canonical]
print(cmp)
cat("\n  dimorphic (non-Concordant) total: dream =",
    sum(cmp[sex_class != "Concordant", dream_canonical]),
    "  lvqw =", sum(cmp[sex_class != "Concordant", lvqw_staging]), "\n\n")

## shared-gene frame ---------------------------------------------------------
setkey(cc, gene); setkey(sc, gene)
shared <- intersect(cc$gene, sc$gene)
cat("shared genes:", length(shared), "\n")
cdt <- cc[J(shared)]; sdt <- sc[J(shared)]
m <- merge(cdt[, .(gene, class_c = sex_class, logFC_M_c = logFC_M, logFC_F_c = logFC_F)],
           sdt[, .(gene, class_s = sex_class, logFC_M_s = logFC_M, logFC_F_s = logFC_F)],
           by = "gene")

## (b) Jaccard of dimorphic set + direction concordance ----------------------
cat("\n---- (b) dimorphic-set Jaccard + direction concordance ----\n")
dim_c <- m[class_c != "Concordant", gene]
dim_s <- m[class_s != "Concordant", gene]
inter <- length(intersect(dim_c, dim_s)); uni <- length(union(dim_c, dim_s))
jacc <- if (uni > 0) inter / uni else NA_real_
cat(sprintf("  dimorphic dream=%d  lvqw=%d  intersect=%d  union=%d  Jaccard=%.4f\n",
            length(dim_c), length(dim_s), inter, uni, jacc))
dir_M <- m[is.finite(logFC_M_c) & is.finite(logFC_M_s),
           mean(sign(logFC_M_c) == sign(logFC_M_s))]
dir_F <- m[is.finite(logFC_F_c) & is.finite(logFC_F_s),
           mean(sign(logFC_F_c) == sign(logFC_F_s))]
cat(sprintf("  direction concordance (all shared): logFC_M=%.4f  logFC_F=%.4f\n", dir_M, dir_F))
sub <- m[gene %in% union(dim_c, dim_s)]
cat(sprintf("  direction concordance (dimorphic-union): logFC_M=%.4f  logFC_F=%.4f\n",
            sub[is.finite(logFC_M_c) & is.finite(logFC_M_s), mean(sign(logFC_M_c)==sign(logFC_M_s))],
            sub[is.finite(logFC_F_c) & is.finite(logFC_F_s), mean(sign(logFC_F_c)==sign(logFC_F_s))]))

## (c) new vs lost dimorphic -------------------------------------------------
cat("\n---- (c) new vs lost dimorphic genes ----\n")
new_dim  <- setdiff(dim_s, dim_c)   # dimorphic in lvqw, not in dream
lost_dim <- setdiff(dim_c, dim_s)   # dimorphic in dream, not in lvqw
cat("  NEW dimorphic (lvqw-only):", length(new_dim), "\n")
cat("  LOST dimorphic (dream-only):", length(lost_dim), "\n")
cat("  retained dimorphic:", length(intersect(dim_c, dim_s)), "\n")
# break new/lost down by the lvqw / dream class they carry
cat("  NEW dimorphic by lvqw class:\n"); print(table(m[gene %in% new_dim, class_s]))
cat("  LOST dimorphic by dream class:\n"); print(table(m[gene %in% lost_dim, class_c]))

## (d) SE ratio (lvqw vs dream interaction SE) -------------------------------
cat("\n---- (d) interaction SE ratio (lvqw / dream) ----\n")
# sex_interaction_dream_v3.csv carries se_int (added by 04). Fall back to
# deriving SE = |logFC / t| from the interaction file if se_int absent.
get_se <- function(dt) {
  if ("se_int" %in% names(dt)) return(dt[, .(gene, se = se_int)])
  dt[, .(gene, se = abs(logFC / t))]
}
se_c <- get_se(ci); se_s <- get_se(si)
se_m <- merge(se_c, se_s, by = "gene", suffixes = c("_c","_s"))
se_m <- se_m[is.finite(se_c) & is.finite(se_s) & se_c > 0]
ratio <- se_m[, se_s / se_c]
cat(sprintf("  n genes with both SE: %d\n", nrow(se_m)))
cat(sprintf("  median se_int ratio (lvqw/dream): %.4f\n", median(ratio, na.rm = TRUE)))
cat(sprintf("  q25 / q75: %.4f / %.4f\n",
            quantile(ratio, .25, na.rm = TRUE), quantile(ratio, .75, na.rm = TRUE)))
cat(sprintf("  PREDICTION CHECK: median ratio < 1.0 (LVQW SEs smaller)?  %s\n",
            if (isTRUE(median(ratio, na.rm = TRUE) < 1.0)) "YES (as predicted)" else "NO"))

## (e) example class flips ---------------------------------------------------
cat("\n---- (e) example class flips (dream -> lvqw) ----\n")
flips <- m[class_c != class_s]
cat("  total class flips on shared genes:", nrow(flips), "\n")
cat("  top flip transitions:\n")
print(head(sort(table(paste(flips$class_c, "->", flips$class_s)), decreasing = TRUE), 12))
# a few concrete examples, prioritising flips INTO a dimorphic class
ex <- flips[class_s != "Concordant"]
ex <- ex[order(-abs(logFC_M_s - logFC_F_s))][1:min(12, .N)]
gm <- if ("gene_symbol" %in% names(sc)) sc[, .(gene, gene_symbol)] else NULL
if (!is.null(gm)) ex <- merge(ex, gm, by = "gene", all.x = TRUE)
cat("  example genes flipping INTO dimorphic (largest |M-F| gap):\n")
print(ex[, .(gene,
             gene_symbol = if ("gene_symbol" %in% names(ex)) gene_symbol else NA,
             dream = class_c, lvqw = class_s,
             lfcF_d = round(logFC_F_c,3), lfcF_l = round(logFC_F_s,3),
             lfcM_d = round(logFC_M_c,3), lfcM_l = round(logFC_M_s,3))])

## write machine-readable summaries -----------------------------------------
fwrite(cmp, file.path(STAGE, "divergence_class_counts.csv"))
fwrite(data.table(metric = c("jaccard_dimorphic","dir_conc_M_all","dir_conc_F_all",
                             "n_new_dimorphic","n_lost_dimorphic","n_retained_dimorphic",
                             "median_se_ratio_lvqw_over_dream","n_class_flips"),
                  value  = c(jacc, dir_M, dir_F, length(new_dim), length(lost_dim),
                             length(intersect(dim_c, dim_s)),
                             median(ratio, na.rm = TRUE), nrow(flips))),
       file.path(STAGE, "divergence_summary.csv"))
fwrite(flips, file.path(STAGE, "divergence_class_flips.csv"))
cat("\nWrote: divergence_{class_counts,summary,class_flips}.csv to lvqw_staging/\n")
cat("Done.\n")
