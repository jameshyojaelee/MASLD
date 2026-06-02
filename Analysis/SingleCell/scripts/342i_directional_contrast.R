#!/usr/bin/env Rscript
# 342i_directional_contrast.R - Directional concordance test per critical-review recommendation.
# Compute a SIGNED per-sample score that uses gene_sign * z-scored expression.
# This avoids the "ssGSEA inflates both UP and DOWN" artifact.
# Score = mean(z-scored expression of UP genes) - mean(z-scored expression of DOWN genes)
#       = mean( sign_g * z_g )
# where sign_g = +1 for genes UP in polyploid, -1 for genes DOWN in polyploid.

suppressPackageStartupMessages({
  library(data.table)
  library(lme4); library(lmerTest)
})
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
setwd(ROOT)
SIG_DIR <- "data/ploidy_signatures"
OUT_DIR <- "Analysis/SingleCell/results_gpu_v2/ploidy"

cat("[", format(Sys.time()), "] Loading bulk log-CPM\n", sep="")
expr <- readRDS("RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/corrected_logcpm.rds")
gm <- fread("data/gencode_v49_gene_metadata.tsv.gz")
ens2sym <- setNames(gm$gene_name, gm$gene_id)
rownames(expr) <- ens2sym[rownames(expr)]
expr <- expr[!is.na(rownames(expr)) & rownames(expr) != "", ]
# Collapse duplicates by max expression
dup <- duplicated(rownames(expr))
if (sum(dup) > 0) {
  ord <- order(rowMeans(expr), decreasing = TRUE)
  expr <- expr[ord, ]
  expr <- expr[!duplicated(rownames(expr)), ]
}
cat("  matrix after symbol mapping:", dim(expr), "\n")
# Z-score per gene across samples
expr_z <- t(scale(t(expr)))
expr_z[is.na(expr_z)] <- 0  # genes with zero variance
cat("  z-scored\n")

# Load signed signatures
load_genes <- function(f, col = "human_symbol") {
  d <- fread(f)
  if (!col %in% names(d)) col <- intersect(c("gene_symbol", "Gene", "human_symbol"), names(d))[1]
  unique(toupper(as.character(d[[col]])))
}
ric_up <- load_genes(file.path(SIG_DIR, "richter2021_polyploid_up.tsv"))
ric_dn <- load_genes(file.path(SIG_DIR, "richter2021_polyploid_down.tsv"))
kat_up <- load_genes(file.path(SIG_DIR, "katsuda2019_polyploid_up.tsv"))
kat_dn <- load_genes(file.path(SIG_DIR, "katsuda2019_polyploid_down.tsv"))
cons_up <- load_genes(file.path(SIG_DIR, "consensus_polyploid_up.tsv"))
cons_dn <- load_genes(file.path(SIG_DIR, "consensus_polyploid_down.tsv"))

# Match to expr genes
expr_genes_upper <- toupper(rownames(expr))
sym_map <- setNames(rownames(expr), expr_genes_upper)
in_expr <- function(g) sym_map[g][!is.na(sym_map[g])]

ric_up_e <- in_expr(ric_up); ric_dn_e <- in_expr(ric_dn)
kat_up_e <- in_expr(kat_up); kat_dn_e <- in_expr(kat_dn)
cons_up_e <- in_expr(cons_up); cons_dn_e <- in_expr(cons_dn)
cat("\nSignature gene coverage in bulk expression matrix:\n")
cat("  Richter UP/DN:", length(ric_up_e), "/", length(ric_dn_e), "\n")
cat("  Katsuda UP/DN:", length(kat_up_e), "/", length(kat_dn_e), "\n")
cat("  Consensus UP/DN:", length(cons_up_e), "/", length(cons_dn_e), "\n")

# Compute signed score per sample
signed_score <- function(up_genes, dn_genes) {
  up_in <- intersect(up_genes, rownames(expr_z))
  dn_in <- intersect(dn_genes, rownames(expr_z))
  if (length(up_in) < 2 && length(dn_in) < 2) return(rep(NA, ncol(expr_z)))
  up_mean <- if (length(up_in) >= 2) colMeans(expr_z[up_in, , drop = FALSE]) else 0
  dn_mean <- if (length(dn_in) >= 2) colMeans(expr_z[dn_in, , drop = FALSE]) else 0
  up_mean - dn_mean
}

scores <- data.table(
  sample_id = colnames(expr_z),
  directional_richter   = signed_score(ric_up_e, ric_dn_e),
  directional_katsuda   = signed_score(kat_up_e, kat_dn_e),
  directional_consensus = signed_score(cons_up_e, cons_dn_e),
  # As negative control, single-direction means too (should both go up if non-specific shift)
  mean_z_richter_up   = colMeans(expr_z[intersect(ric_up_e, rownames(expr_z)), , drop=FALSE]),
  mean_z_richter_down = colMeans(expr_z[intersect(ric_dn_e, rownames(expr_z)), , drop=FALSE]),
  mean_z_katsuda_up   = colMeans(expr_z[intersect(kat_up_e, rownames(expr_z)), , drop=FALSE]),
  mean_z_katsuda_down = colMeans(expr_z[intersect(kat_dn_e, rownames(expr_z)), , drop=FALSE])
)

# Merge with metadata + hep fraction
meta <- fread("RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
ct <- fread("RNA-seq/results/celltype_attribution/persample_celltype_proportions.csv")
ct <- ct[, .(sample_id, hepatocyte_fraction = Hepatocytes)]
m <- merge(scores, meta, by = "sample_id")
m <- merge(m, ct, by = "sample_id", all.x = TRUE)

# Save per-sample directional scores (used by figure 03 to plot jittered samples per stage)
fwrite(m[, c("sample_id", "dataset", "fibrosis_stage", "nas_score", "hepatocyte_fraction",
             grep("^directional|^mean_z", names(m), value = TRUE)), with = FALSE],
       file.path(OUT_DIR, "directional_per_sample.csv"))

cat("\n[", format(Sys.time()), "] Disease tests on DIRECTIONAL scores\n", sep="")
score_cols <- grep("^directional|^mean_z", names(m), value = TRUE)
results <- list()
for (sc in score_cols) {
  sub <- m[!is.na(fibrosis_stage) & !is.na(get(sc))]
  # With cohort RE + hep fraction
  s <- if (sum(!is.na(sub$hepatocyte_fraction)) > 0.8 * nrow(sub)) {
    sub_full <- sub[!is.na(hepatocyte_fraction)]
    mm <- tryCatch(lmer(sub_full[[sc]] ~ fibrosis_stage + hepatocyte_fraction + (1|dataset),
                        data = sub_full), error = function(e) NULL)
    if (!is.null(mm)) summary(mm)$coefficients["fibrosis_stage", ] else NA
  } else NA
  # Unadjusted
  m_u <- lm(sub[[sc]] ~ sub$fibrosis_stage)
  cf <- summary(m_u)$coefficients["sub$fibrosis_stage", ]
  results[[sc]] <- data.table(
    signature = sc,
    n = nrow(sub),
    est_unadj = cf["Estimate"],
    p_unadj = cf["Pr(>|t|)"],
    est_mixed_hep = if (length(s) > 1) s["Estimate"] else NA,
    p_mixed_hep   = if (length(s) > 1) s["Pr(>|t|)"]  else NA
  )
}
res <- rbindlist(results)
res[, padj_mixed := p.adjust(p_mixed_hep, method = "BH")]
res <- res[order(p_mixed_hep)]
fwrite(res, file.path(OUT_DIR, "directional_concordance_tests.csv"))

cat("\n=== DIRECTIONAL CONCORDANCE TESTS (mixed-effects, hep-adjusted) ===\n")
print(res)

# Per-stage means for directional scores
m[, fibrosis_stage_f := factor(fibrosis_stage, levels = 0:4)]
per_stage <- m[!is.na(fibrosis_stage_f), lapply(.SD, mean, na.rm=TRUE),
              by = fibrosis_stage_f, .SDcols = score_cols][order(fibrosis_stage_f)]
fwrite(per_stage, file.path(OUT_DIR, "directional_per_stage.csv"))
cat("\n=== Directional scores by fibrosis stage ===\n")
print(per_stage)

# CRITICAL INTERPRETATION:
# If directional_* is null AND mean_z_*_up + mean_z_*_down both go UP -> non-specific shift
# If directional_* significant in correct direction -> genuine polyploid signal
cat("\n=== INTERPRETATION ===\n")
for (sig in c("directional_richter", "directional_katsuda", "directional_consensus")) {
  if (sig %in% res$signature) {
    row <- res[signature == sig]
    sign <- ifelse(row$est_mixed_hep > 0, "+", "-")
    cat(sig, ": est=", round(row$est_mixed_hep, 4), " p_mixed_hep=", signif(row$p_mixed_hep,3),
        " padj=", signif(row$padj_mixed, 3), " direction=", sign, "\n", sep = "")
  }
}

cat("\n=== Magnitude comparison: up-mean vs down-mean per stage ===\n")
print(per_stage[, .(fibrosis_stage_f, mean_z_richter_up, mean_z_richter_down,
                   mean_z_katsuda_up, mean_z_katsuda_down,
                   directional_richter, directional_katsuda)])

cat("\n[", format(Sys.time()), "] DONE directional contrast tests\n", sep="")
