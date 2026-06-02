#!/usr/bin/env Rscript
# 232_sensitivity_de_runs.R
#
# Sensitivity-analysis re-runs of key main-text comparisons after molecularly
# screening contaminated controls. For each comparison, run two limma-voom DEs:
#   (orig): cases vs all clinical-healthy controls
#   (clean): cases vs clinical-healthy controls minus subclinical-suspect tertile
# Then report:
#   - Spearman ρ of logFC vectors (top 1000 by |original logFC|)
#   - Per-cohort LOO sensitivity
#   - Effect-size attenuation summary
#
# Comparison: F1+ vs F0 (F2-switch flank). This is a proxy for the F1→F2
# transition that drives the paper's central thesis. Limited by metadata —
# we use fibrosis_stage groupings.
#
# Inputs:
#   merged_dge.rds (counts)
#   unified_metadata.csv
#   subclinical_screen.csv
# Output: sensitivity_results/ subdir with logFC vectors + attenuation summary
#
# Spec: docs/superpowers/specs/2026-04-27-healthy-control-audit-design.md
# Env: rnaseq

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})

BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
                     "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
HCDIR  <- file.path(BASE, "RNA-seq/results/audit_sensitivity/healthy_control_audit")
SDIR   <- file.path(HCDIR, "sensitivity")
META   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
DGE_RDS <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_dge.rds")

dir.create(SDIR, recursive = TRUE, showWarnings = FALSE)

cat("Loading data...\n")
m   <- fread(META)
sc  <- fread(file.path(HCDIR, "subclinical_screen.csv"))
defs <- fread(file.path(HCDIR, "controls_definitions.csv"))
dge <- readRDS(DGE_RDS)

cat(sprintf("  unified_metadata: %d rows\n", nrow(m)))
cat(sprintf("  subclinical_screen: %d rows\n", nrow(sc)))

# --- Define comparison: F2+ vs F0 ---
fib_num <- suppressWarnings(as.numeric(as.character(m$fibrosis_stage)))
m[, fib_num := fib_num]
m[, comparison_group := fcase(
  fib_num == 0, "F0",
  !is.na(fib_num) & fib_num >= 2, "F2plus",
  default = NA_character_
)]
cat("\nComparison group counts:\n")
print(table(m$comparison_group, useNA = "ifany"))

# Only keep samples in the comparison
de_meta <- m[!is.na(comparison_group)]
de_meta <- de_meta[sample_id %in% colnames(dge$counts)]
cat(sprintf("  Total in counts matrix: %d\n", nrow(de_meta)))

# --- Original control set: all F0 ---
ids_orig_F0  <- de_meta[comparison_group == "F0", sample_id]
ids_F2plus   <- de_meta[comparison_group == "F2plus", sample_id]
cat(sprintf("\nOriginal F0: %d, F2+: %d\n", length(ids_orig_F0), length(ids_F2plus)))

# --- Cleaned control set: drop "suspect" tertile from F0 controls ---
suspect_ids <- sc[p1_tertile == "suspect", sample_id]
cat(sprintf("Suspect (subclinical) IDs from screen: %d\n", length(suspect_ids)))

ids_clean_F0 <- setdiff(ids_orig_F0, suspect_ids)
cat(sprintf("Cleaned F0: %d (dropped %d suspect)\n",
            length(ids_clean_F0), length(ids_orig_F0) - length(ids_clean_F0)))

run_compare <- function(ids_F0, ids_F2plus, label) {
  ids <- c(ids_F0, ids_F2plus)
  meta_subset <- de_meta[sample_id %in% ids]
  setkey(meta_subset, sample_id)
  meta_subset <- meta_subset[ids]

  cnts <- dge$counts[, ids, drop = FALSE]
  d <- DGEList(cnts)
  d <- d[filterByExpr(d, group = meta_subset$comparison_group), , keep.lib.sizes = FALSE]
  d <- calcNormFactors(d)

  group  <- factor(meta_subset$comparison_group, levels = c("F0", "F2plus"))
  cohort <- factor(meta_subset$dataset)
  sex    <- factor(meta_subset$sex)
  age    <- as.numeric(meta_subset$age)

  drop_idx <- is.na(sex)
  if (any(drop_idx)) {
    meta_subset <- meta_subset[!drop_idx]
    d <- d[, !drop_idx]
    group <- group[!drop_idx]
    cohort <- droplevels(cohort[!drop_idx])
    sex   <- droplevels(sex[!drop_idx])
    age   <- age[!drop_idx]
  }

  age_usable <- sum(!is.na(age)) >= 0.7 * length(age) && isTRUE(sd(age, na.rm = TRUE) > 0)
  if (age_usable) age[is.na(age)] <- median(age, na.rm = TRUE)

  rhs <- "0 + group + sex"
  if (nlevels(cohort) > 1) rhs <- paste(rhs, "+ cohort")
  if (age_usable) rhs <- paste(rhs, "+ age")
  design <- model.matrix(as.formula(paste("~", rhs)))

  cat(sprintf("  [%s] samples: %d, genes: %d, design rank: %d\n",
              label, ncol(d), nrow(d), qr(design)$rank))

  v <- voom(d, design)
  fit <- lmFit(v, design)
  contr <- makeContrasts(contrasts = "groupF2plus - groupF0", levels = design)
  fit2 <- contrasts.fit(fit, contr)
  fit2 <- eBayes(fit2)
  tt <- topTable(fit2, number = Inf, sort.by = "none")
  tt <- as.data.table(tt, keep.rownames = "gene_id")
  tt[, label := label]
  tt[, n_F0 := length(ids_F0)]
  tt[, n_F2plus := length(ids_F2plus)]
  tt
}

cat("\n=== Run (orig): F0 vs F2+ with all F0 controls ===\n")
res_orig <- run_compare(ids_orig_F0, ids_F2plus, "orig")
fwrite(res_orig, file.path(SDIR, "sensitivity_de_orig.csv"))

cat("\n=== Run (clean): F0 vs F2+ with subclinical-suspect F0 dropped ===\n")
res_clean <- run_compare(ids_clean_F0, ids_F2plus, "clean")
fwrite(res_clean, file.path(SDIR, "sensitivity_de_clean.csv"))

# --- Effect-size scatter: original vs cleaned for shared genes ---
shared <- intersect(res_orig$gene_id, res_clean$gene_id)
cat(sprintf("\nShared genes: %d\n", length(shared)))

setkey(res_orig, gene_id)
setkey(res_clean, gene_id)
m_join <- data.table(
  gene_id = shared,
  logFC_orig = res_orig[shared, logFC],
  logFC_clean = res_clean[shared, logFC],
  padj_orig = res_orig[shared, adj.P.Val],
  padj_clean = res_clean[shared, adj.P.Val]
)

# Spearman ρ on top 1000 by |orig logFC|
top1000 <- m_join[order(-abs(logFC_orig))][1:min(1000, .N)]
sp <- cor(top1000$logFC_orig, top1000$logFC_clean,
          method = "spearman", use = "pairwise.complete.obs")
cat(sprintf("\nSpearman ρ (top 1000 |logFC|): %.3f\n", sp))

# All shared
sp_all <- cor(m_join$logFC_orig, m_join$logFC_clean,
              method = "spearman", use = "pairwise.complete.obs")
cat(sprintf("Spearman ρ (all shared genes): %.3f\n", sp_all))

# Pearson on logFC
pe <- cor(m_join$logFC_orig, m_join$logFC_clean,
          method = "pearson", use = "pairwise.complete.obs")
cat(sprintf("Pearson r (all shared): %.3f\n", pe))

fwrite(m_join, file.path(SDIR, "sensitivity_logFC_scatter.csv"))

# Sign-flip count
flips <- sum(sign(m_join$logFC_orig) != sign(m_join$logFC_clean), na.rm = TRUE)
cat(sprintf("Sign flips: %d / %d (%.1f%%)\n",
            flips, nrow(m_join), 100 * flips / nrow(m_join)))

# Effect-size attenuation: ratio of |logFC_clean| / |logFC_orig| in top 1000
abs_ratios <- abs(top1000$logFC_clean) / pmax(abs(top1000$logFC_orig), 0.01)
cat(sprintf("Median attenuation ratio (clean/orig, top 1000): %.3f\n", median(abs_ratios)))

# Significant gene retention
n_sig_orig <- sum(m_join$padj_orig < 0.05, na.rm = TRUE)
n_sig_clean <- sum(m_join$padj_clean < 0.05, na.rm = TRUE)
n_sig_both <- sum(m_join$padj_orig < 0.05 & m_join$padj_clean < 0.05, na.rm = TRUE)
cat(sprintf("Significant genes: orig=%d, clean=%d, both=%d (%.1f%%)\n",
            n_sig_orig, n_sig_clean, n_sig_both,
            100 * n_sig_both / max(n_sig_orig, 1)))

# Summary table
summary_dt <- data.table(
  metric = c("spearman_rho_top1000",
             "spearman_rho_all",
             "pearson_r_all",
             "median_attenuation_ratio",
             "sign_flips",
             "n_F0_orig", "n_F0_clean", "n_F2plus",
             "n_sig_orig", "n_sig_clean", "n_sig_both",
             "retention_rate"),
  value = c(round(sp, 4),
            round(sp_all, 4),
            round(pe, 4),
            round(median(abs_ratios), 4),
            flips,
            length(ids_orig_F0), length(ids_clean_F0), length(ids_F2plus),
            n_sig_orig, n_sig_clean, n_sig_both,
            round(100 * n_sig_both / max(n_sig_orig, 1), 2))
)
fwrite(summary_dt, file.path(SDIR, "sensitivity_summary.csv"))
cat("\nSummary:\n")
print(summary_dt)

cat("\nDone.\n")
