#!/usr/bin/env Rscript
# sex_v3/07b2_v6_loco_magnitude.R
# ---------------------------------------------------------------------------
# Pillar 8B — Magnitude-aware LOCO with bootstrap CIs.
#
# Inputs: intermediates/loco_v5/fold_{1..5}.rds  (each carries loco_dt with
# columns gene, beta_int_pred, se_int_pred, p_int_pred, beta_int_obs,
# se_int_obs, p_int_obs, sign_match, cohort_held, n_held, fold, ...)
#
# Per fold × class (joining the v5 interaction classifier):
#   - Spearman rho(beta_int_pred, beta_int_obs) with 1000-resample bootstrap CI
#   - sign_concord = mean(sign_match)
#   - magnitude_concord = fraction(|beta_obs|/|beta_pred| in [0.5, 2.0])
#   - holdout_padj_lt_0.20 = fraction(p.adjust(p_int_obs, "BH") < 0.20)
#
# Output: loco_concordance_per_class_v6.csv
# Cols: class, fold, n_genes, spearman_rho, rho_ci_lo, rho_ci_hi,
#       sign_concord, magnitude_concord, holdout_padj_lt_0.20
#
# Compute: ~1k bootstrap iterations × 5 folds × ~10 classes; <30 min single-core
# but uses 4 CPU for parallel rho-bootstrap via mclapply. cpu / 32G / 2h.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SEXV3 <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                   "results/integration/sex_v3")

# Contrast routing (sex_v3_utils.R::contrast_paths) — overrides SEXV3 / IDIR
# when CONTRAST_NAME != "disease_vs_ctrl"; default preserves legacy layout.
if (!exists("contrast_paths")) source(file.path(
  Sys.getenv("MASLD_PROJECT_ROOT",
             "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3/sex_v3_utils.R"))
.cpaths <- contrast_paths()
SEXV3 <- .cpaths$sexv3
IDIR  <- .cpaths$idir
dir.create(IDIR, recursive = TRUE, showWarnings = FALSE)
LOCO  <- file.path(SEXV3, "intermediates/loco_v5")
UTILS <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                   "scripts/sex_v3/sex_v3_utils.R")
if (file.exists(UTILS)) source(UTILS)

OUT_CSV <- file.path(SEXV3, "loco_concordance_per_class_v6.csv")
LOG_TXT <- file.path(SEXV3, "logs/loco_magnitude_v6.log")
dir.create(dirname(LOG_TXT), showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) {
  msg <- sprintf("[%s] %s",
                 format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                 paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_TXT, append = TRUE)
}

N_CPU  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
N_BOOT <- as.integer(Sys.getenv("LOCO_N_BOOT", "1000"))
SEED   <- 42

# ---------------------------------------------------------------------------
# 1. Load LOCO folds
# ---------------------------------------------------------------------------
rds_files <- list.files(LOCO, pattern = "^fold_.*\\.rds$", full.names = TRUE)
log_msg("LOCO folds found: ", length(rds_files))
if (length(rds_files) == 0) stop("No LOCO fold files under ", LOCO)

combined <- rbindlist(lapply(rds_files, function(f) {
  o <- readRDS(f)
  d <- o$loco_dt
  d[, cohort_held := o$cohort_held]
  d[, n_held := o$n_held]
  d[, fold := o$fold]
  d
}), fill = TRUE)
log_msg("Combined LOCO rows: ", nrow(combined),
        "  unique genes: ", uniqueN(combined$gene),
        "  folds: ", uniqueN(combined$fold))

# ---------------------------------------------------------------------------
# 2. Load v5 classifier for per-class assignment
# ---------------------------------------------------------------------------
v5_csv <- file.path(SEXV3, "interaction_classifier_v5.csv")
class_dt <- NULL
class_src <- NA_character_
if (file.exists(v5_csv)) {
  v5 <- fread(v5_csv)
  if ("class_v5_interaction" %in% names(v5)) {
    class_dt <- v5[, .(gene, class = class_v5_interaction)]
    class_src <- "v5:class_v5_interaction"
  }
}
if (is.null(class_dt)) {
  log_msg("WARNING: no v5 classifier — will report ALL-class concordance only")
  class_dt <- data.table(gene = unique(combined$gene), class = "ALL")
  class_src <- "fallback:ALL"
}
log_msg("Class source: ", class_src, " (", nrow(class_dt), " gene assignments)")

m <- merge(combined, class_dt, by = "gene", all.x = TRUE)
m[is.na(class), class := "Unclassified"]

# ---------------------------------------------------------------------------
# 3. Per fold × class metric computation
# ---------------------------------------------------------------------------
spearman_boot_ci <- function(x, y, n_boot = N_BOOT, seed = SEED,
                              alpha = 0.05, ncpu = N_CPU) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 5) return(c(rho = NA_real_, lo = NA_real_, hi = NA_real_))
  rho <- suppressWarnings(cor(x, y, method = "spearman"))
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, n_boot)
  boot_rhos <- unlist(mclapply(seq_len(n_boot), function(b) {
    set.seed(seeds[b])
    idx <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  }, mc.cores = min(ncpu, n_boot)))
  ci <- quantile(boot_rhos, c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  c(rho = rho, lo = unname(ci[1]), hi = unname(ci[2]))
}

magnitude_concord <- function(beta_pred, beta_obs) {
  ok <- is.finite(beta_pred) & is.finite(beta_obs) & abs(beta_pred) > 1e-8
  if (sum(ok) < 5) return(NA_real_)
  ratio <- abs(beta_obs[ok]) / abs(beta_pred[ok])
  mean(ratio >= 0.5 & ratio <= 2.0, na.rm = TRUE)
}

holdout_padj_lt_thr <- function(p_obs, thr = 0.20) {
  ok <- is.finite(p_obs)
  if (sum(ok) < 5) return(NA_real_)
  mean(p.adjust(p_obs[ok], "BH") < thr, na.rm = TRUE)
}

classes <- sort(unique(m$class))
folds   <- sort(unique(m$fold))
log_msg("Classes (", length(classes), "): ", paste(classes, collapse = ", "))
log_msg("Folds: ", paste(folds, collapse = ", "))

rows <- list()
for (cl in classes) {
  for (fd in folds) {
    sub <- m[class == cl & fold == fd]
    n_genes <- nrow(sub)
    if (n_genes < 5) {
      rows[[length(rows) + 1]] <- data.table(
        class = cl, fold = fd, n_genes = n_genes,
        spearman_rho = NA_real_, rho_ci_lo = NA_real_, rho_ci_hi = NA_real_,
        sign_concord = NA_real_,
        magnitude_concord = NA_real_,
        holdout_padj_lt_0.20 = NA_real_
      )
      next
    }
    sp <- spearman_boot_ci(sub$beta_int_pred, sub$beta_int_obs)
    sc <- mean(sub$sign_match, na.rm = TRUE)
    mc <- magnitude_concord(sub$beta_int_pred, sub$beta_int_obs)
    hp <- holdout_padj_lt_thr(sub$p_int_obs, 0.20)
    rows[[length(rows) + 1]] <- data.table(
      class = cl, fold = fd, n_genes = n_genes,
      spearman_rho = sp[["rho"]],
      rho_ci_lo = sp[["lo"]],
      rho_ci_hi = sp[["hi"]],
      sign_concord = sc,
      magnitude_concord = mc,
      holdout_padj_lt_0.20 = hp
    )
  }
}
out_dt <- rbindlist(rows, fill = TRUE)
log_msg("Output rows: ", nrow(out_dt))

# Per-class pooled (across folds, weighted by n_genes)
per_class_summary <- out_dt[, .(
  n_folds = sum(!is.na(spearman_rho)),
  total_genes = sum(n_genes, na.rm = TRUE),
  mean_rho = weighted.mean(spearman_rho, n_genes, na.rm = TRUE),
  mean_sign = weighted.mean(sign_concord, n_genes, na.rm = TRUE),
  mean_magnitude = weighted.mean(magnitude_concord, n_genes, na.rm = TRUE),
  mean_holdout_padj = weighted.mean(holdout_padj_lt_0.20, n_genes, na.rm = TRUE)
), by = class][order(-total_genes)]

log_msg("Per-class pooled summary:")
print(per_class_summary)

if (exists("write_atomic_csv")) {
  write_atomic_csv(out_dt, OUT_CSV)
} else {
  fwrite(out_dt, OUT_CSV)
}
log_msg("Wrote: ", OUT_CSV)

# Acceptance gate print
for (cls in c("Female_biased", "Male_biased", "Male_biased_F_underpowered")) {
  row <- per_class_summary[class == cls]
  if (nrow(row) > 0) {
    v <- row$mean_rho
    verdict <- if (is.finite(v) && v >= 0.5) "PASS" else "FAIL/n/a"
    log_msg(sprintf("  [%s]  mean_spearman_rho = %.4f   gate>=0.5 : %s",
                    cls, v, verdict))
  }
}

cat("\nDone:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
