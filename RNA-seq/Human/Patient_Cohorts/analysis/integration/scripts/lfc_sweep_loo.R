#!/usr/bin/env Rscript
# lfc_sweep_loo.R
# ---------------------------------------------------------------------------
# Lean LFC threshold sweep for new contrasts (MASH-vs-MASL, MASH-vs-Healthy),
# parallel to the CV-of-DEG-count plateau analysis that originally chose the
# Tier 1 |logFC|>0.5 cutoff for the disease_vs_control mega-analysis.
#
# At each |logFC| threshold in LFC_GRID and a fixed padj cut, we compute:
#   - n_full         : DEG count from the full dream mega CSV
#   - n_fold_mean/sd : DEG count across leave-one-cohort-out folds
#   - cv             : sd(n_fold) / mean(n_fold)              <- plateau metric
#   - jaccard_mean   : mean Jaccard(fold DEGs, full DEGs)
#   - dir_concord_mean: mean sign concordance on the intersection
#
# Intentionally simpler than compare_lfc_threshold_sweep.R (no CPSS, no
# bootstrap, no k-fold — just LOO).
#
# Env vars:
#   CONTRAST  : disease_vs_control | mash_vs_masl | mash_vs_healthy   (required)
#   MASH_DEF  : borderline_grouped | strict   (default borderline_grouped;
#                                              ignored for disease_vs_control)
#   LFC_GRID  : comma-separated thresholds (default 0..1.5 step 0.1)
#   PADJ_CUT  : padj cutoff (default 0.05)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

BASE_ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
BASE      <- file.path(BASE_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
source(file.path(BASE_ROOT, "scripts/figures/publication_theme.R"))

# ---------------------------------------------------------------------------
# Env-var contract
# ---------------------------------------------------------------------------
CONTRAST <- Sys.getenv("CONTRAST", "")
MASH_DEF <- Sys.getenv("MASH_DEF", "borderline_grouped")
LFC_GRID <- Sys.getenv("LFC_GRID",
                       "0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5")
PADJ_CUT <- as.numeric(Sys.getenv("PADJ_CUT", "0.05"))

if (!CONTRAST %in% c("disease_vs_control", "mash_vs_masl", "mash_vs_healthy", "masl_vs_healthy")) {
  stop("CONTRAST must be one of: disease_vs_control, mash_vs_masl, mash_vs_healthy, masl_vs_healthy. Got: '",
       CONTRAST, "'")
}
if (CONTRAST %in% c("mash_vs_masl", "mash_vs_healthy") && !MASH_DEF %in% c("borderline_grouped", "strict")) {
  stop("MASH_DEF must be 'borderline_grouped' or 'strict'. Got: '", MASH_DEF, "'")
}

lfc_grid <- as.numeric(strsplit(LFC_GRID, ",", fixed = TRUE)[[1]])
if (any(is.na(lfc_grid))) stop("LFC_GRID parse error: ", LFC_GRID)

# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------
contrast_tag <- if (CONTRAST %in% c("disease_vs_control", "masl_vs_healthy")) {
  CONTRAST
} else if (MASH_DEF == "strict") {
  paste0(CONTRAST, "_strict")
} else {
  CONTRAST
}

resolve_paths <- function(contrast, mash_def) {
  if (contrast == "disease_vs_control") {
    full <- file.path(BASE, "results/integration/dream_results.csv")
    loo  <- file.path(BASE, "results/integration/loo_cv/disease_vs_control")
    if (!dir.exists(loo)) {
      # backward-compat with legacy dream_loo_cv.R output
      legacy <- file.path(BASE, "results/integration/loo_cv")
      if (dir.exists(legacy)) loo <- legacy
    }
    return(list(full = full, loo = loo))
  }
  if (contrast == "mash_vs_masl" && mash_def == "borderline_grouped") {
    return(list(
      full = file.path(BASE, "results/disease_signatures/mash_vs_masl_dream.csv"),
      loo  = file.path(BASE, "results/integration/loo_cv/mash_vs_masl")))
  }
  if (contrast == "mash_vs_masl" && mash_def == "strict") {
    return(list(
      full = file.path(BASE, "results/disease_signatures/mash_vs_masl_dream_strict.csv"),
      loo  = file.path(BASE, "results/integration/loo_cv/mash_vs_masl_strict")))
  }
  if (contrast == "mash_vs_healthy" && mash_def == "borderline_grouped") {
    return(list(
      full = file.path(BASE, "results/disease_signatures/mash_vs_healthy_dream.csv"),
      loo  = file.path(BASE, "results/integration/loo_cv/mash_vs_healthy")))
  }
  if (contrast == "mash_vs_healthy" && mash_def == "strict") {
    return(list(
      full = file.path(BASE, "results/disease_signatures/mash_vs_healthy_dream_strict.csv"),
      loo  = file.path(BASE, "results/integration/loo_cv/mash_vs_healthy_strict")))
  }
  if (contrast == "masl_vs_healthy") {
    return(list(
      full = file.path(BASE, "results/disease_signatures/masl_vs_healthy_dream.csv"),
      loo  = file.path(BASE, "results/integration/loo_cv/masl_vs_healthy")))
  }
  stop("Unhandled contrast/mash_def combination: ", contrast, " / ", mash_def)
}

paths <- resolve_paths(CONTRAST, MASH_DEF)
if (!file.exists(paths$full)) stop("Missing full dream CSV: ", paths$full)
if (!dir.exists(paths$loo))   stop("Missing LOO folds dir: ", paths$loo)

fold_files <- list.files(paths$loo, pattern = "^dream_loo_.*\\.csv$", full.names = TRUE)
if (length(fold_files) == 0) {
  stop("No LOO fold files (dream_loo_*.csv) under: ", paths$loo)
}

cat(sprintf("CONTRAST=%s  MASH_DEF=%s  PADJ_CUT=%g\n", CONTRAST, MASH_DEF, PADJ_CUT))
cat(sprintf("Full dream : %s\n", paths$full))
cat(sprintf("LOO dir    : %s  (%d folds)\n", paths$loo, length(fold_files)))

# ---------------------------------------------------------------------------
# Load + column normalization
# ---------------------------------------------------------------------------
read_dream <- function(path) {
  dt <- fread(path)
  if (!"padj" %in% names(dt) && "adj.P.Val" %in% names(dt)) {
    setnames(dt, "adj.P.Val", "padj")
  }
  if (!all(c("gene", "logFC", "padj") %in% names(dt))) {
    stop("Required columns gene/logFC/padj missing in: ", path,
         "  (have: ", paste(names(dt), collapse = ", "), ")")
  }
  dt[, .(gene, logFC, padj)]
}

full <- read_dream(paths$full)
folds <- lapply(fold_files, read_dream)
names(folds) <- sub("^dream_loo_", "", sub("\\.csv$", "", basename(fold_files)))

# ---------------------------------------------------------------------------
# Sweep
# ---------------------------------------------------------------------------
sweep_one <- function(x) {
  full_sig  <- full[!is.na(padj) & padj < PADJ_CUT & abs(logFC) > x]
  full_set  <- unique(full_sig$gene)
  n_full    <- length(full_set)
  full_sign <- setNames(sign(full_sig$logFC), full_sig$gene)

  per_fold <- lapply(folds, function(dt) {
    sig <- dt[!is.na(padj) & padj < PADJ_CUT & abs(logFC) > x]
    set <- unique(sig$gene)
    n_f <- length(set)
    inter <- intersect(set, full_set)
    uni   <- length(union(set, full_set))
    jacc  <- if (uni == 0) NA_real_ else length(inter) / uni
    if (length(inter) == 0) {
      sign_f <- NA_real_
    } else {
      fold_sign <- setNames(sign(sig$logFC), sig$gene)[inter]
      sign_f    <- mean(fold_sign == full_sign[inter])
    }
    list(n = n_f, jacc = jacc, sign = sign_f)
  })

  n_vec    <- vapply(per_fold, function(z) z$n,    numeric(1))
  jacc_vec <- vapply(per_fold, function(z) z$jacc, numeric(1))
  sign_vec <- vapply(per_fold, function(z) z$sign, numeric(1))

  n_mean <- mean(n_vec)
  n_sd   <- sd(n_vec)
  cv     <- if (n_full == 0 || n_mean == 0) NA_real_ else n_sd / n_mean

  data.table(
    contrast         = CONTRAST,
    mash_def         = if (CONTRAST == "disease_vs_control") NA_character_ else MASH_DEF,
    lfc_thr          = x,
    padj_cut         = PADJ_CUT,
    n_full           = n_full,
    n_fold_mean      = n_mean,
    n_fold_sd        = n_sd,
    cv               = cv,
    jaccard_mean     = mean(jacc_vec, na.rm = TRUE),
    dir_concord_mean = mean(sign_vec, na.rm = TRUE),
    n_folds          = length(folds))
}

res <- rbindlist(lapply(lfc_grid, sweep_one))
cat("\n=== Sweep results ===\n"); print(res)

# ---------------------------------------------------------------------------
# Recommended cutoff picker
# ---------------------------------------------------------------------------
cand <- res[!is.na(cv) &
              dir_concord_mean >= 0.9 &
              jaccard_mean    >= 0.7 &
              n_full          >= 100]

rec_thr <- NA_real_
plateau_thrs <- numeric(0)

if (nrow(cand) == 0) {
  warning("No sweep rows passed filters (dir_concord>=0.9, jaccard>=0.7, n_full>=100). ",
          "Skipping recommended-cutoff pick.")
} else {
  cv_min       <- min(cand$cv)
  plateau_thrs <- cand[cv <= cv_min * 1.10, lfc_thr]
  rec_thr      <- min(plateau_thrs)
  rec_row      <- cand[lfc_thr == rec_thr][1]
  tag_disp     <- contrast_tag
  cat(sprintf(
    "\nRECOMMENDED CUTOFF for %s: padj < %g & |logFC| > %g (CV = %.4f, n_full = %d, jaccard = %.3f)\n",
    tag_disp, PADJ_CUT, rec_thr, rec_row$cv, rec_row$n_full, rec_row$jaccard_mean))
}

# ---------------------------------------------------------------------------
# Output CSV
# ---------------------------------------------------------------------------
audit_dir <- file.path(BASE_ROOT, "RNA-seq/results/audit_sensitivity")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
csv_path <- file.path(audit_dir, sprintf("lfc_sweep_%s.csv", contrast_tag))
fwrite(res, csv_path)
cat(sprintf("Saved: %s\n", csv_path))

# ---------------------------------------------------------------------------
# 4-panel figure
# ---------------------------------------------------------------------------
fig_dir <- file.path(BASE_ROOT, "figures/supplementary/lfc_sweep")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

up_col   <- masld_colors$up
down_col <- masld_colors$down
ribb_col <- masld_colors$up

# Panel A — DEG counts
plot_dt_A <- rbind(
  data.table(lfc_thr = res$lfc_thr, value = res$n_full,
             lo = res$n_full, hi = res$n_full, series = "Full dream"),
  data.table(lfc_thr = res$lfc_thr, value = res$n_fold_mean,
             lo = pmax(res$n_fold_mean - res$n_fold_sd, 0),
             hi = res$n_fold_mean + res$n_fold_sd,
             series = "LOO fold mean ± sd"))

pA <- ggplot(plot_dt_A, aes(x = lfc_thr, y = value, colour = series, fill = series)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.1) +
  scale_colour_manual(values = c("Full dream" = down_col,
                                  "LOO fold mean ± sd" = up_col)) +
  scale_fill_manual(values = c("Full dream" = down_col,
                                "LOO fold mean ± sd" = up_col)) +
  scale_y_continuous(labels = comma) +
  labs(title = "A — DEG count vs |logFC| threshold",
       x = "|logFC| threshold", y = "DEGs", colour = NULL, fill = NULL) +
  theme_masld() + theme_pub() + theme(legend.position = "bottom")

# Panel B — CV with plateau shading + recommended cutoff
pB <- ggplot(res, aes(x = lfc_thr, y = cv))
if (length(plateau_thrs) > 0) {
  pB <- pB + annotate("rect",
                       xmin = min(plateau_thrs) - 0.025,
                       xmax = max(plateau_thrs) + 0.025,
                       ymin = -Inf, ymax = Inf,
                       alpha = 0.12, fill = up_col)
}
pB <- pB +
  geom_line(linewidth = 0.5, colour = up_col) +
  geom_point(size = 1.1, colour = up_col)
# Mark canonical uniform Tier 1 cutoff |LFC|>0.5 (2026-05-13 decision).
TIER1_LFC <- 0.5
pB <- pB + geom_vline(xintercept = TIER1_LFC, linetype = "dashed",
                       colour = "#FFB300", linewidth = 0.5)
cv_at_t1 <- res[lfc_thr == TIER1_LFC, cv]
n_at_t1  <- res[lfc_thr == TIER1_LFC, n_full]
pB <- pB +
  labs(title = "B — CV of DEG count across LOO folds",
       x = "|logFC| threshold",
       y = "CV = sd(n_fold) / mean(n_fold)",
       caption = sprintf("Uniform Tier 1: |logFC|>%g (gold dashed). CV at Tier 1 = %.3f, n_full = %s.",
                         TIER1_LFC,
                         ifelse(length(cv_at_t1), cv_at_t1, NA_real_),
                         format(ifelse(length(n_at_t1), n_at_t1, NA_integer_), big.mark = ","))) +
  theme_masld() + theme_pub()

# Panel C — Jaccard
pC <- ggplot(res, aes(x = lfc_thr, y = jaccard_mean)) +
  geom_hline(yintercept = 0.7, linetype = "dashed", colour = "grey25", linewidth = 0.4) +
  geom_line(linewidth = 0.5, colour = down_col) +
  geom_point(size = 1.1, colour = down_col) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "C — Jaccard(fold DEGs, full DEGs)",
       x = "|logFC| threshold", y = "Mean Jaccard across folds") +
  theme_masld() + theme_pub()

# Panel D — Direction concordance
pD <- ggplot(res, aes(x = lfc_thr, y = dir_concord_mean)) +
  geom_hline(yintercept = 0.9, linetype = "dashed", colour = "grey25", linewidth = 0.4) +
  geom_line(linewidth = 0.5, colour = masld_colors$conserved) +
  geom_point(size = 1.1, colour = masld_colors$conserved) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "D — Direction concordance on intersection",
       x = "|logFC| threshold", y = "Mean sign concordance") +
  theme_masld() + theme_pub()

# Strip the "_strict" / "_borderline_grouped" suffix from the display title;
# keep contrast_tag for the filename only.
display_title <- sub("_strict$", "", contrast_tag)
display_title <- gsub("_", " ", display_title)
display_title <- gsub("mash vs masl",    "MASH vs MASL",    display_title, ignore.case = TRUE)
display_title <- gsub("mash vs healthy", "MASH vs Healthy", display_title, ignore.case = TRUE)
display_title <- gsub("masl vs healthy", "MASL vs Healthy", display_title, ignore.case = TRUE)
display_title <- gsub("disease vs control", "Disease vs Control", display_title, ignore.case = TRUE)

combined <- (pA | pB) / (pC | pD) + plot_layout(ncol = 1) +
  plot_annotation(
    title = sprintf("LFC sweep — %s (padj < %g)", display_title, PADJ_CUT),
    subtitle = sprintf("LOO across %d cohorts; %d thresholds; uniform Tier 1 = |LFC|>%g",
                       length(folds), nrow(res), TIER1_LFC),
    theme = theme(plot.title    = element_text(size = PUB_TITLE, face = "bold"),
                  plot.subtitle = element_text(size = PUB_SUBTITLE, colour = "grey30")))

pdf_path <- file.path(fig_dir, sprintf("%s_lfc_sweep.pdf", contrast_tag))
ggsave(pdf_path, combined, width = 9, height = 8, device = cairo_pdf)
cat(sprintf("Saved: %s\n", pdf_path))
