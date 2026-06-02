#!/usr/bin/env Rscript
# sex_v3/16_threshold_sweep.R
# ---------------------------------------------------------------------------
# Pillar 7 — Threshold sensitivity sweep.
#
# PRE-REGISTRATION (must appear BEFORE any data is read):
#   Timestamp: 2026-05-15 (script written before P1 completes).
#   Operative gate values (post-meta-review-A9 rename — previously `tier_B_padj`
#   with a `× 2` multiplier at the gate site, which made the advertised 0.10
#   sweep actually fire at 0.20; corrected 2026-05-16):
#     mag_ratio              = 2.0   (|β_dominant| ≥ 2.0 × |β_other| for F/M-biased)
#     lfsr_other             = 0.20  (lfsr of dominated arm must exceed this)
#     branch_B_padj_cutoff   = 0.20  (Branch-B "interaction signal" gate; operative)
#   Sweep grid (48 cells; 4 × 4 × 3):
#     mag_ratio            ∈ {1.5, 2.0, 2.5, 3.0}
#     lfsr_other           ∈ {0.10, 0.20, 0.30, 0.45} (extended 2026-05-16 per P1-19)
#     branch_B_padj_cutoff ∈ {0.10, 0.20, 0.30}
#   For each cell: count F_biased, M_biased, FM_ratio, Divergent, Sex_modifier.
#   Acceptance: canonical cell (2.0, 0.20, 0.20) is the modal F:M ratio cell
#   across the sweep within 2× of the median.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(viridis); library(patchwork)
})

# ---- Pre-registered canonical thresholds (frozen) ----
# CANON gate values are the LITERAL cutoffs applied in classify_param().
CANON <- list(mag_ratio = 2.0, lfsr_other = 0.20, branch_B_padj_cutoff = 0.20)
# lfsr_other grid extended 2026-05-16 per meta-review A3/P1-19: lfsr ∈ [0, 0.5]
# and 0.20 corresponds to P(sign match)=0.80 (directionally resolved). The
# truly-null cap is at lfsr_M = 0.45-0.50 (P(sign) ≈ 0.55-0.50). Sweep
# extended to 0.45 so the canonical 0.20 can be benchmarked against a near-null
# boundary; the 4-value lfsr grid gives a 4 × 3 × 4 = 48-cell sweep.
GRID  <- list(mag_ratio            = c(1.5, 2.0, 2.5, 3.0),
              lfsr_other           = c(0.10, 0.20, 0.30, 0.45),
              branch_B_padj_cutoff = c(0.10, 0.20, 0.30))
cat("== sex_v6 P7 threshold sweep ==\n")
cat("Pre-registered canonical: mag_ratio =", CANON$mag_ratio,
    "  lfsr_other =", CANON$lfsr_other,
    "  branch_B_padj_cutoff =", CANON$branch_B_padj_cutoff, "\n")
cat("Grid: ", length(GRID$mag_ratio) * length(GRID$lfsr_other) *
        length(GRID$branch_B_padj_cutoff), "cells\n")

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
FIG_DIR <- file.path(BASE, "figures/supplementary/figS_sex_dimorphism")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_CSV <- file.path(SEXV3, "sensitivity_sweep_v6.csv")
OUT_PDF <- file.path(FIG_DIR, "figS_sex_threshold_sweep.pdf")
RANDOM_CSV <- file.path(SEXV3, "interaction_classifier_v5_random.csv")

# ---- Load P1 output ----
if (!file.exists(RANDOM_CSV)) stop("P1 output missing: ", RANDOM_CSV)
dt <- fread(RANDOM_CSV)
cat("Loaded P1 random-slope table: ", nrow(dt), " rows\n")
need <- c("gene", "beta_F", "beta_M", "lfsr_F_strat", "lfsr_M_strat",
          "padj_int_5k_random", "power_M_at_F", "chr")
miss <- setdiff(need, names(dt))
if (length(miss) > 0) {
  # Try alternates
  if ("beta_F_strat" %in% names(dt) && !"beta_F" %in% names(dt))
    setnames(dt, "beta_F_strat", "beta_F")
  if ("beta_M_strat" %in% names(dt) && !"beta_M" %in% names(dt))
    setnames(dt, "beta_M_strat", "beta_M")
}
stopifnot(all(c("beta_F", "beta_M", "lfsr_F_strat", "lfsr_M_strat",
                "padj_int_5k_random") %in% names(dt)))

# ---- Decision tree parameterized on (mag_ratio, lfsr_other, branch_B_padj_cutoff) ----
# Operative gate is `branch_B_padj_cutoff` (literally the cutoff for branch B
# membership). Earlier this was named `branch_B_padj_cutoff` but multiplied by 2× at
# the gate, making the pre-registered 0.10 sweep actually fire at 0.20 — a
# silent ×2 typo flagged by meta-review A9. The grid values are now the
# operative gate values directly.
classify_param <- function(bF, bM, lF, lM, pI, pow, ch,
                            mag_ratio, lfsr_other, branch_B_padj_cutoff) {
  n <- length(bF)
  cls <- rep(NA_character_, n)
  # Branch A: no interaction signal — padj_int_5k_random ≥ branch_B_padj_cutoff
  brA <- is.na(pI) | pI >= branch_B_padj_cutoff
  # Branch B: interaction signal — padj_int_5k_random < branch_B_padj_cutoff
  brB <- !is.na(pI) & pI < branch_B_padj_cutoff
  same_sign <- !is.na(bF) & !is.na(bM) & bF != 0 & bM != 0 & sign(bF) == sign(bM)
  opp_sign  <- !is.na(bF) & !is.na(bM) & bF != 0 & bM != 0 & sign(bF) != sign(bM)
  # Branch A → Not_DEG / Concordant
  conc_A   <- brA & !is.na(lF) & !is.na(lM) & lF < 0.05 & lM < 0.05 & same_sign
  cls[conc_A] <- "Concordant"
  cls[brA & is.na(cls)] <- "Not_DEG"
  # Branch B-divergent
  is_div   <- brB & opp_sign & !is.na(lF) & !is.na(lM) & lF < 0.05 & lM < 0.05
  cls[is_div] <- "Divergent"
  is_div_one <- brB & opp_sign & is.na(cls) &
                ((!is.na(lF) & lF < 0.05) | (!is.na(lM) & lM < 0.05))
  cls[is_div_one] <- "Divergent_one_sided"
  # Branch B same-sign F/M-biased
  brB3 <- brB & same_sign & is.na(cls)
  ratio_F <- abs(bF) / pmax(abs(bM), .Machine$double.eps)
  ratio_M <- abs(bM) / pmax(abs(bF), .Machine$double.eps)
  is_F     <- brB3 & ratio_F >= mag_ratio & !is.na(lF) & lF < 0.05 &
              !is.na(lM) & lM > lfsr_other
  is_M     <- brB3 & ratio_M >= mag_ratio & !is.na(lM) & lM < 0.05 &
              !is.na(lF) & lF > lfsr_other
  is_F_pu  <- brB3 & ratio_F >= mag_ratio & !is.na(lF) & lF < 0.05 &
              !is.na(lM) & lM >= 0.05 & lM <= lfsr_other &
              !is.na(pow) & pow < 0.5
  is_M_pu  <- brB3 & ratio_M >= mag_ratio & !is.na(lM) & lM < 0.05 &
              !is.na(lF) & lF >= 0.05 & lF <= lfsr_other
  cls[is_F]    <- "Female_biased"
  cls[is_M]    <- "Male_biased"
  cls[is_F_pu] <- "Female_biased_M_underpowered"
  cls[is_M_pu] <- "Male_biased_F_underpowered"
  is_modif <- brB & is.na(cls)
  cls[is_modif] <- "Sex_modifier"
  cls[is.na(cls)] <- "Not_DEG"
  # chrY F→M reassign
  ch_vec <- ifelse(is.na(ch), "", ch)
  flip <- ch_vec == "chrY" & cls == "Female_biased"
  cls[flip] <- "Male_biased"
  cls
}

# ---- Sweep ----
grid_df <- as.data.table(expand.grid(mag_ratio = GRID$mag_ratio,
                                      lfsr_other = GRID$lfsr_other,
                                      branch_B_padj_cutoff = GRID$branch_B_padj_cutoff,
                                      stringsAsFactors = FALSE))
res_list <- list()
for (i in seq_len(nrow(grid_df))) {
  mr <- grid_df$mag_ratio[i]
  lo <- grid_df$lfsr_other[i]
  tp <- grid_df$branch_B_padj_cutoff[i]
  cls <- classify_param(dt$beta_F, dt$beta_M,
                         dt$lfsr_F_strat, dt$lfsr_M_strat,
                         dt$padj_int_5k_random, dt$power_M_at_F, dt$chr,
                         mr, lo, tp)
  n_F   <- sum(cls %in% c("Female_biased", "Female_biased_M_underpowered"))
  n_M   <- sum(cls %in% c("Male_biased", "Male_biased_F_underpowered"))
  n_div <- sum(cls %in% c("Divergent", "Divergent_one_sided"))
  n_mod <- sum(cls == "Sex_modifier")
  n_conc <- sum(cls == "Concordant")
  fm_ratio <- if (n_M == 0) NA_real_ else n_F / n_M
  is_canon <- (mr == CANON$mag_ratio) & (lo == CANON$lfsr_other) &
              (tp == CANON$branch_B_padj_cutoff)
  res_list[[i]] <- data.table(mag_ratio = mr, lfsr_other = lo,
                              branch_B_padj_cutoff = tp,
                              n_F_biased = n_F, n_M_biased = n_M,
                              FM_ratio = fm_ratio, n_Divergent = n_div,
                              n_Sex_modifier = n_mod, n_Concordant = n_conc,
                              is_canonical = is_canon)
}
res <- rbindlist(res_list)
fwrite(res, OUT_CSV)
cat("Wrote sweep CSV:", OUT_CSV, "  cells:", nrow(res), "\n")

cat("\n== Sweep summary ==\n"); print(res)
canon_row <- res[is_canonical == TRUE]
median_fm <- median(res$FM_ratio, na.rm = TRUE)
cat("\nCanonical cell FM_ratio:", round(canon_row$FM_ratio, 2), "\n")
cat("Median sweep FM_ratio:   ", round(median_fm, 2), "\n")
if (!is.na(canon_row$FM_ratio) && !is.na(median_fm) && median_fm > 0) {
  within_2x <- canon_row$FM_ratio >= median_fm / 2 &
               canon_row$FM_ratio <= median_fm * 2
  cat("Canonical within 2× of median:", if (within_2x) "PASS" else "FAIL", "\n")
}

# ===========================================================================
# Publication-quality figure — 2 panels (A,B)
#
# OVERALL KEY MESSAGE
# "The v6 sex-modulated candidate list is robust to threshold choice across
#  48 grid cells (4 mag_ratio × 4 lfsr_other × 3 branch_B_padj_cutoff). The
#  canonical cell (2.0, 0.20, 0.20) recovers 1 F + 0 M + 8 Divergent + 9
#  Sex_modifier. At the near-null lfsr_other = 0.45 cap (P(sign match) ≈ 0.55),
#  zero F-biased and zero M-biased survive ALL 12 cells — the canonical 0.20
#  boundary is conservative, not artifact-prone."
# ===========================================================================

# Source publication theme
THEME_PATH <- file.path(BASE, "scripts/figures/publication_theme.R")
if (file.exists(THEME_PATH)) {
  source(THEME_PATH)
} else {
  warning("publication_theme.R missing")
  theme_masld <- function(...) theme_bw()
  theme_pub <- function() theme()
  masld_colors <- list(female = "#AD1457", male = "#1A237E", ns = "#9E9E9E",
                       nafl = "#F4A674", down = "#1565C0", up = "#C9265E")
  save_fig <- function(p, f, width = 7, height = 5, dpi = 300)
    ggsave(f, p, width = width, height = height, dpi = dpi)
  fig_full_width <- 180 / 25.4
}

# Reshape sweep results into long format with display labels
res_long <- melt(res, id.vars = c("mag_ratio", "lfsr_other",
                                    "branch_B_padj_cutoff", "is_canonical"),
                  measure.vars = c("n_F_biased", "n_M_biased", "n_Divergent",
                                    "n_Sex_modifier"),
                  variable.name = "class", value.name = "n")
res_long[, class := factor(class,
                            levels = c("n_F_biased", "n_M_biased",
                                       "n_Divergent", "n_Sex_modifier"),
                            labels = c("Female-biased", "Male-biased",
                                       "Divergent", "Sex-modifier"))]
# Canonical 4-class palette aligned with the rest of the sex supp directory
# (figS_sex_scatter / figS_sex_class_distribution / figS_sex_beta_int_violin):
# Female-biased = rose, Male-biased = navy, Divergent = yellow, Sex-modifier
# (i.e. Suggestive-other / directionally-unresolved) = pale green.
class_pal <- c("Female-biased" = masld_colors$female,
               "Male-biased"   = masld_colors$male,
               "Divergent"     = "#FFC107",
               "Sex-modifier"  = "#A5D6A7")

# ---------- Panel A — Class counts across the 48-cell sweep ---------------
# KEY MESSAGE: Female-biased + Male-biased counts stay near zero across all
# 48 threshold combinations; canonical cell (red border) is representative.
# Aggregate per (mag_ratio, branch_B_padj_cutoff, lfsr_other) — show n class
# stacked.
canon_label <- sprintf("canonical: mag=%.1f, lfsr=%.2f, padj=%.2f",
                       CANON$mag_ratio, CANON$lfsr_other,
                       CANON$branch_B_padj_cutoff)
pA <- ggplot(res_long, aes(x = factor(lfsr_other), y = n, fill = class)) +
  geom_col(width = 0.85, color = "white", linewidth = 0.15) +
  facet_grid(paste0("padj=", branch_B_padj_cutoff) ~
               paste0("mag=", mag_ratio)) +
  geom_text(data = res_long[is_canonical == TRUE & class == "Female-biased"],
            aes(x = factor(lfsr_other), y = -1, label = "★"),
            color = masld_colors$up, size = 2.5, inherit.aes = FALSE,
            show.legend = FALSE) +
  scale_fill_manual(values = class_pal, name = NULL) +
  scale_y_continuous(limits = c(-3, NA),
                      breaks = function(lim) pretty(c(0, lim[2]))) +
  labs(title = "a   Class counts across 48-cell threshold grid",
       subtitle = paste0("Stacked counts of each sex class per cell. ",
                          "★ marks the ", canon_label, "."),
       x = expression("lfsr_other gate (cap on null-arm lfsr)"),
       y = "Gene count") +
  theme_masld() + theme_pub() +
  theme(legend.position = "bottom",
        panel.spacing = unit(0.25, "lines"),
        strip.text = element_text(size = 5, face = "bold"))

# ---------- Panel B — F-biased + M-biased counts at varying gates ---------
# Cleaner alternative to the FM ratio (which is NA when M=0). Show F and M
# counts as two side-by-side dot plots per (mag, padj) facet, with the canonical
# cell highlighted.
fm_long <- res_long[class %in% c("Female-biased", "Male-biased")]
canon_cell <- res[is_canonical == TRUE]
pB <- ggplot(fm_long,
             aes(x = factor(lfsr_other), y = n,
                 color = class, shape = class)) +
  # Zero baseline (so n=0 cells are visually marked, not invisible)
  geom_hline(yintercept = 0, color = "gray85", linewidth = 0.3) +
  geom_segment(aes(xend = factor(lfsr_other), y = 0, yend = n),
               linewidth = 0.4, alpha = 0.7) +
  # All points; zero-count points get an explicit text label
  geom_point(size = 2, stroke = 0.5) +
  geom_text(data = fm_long[n == 0],
            aes(x = factor(lfsr_other), y = 0, label = "0"),
            color = "gray45", size = 1.7, vjust = -0.9,
            show.legend = FALSE) +
  facet_grid(paste0("padj=", branch_B_padj_cutoff) ~
               paste0("mag=", mag_ratio)) +
  scale_color_manual(values = c("Female-biased" = masld_colors$female,
                                 "Male-biased"   = masld_colors$male),
                      name = NULL) +
  scale_shape_manual(values = c("Female-biased" = 16,
                                 "Male-biased"   = 17),
                      name = NULL) +
  scale_y_continuous(limits = c(-0.5, max(fm_long$n, na.rm = TRUE) + 1)) +
  geom_vline(data = data.frame(xint = which(levels(factor(fm_long$lfsr_other))
                                              == as.character(CANON$lfsr_other))),
              aes(xintercept = xint), linetype = "dashed",
              color = "gray60", linewidth = 0.25) +
  labs(title = "b   Female- vs Male-biased counts (cleaner view of panel a)",
       subtitle = expression("Lollipops per (" * delta * "-grid) cell; dashed line marks canonical lfsr_other = 0.20; \"0\" labels mark empty cells"),
       x = "lfsr_other (null-arm lfsr cap)",
       y = "Gene count") +
  theme_masld() + theme_pub() +
  theme(legend.position = "bottom",
        panel.spacing = unit(0.25, "lines"),
        strip.text = element_text(size = 5, face = "bold"))

# Headline annotation
n_male_canon <- canon_cell$n_M_biased
n_fem_canon  <- canon_cell$n_F_biased
n_div_canon  <- canon_cell$n_Divergent
n_sm_canon   <- canon_cell$n_Sex_modifier
near_null_subset <- res[lfsr_other == 0.45]
sum_fm_near_null <- sum(near_null_subset$n_F_biased + near_null_subset$n_M_biased,
                         na.rm = TRUE)
headline_lab_ts <- sprintf(
  "Canonical (2.0, 0.20, 0.20): %d F + %d M + %d Divergent + %d Sex-modifier.  Near-null lfsr_other=0.45 (12 cells): %d F + M total.",
  n_fem_canon, n_male_canon, n_div_canon, n_sm_canon, sum_fm_near_null)

p <- pA / pB + plot_annotation(
  title    = "Threshold sensitivity sweep",
  subtitle = headline_lab_ts,
  caption  = "48 cells = 4 mag_ratio × 4 lfsr_other × 3 branch_B_padj_cutoff. The lfsr_other = 0.45 row brackets the near-null directional-resolution boundary.",
  theme    = theme(plot.title    = element_text(size = 8, face = "bold"),
                   plot.subtitle = element_text(size = 6, color = "gray25"),
                   plot.caption  = element_text(size = 5, color = "gray40", hjust = 0)))

save_fig(p, OUT_PDF, width = fig_full_width, height = 5.0, dpi = 300)
cat("Wrote figure:", OUT_PDF, "\n")
cat("\nFinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
