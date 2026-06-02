#!/usr/bin/env Rscript
# NOTE (Phase 10 swap, 2026-05-06): This script intentionally loads the FROZEN
# UKBB v1 sghatan snapshot as a sensitivity-comparison panel. v1 is no longer
# the production EUR LD reference — PolyFun replaced it (concordance r=0.986).
# v1 archive: archive/ld_panel_v1_sghatan_2026-05-06/
#
# figS_polyfun_fm_concordance.R
# Compares per-variant SuSiE PIP across LD panels for the standalone
# fine-mapping (Phase 9e):
#
#   panels:  EUR_0.5Mb (sghatan v1), 1kg_eur_0.5Mb, polyfun_eur_0.5Mb
#   metric:  per-variant PIP from each .tsv produced by 03_run_fm_per_locus.R
#   loci:    140 PolyFun TSVs (one per per-GWAS lead SNP)
#
# Output panels:
#   A: pairwise PIP scatter (PolyFun vs sghatan v1 + 1KG vs sghatan), all variants
#   B: per-locus Pearson r distribution (histogram), faceted by convergence
#   C: top-PIP variant per locus (PolyFun vs sghatan v1 PIP scatter)
#   D: agreement breakdown by convergence status (stacked bar)
#   E: example "agree" vs "disagree" locus PIP traces
#
# Fig: figures/supplementary/figS_polyfun_fm_concordance.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(ggrastr)
})

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(PROJ, "scripts/figures/publication_theme.R"))
source(file.path(PROJ, "scripts/figures/load_figure_data.R"))
FM_DIR <- file.path(PROJ, "GWAS/finemapping")
# Goes alongside the other 4-way / 3-way EUR panel comparisons.
OUT_DIR <- file.path(PROJ, "figures/supplementary/figS04_coloc/ld_panel_4way_eur/fm")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PANEL_COLORS <- c("UKBB v1 (sghatan)" = "#1F77B4", "1KG EUR" = "#2CA02C",
                  "TOP-LD EUR" = "#D62728", "PolyFun" = "#9467BD")

# ----------------------------------------------------------------------------
# 1. Load all PolyFun TSVs and find matching sghatan v1 + 1KG
# ----------------------------------------------------------------------------
load_panel_tsv <- function(path) {
  if (!file.exists(path)) return(NULL)
  d <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(d) || !"PIP" %in% names(d) || !"position" %in% names(d)) return(NULL)
  d[, .(position, PIP, converged = !grepl("notconverged", basename(path)))]
}

cat("Scanning PolyFun TSVs...\n")
pf_files <- Sys.glob(file.path(FM_DIR, "output/*/polyfun_eur_0.5Mb/susie/*.tsv"))
cat(sprintf("  %d PolyFun FM TSVs\n", length(pf_files)))

per_locus <- list()
all_var <- list()
for (pf in pf_files) {
  parts <- strsplit(pf, "/")[[1]]
  study <- parts[length(parts) - 3]
  bn <- sub("_cov.*", "", basename(pf))

  v1_glob <- Sys.glob(file.path(FM_DIR, "output", study, "EUR_0.5Mb/susie",
                                 paste0(bn, "_cov*.tsv")))
  kg_glob <- Sys.glob(file.path(FM_DIR, "output", study, "1kg_eur_0.5Mb/susie",
                                 paste0(bn, "_cov*.tsv")))

  d_pf <- load_panel_tsv(pf)
  d_v1 <- if (length(v1_glob) > 0) load_panel_tsv(v1_glob[1]) else NULL
  d_kg <- if (length(kg_glob) > 0) load_panel_tsv(kg_glob[1]) else NULL
  if (is.null(d_pf)) next

  conv_pf <- d_pf$converged[1]
  conv_v1 <- if (!is.null(d_v1)) d_v1$converged[1] else NA
  conv_kg <- if (!is.null(d_kg)) d_kg$converged[1] else NA

  # Per-locus Pearson r (PolyFun vs each)
  r_v1 <- NA_real_
  r_kg <- NA_real_
  if (!is.null(d_v1)) {
    m <- merge(d_pf[, .(position, PIP_pf = PIP)], d_v1[, .(position, PIP_v1 = PIP)], by = "position")
    if (nrow(m) > 50 && sd(m$PIP_pf) > 0 && sd(m$PIP_v1) > 0) r_v1 <- cor(m$PIP_pf, m$PIP_v1)
  }
  if (!is.null(d_kg)) {
    m2 <- merge(d_pf[, .(position, PIP_pf = PIP)], d_kg[, .(position, PIP_kg = PIP)], by = "position")
    if (nrow(m2) > 50 && sd(m2$PIP_pf) > 0 && sd(m2$PIP_kg) > 0) r_kg <- cor(m2$PIP_pf, m2$PIP_kg)
  }

  per_locus[[length(per_locus) + 1]] <- data.table(
    study = study, locus = bn,
    n_pf = nrow(d_pf),
    n_v1 = if (!is.null(d_v1)) nrow(d_v1) else NA_integer_,
    n_kg = if (!is.null(d_kg)) nrow(d_kg) else NA_integer_,
    r_pf_v1 = r_v1, r_pf_kg = r_kg,
    conv_pf = conv_pf, conv_v1 = conv_v1, conv_kg = conv_kg,
    top_pf_pip = max(d_pf$PIP, na.rm = TRUE),
    top_v1_pip = if (!is.null(d_v1)) max(d_v1$PIP, na.rm = TRUE) else NA_real_,
    top_kg_pip = if (!is.null(d_kg)) max(d_kg$PIP, na.rm = TRUE) else NA_real_
  )

  # Per-variant comparison rows (all panels matched)
  if (!is.null(d_v1)) {
    m <- merge(d_pf[, .(position, PIP_pf = PIP)], d_v1[, .(position, PIP_v1 = PIP)], by = "position")
    if (nrow(m) > 0) {
      m[, study := study]; m[, locus := bn]
      m[, conv_pf := conv_pf]; m[, conv_v1 := conv_v1]
      all_var[[length(all_var) + 1]] <- m
    }
  }
}

dt <- rbindlist(per_locus, fill = TRUE)
var_dt <- rbindlist(all_var, fill = TRUE)
cat(sprintf("  Per-locus rows: %d\n", nrow(dt)))
cat(sprintf("  Per-variant rows (PolyFun vs v1 paired): %d\n", nrow(var_dt)))

# Drop loci where r is NA
dt_v1 <- dt[!is.na(r_pf_v1)]
dt_kg <- dt[!is.na(r_pf_kg)]
cat(sprintf("  Loci with v1 comparison: %d\n", nrow(dt_v1)))
cat(sprintf("  Loci with 1KG comparison: %d\n", nrow(dt_kg)))

# ----------------------------------------------------------------------------
# Panel A: per-variant PIP scatter (PolyFun vs sghatan v1)
# ----------------------------------------------------------------------------
var_dt[, conv_status := ifelse(conv_pf & conv_v1, "Both converged",
                       ifelse(!conv_pf & !conv_v1, "Both not-converged",
                              "One not-converged"))]
var_dt[, conv_status := factor(conv_status,
        levels = c("Both converged", "One not-converged", "Both not-converged"))]

p_A <- ggplot(var_dt, aes(x = PIP_v1, y = PIP_pf, color = conv_status)) +
  geom_abline(slope = 1, intercept = 0, color = "grey40", linetype = "dashed") +
  rasterise(geom_point(alpha = 0.25, size = 0.6), dpi = 300) +
  scale_color_manual(values = c("Both converged" = "#2C7FB8",
                                 "One not-converged" = "#FFA500",
                                 "Both not-converged" = "#D62728")) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "PIP (sghatan UKBB v1)", y = "PIP (PolyFun)",
       title = sprintf("PolyFun vs sghatan v1 (n=%s vars)",
                       format(nrow(var_dt), big.mark=",")),
       color = NULL) +
  theme_masld(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

# ----------------------------------------------------------------------------
# Panel B: per-locus Pearson r distribution
# ----------------------------------------------------------------------------
dt_v1[, conv_status := ifelse(conv_pf & conv_v1, "Both converged",
                      ifelse(!conv_pf & !conv_v1, "Both not-converged",
                             "One not-converged"))]
dt_v1[, conv_status := factor(conv_status,
       levels = c("Both converged", "One not-converged", "Both not-converged"))]

p_B <- ggplot(dt_v1, aes(x = r_pf_v1, fill = conv_status)) +
  geom_histogram(binwidth = 0.05, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = c("Both converged" = "#2C7FB8",
                                "One not-converged" = "#FFA500",
                                "Both not-converged" = "#D62728")) +
  geom_vline(xintercept = 0.95, color = "darkgreen", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = 0.93, y = max(table(cut(dt_v1$r_pf_v1, seq(-0.05, 1.05, 0.05))))*0.9,
           label = "r = 0.95", angle = 90, hjust = 1, vjust = -0.5, size = 3, color = "darkgreen") +
  labs(x = "Per-locus Pearson r", y = "# loci",
       title = sprintf("Per-locus PIP r (n=%d)", nrow(dt_v1)),
       fill = NULL) +
  theme_masld(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

# ----------------------------------------------------------------------------
# Panel C: top-PIP per locus scatter
# ----------------------------------------------------------------------------
p_C <- ggplot(dt_v1, aes(x = top_v1_pip, y = top_pf_pip, color = conv_status)) +
  geom_abline(slope = 1, intercept = 0, color = "grey40", linetype = "dashed") +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Both converged" = "#2C7FB8",
                                 "One not-converged" = "#FFA500",
                                 "Both not-converged" = "#D62728")) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Max PIP (sghatan v1)", y = "Max PIP (PolyFun)",
       title = "Top-PIP per locus", color = NULL) +
  theme_masld(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

# ----------------------------------------------------------------------------
# Panel D: convergence-stratified r summary
# ----------------------------------------------------------------------------
summary_dt <- dt_v1[, .(
  n = .N,
  median_r = median(r_pf_v1),
  pct_r_above_95 = 100 * mean(r_pf_v1 > 0.95),
  pct_r_above_50 = 100 * mean(r_pf_v1 > 0.50)
), by = conv_status][order(conv_status)]

cat("\n## Per-convergence-class summary\n")
print(summary_dt)

p_D <- ggplot(dt_v1, aes(x = conv_status, y = r_pf_v1, fill = conv_status)) +
  geom_boxplot(outlier.alpha = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 0.8, color = "grey30") +
  scale_fill_manual(values = c("Both converged" = "#2C7FB8",
                                "One not-converged" = "#FFA500",
                                "Both not-converged" = "#D62728"),
                    guide = "none") +
  geom_hline(yintercept = 0.95, color = "darkgreen", linetype = "dashed") +
  labs(x = NULL, y = "Per-locus PIP r",
       title = "r vs convergence status") +
  theme_masld(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 20, hjust = 1))

# ----------------------------------------------------------------------------
# Panel E: example loci — one good (PNPLA3), one mediocre, one bad
# ----------------------------------------------------------------------------
example_loci <- c(
  "good"  = "2019_31311600_NAFLD_EUR/UKBB_22.44324727",
  "mid"   = NA,  # filled below
  "bad"   = NA
)
mid_locus <- dt_v1[r_pf_v1 > 0.7 & r_pf_v1 < 0.9 & conv_pf & conv_v1][1]
bad_locus <- dt_v1[r_pf_v1 < 0.5 & !conv_pf | !conv_v1][order(r_pf_v1)][1]
if (nrow(mid_locus) > 0) example_loci["mid"] <- paste(mid_locus$study, mid_locus$locus, sep = "/")
if (nrow(bad_locus) > 0) example_loci["bad"] <- paste(bad_locus$study, bad_locus$locus, sep = "/")

ex_data <- list()
for (label in names(example_loci)) {
  if (is.na(example_loci[label])) next
  parts <- strsplit(example_loci[label], "/")[[1]]
  study <- parts[1]; locus <- parts[2]
  pf <- Sys.glob(file.path(FM_DIR, "output", study, "polyfun_eur_0.5Mb/susie",
                           paste0(locus, "_cov*.tsv")))[1]
  v1 <- Sys.glob(file.path(FM_DIR, "output", study, "EUR_0.5Mb/susie",
                           paste0(locus, "_cov*.tsv")))[1]
  d_pf <- load_panel_tsv(pf); d_v1 <- load_panel_tsv(v1)
  if (is.null(d_pf) || is.null(d_v1)) next
  m <- merge(d_pf[, .(position, PolyFun = PIP)], d_v1[, .(position, sghatan_v1 = PIP)], by = "position")
  m_long <- melt(m, id.vars = "position", variable.name = "panel", value.name = "PIP")
  m_long[, locus_label := sprintf("%s — %s\n%s",
    label, paste0(study, " ", locus),
    sprintf("r = %.3f", cor(m$PolyFun, m$sghatan_v1)))]
  ex_data[[label]] <- m_long
}
ex_dt <- rbindlist(ex_data)

p_E <- ggplot(ex_dt, aes(x = position, y = PIP, color = panel)) +
  rasterise(geom_point(alpha = 0.6, size = 1), dpi = 300) +
  scale_color_manual(values = c("PolyFun" = "#9467BD", "sghatan_v1" = "#1F77B4")) +
  facet_wrap(~ locus_label, ncol = 1, scales = "free") +
  labs(x = "Genomic position", y = "PIP", title = "Example loci") +
  theme_masld(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(size = 10, face = "bold"),
        legend.position = "bottom")

# ----------------------------------------------------------------------------
# Compose figure
# ----------------------------------------------------------------------------
top_row <- p_A + p_B + plot_layout(widths = c(1, 1.1))
mid_row <- p_C + p_D + plot_layout(widths = c(1, 1))
final <- top_row / mid_row / p_E +
  plot_layout(heights = c(1, 1, 1.5)) +
  plot_annotation(
    title = "PolyFun vs sghatan v1 — fine-mapping concordance",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 16))
  )

OUTFILE <- file.path(OUT_DIR, "concordance_PolyFun_vs_sghatan_v1.pdf")
ggsave(OUTFILE, final, width = 15, height = 18, device = cairo_pdf, limitsize = FALSE)
cat(sprintf("\nWrote: %s\n", OUTFILE))
cat(sprintf("Size: %s\n", system(paste("du -h", shQuote(OUTFILE), "| cut -f1"), intern = TRUE)))

fwrite(dt, file.path(OUT_DIR, "concordance_per_locus.csv"))
cat(sprintf("Wrote: %s\n",
            file.path(OUT_DIR, "concordance_per_locus.csv")))
