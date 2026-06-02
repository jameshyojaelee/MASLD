##############################################################################
# Supplementary Figure: COLOC Reviewer Defense Analyses
#
# 6 panels (3 rows x 2 cols):
#   (a) Prior sensitivity — AFR genes (entirely prior-dependent)
#   (b) Prior sensitivity — EUR genes (robust vs fragile)
#   (c) Power analysis — minimum detectable effect size
#   (d) EAF concordance — EUR eQTL vs cross-ancestry GWAS
#   (e) EAF divergence — detailed cross-population metrics
#   (f) Anchor gene COLOC status tile/heatmap
#
# Demonstrates that AFR/CSA COLOC failures reflect LD/power issues,
# not biological absence.
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT       <- file.path(FIGS04_DIR, "figS_coloc_sensitivity.pdf")
PANEL_DIR <- file.path(FIGS04_DIR, "panels")
dir.create(FIGS04_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

# Pre-initialize all panels as placeholders
p_a <- placeholder("(a) Prior sensitivity: AFR")
p_b <- placeholder("(b) Prior sensitivity: EUR")
p_c <- placeholder("(c) Min detectable effect")
p_d <- placeholder("(d) EAF concordance")
p_e <- placeholder("(e) EAF divergence")
p_f <- placeholder("(f) Anchor gene COLOC")

# ═══════════════════════════════════════════════════════════════════════════
# Panel (a): Prior Sensitivity — AFR Genes
# ═══════════════════════════════════════════════════════════════════════════
tryCatch({
  message("Panel (a): Prior sensitivity — AFR genes...")
  prior_f <- file.path(CAUSAL, "prior_sensitivity", "prior_sensitivity_results.csv")
  if (!file.exists(prior_f)) stop("prior_sensitivity_results.csv not found")

  prior <- fread(prior_f)
  afr <- prior[source == "PanUKBB_AFR_ALT"]
  if (nrow(afr) == 0) stop("No AFR data in prior sensitivity results")

  afr <- add_symbols(afr, "gene")

  # Compute max PP.H4 per gene (for color gradient)
  afr[, max_pp4 := max(PP.H4, na.rm = TRUE), by = symbol]

  # Count how many are entirely prior-dependent (cross 0.5 threshold)
  gene_summary <- afr[, .(min_pp4 = min(PP.H4), max_pp4 = max(PP.H4),
                           pp4_default = PP.H4[which.min(abs(p12 - 5e-6))]),
                       by = symbol]
  n_fragile <- sum(gene_summary$min_pp4 < 0.9 & gene_summary$max_pp4 > 0.9)
  n_total   <- uniqueN(afr$symbol)

  p_a <- ggplot(afr, aes(x = p12, y = PP.H4, group = symbol, color = max_pp4)) +
    geom_line(linewidth = 0.3, alpha = 0.7) +
    geom_point(size = 0.5, alpha = 0.7) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red",
               linewidth = 0.3) +
    scale_x_log10(breaks = c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4),
                  labels = function(x) format(x, scientific = TRUE)) +
    scale_color_gradient(low = "#1565C0", high = "#C2185B",
                         name = "Max PP.H4") +
    annotate("text", x = 1e-6, y = 0.95,
             label = paste0(n_fragile, "/", n_total,
                            " genes entirely\nprior-dependent"),
             size = 2.2, hjust = 0, fontface = "bold", color = "gray30") +
    labs(x = expression(p[12]~"prior"), y = "PP.H4",
         title = "Prior sensitivity: AFR (N=6,636)") +
    theme_masld() +
    theme(legend.position = "right",
          legend.key.height = unit(0.3, "cm"),
          legend.key.width  = unit(0.15, "cm"),
          legend.text = element_text(size = 5),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 5))

  message("  AFR: ", n_fragile, "/", n_total, " genes prior-dependent")
}, error = function(e) message("Panel (a) error: ", e$message))

# ═══════════════════════════════════════════════════════════════════════════
# Panel (b): Prior Sensitivity — EUR Genes
# ═══════════════════════════════════════════════════════════════════════════
tryCatch({
  message("Panel (b): Prior sensitivity — EUR genes...")
  prior_f <- file.path(CAUSAL, "prior_sensitivity", "prior_sensitivity_results.csv")
  if (!file.exists(prior_f)) stop("prior_sensitivity_results.csv not found")

  prior <- fread(prior_f)
  eur <- prior[source == "UKBB_ALT"]
  if (nrow(eur) == 0) stop("No EUR data in prior sensitivity results")

  eur <- add_symbols(eur, "gene")

  # Get top 30 genes by PP.H4 at default p12=5e-6
  eur_default <- eur[, .(pp4_default = PP.H4[which.min(abs(p12 - 5e-6))]),
                     by = symbol]
  setorder(eur_default, -pp4_default)
  top30 <- head(eur_default, 30)$symbol

  eur_top <- eur[symbol %in% top30]

  # Classify: robust = PP.H4 > 0.9 at most stringent prior (p12=1e-6)
  eur_robust <- eur[p12 == min(p12), .(symbol, pp4_strict = PP.H4)]
  eur_top <- merge(eur_top, eur_robust, by = "symbol", all.x = TRUE)
  eur_top[, status := ifelse(pp4_strict > 0.9, "Robust", "Fragile")]

  # Count overall robust/total at PP.H4 > 0.9 and default p12
  n_default <- nrow(eur_default[pp4_default > 0.9])
  n_robust  <- sum(eur_robust$pp4_strict > 0.9, na.rm = TRUE)

  p_b <- ggplot(eur_top, aes(x = p12, y = PP.H4, group = symbol, color = status)) +
    geom_line(linewidth = 0.3, alpha = 0.6) +
    geom_point(size = 0.5, alpha = 0.6) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red",
               linewidth = 0.3) +
    scale_x_log10(breaks = c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4),
                  labels = function(x) format(x, scientific = TRUE)) +
    scale_color_manual(values = c("Robust" = masld_colors$conserved,
                                  "Fragile" = masld_colors$discordant),
                       name = NULL) +
    annotate("text", x = 1e-6, y = 0.95,
             label = paste0(n_robust, "/", n_default,
                            " genes robust\nacross all priors"),
             size = 2.2, hjust = 0, fontface = "bold", color = "gray30") +
    labs(x = expression(p[12]~"prior"), y = "PP.H4",
         title = paste0("Prior sensitivity: EUR (N=361,194)")) +
    theme_masld() +
    theme(legend.position = "right",
          legend.key.height = unit(0.2, "cm"),
          legend.text = element_text(size = 5),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 5))

  message("  EUR: ", n_robust, "/", n_default, " robust at p12=1e-6")
}, error = function(e) message("Panel (b) error: ", e$message))

# ═══════════════════════════════════════════════════════════════════════════
# Panel (c): Power Analysis — Minimum Detectable Effect Size
# ═══════════════════════════════════════════════════════════════════════════
tryCatch({
  message("Panel (c): Power analysis — minimum detectable effect...")
  power_f <- file.path(CAUSAL, "power_analysis", "power_curves.csv")
  if (!file.exists(power_f)) stop("power_curves.csv not found")

  power <- fread(power_f)

  # Map source names to display labels and colors
  source_map <- data.table(
    source = c("UKBB EUR", "BBJ", "PanUKBB CSA", "PanUKBB AFR"),
    label  = c("UKBB EUR (N=361,194)", "BBJ (N=165,084)",
               "PanUKBB CSA (N=8,876)", "PanUKBB AFR (N=6,636)"),
    color  = c("#1565C0", "#7B1FA2", "#E91E63", "#FF8F00")
  )
  power <- merge(power, source_map, by = "source", all.x = TRUE)
  power <- power[!is.na(label)]

  # Order labels by N (largest first)
  power[, label := factor(label, levels = source_map$label)]

  # Annotate values at MAF=0.05
  maf05 <- power[abs(MAF - 0.05) < 0.001]
  maf05_labels <- maf05[, .(MAF = MAF[1], min_beta_80pct = min_beta_80pct[1],
                             label_text = sprintf("%.2f", min_beta_80pct[1])),
                         by = label]

  source_colors <- setNames(source_map$color, source_map$label)

  p_c <- ggplot(power, aes(x = MAF, y = min_beta_80pct,
                            color = label, group = label)) +
    geom_line(linewidth = 0.5) +
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "gray50",
               linewidth = 0.3) +
    scale_color_manual(values = source_colors, name = NULL) +
    scale_x_continuous(breaks = seq(0, 0.5, by = 0.1)) +
    labs(x = "Minor allele frequency (MAF)",
         y = expression("Min |"*beta*"| at 80% power"),
         title = "Minimum detectable effect size (80% power)") +
    theme_masld() +
    theme(legend.position = "right",
          legend.key.height = unit(0.2, "cm"),
          legend.text = element_text(size = 5))

  # Add MAF=0.05 annotations if data exists
  if (nrow(maf05_labels) > 0) {
    p_c <- p_c +
      geom_point(data = maf05, aes(x = MAF, y = min_beta_80pct),
                 size = 1.5, show.legend = FALSE) +
      geom_text_repel(data = maf05_labels,
                      aes(x = MAF, y = min_beta_80pct, label = label_text),
                      size = 2, nudge_x = 0.03, show.legend = FALSE,
                      min.segment.length = 0, segment.size = 0.15)
  }

  message("  Power curves: ", uniqueN(power$source), " sources plotted")
}, error = function(e) message("Panel (c) error: ", e$message))

# ═══════════════════════════════════════════════════════════════════════════
# Panel (d): EAF Concordance — EUR eQTL vs Cross-Ancestry GWAS
# ═══════════════════════════════════════════════════════════════════════════
tryCatch({
  message("Panel (d): EAF concordance — bar chart...")
  eaf_f <- file.path(CAUSAL, "eaf_concordance", "eaf_concordance_summary.csv")
  if (!file.exists(eaf_f)) stop("eaf_concordance_summary.csv not found")

  eaf <- fread(eaf_f)

  # Clean ancestry labels
  eaf[, ancestry_label := fifelse(ancestry == "AFR", "AFR (N=6,636)",
                           fifelse(ancestry == "CSA", "CSA (N=8,876)",
                                   ancestry))]

  ancestry_colors <- c("AFR (N=6,636)" = "#FF8F00", "CSA (N=8,876)" = "#E91E63")

  # Show Pearson r as grouped bars with reference line at r=1 (same ancestry)
  eaf_r <- eaf[, .(ancestry_label, pearson_r, n_variants)]

  p_d <- ggplot(eaf_r, aes(x = ancestry_label, y = pearson_r,
                            fill = ancestry_label)) +
    geom_col(width = 0.6, show.legend = FALSE) +
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray50",
               linewidth = 0.3) +
    geom_text(aes(label = sprintf("r = %.3f\nn = %s", pearson_r,
                                  format(n_variants, big.mark = ","))),
              vjust = -0.3, size = 2.2) +
    scale_fill_manual(values = ancestry_colors) +
    scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25)) +
    annotate("text", x = 1.5, y = 1.05,
             label = "Expected (same ancestry)", size = 2,
             color = "gray50", fontface = "italic") +
    labs(x = NULL, y = "Pearson r (EUR eQTL EAF vs GWAS EAF)",
         title = "EAF concordance: EUR eQTL vs GWAS") +
    theme_masld()

  message("  EAF concordance: AFR r=", eaf[ancestry == "AFR", pearson_r],
          ", CSA r=", eaf[ancestry == "CSA", pearson_r])
}, error = function(e) message("Panel (d) error: ", e$message))

# ═══════════════════════════════════════════════════════════════════════════
# Panel (e): EAF Divergence — Detailed Multi-Metric Comparison
# ═══════════════════════════════════════════════════════════════════════════
tryCatch({
  message("Panel (e): EAF divergence — multi-metric...")
  eaf_f <- file.path(CAUSAL, "eaf_concordance", "eaf_concordance_summary.csv")
  if (!file.exists(eaf_f)) stop("eaf_concordance_summary.csv not found")

  eaf <- fread(eaf_f)

  # Reshape to long for faceted comparison
  metrics <- rbindlist(list(
    eaf[, .(ancestry, metric = "Pearson r", value = pearson_r)],
    eaf[, .(ancestry, metric = "Median |EAF diff|", value = median_eaf_diff)],
    eaf[, .(ancestry, metric = "% variants |diff| > 0.3", value = pct_large_diff_0.3)]
  ))

  # Ensure consistent ordering
  metrics[, metric := factor(metric,
    levels = c("Pearson r", "Median |EAF diff|", "% variants |diff| > 0.3"))]
  metrics[, ancestry_label := fifelse(ancestry == "AFR", "AFR", "CSA")]

  ancestry_colors_e <- c("AFR" = "#FF8F00", "CSA" = "#E91E63")

  p_e <- ggplot(metrics, aes(x = ancestry_label, y = value,
                              fill = ancestry_label)) +
    geom_col(width = 0.6, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.2f", value)), vjust = -0.3, size = 2.2) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = ancestry_colors_e) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(x = NULL, y = NULL,
         title = "Cross-population allele frequency divergence") +
    theme_masld() +
    theme(strip.text = element_text(size = 6))

  message("  EAF metrics plotted for AFR and CSA")
}, error = function(e) message("Panel (e) error: ", e$message))

# ═══════════════════════════════════════════════════════════════════════════
# Panel (f): Anchor Gene COLOC Status — Tile Heatmap
# ═══════════════════════════════════════════════════════════════════════════
tryCatch({
  message("Panel (f): Anchor gene COLOC status tile...")

  # 10 canonical MASLD risk genes
  anchor_genes <- c("PNPLA3", "TM6SF2", "HSD17B13", "MBOAT7", "GCKR",
                     "RORA", "DGAT2", "THRB", "SLC39A8", "SORT1")

  # 4 Broadaway European GWAS
  gwas_dirs <- c(
    "UKBB ALT" = "broadaway_ukbb",
    "UKBB AST" = "broadaway_ukbb_ast",
    "UKBB GGT" = "broadaway_ukbb_ggt",
    "PDFF"     = "broadaway_pdff"
  )

  # Load COLOC results and extract PP.H4 for each anchor gene
  tile_data <- rbindlist(lapply(names(gwas_dirs), function(gwas_label) {
    f <- file.path(CAUSAL, gwas_dirs[[gwas_label]], "coloc_results.csv")
    if (!file.exists(f)) return(NULL)

    dt <- fread(f)
    dt <- add_symbols(dt, "gene")

    # Get best PP.H4 per symbol
    dt_best <- dt[, .(PP.H4 = max(PP.H4, na.rm = TRUE)), by = symbol]

    rbindlist(lapply(anchor_genes, function(g) {
      row <- dt_best[symbol == g]
      if (nrow(row) == 0) {
        data.table(gene = g, gwas = gwas_label, PP.H4 = NA_real_,
                   status = "Untested")
      } else {
        pp4 <- row$PP.H4[1]
        data.table(gene = g, gwas = gwas_label, PP.H4 = pp4,
                   status = fifelse(pp4 > 0.9, "Colocalized", "Not colocalized"))
      }
    }))
  }))

  if (nrow(tile_data) == 0) stop("No COLOC data loaded for anchor genes")

  # Label for tiles
  tile_data[, pp4_label := fifelse(is.na(PP.H4), "N/A",
                                   fifelse(PP.H4 < 0.001,
                                           sprintf("%.0e", PP.H4),
                                           sprintf("%.3f", PP.H4)))]

  # Order genes: successes first, then failures, then untested
  gene_max <- tile_data[, .(max_pp4 = max(PP.H4, na.rm = TRUE),
                             any_coloc = any(status == "Colocalized")),
                         by = gene]
  gene_max[is.na(any_coloc), any_coloc := FALSE]
  # Separate: colocalized (desc by max PP.H4), not colocalized, untested
  gene_order <- c(
    gene_max[any_coloc == TRUE][order(-max_pp4), gene],
    gene_max[any_coloc == FALSE & is.finite(max_pp4)][order(-max_pp4), gene],
    gene_max[!is.finite(max_pp4), gene]
  )
  tile_data[, gene := factor(gene, levels = rev(gene_order))]
  tile_data[, gwas := factor(gwas, levels = names(gwas_dirs))]

  # Count successes/failures at the gene level
  n_success  <- sum(gene_max$any_coloc == TRUE)
  n_untested <- sum(!is.finite(gene_max$max_pp4))
  n_fail     <- length(anchor_genes) - n_success - n_untested

  status_colors <- c(
    "Colocalized"      = masld_colors$conserved,  # "#00695C"
    "Not colocalized"  = masld_colors$discordant,       # "#E91E63"
    "Untested"         = masld_colors$ns                 # "#BDBDBD"
  )

  p_f <- ggplot(tile_data, aes(x = gwas, y = gene, fill = status)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = pp4_label), size = 1.8, color = "white",
              fontface = "bold") +
    scale_fill_manual(values = status_colors, name = NULL) +
    labs(x = NULL, y = NULL,
         title = "Anchor gene COLOC status",
         subtitle = paste0(n_fail, "/", length(anchor_genes),
                           " fail, ", n_untested, " untested \u2014 coding variant architecture")) +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          axis.text.y = element_text(face = "italic", size = 6),
          legend.position = "bottom",
          legend.key.size = unit(0.25, "cm"),
          plot.subtitle = element_text(size = 6, color = "gray40"))

  message("  Anchor genes: ", n_success, " colocalized, ", n_fail,
          " fail, ", n_untested, " untested across ", length(gwas_dirs), " GWAS")
}, error = function(e) message("Panel (f) error: ", e$message))

# ═══════════════════════════════════════════════════════════════════════════
# Save individual panels
# ═══════════════════════════════════════════════════════════════════════════
message("Saving individual panels...")
panel_list <- list(
  a_prior_afr         = p_a,
  b_prior_eur         = p_b,
  c_power_analysis    = p_c,
  d_eaf_concordance   = p_d,
  e_eaf_divergence    = p_e,
  f_anchor_gene_coloc = p_f
)

for (nm in names(panel_list)) {
  panel_path <- file.path(PANEL_DIR, paste0("panel_", nm, ".pdf"))
  save_fig(panel_list[[nm]], panel_path, width = fig_half_width, height = 2.8)
  message("  Saved: ", panel_path)
}

# ═══════════════════════════════════════════════════════════════════════════
# Assemble composite figure (3 rows x 2 cols)
# ═══════════════════════════════════════════════════════════════════════════
message("Assembling 3x2 composite figure...")
fig <- (p_a | p_b) / (p_c | p_d) / (p_e | p_f) +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_fig_tall(fig, OUT, height = 10)
message("Saved composite: ", OUT)

message("=== figS_coloc_sensitivity.R complete ===")
