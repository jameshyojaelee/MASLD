##############################################################################
# Supplementary Figure: SuSiE-COLOC vs ABF Comparison
#
# 3 panels:
#   (a) ABF vs SuSiE PP.H4 scatter (UKBB ALT, primary GWAS)
#   (b) Concordance bar chart across 6 EUR GWAS
#   (c) Multi-signal gene bars: PP.H4 per GWAS (genes with n_cs_pairs > 1)
#
# Data: RNA-seq/results/causal_inference/susie_comparison/
#   - susie_abf_comparison.csv  (per-gene, per-GWAS)
#   - susie_abf_summary.csv     (per-GWAS summary)
#   - susie_multi_signal_genes.csv (multi-signal genes)
##############################################################################

set.seed(42)

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

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SUSIE_DIR <- file.path(CAUSAL, "susie_comparison")
OUT       <- file.path(FIGS04_DIR, "figS_susie_comparison.pdf")
PANEL_DIR <- file.path(FIGS04_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

# Pre-initialize all panels as placeholders
p_a <- placeholder("(a) ABF vs SuSiE PP.H4 scatter")
p_b <- placeholder("(b) Concordance across EUR GWAS")
p_c <- placeholder("(c) Multi-signal genes")

# ═══════════════════════════════════════════════════════════════════════════
# Panel (a): ABF vs SuSiE PP.H4 Scatter — UKBB ALT
# ═══════════════════════════════════════════════════════════════════════════
tryCatch({
  message("Panel (a): ABF vs SuSiE PP.H4 scatter (UKBB ALT)...")
  comp_f <- file.path(SUSIE_DIR, "susie_abf_comparison.csv")
  if (!file.exists(comp_f)) stop("susie_abf_comparison.csv not found")

  comp <- fread(comp_f)

  # Filter to UKBB ALT (primary, largest gene set)
  # Try common naming patterns for the GWAS column
  gwas_col <- intersect(names(comp), c("gwas", "GWAS", "source"))[1]
  if (is.na(gwas_col)) stop("No GWAS identifier column found")

  alt_patterns <- c("UKBB_ALT", "UKBB ALT", "ukbb_alt", "broadaway_ukbb")
  alt_hits <- comp[get(gwas_col) %in% alt_patterns]
  if (nrow(alt_hits) == 0) {
    # Fallback: take the GWAS with the most genes
    gwas_counts <- comp[, .N, by = gwas_col]
    setorder(gwas_counts, -N)
    primary_gwas <- gwas_counts[[gwas_col]][1]
    alt_hits <- comp[get(gwas_col) == primary_gwas]
    message("  Using GWAS with most genes: ", primary_gwas, " (n=", nrow(alt_hits), ")")
  }

  dt <- alt_hits
  if (nrow(dt) == 0) stop("No data for primary GWAS")

  # Ensure required columns
  stopifnot(all(c("PP.H4.abf", "PP.H4.susie") %in% names(dt)))

  # Add gene symbols if gene column contains Ensembl IDs
  gene_col_name <- intersect(names(dt), c("gene", "gene_id", "ensembl_id"))[1]
  if (!is.na(gene_col_name) && any(grepl("^ENSG", dt[[gene_col_name]]))) {
    dt <- add_symbols(dt, gene_col_name)
  } else if ("symbol" %in% names(dt)) {
    # Already has symbol
  } else if (!is.na(gene_col_name)) {
    dt[, symbol := get(gene_col_name)]
  }

  # Ensure n_cs_pairs column exists (default to 1 if missing)
  if (!"n_cs_pairs" %in% names(dt)) dt[, n_cs_pairs := 1L]
  dt[is.na(n_cs_pairs), n_cs_pairs := 1L]

  # Compute concordance: both agree on threshold direction
  threshold <- 0.5
  dt[, abf_sig   := PP.H4.abf   > threshold]
  dt[, susie_sig := PP.H4.susie > threshold]
  dt[, concordant := (abf_sig == susie_sig)]
  concordance_pct <- round(100 * mean(dt$concordant, na.rm = TRUE), 1)

  # Spearman correlation
  rho <- cor(dt$PP.H4.abf, dt$PP.H4.susie, method = "spearman", use = "complete.obs")
  rho_label <- sprintf("rho == %.3f", rho)
  conc_label <- paste0(concordance_pct, "% concordant")

  # Identify discordant genes for labeling (large PP.H4 difference)
  dt[, pp4_diff := abs(PP.H4.abf - PP.H4.susie)]
  dt[, discordant := (abf_sig != susie_sig) & pp4_diff > 0.3]
  # Top 10 most discordant
  label_dt <- dt[discordant == TRUE][order(-pp4_diff)][1:min(.N, 10)]

  # Cap n_cs_pairs for color scale
  dt[, cs_pairs_capped := pmin(n_cs_pairs, 5)]

  p_a <- ggplot(dt, aes(x = PP.H4.abf, y = PP.H4.susie)) +
    rasterize_layer(
      geom_point(aes(color = cs_pairs_capped), size = 0.6, alpha = 0.6)
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "gray50", linewidth = 0.3) +
    scale_color_gradient(low = masld_colors$ns, high = masld_colors$up,
                         name = "CS pairs",
                         breaks = c(1, 2, 3, 4, 5),
                         labels = c("1", "2", "3", "4", "5+")) +
    annotate("text", x = 0.02, y = 0.95, label = rho_label,
             parse = TRUE, size = 2.2, hjust = 0, color = "gray20") +
    annotate("text", x = 0.02, y = 0.87, label = conc_label,
             size = 2.2, hjust = 0, color = "gray20") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "PP.H4 (ABF)", y = "PP.H4 (SuSiE)",
         title = "ABF vs SuSiE colocalization (UKBB ALT)") +
    theme_masld() +
    theme(legend.position = "right",
          legend.key.height = unit(0.3, "cm"),
          legend.key.width  = unit(0.15, "cm"),
          legend.text = element_text(size = 5))

  # Add labels for discordant genes

  if (nrow(label_dt) > 0 && "symbol" %in% names(label_dt)) {
    p_a <- p_a +
      geom_text_repel(data = label_dt,
                      aes(x = PP.H4.abf, y = PP.H4.susie, label = symbol),
                      size = 1.8, color = masld_colors$up,
                      max.overlaps = 15, min.segment.length = 0,
                      segment.size = 0.15, segment.color = "gray60",
                      fontface = "italic", seed = 42)
  }

  message("  n=", nrow(dt), " genes; Spearman rho=", round(rho, 3),
          "; concordance=", concordance_pct, "%")
}, error = function(e) message("Panel (a) error: ", e$message))

# ═══════════════════════════════════════════════════════════════════════════
# Panel (b): Concordance Bar Chart Across 6 EUR GWAS
# ═══════════════════════════════════════════════════════════════════════════
tryCatch({
  message("Panel (b): Concordance across EUR GWAS...")
  summ_f <- file.path(SUSIE_DIR, "susie_abf_summary.csv")
  if (!file.exists(summ_f)) stop("susie_abf_summary.csv not found")

  summ <- fread(summ_f)
  if (nrow(summ) == 0) stop("Summary table is empty")

  # Identify GWAS name column
  gwas_col_s <- intersect(names(summ), c("gwas", "GWAS", "source"))[1]
  if (is.na(gwas_col_s)) stop("No GWAS identifier column found in summary")

  # Try to find concordance and fallback columns
  conc_col <- intersect(names(summ), c("concordance_pct", "pct_concordant",
                                        "concordance", "concordant_pct"))[1]
  fall_col <- intersect(names(summ), c("fallback_rate", "susie_fallback_rate",
                                        "pct_fallback", "fallback_pct"))[1]

  if (is.na(conc_col)) stop("No concordance column found in summary")

  # Build long-form data for plotting
  bar_dt <- summ[, .(gwas = get(gwas_col_s), concordance = get(conc_col))]

  # Add fallback rate if available
  has_fallback <- !is.na(fall_col)
  if (has_fallback) {
    bar_dt[, fallback := summ[[fall_col]]]
  }

  # Shorten GWAS labels for display
  bar_dt[, gwas_short := gsub("broadaway_|UKBB_|UKBB |Broadaway_", "", gwas)]
  bar_dt[, gwas_short := gsub("ukbb_", "", gwas_short)]

  # Order by concordance descending
  setorder(bar_dt, -concordance)
  bar_dt[, gwas_short := factor(gwas_short, levels = rev(gwas_short))]

  # Melt for grouped bars if fallback available
  if (has_fallback) {
    bar_long <- melt(bar_dt, id.vars = "gwas_short",
                     measure.vars = c("concordance", "fallback"),
                     variable.name = "metric", value.name = "pct")
    bar_long[, metric_label := fifelse(metric == "concordance",
                                       "Concordant (PP.H4 > 0.5)",
                                       "SuSiE fallback to ABF")]

    bar_colors <- c("Concordant (PP.H4 > 0.5)" = masld_colors$down,
                     "SuSiE fallback to ABF"     = masld_colors$ns)

    p_b <- ggplot(bar_long, aes(x = gwas_short, y = pct, fill = metric_label)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6) +
      geom_text(aes(label = sprintf("%.0f%%", pct)),
                position = position_dodge(width = 0.7),
                hjust = -0.1, size = 1.8) +
      scale_fill_manual(values = bar_colors, name = NULL) +
      scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25),
                         labels = function(x) paste0(x, "%")) +
      coord_flip() +
      labs(x = NULL, y = "Percentage",
           title = "SuSiE vs ABF concordance by GWAS") +
      theme_masld() +
      theme(legend.position = "bottom",
            legend.key.size = unit(0.25, "cm"),
            legend.text = element_text(size = 5))
  } else {
    # Simple bar chart without fallback
    p_b <- ggplot(bar_dt, aes(x = gwas_short, y = concordance)) +
      geom_col(fill = masld_colors$down, width = 0.6) +
      geom_text(aes(label = sprintf("%.0f%%", concordance)),
                hjust = -0.1, size = 1.8) +
      scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25),
                         labels = function(x) paste0(x, "%")) +
      coord_flip() +
      labs(x = NULL, y = "Concordance (%)",
           title = "SuSiE vs ABF concordance by GWAS") +
      theme_masld()
  }

  message("  ", nrow(bar_dt), " GWAS plotted; mean concordance=",
          round(mean(bar_dt$concordance, na.rm = TRUE), 1), "%")
}, error = function(e) message("Panel (b) error: ", e$message))

# ═══════════════════════════════════════════════════════════════════════════
# Panel (c): Multi-Signal Gene Lollipop (n_cs_pairs > 1)
# ═══════════════════════════════════════════════════════════════════════════
tryCatch({
  message("Panel (c): Multi-signal gene bars (PP.H4 per GWAS)...")
  multi_f <- file.path(SUSIE_DIR, "susie_multi_signal_genes.csv")
  if (!file.exists(multi_f)) stop("susie_multi_signal_genes.csv not found")

  multi <- fread(multi_f)
  if (nrow(multi) == 0) stop("No multi-signal genes found")

  # Identify columns
  gwas_col_m <- intersect(names(multi), c("gwas", "GWAS", "source"))[1]
  if (is.na(gwas_col_m)) stop("No GWAS identifier column in multi-signal data")

  # Ensure PP.H4.susie exists
  pp4_col <- intersect(names(multi), c("PP.H4.susie", "PP.H4_susie", "PP.H4"))[1]
  if (is.na(pp4_col)) stop("No SuSiE PP.H4 column found")
  if (pp4_col != "PP.H4.susie") setnames(multi, pp4_col, "PP.H4.susie")

  # Ensure n_cs_pairs column exists
  cs_col <- intersect(names(multi), c("n_cs_pairs", "n_cs", "cs_pairs"))[1]
  if (is.na(cs_col)) stop("No credible set pairs column found")
  if (cs_col != "n_cs_pairs") setnames(multi, cs_col, "n_cs_pairs")

  # Add gene symbols
  gene_col_m <- intersect(names(multi), c("gene", "gene_id", "ensembl_id"))[1]
  if (!is.na(gene_col_m) && any(grepl("^ENSG", multi[[gene_col_m]]))) {
    multi <- add_symbols(multi, gene_col_m)
  } else if ("symbol" %in% names(multi)) {
    # Already has symbol
  } else if (!is.na(gene_col_m)) {
    multi[, symbol := get(gene_col_m)]
  }

  if (!"symbol" %in% names(multi)) stop("Cannot resolve gene symbols")

  # Keep only multi-signal (n_cs_pairs > 1)
  multi <- multi[n_cs_pairs > 1]
  if (nrow(multi) == 0) stop("No genes with n_cs_pairs > 1")

  # Sort genes by max PP.H4.susie across GWAS
  gene_order_m <- multi[, .(max_pp4 = max(PP.H4.susie, na.rm = TRUE)),
                         by = symbol]
  setorder(gene_order_m, max_pp4)
  multi[, symbol := factor(symbol, levels = gene_order_m$symbol)]

  # Shorten GWAS labels
  multi[, gwas_label := gsub("broadaway_|UKBB_|UKBB |Broadaway_", "",
                              get(gwas_col_m))]
  multi[, gwas_label := gsub("ukbb_", "", gwas_label)]

  # GWAS color palette (up to 6 EUR GWAS)
  gwas_levels <- sort(unique(multi$gwas_label))
  gwas_pal <- c("#1565C0", "#C2185B", "#7B1FA2", "#E91E63", "#42A5F5", "#00695C")
  names(gwas_pal) <- gwas_levels[seq_len(min(length(gwas_levels), 6))]

  # De-lollipop (2026-06-19): ONE horizontal BAR per gene = best (max) PP.H4.susie
  # across GWAS, coloured by the supporting GWAS; the credible-set-pair count
  # (the "multi-signal" metric) is annotated at the bar end. (Dodged per-GWAS bars
  # render as unreadable hairlines at 30 genes x 4 GWAS — single best bar is the
  # readable, faithful bar form.) Cap at top 24 genes by PP.H4.
  best <- multi[multi[, .I[which.max(PP.H4.susie)], by = symbol]$V1]
  ncs  <- multi[, .(n_cs = max(n_cs_pairs, na.rm = TRUE)), by = symbol]
  best[, n_cs := ncs$n_cs[match(symbol, ncs$symbol)]]
  setorder(best, PP.H4.susie)
  best <- tail(best, 24)
  best[, symbol := factor(symbol, levels = symbol)]

  p_c <- ggplot(best, aes(x = PP.H4.susie, y = symbol, fill = gwas_label)) +
    geom_col(width = 0.74) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50",
               linewidth = 0.3) +
    geom_text(aes(label = n_cs), hjust = -0.45, size = 1.7, color = "black") +
    scale_fill_manual(values = gwas_pal, name = "Best GWAS") +
    scale_x_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25),
                       expand = expansion(mult = c(0, 0.04))) +
    labs(x = "PP.H4 (SuSiE)  ·  number = credible-set pairs", y = NULL,
         title = "Multi-signal genes (>1 credible-set pair)") +
    theme_masld() +
    theme(axis.text.y = element_text(face = "italic", size = 5, color = "black"),
          legend.position = "right",
          legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = 5),
          panel.grid.major.y = element_blank())

  message("  ", uniqueN(best$symbol), " multi-signal genes shown; max n_cs_pairs=",
          max(best$n_cs))
}, error = function(e) message("Panel (c) error: ", e$message))

# ═══════════════════════════════════════════════════════════════════════════
# Save individual panels
# ═══════════════════════════════════════════════════════════════════════════
message("Saving individual panels...")
panel_list <- list(
  a_scatter     = p_a,
  b_concordance = p_b,
  c_bars        = p_c
)

for (nm in names(panel_list)) {
  panel_path <- file.path(PANEL_DIR, paste0("panel_", nm, ".pdf"))
  save_fig(panel_list[[nm]], panel_path, width = fig_half_width, height = 3)
  message("  Saved: ", panel_path)
}

# ═══════════════════════════════════════════════════════════════════════════
# Assemble composite figure
# ═══════════════════════════════════════════════════════════════════════════
message("Assembling composite figure...")
fig <- (p_a | p_b) / p_c +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "plain"))

save_fig_tall(fig, OUT, height = 7)
message("Saved composite: ", OUT)

message("=== figS_susie_comparison.R complete ===")
