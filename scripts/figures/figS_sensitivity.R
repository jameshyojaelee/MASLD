##############################################################################
# Supplementary Figure: Threshold Sensitivity Analysis
# 6 panels (2x3): (a) Dream DEG thresholds, (b) COLOC PP.H4,
#   (c) MR/TWAS p-value, (d) Multi-evidence sources_active,
#   (e) Cross-species concordance, (f) Deconvolution attribution
# Shows all key findings are stable across reasonable threshold ranges.
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
})


BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT     <- file.path(FIGS_SENS_DIR, "figS_sensitivity.pdf")
CSV_OUT <- file.path(FIGS_SENS_DIR, "figS_sensitivity_data.csv")

# Collect all sweep data for companion CSV
all_sweep_parts <- list()

# Pre-initialize all panels as placeholders
p_a <- placeholder("(a) Dream DEG thresholds")
p_b <- placeholder("(b) COLOC PP.H4")
p_c <- placeholder("(c) MR/TWAS p-value")
p_d <- placeholder("(d) Multi-evidence sources")
p_e <- placeholder("(e) Cross-species concordance")
p_f <- placeholder("(f) Deconvolution attribution")

# ═══════════════════════════════════════════════════════════════════════════
# Panel (a): Dream DEG Threshold Sensitivity
# ═══════════════════════════════════════════════════════════════════════════
tryCatch({
  message("Panel (a): Dream DEG threshold sweep (padj + |logFC|)...")
  dream <- load_dream_results()
  if (!is.null(dream)) {
    sig_vals <- c(0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2)
    lfc_vals <- c(0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0)
    sig_label <- "padj threshold"
    primary_sig <- 0.05; primary_lfc <- 0.3

    sweep_a <- rbindlist(lapply(sig_vals, function(sv) {
      rbindlist(lapply(lfc_vals, function(lf) {
        n <- sum(dream$dream_padj < sv & abs(dream$dream_logFC) >= lf, na.rm = TRUE)
        data.table(sig_threshold = sv, lfc_threshold = lf, n_degs = n)
      }))
    }))
    sweep_a[, panel := "a_dream"]
    all_sweep_parts[["a"]] <- sweep_a

    # Primary threshold count
    primary_n <- sweep_a[sig_threshold == primary_sig & lfc_threshold == primary_lfc, n_degs]

    # LFC as factor for color mapping
    sweep_a[, lfc_label := factor(sprintf("|LFC| >= %.2f", lfc_threshold),
                                  levels = sprintf("|LFC| >= %.2f", lfc_vals))]

    blue_grad <- colorRampPalette(c("#90CAF9", "#0D47A1"))(length(lfc_vals))

    p_a <- ggplot(sweep_a, aes(x = sig_threshold, y = n_degs / 1000,
                                color = lfc_label, group = lfc_label)) +
      geom_line(linewidth = 0.5) +
      geom_point(size = 0.8) +
      # Crosshair at primary threshold
      geom_vline(xintercept = primary_sig, linetype = "dashed", color = "red", linewidth = 0.3) +
      geom_hline(yintercept = primary_n / 1000, linetype = "dashed", color = "red", linewidth = 0.3) +
      annotate("text", x = primary_sig + 0.02, y = primary_n / 1000,
               label = paste0(format(primary_n, big.mark = ","), " DEGs"),
               size = 2, color = "red", hjust = 0, vjust = -0.5) +
      scale_color_manual(values = setNames(blue_grad, levels(sweep_a$lfc_label))) +
      scale_x_continuous(breaks = sig_vals) +
      labs(x = sig_label, y = "DEGs (thousands)", color = NULL,
           title = "Integrated DEG thresholds") +
      theme_masld() +
      theme(legend.position = "right",
            legend.key.height = unit(0.2, "cm"),
            legend.text = element_text(size = 5))

    message("  Primary: ", format(primary_n, big.mark = ","),
            " DEGs at ", sig_label, "=", primary_sig, ", |LFC|>=", primary_lfc)
  }
}, error = function(e) message("Panel (a) error: ", e$message))

# ═══════════════════════════════════════════════════════════════════════════
# Panel (b): COLOC PP.H4 Threshold Sensitivity
# ═══════════════════════════════════════════════════════════════════════════
tryCatch({
  message("Panel (b): COLOC PP.H4 sweep...")

  # Load COLOC files (Ghodsian excluded: 0 genes at all thresholds)
  # 2026 portfolio refactor: cirrhosis/HCC GWAS dropped (FinnGen_HCC removed)
  coloc_sources <- list(
    list(dir = "broadaway_ukbb",    label = "UKBB ALT",    ancestry = "European", lty = "solid"),
    list(dir = "broadaway_ukbb_ast", label = "UKBB AST",   ancestry = "European", lty = "solid"),
    list(dir = "broadaway_ukbb_ggt", label = "UKBB GGT",   ancestry = "European", lty = "solid"),
    list(dir = "broadaway_pdff",    label = "PDFF",         ancestry = "European", lty = "solid"),
    list(dir = "finngen_nafld",     label = "FinnGen NAFLD", ancestry = "FinnGen",  lty = "dashed"),
    list(dir = "finngen_nash",      label = "FinnGen NASH",  ancestry = "FinnGen",  lty = "dashed"),
    list(dir = "bbj_alt",           label = "BBJ ALT",      ancestry = "BBJ",      lty = "dotted"),
    list(dir = "bbj_ast",           label = "BBJ AST",      ancestry = "BBJ",      lty = "dotted"),
    list(dir = "bbj_ggt",           label = "BBJ GGT",      ancestry = "BBJ",      lty = "dotted")
  )

  pp_cutoffs <- seq(0.1, 0.9, by = 0.05)

  sweep_b <- rbindlist(lapply(coloc_sources, function(src) {
    f <- file.path(CAUSAL, src$dir, "coloc_results.csv")
    if (!file.exists(f)) return(NULL)
    dt <- fread(f, select = c("gene", "PP.H4"))
    rbindlist(lapply(pp_cutoffs, function(pp) {
      n <- sum(dt$PP.H4 > pp, na.rm = TRUE)
      data.table(pp_threshold = pp, n_genes = n,
                 source = src$label, ancestry = src$ancestry, lty = src$lty)
    }))
  }))
  sweep_b[, panel := "b_coloc"]
  all_sweep_parts[["b"]] <- sweep_b[, .(pp_threshold, n_genes, source, ancestry, panel, lty)]

  # Color by ancestry group
  # RESTORED 2026-04-09: FinnGen color restored
  ancestry_colors <- c(
    "European" = "#1565C0",
    "FinnGen"  = "#2E7D32",
    "BBJ"      = "#7B1FA2"
  )

  p_b <- ggplot(sweep_b, aes(x = pp_threshold, y = n_genes,
                              color = ancestry, group = source, linetype = lty)) +
    geom_line(linewidth = 0.4) +
    geom_point(size = 0.5) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "red", linewidth = 0.3) +
    scale_color_manual(values = ancestry_colors) +
    scale_linetype_identity() +
    labs(x = "PP.H4 threshold", y = "Colocalized genes", color = NULL,
         title = "COLOC PP.H4 thresholds") +
    theme_masld() +
    theme(legend.position = "right",
          legend.key.height = unit(0.2, "cm"),
          legend.text = element_text(size = 5))

  # Report primary counts
  primary_b <- sweep_b[pp_threshold == 0.5, .(source, n_genes)]
  message("  PP.H4 > 0.5: ", paste(primary_b$source, primary_b$n_genes, sep = "=", collapse = ", "))

}, error = function(e) message("Panel (b) error: ", e$message))

# ═══════════════════════════════════════════════════════════════════════════
# Panel (c): TWAS / sc-TWAS p-value Sensitivity (MR trace dropped 2026-04-22)
# ═══════════════════════════════════════════════════════════════════════════
tryCatch({
  message("Panel (c): TWAS p-value sweep (MR ditched)...")

  neglog_cuts <- seq(1, 5, by = 0.25)
  p_cuts <- 10^(-neglog_cuts)

  # MR trace removed 2026-04-22 — MR ditched from paper.
  mr_sweep <- NULL

  # TWAS
  twas <- load_twas_results("ghodsian")
  twas_sweep <- NULL
  if (!is.null(twas)) {
    twas_sweep <- rbindlist(lapply(p_cuts, function(pc) {
      data.table(neglog10p = -log10(pc),
                 n_genes = sum(twas$pvalue < pc, na.rm = TRUE),
                 analysis = "TWAS (Ghodsian)")
    }))
  }

  # sc-TWAS
  sceqtl_f <- file.path(CAUSAL, "sceqtl_twas", "sceqtl_twas_all_results.csv")
  sctwas_sweep <- NULL
  if (file.exists(sceqtl_f)) {
    sctwas <- fread(sceqtl_f)
    # Best p-value per gene across cell types/GWAS/conditions
    sctwas_best <- sctwas[, .(best_p = min(pval, na.rm = TRUE)), by = exposure]
    sctwas_sweep <- rbindlist(lapply(p_cuts, function(pc) {
      data.table(neglog10p = -log10(pc),
                 n_genes = sum(sctwas_best$best_p < pc, na.rm = TRUE),
                 analysis = "sc-TWAS")
    }))
  }

  sweep_c <- rbindlist(list(mr_sweep, twas_sweep, sctwas_sweep), use.names = TRUE)
  if (nrow(sweep_c) > 0) {
    sweep_c[, panel := "c_mr_twas"]
    all_sweep_parts[["c"]] <- sweep_c

    causal_colors <- c(
      # MR color removed 2026-04-22 — MR ditched from paper.
      "TWAS (Ghodsian)" = "#E91E63",
      "sc-TWAS"          = "#7B1FA2"
    )

    p_c <- ggplot(sweep_c, aes(x = neglog10p, y = n_genes,
                                color = analysis, group = analysis)) +
      geom_line(linewidth = 0.5) +
      geom_point(size = 0.8) +
      geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.3) +
      annotate("text", x = -log10(0.05), y = Inf, label = "p=0.05",
               size = 2, color = "red", vjust = 1.5, hjust = -0.1) +
      scale_color_manual(values = causal_colors) +
      labs(x = expression(-log[10](p)), y = "Significant genes", color = NULL,
           title = "MR/TWAS p-value thresholds") +
      theme_masld() +
      theme(legend.position = "right",
            legend.key.height = unit(0.2, "cm"),
            legend.text = element_text(size = 5))

    # Report counts near p<0.05
    near_005 <- sweep_c[abs(neglog10p - (-log10(0.05))) < 0.15]
    if (nrow(near_005) > 0) {
      message("  Near p~0.05: ",
              paste(near_005[, paste(analysis, n_genes, sep = "=")], collapse = ", "))
    }
  }
}, error = function(e) message("Panel (c) error: ", e$message))

# ═══════════════════════════════════════════════════════════════════════════
# Panel (d): Multi-Evidence sources_active Sensitivity
# ═══════════════════════════════════════════════════════════════════════════
tryCatch({
  message("Panel (d): Multi-evidence sources_active sweep...")

  atlas <- load_multi_evidence()
  posctrl <- load_positive_controls()

  if (!is.null(atlas) && "sources_active" %in% names(atlas)) {
    src_vals <- 1:7

    # Match positive controls to atlas
    posctrl_genes <- character(0)
    if (!is.null(posctrl)) {
      # Filter to expression-driven controls
      if ("control_type" %in% names(posctrl)) {
        posctrl <- posctrl[control_type == "Expression_driven"]
      }
      ctrl_sym <- unique(posctrl$gene)
      posctrl_genes <- ctrl_sym[ctrl_sym %in% atlas$human_symbol]
      message("  Matched ", length(posctrl_genes), " expression-driven positive controls to atlas")
    }

    sweep_d <- rbindlist(lapply(src_vals, function(k) {
      genes_at_k <- atlas[sources_active >= k]
      n_genes <- nrow(genes_at_k)
      n_ctrl <- sum(genes_at_k$human_symbol %in% posctrl_genes)
      pct_ctrl <- if (length(posctrl_genes) > 0) 100 * n_ctrl / length(posctrl_genes) else NA_real_
      data.table(min_sources = k, n_genes = n_genes, n_controls = n_ctrl, pct_controls = pct_ctrl)
    }))
    sweep_d[, panel := "d_multi_evidence"]
    all_sweep_parts[["d"]] <- sweep_d

    # Scaling factor for dual y-axis
    max_genes <- max(sweep_d$n_genes, na.rm = TRUE)
    scale_factor <- max_genes / 100

    p_d <- ggplot(sweep_d, aes(x = min_sources)) +
      geom_col(aes(y = n_genes), fill = masld_colors$up, alpha = 0.7, width = 0.7) +
      geom_line(aes(y = pct_controls * scale_factor, group = 1),
                color = masld_colors$conserved, linewidth = 0.7) +
      geom_point(aes(y = pct_controls * scale_factor),
                 color = masld_colors$conserved, size = 1.5) +
      scale_y_continuous(
        name = "Genes at threshold",
        labels = comma,
        sec.axis = sec_axis(~ . / scale_factor, name = "% positive controls recovered")
      ) +
      scale_x_continuous(breaks = src_vals) +
      labs(x = "Minimum sources active", title = "Multi-evidence convergence") +
      theme_masld() +
      theme(axis.title.y.right = element_text(color = masld_colors$conserved),
            axis.text.y.right = element_text(color = masld_colors$conserved))

    message("  At sources>=1: ", sweep_d[min_sources == 1, n_genes], " genes, ",
            round(sweep_d[min_sources == 1, pct_controls], 1), "% controls")
  }
}, error = function(e) message("Panel (d) error: ", e$message))

# ═══════════════════════════════════════════════════════════════════════════
# Panel (e): Cross-Species Conserved Sensitivity
# ═══════════════════════════════════════════════════════════════════════════
tryCatch({
  message("Panel (e): Cross-species concordance sweep...")

  conc <- load_concordance_atlas()
  if (!is.null(conc) && "n_concordant" %in% names(conc)) {
    k_vals <- 1:5

    sweep_e <- rbindlist(lapply(k_vals, function(k) {
      data.table(min_concordant = k,
                 n_genes = sum(conc$n_concordant >= k, na.rm = TRUE))
    }))
    sweep_e[, panel := "e_concordance"]
    all_sweep_parts[["e"]] <- sweep_e

    cc_n <- sweep_e[min_concordant == 3, n_genes]

    p_e <- ggplot(sweep_e, aes(x = min_concordant, y = n_genes)) +
      geom_step(color = masld_colors$conserved, linewidth = 0.7, direction = "mid") +
      geom_point(color = masld_colors$conserved, size = 1.5) +
      # Reference lines at primary threshold (k=3)
      geom_hline(yintercept = cc_n, linetype = "dashed", color = "red", linewidth = 0.3) +
      geom_vline(xintercept = 3, linetype = "dashed", color = "red", linewidth = 0.3) +
      annotate("text", x = 3.2, y = cc_n + 200,
               label = paste0(format(cc_n, big.mark = ","), " Conserved"),
               size = 2, color = "red", hjust = 0) +
      scale_x_continuous(breaks = k_vals) +
      scale_y_continuous(labels = comma) +
      labs(x = "Min concordant mouse models", y = "Genes",
           title = "Cross-species concordance") +
      theme_masld()

    message("  At k>=3: ", sweep_e[min_concordant == 3, n_genes], " genes")
  }
}, error = function(e) message("Panel (e) error: ", e$message))

# ═══════════════════════════════════════════════════════════════════════════
# Panel (f): Deconvolution Attribution Sensitivity
# ═══════════════════════════════════════════════════════════════════════════
tryCatch({
  message("Panel (f): Deconvolution attribution sweep...")

  deconv <- load_deconv_attribution()
  if (!is.null(deconv)) {
    padj_vals_f <- c(0.01, 0.025, 0.05, 0.075, 0.1)
    lfc_vals_f  <- c(0.25, 0.5, 0.75, 1.0)

    sweep_f <- rbindlist(lapply(padj_vals_f, function(pa) {
      rbindlist(lapply(lfc_vals_f, function(lf) {
        sig_adj   <- deconv$padj_adj < pa & abs(deconv$logFC_adj) >= lf
        sig_unadj <- deconv$padj_unadj < pa & abs(deconv$logFC_unadj) >= lf
        hep_intrinsic <- sig_adj & sig_unadj
        data.table(padj_threshold = pa, lfc_threshold = lf,
                   n_hep_intrinsic = sum(hep_intrinsic, na.rm = TRUE))
      }))
    }))
    sweep_f[, panel := "f_deconv"]
    all_sweep_parts[["f"]] <- sweep_f

    primary_f <- sweep_f[padj_threshold == 0.05 & lfc_threshold == 0.5, n_hep_intrinsic]

    sweep_f[, lfc_label := factor(sprintf("|LFC| >= %.2f", lfc_threshold),
                                  levels = sprintf("|LFC| >= %.2f", lfc_vals_f))]

    blue_grad_f <- colorRampPalette(c("#90CAF9", "#0D47A1"))(length(lfc_vals_f))

    p_f <- ggplot(sweep_f, aes(x = padj_threshold, y = n_hep_intrinsic,
                                color = lfc_label, group = lfc_label)) +
      geom_line(linewidth = 0.5) +
      geom_point(size = 0.8) +
      geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.3) +
      geom_hline(yintercept = primary_f, linetype = "dashed", color = "red", linewidth = 0.3) +
      annotate("text", x = 0.05, y = primary_f + 200,
               label = paste0(format(primary_f, big.mark = ","), " genes"),
               size = 2, color = "red", hjust = 0.5) +
      scale_color_manual(values = setNames(blue_grad_f, levels(sweep_f$lfc_label))) +
      scale_x_continuous(breaks = padj_vals_f) +
      labs(x = "padj threshold", y = "Hepatocyte-intrinsic genes", color = NULL,
           title = "Deconvolution attribution") +
      theme_masld() +
      theme(legend.position = "right",
            legend.key.height = unit(0.2, "cm"),
            legend.text = element_text(size = 5))

    message("  Primary: ", format(primary_f, big.mark = ","),
            " hep-intrinsic at padj=0.05, |LFC|>=0.5")
  }
}, error = function(e) message("Panel (f) error: ", e$message))

# ═══════════════════════════════════════════════════════════════════════════
# Save individual panels
# ═══════════════════════════════════════════════════════════════════════════
PANEL_DIR <- file.path(FIGS_SENS_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

panel_list <- list(a = p_a, b = p_b, c = p_c, d = p_d, e = p_e, f = p_f)
panel_w <- fig_half_width
panel_h <- 2.5

for (nm in names(panel_list)) {
  panel_path <- file.path(PANEL_DIR, paste0("panel_", nm, ".pdf"))
  save_fig(panel_list[[nm]], panel_path, width = panel_w, height = panel_h)
  message("  Saved panel ", nm, ": ", panel_path)
}

# ═══════════════════════════════════════════════════════════════════════════
# Assemble composite figure
# ═══════════════════════════════════════════════════════════════════════════
message("Assembling 2x3 figure...")
fig <- (p_a | p_b | p_c) / (p_d | p_e | p_f)
fig <- auto_tag(fig)
save_fig(fig, OUT, width = fig_full_width, height = 5)
message("Saved: ", OUT)

# Write companion CSV
all_sweep <- rbindlist(all_sweep_parts, fill = TRUE)
fwrite(all_sweep, CSV_OUT)
message("Saved companion CSV: ", CSV_OUT, " (", nrow(all_sweep), " rows)")

message("=== figS_sensitivity.R complete ===")
