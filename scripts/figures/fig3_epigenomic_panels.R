##############################################################################
# Figure 3 panels (f)-(i): Epigenomic Regulatory Evidence
#   (f) SCENIC+ disease regulon activity bar
#   (g) SCENIC+ TF-target regulatory network
#   (h) chromVAR TF motif enrichment
#   (i) Epigenomic-transcriptomic convergence scatter
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(ggrepel)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

dir.create(file.path(FIG3_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)
OUT <- file.path(FIG3_DIR, "panels", "fig3_epigenomic_panels.pdf")

# Pre-initialize all panels with placeholders
p_f <- placeholder("Panel f: SCENIC+ regulon activity")
p_g <- placeholder("Panel g: TF-target network")
p_h <- placeholder("Panel h: chromVAR motif enrichment")
p_i <- placeholder("Panel i: Epigenomic-transcriptomic convergence")

# Drug target TFs to highlight across panels
DRUG_TARGET_TFS <- c("THRB", "NR1H4", "PPARA", "PPARG")

# ==========================================================================
# Panel (f): SCENIC+ disease regulon activity bar (horizontal)
#   24 TFs ordered by regulon_activity_diff (most negative at top)
#   Blue = downregulated in MASLD, magenta = upregulated
#   Stars for significance; drug targets annotated
# ==========================================================================

regulons <- load_disease_regulons()

if (!is.null(regulons) && nrow(regulons) > 0) {
  reg_f <- copy(regulons)

  # Ensure required columns
  if (all(c("tf_name", "regulon_activity_diff", "activity_padj") %in% names(reg_f))) {

    # Order by regulon_activity_diff (most negative at top of bar chart)
    setorder(reg_f, regulon_activity_diff)
    reg_f[, tf_name := factor(tf_name, levels = rev(reg_f$tf_name))]

    # Significance stars
    reg_f[, stars := fifelse(
      !is.na(activity_padj) & activity_padj < 0.001, "***",
      fifelse(activity_padj < 0.01, "**",
        fifelse(activity_padj < 0.05, "*", ""))
    )]

    # Color by direction
    reg_f[, direction := fifelse(regulon_activity_diff >= 0, "up", "down")]

    # Determine which labels to bold (drug targets)
    tf_levels <- levels(reg_f$tf_name)
    y_faces <- ifelse(tf_levels %in% DRUG_TARGET_TFS, "bold", "plain")
    names(y_faces) <- tf_levels

    # Position for stars: slightly beyond the bar end
    bar_range <- max(abs(reg_f$regulon_activity_diff), na.rm = TRUE)
    nudge <- bar_range * 0.06
    reg_f[, star_x := fifelse(
      regulon_activity_diff >= 0,
      regulon_activity_diff + nudge,
      regulon_activity_diff - nudge
    )]

    p_f <- ggplot(reg_f, aes(x = regulon_activity_diff, y = tf_name,
                              fill = direction)) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_text(aes(x = star_x, label = stars), size = 2, vjust = 0.5,
                show.legend = FALSE) +
      # Annotate drug targets with a small triangle marker
      geom_point(data = reg_f[tf_name %in% DRUG_TARGET_TFS],
                 aes(x = 0, y = tf_name),
                 shape = 18, size = 1.5, color = "#00695C",
                 inherit.aes = FALSE) +
      scale_fill_manual(
        values = c("down" = epigenomic_colors[["regulon_down"]],
                   "up"   = epigenomic_colors[["regulon_up"]]),
        labels = c("down" = "Suppressed in MASLD",
                   "up"   = "Activated in MASLD"),
        name = NULL
      ) +
      geom_vline(xintercept = 0, linewidth = 0.3, color = "gray40") +
      scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
      labs(x = "Regulon activity change\n(Normal \u2192 MASLD)",
           y = NULL,
           title = paste0("SCENIC+ disease regulons (n = ",
                          nrow(reg_f), ")")) +
      theme_masld() +
      theme(axis.text.y = element_text(size = 5.5, face = y_faces),
            legend.position = "inside",
            legend.position.inside = c(0.75, 0.15),
            legend.background = element_blank(),
            legend.key = element_blank(),
            legend.key.size = unit(0.25, "cm"))

    message("Panel f: ", nrow(reg_f), " disease regulons plotted")
  } else {
    message("Panel f: missing required columns in disease_regulons.csv")
  }
}

# ==========================================================================
# Panel (g): SCENIC+ TF-target regulatory network (bipartite dot plot)
#   Top 10 most significant TFs connected to their target genes
#   TF color by direction; dot size by n_enhancers
# ==========================================================================

hep_regulons <- load_hepatocyte_regulons()

if (!is.null(hep_regulons) && nrow(hep_regulons) > 0 &&
    !is.null(regulons) && nrow(regulons) > 0) {

  # Identify top 10 disease-significant TFs by activity_padj
  sig_tfs <- regulons[!is.na(activity_padj) & activity_padj < 0.05]
  setorder(sig_tfs, activity_padj)
  top10_tfs <- head(sig_tfs$tf_name, 10)

  if (length(top10_tfs) > 0) {
    # Filter hepatocyte regulon links to top TFs
    links <- hep_regulons[tf_name %in% top10_tfs]

    if (nrow(links) > 0) {
      # Get direction for each TF
      tf_dir <- regulons[tf_name %in% top10_tfs,
                          .(tf_name, regulon_activity_diff)]
      tf_dir[, direction := fifelse(regulon_activity_diff >= 0, "up", "down")]
      links <- merge(links, tf_dir[, .(tf_name, direction)],
                     by = "tf_name", all.x = TRUE)

      # Assign y-coordinates: TFs on left, targets on right
      # TFs ordered by significance (top = most significant)
      tf_order <- top10_tfs
      tf_y <- data.table(tf_name = tf_order,
                          tf_y = seq(length(tf_order), 1, by = -1))

      # For each TF, spread its targets evenly
      links <- merge(links, tf_y, by = "tf_name", all.x = TRUE)

      # Collect unique targets and assign y positions
      unique_targets <- unique(links$target_gene)
      # Space targets evenly across the y range
      n_targets <- length(unique_targets)
      target_y_vals <- seq(length(tf_order), 1,
                           length.out = max(n_targets, 1))
      target_y <- data.table(target_gene = unique_targets,
                              tgt_y = target_y_vals[seq_len(n_targets)])
      links <- merge(links, target_y, by = "target_gene", all.x = TRUE)

      # Ensure n_enhancers is numeric
      if ("n_enhancers" %in% names(links)) {
        links[, n_enhancers := as.numeric(n_enhancers)]
        links[is.na(n_enhancers), n_enhancers := 1]
      } else {
        links[, n_enhancers := 1]
      }

      # Build plot: segments connecting TF to targets, dots at each end
      x_tf  <- 0
      x_tgt <- 1

      p_g <- ggplot() +
        # Segments connecting TF to target
        geom_segment(data = links,
                     aes(x = x_tf, xend = x_tgt, y = tf_y, yend = tgt_y,
                         color = direction),
                     linewidth = 0.3, alpha = 0.4) +
        # TF dots (left)
        geom_point(data = unique(links[, .(tf_name, tf_y, direction)]),
                   aes(x = x_tf, y = tf_y, color = direction),
                   size = 2.5, shape = 16) +
        # Target dots (right), sized by n_enhancers
        geom_point(data = links[, .(n_enhancers = max(n_enhancers)),
                                 by = .(target_gene, tgt_y)],
                   aes(x = x_tgt, y = tgt_y, size = n_enhancers),
                   color = "gray40", shape = 16, alpha = 0.7) +
        # TF labels (left side)
        geom_text(data = unique(links[, .(tf_name, tf_y, direction)]),
                  aes(x = x_tf - 0.05, y = tf_y, label = tf_name,
                      fontface = ifelse(tf_name %in% DRUG_TARGET_TFS,
                                        "bold.italic", "italic")),
                  hjust = 1, size = 2, color = "black") +
        # Target labels (right side)
        geom_text(data = links[, .(target_gene, tgt_y)][!duplicated(target_gene)],
                  aes(x = x_tgt + 0.05, y = tgt_y, label = target_gene),
                  hjust = 0, size = 1.8, fontface = "italic", color = "gray30") +
        scale_color_manual(
          values = c("down" = epigenomic_colors[["regulon_down"]],
                     "up"   = epigenomic_colors[["regulon_up"]]),
          labels = c("down" = "Suppressed", "up" = "Activated"),
          name = "Regulon direction"
        ) +
        scale_size_continuous(range = c(1, 3.5), name = "Enhancers",
                              breaks = pretty_breaks(3)) +
        coord_cartesian(xlim = c(-0.45, 1.55), clip = "off") +
        labs(title = paste0("TF \u2192 target links (top ",
                            length(top10_tfs), " TFs)")) +
        theme_masld() +
        theme(axis.text  = element_blank(),
              axis.title = element_blank(),
              axis.line  = element_blank(),
              axis.ticks = element_blank(),
              legend.position = "inside",
              legend.position.inside = c(0.50, 0.10),
              legend.background = element_blank(),
              legend.key = element_blank(),
              legend.key.size = unit(0.25, "cm"),
              legend.direction = "horizontal")

      message("Panel g: ", length(top10_tfs), " TFs, ",
              nrow(links), " links plotted")
    } else {
      message("Panel g: no TF-target links for top TFs")
    }
  } else {
    message("Panel g: no significant disease regulons (padj < 0.05)")
  }
}

# ==========================================================================
# Panel (h): chromVAR TF motif enrichment
#   Filter corrupted motifs (|mean_deviation| > 5), show balanced top 15 up + 15 down
#   Color: direction; significance markers for padj < 0.05
# ==========================================================================

chromvar <- load_chromvar_hepatocyte()

if (!is.null(chromvar) && nrow(chromvar) > 0) {
  cv <- copy(chromvar)

  if (all(c("tf_name", "logFC_deviation") %in% names(cv))) {

    cv <- cv[!is.na(logFC_deviation)]

    # Filter corrupted motifs (numeric overflow from rare motifs)
    n_before <- nrow(cv)
    cv <- cv[abs(mean_deviation_masld) < 5 & abs(mean_deviation_normal) < 5]
    n_dropped <- n_before - nrow(cv)
    if (n_dropped > 0)
      message("  chromVAR: filtered ", n_dropped, " corrupted motifs (|dev| > 5)")

    # Significant only, balanced top 15 up + 15 down
    cv_sig <- cv[!is.na(padj) & padj < 0.05]
    cv_up   <- head(cv_sig[logFC_deviation > 0][order(-logFC_deviation)], 15)
    cv_down <- head(cv_sig[logFC_deviation < 0][order(logFC_deviation)], 15)
    cv_top  <- rbind(cv_up, cv_down)

    setorder(cv_top, logFC_deviation)
    cv_top[, tf_name := factor(tf_name, levels = cv_top$tf_name)]
    cv_top[, direction := fifelse(logFC_deviation >= 0, "up", "down")]

    cv_top[, stars := ""]
    if ("padj" %in% names(cv_top)) {
      cv_top[!is.na(padj) & padj < 0.001, stars := "***"]
      cv_top[stars == "" & !is.na(padj) & padj < 0.01, stars := "**"]
      cv_top[stars == "" & !is.na(padj) & padj < 0.05, stars := "*"]
    }

    bar_range_h <- max(abs(cv_top$logFC_deviation), na.rm = TRUE)
    nudge_h <- bar_range_h * 0.06
    cv_top[, star_x := fifelse(
      logFC_deviation >= 0,
      logFC_deviation + nudge_h,
      logFC_deviation - nudge_h
    )]

    # Bold known MASLD TFs
    known_tfs <- c("HNF4A", "PPARA", "NR1H4", "Nr1H4", "Nr1h3", "NR1H3",
                   "THRB", "HNF1A", "FOXA1", "FOXA2", "ETS1", "CEBPB")
    tf_levels <- levels(cv_top$tf_name)
    y_faces <- ifelse(toupper(tf_levels) %in% toupper(known_tfs), "bold", "plain")
    names(y_faces) <- tf_levels

    n_sig_up <- nrow(cv_sig[logFC_deviation > 0])
    n_sig_dn <- nrow(cv_sig[logFC_deviation < 0])

    p_h <- ggplot(cv_top, aes(x = logFC_deviation, y = tf_name, fill = direction)) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_text(aes(x = star_x, label = stars),
                size = 2, vjust = 0.5, show.legend = FALSE) +
      scale_fill_manual(
        values = c("down" = epigenomic_colors[["regulon_down"]],
                   "up"   = epigenomic_colors[["regulon_up"]]),
        labels = c("down" = paste0("Lower (", n_sig_dn, " sig)"),
                   "up"   = paste0("Higher (", n_sig_up, " sig)")),
        name = NULL
      ) +
      geom_vline(xintercept = 0, linewidth = 0.3, color = "gray40") +
      scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
      labs(x = "Motif accessibility change",
           y = NULL,
           title = paste0("chromVAR hepatocyte motifs (top 15\u2191 + 15\u2193)")) +
      theme_masld() +
      theme(axis.text.y = element_text(size = 5, face = y_faces),
            legend.position = "inside",
            legend.position.inside = c(0.80, 0.15),
            legend.background = element_blank(),
            legend.key = element_blank(),
            legend.key.size = unit(0.25, "cm"))

    message("Panel h: chromVAR done (", nrow(cv_top), " motifs, ",
            n_sig_up, " up + ", n_sig_dn, " down significant)")
  } else {
    message("Panel h: missing required columns in chromvar data")
  }
}

# ==========================================================================
# Panel (i): Epigenomic-transcriptomic convergence scatter
#   X-axis: dream_logFC (from dream results)
#   Y-axis: regulon_activity_diff (from SCENIC+, for genes in hepatocyte
#           regulon targets)
#   Spearman correlation annotation; highlight key genes
# ==========================================================================

dream <- load_dream_results()

if (!is.null(dream) && nrow(dream) > 0 &&
    !is.null(hep_regulons) && nrow(hep_regulons) > 0 &&
    !is.null(regulons) && nrow(regulons) > 0) {

  # Build a table: for each target_gene in hepatocyte_regulons,
  # get its associated TF's regulon_activity_diff and the gene's dream_logFC
  hr <- copy(hep_regulons)

  # Merge TF-level regulon_activity_diff
  tf_info <- regulons[, .(tf_name, regulon_activity_diff, activity_padj)]
  hr <- merge(hr, tf_info, by = "tf_name", all.x = TRUE,
              suffixes = c("", ".tf"))

  # Merge dream results for target genes
  dream_slim <- dream[, .(symbol, dream_logFC, dream_padj)]
  dream_slim <- dream_slim[!duplicated(symbol)]
  conv <- merge(hr, dream_slim, by.x = "target_gene", by.y = "symbol",
                all.x = TRUE)

  # Keep only genes with both values
  conv <- conv[!is.na(dream_logFC) & !is.na(regulon_activity_diff)]

  if (nrow(conv) > 0) {
    # If a gene appears under multiple TFs, take the one with most
    # significant TF
    conv <- conv[order(activity_padj)]
    conv <- conv[!duplicated(target_gene)]

    # Significance classification
    conv[, sig_class := fifelse(
      !is.na(dream_padj) & dream_padj < 0.1 &
        !is.na(activity_padj) & activity_padj < 0.05,
      "Both significant",
      "Other"
    )]

    # Spearman correlation
    cor_test <- tryCatch(
      cor.test(conv$dream_logFC, conv$regulon_activity_diff,
               method = "spearman"),
      error = function(e) NULL
    )
    if (!is.null(cor_test)) {
      cor_label <- sprintf("rho = %.3f, p = %.2e",
                           cor_test$estimate, cor_test$p.value)
    } else {
      cor_label <- "Spearman: N/A"
    }

    # Genes to highlight
    highlight_genes <- c("THRB", "NR1H4", "PPARA", "PPARG", "HNF4A")
    label_dt <- conv[target_gene %in% highlight_genes |
                     tf_name %in% highlight_genes]
    # Also label targets of highlighted TFs if they themselves are notable
    # Use target_gene as label
    label_dt[, label := target_gene]

    # Colors
    sig_colors <- c("Both significant" = epigenomic_colors[["convergent"]],
                    "Other"            = "gray70")

    p_i <- ggplot(conv, aes(x = dream_logFC, y = regulon_activity_diff)) +
      rasterize_layer(
        geom_point(aes(color = sig_class), size = 0.8, alpha = 0.6, shape = 16)
      ) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 0.5,
                  color = "gray30", linetype = "dashed") +
      geom_hline(yintercept = 0, linewidth = 0.2, color = "gray50") +
      geom_vline(xintercept = 0, linewidth = 0.2, color = "gray50") +
      scale_color_manual(values = sig_colors, name = NULL) +
      labs(x = expression("Transcriptomic log"[2]*"FC"),
           y = "Regulon activity change\n(SCENIC+)",
           title = "SCENIC+ transcriptomic convergence") +
      annotate("text", x = Inf, y = Inf, label = cor_label,
               hjust = 1.05, vjust = 1.5, size = 2.2, color = "gray20") +
      annotate("text", x = Inf, y = Inf,
               label = paste0("n = ", nrow(conv), " genes"),
               hjust = 1.05, vjust = 3.0, size = 2, color = "gray40") +
      theme_masld() +
      theme(legend.position = "inside",
            legend.position.inside = c(0.25, 0.90),
            legend.background = element_blank(),
            legend.key = element_blank(),
            legend.key.size = unit(0.25, "cm"))

    # Add repel labels if any highlight genes found
    if (nrow(label_dt) > 0) {
      p_i <- p_i +
        geom_label_repel(
          data = label_dt,
          aes(label = label),
          size = 1.8, max.overlaps = 15,
          label.padding = 0.1, segment.size = 0.15,
          min.segment.length = 0, fontface = "italic",
          color = "black", fill = "white", alpha = 0.85,
          show.legend = FALSE
        )
    }

    message("Panel i: ", nrow(conv), " convergent genes plotted, ", cor_label)
  } else {
    message("Panel i: no overlapping genes between dream and SCENIC+ targets")
  }
}

# ==========================================================================
# Panel (j): chromVAR-transcriptomic convergence scatter
#   X-axis: dream logFC of the TF gene itself
#   Y-axis: chromVAR motif accessibility change (logFC_deviation)
#   Shows whether TFs with altered expression also have altered chromatin
# ==========================================================================

p_j <- placeholder("Panel j: chromVAR-transcriptomic convergence")

if (!is.null(chromvar) && nrow(chromvar) > 0 &&
    !is.null(dream) && nrow(dream) > 0) {

  cv_conv <- copy(chromvar)
  cv_conv <- cv_conv[!is.na(logFC_deviation)]

  # Filter corrupted motifs (same as panel h)
  cv_conv <- cv_conv[abs(mean_deviation_masld) < 5 & abs(mean_deviation_normal) < 5]

  # Restrict to significant motifs (padj < 0.05)
  cv_conv <- cv_conv[!is.na(padj) & padj < 0.05]
  message("  chromVAR convergence: ", nrow(cv_conv), " significant motifs (padj<0.05)")

  # Deduplicate case variants: normalize tf_name to uppercase for matching,
  # then keep the entry with smallest padj per unique TF
  cv_conv[, tf_upper := toupper(gsub("::.*", "", tf_name))]
  setorder(cv_conv, padj)
  cv_conv <- cv_conv[!duplicated(tf_upper)]
  message("  After deduplication: ", nrow(cv_conv), " unique TFs")

  # Merge dream logFC for the TF itself
  dream_slim2 <- dream[, .(symbol, dream_logFC, dream_padj)]
  dream_slim2 <- dream_slim2[!duplicated(symbol)]
  dream_slim2[, symbol_upper := toupper(symbol)]

  cv_merge <- merge(cv_conv, dream_slim2,
                    by.x = "tf_upper", by.y = "symbol_upper",
                    all.x = FALSE)

  cv_merge <- cv_merge[!is.na(dream_logFC) & !is.na(logFC_deviation)]
  message("  Merged with dream: ", nrow(cv_merge), " TFs with both data")

  if (nrow(cv_merge) > 0) {
    # All chromVAR motifs already significant; classify by dream significance
    cv_merge[, sig_class := fifelse(
      !is.na(dream_padj) & dream_padj < 0.1,
      "Both significant",
      "chromVAR only"
    )]

    # Spearman correlation
    cor_test_cv <- tryCatch(
      cor.test(cv_merge$dream_logFC, cv_merge$logFC_deviation,
               method = "spearman"),
      error = function(e) NULL
    )
    if (!is.null(cor_test_cv)) {
      cor_label_cv <- sprintf("rho = %.3f, p = %.2e",
                               cor_test_cv$estimate, cor_test_cv$p.value)
    } else {
      cor_label_cv <- "Spearman: N/A"
    }

    # ~12 labels spanning all quadrants: drug targets + top chromVAR hits
    highlight_cv <- c(
      "THRB", "NR1H4", "PPARA", "HNF4A",    # drug targets (bottom-left)
      "GMEB1", "KLF7", "FOSL1",              # top (increased accessibility)
      "ZEB1", "NR1H3", "NFIA",               # bottom (decreased accessibility)
      "ETS1", "CEBPB"                         # central, biologically key
    )
    label_cv <- cv_merge[tf_upper %in% highlight_cv]
    label_cv <- label_cv[!duplicated(tf_upper)]
    # Use uppercase canonical name for labels
    label_cv[, label := tf_upper]

    sig_colors_cv <- c("Both significant" = epigenomic_colors[["convergent"]],
                       "chromVAR only"    = "gray70")

    p_j <- ggplot(cv_merge, aes(x = dream_logFC, y = logFC_deviation)) +
      rasterize_layer(
        geom_point(aes(color = sig_class), size = 0.8, alpha = 0.6, shape = 16)
      ) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 0.5,
                  color = "gray30", linetype = "dashed") +
      geom_hline(yintercept = 0, linewidth = 0.2, color = "gray50") +
      geom_vline(xintercept = 0, linewidth = 0.2, color = "gray50") +
      scale_color_manual(values = sig_colors_cv, name = NULL) +
      labs(x = expression("Transcriptomic log"[2]*"FC"),
           y = "Motif accessibility change\n(chromVAR)",
           title = "chromVAR transcriptomic convergence") +
      annotate("text", x = Inf, y = Inf, label = cor_label_cv,
               hjust = 1.05, vjust = 1.5, size = 2.2, color = "gray20") +
      annotate("text", x = Inf, y = Inf,
               label = paste0("n = ", nrow(cv_merge), " TFs"),
               hjust = 1.05, vjust = 3.0, size = 2, color = "gray40") +
      theme_masld() +
      theme(legend.position = "inside",
            legend.position.inside = c(0.25, 0.90),
            legend.background = element_blank(),
            legend.key = element_blank(),
            legend.key.size = unit(0.25, "cm"))

    # Add repel labels
    if (nrow(label_cv) > 0) {
      p_j <- p_j +
        geom_label_repel(
          data = label_cv,
          aes(label = label),
          size = 1.8, max.overlaps = 30,
          label.padding = 0.1, segment.size = 0.15,
          min.segment.length = 0, fontface = "italic",
          color = "black", fill = "white", alpha = 0.85,
          show.legend = FALSE
        )
    }

    message("Panel j: ", nrow(cv_merge), " TFs plotted (", nrow(label_cv),
            " labeled), ", cor_label_cv)
  } else {
    message("Panel j: no TFs with both chromVAR and dream data")
  }
}

# ==========================================================================
# Save individual panels
# ==========================================================================
EPIG_PANELS <- file.path(FIG3_DIR, "panels")

save_fig(p_f, file.path(EPIG_PANELS, "fig3a_scenic_regulons.pdf"),
         width = fig_half_width, height = 3.5, dpi = 300)
save_fig(p_g, file.path(EPIG_PANELS, "fig3b_tf_target_network.pdf"),
         width = fig_half_width, height = 3.5, dpi = 300)
save_fig(p_h, file.path(EPIG_PANELS, "fig3c_chromvar_motifs.pdf"),
         width = fig_half_width, height = 3.5, dpi = 300)
save_fig(p_i, file.path(EPIG_PANELS, "fig3d_scenic_convergence.pdf"),
         width = fig_half_width, height = 3.5, dpi = 300)
save_fig(p_j, file.path(EPIG_PANELS, "fig3e_chromvar_convergence.pdf"),
         width = fig_half_width, height = 3.5, dpi = 300)

message("Individual epigenomic panels saved to ", EPIG_PANELS)

# ==========================================================================
# Assemble: 2-row x 2-column layout (panels f-i)
# ==========================================================================
fig <- (p_f | p_g) / (p_h | p_i) +
  plot_layout(heights = c(1, 1)) +
  patchwork::plot_annotation(tag_levels = "a",
                             tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_fig(fig, OUT, width = fig_full_width, height = 7)
message("Fig 3 epigenomic panels saved to ", OUT)
