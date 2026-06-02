##############################################################################
# Figure 5: Translational Convergence and Validation
# 6 panels: (a) Clinical drug validation matrix,
#   (b) LINCS L1000 top compounds + network proximity,
#   (c) Spatial SVG validation, (d) Proteomics direction concordance,
#   (e) Network convergence, (f) Atlas-guided target prioritization
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

OUT <- file.path(FIGS_THERA_DIR, "fig5_translation.pdf")

# Pre-initialize all panels with placeholders
p_a <- placeholder("(a) Clinical drug validation")
p_b <- placeholder("(b) LINCS + network proximity")
p_c <- placeholder("(c) Spatial SVG validation")
p_d <- placeholder("(d) Proteomics concordance")
p_e <- placeholder("(e) Network convergence")
p_f <- placeholder("(f) Target prioritization")

# ==========================================================================
# (a) Clinical drug target × evidence source heatmap (magnitude)
# ==========================================================================
drug_val <- load_drug_validation_table()
atlas    <- load_multi_evidence()
if (!is.null(drug_val) && !is.null(atlas)) {

  # --- 1. Deduplicate to unique target genes, concatenate drug names ------
  gene_drugs <- drug_val[, .(
    drugs = paste(unique(sub("\\s*\\(.*", "", drug)), collapse = ", "),
    atlas_support = atlas_support[1],
    # S4 epigenomic flags (from drug_val, not atlas)
    scenic_grn      = any(scenic_grn_target == TRUE, na.rm = TRUE),
    disease_regulon = any(is_disease_regulon_target == TRUE, na.rm = TRUE),
    promoter_acc    = any(mouse_promoter_accessible == TRUE, na.rm = TRUE)
  ), by = target_gene]

  # --- 2. Merge atlas columns --------------------------------------------
  # mr_pval column removed 2026-04-22 — MR ditched from paper.
  atlas_sub <- atlas[human_symbol %in% gene_drugs$target_gene,
    .(human_symbol, dream_padj, dream_logFC,
      best_liver_enzyme_pp4, broadaway_coloc_pp4, pdff_coloc_pp4,
      sceqtl_coloc_pp4_hep, twas_pval,
      essentiality_chronos)]
  # Keep first match per gene (should be unique)
  atlas_sub <- atlas_sub[!duplicated(human_symbol)]

  gene_dt <- merge(gene_drugs, atlas_sub,
                   by.x = "target_gene", by.y = "human_symbol", all.x = TRUE)

  # --- 3. Compute normalised evidence strength per source (0–1) ----------

  # S1: Transcriptomic — padj < 0.05 + |logFC| > 0.5 gate, then score by -log10
  gene_dt[, S1 := ifelse(!is.na(dream_padj) & dream_padj < 0.05 & abs(dream_logFC) > 0.5,
                          pmin(-log10(pmax(dream_padj, 1e-300)) / 30, 1), 0)]

  # S2: Genetic/Causal — max across COLOC PP4s + TWAS score
  # (mr_score removed 2026-04-22 — MR ditched from paper)
  gene_dt[, twas_score := ifelse(!is.na(twas_pval) & twas_pval < 0.05,
                                  1 - twas_pval, 0)]
  gene_dt[, S2 := pmax(
    fifelse(is.na(best_liver_enzyme_pp4), 0, best_liver_enzyme_pp4),
    fifelse(is.na(broadaway_coloc_pp4),   0, broadaway_coloc_pp4),
    fifelse(is.na(pdff_coloc_pp4),        0, pdff_coloc_pp4),
    fifelse(is.na(sceqtl_coloc_pp4_hep),  0, sceqtl_coloc_pp4_hep),
    twas_score)]

  # S3: Essentiality — |chronos| / 2, capped at 1; NA if missing
  gene_dt[, S3 := fifelse(!is.na(essentiality_chronos),
                           pmin(abs(essentiality_chronos) / 2, 1),
                           NA_real_)]

  # S4: Epigenomic — fraction of 3 flags TRUE
  gene_dt[, S4 := (as.numeric(scenic_grn) +
                    as.numeric(disease_regulon) +
                    as.numeric(promoter_acc)) / 3]

  # S5: Spatial — check atlas columns
  gene_dt[, S5 := NA_real_]
  if ("spatial_morans_i" %in% names(atlas)) {
    sp_sub <- atlas[human_symbol %in% gene_dt$target_gene,
      .(human_symbol, spatial_morans_i)]
    sp_sub <- sp_sub[!duplicated(human_symbol)]
    gene_dt <- merge(gene_dt, sp_sub, by.x = "target_gene",
                     by.y = "human_symbol", all.x = TRUE)
    gene_dt[, S5 := fifelse(!is.na(spatial_morans_i),
                             pmin(spatial_morans_i, 1), NA_real_)]
    gene_dt[, spatial_morans_i := NULL]
  }

  # S6: Single-cell — check atlas columns
  gene_dt[, S6 := NA_real_]
  if ("sc_n_celltypes_sig" %in% names(atlas)) {
    sc_sub <- atlas[human_symbol %in% gene_dt$target_gene,
      .(human_symbol, sc_n_celltypes_sig)]
    sc_sub <- sc_sub[!duplicated(human_symbol)]
    gene_dt <- merge(gene_dt, sc_sub, by.x = "target_gene",
                     by.y = "human_symbol", all.x = TRUE)
    gene_dt[, S6 := fifelse(!is.na(sc_n_celltypes_sig) & sc_n_celltypes_sig > 0,
                             pmin(sc_n_celltypes_sig / 5, 1), NA_real_)]
    gene_dt[, sc_n_celltypes_sig := NULL]
  }

  # --- 4. Row ordering: tier then total evidence --------------------------
  support_order <- c("Strong", "Moderate", "Absent")
  gene_dt[, tier_rank := match(atlas_support, support_order)]
  gene_dt[, total_ev  := rowSums(.SD, na.rm = TRUE),
          .SDcols = c("S1", "S2", "S3", "S4", "S5", "S6")]
  setorder(gene_dt, tier_rank, -total_ev)

  gene_dt[, gene_label := paste0(target_gene, " (", drugs, ")")]
  gene_dt[, gene_label := factor(gene_label, levels = rev(gene_label))]

  # --- 5. Melt to long format --------------------------------------------
  # P0-E fix 2026-05-28: source numbering aligned to bayesian_posterior.csv
  # columns (S1 human bulk, S2 genetic, S3 essentiality, S4 epigenomic,
  # S5 spatial, S6 single-cell). Labels kept consistent with panel (f).
  source_labels <- c(S1 = "Human bulk",   S2 = "Genetic",
                     S3 = "Essentiality", S4 = "Epigenomic",
                     S5 = "Spatial",      S6 = "Single-cell")
  tile_long <- melt(gene_dt[, .(gene_label, atlas_support,
                                 S1, S2, S3, S4, S5, S6)],
                    id.vars = c("gene_label", "atlas_support"),
                    variable.name = "source", value.name = "strength")
  tile_long[, source_label := factor(source_labels[as.character(source)],
                                      levels = source_labels)]

  # --- 6. Build heatmap ---------------------------------------------------
  p_heat <- ggplot(tile_long,
                   aes(x = source_label, y = gene_label, fill = strength)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(is.na(strength), "",
                                  sprintf("%.2f", strength)),
                  color = ifelse(!is.na(strength) & strength > 0.5,
                                  "light", "dark")),
              size = 2, show.legend = FALSE) +
    scale_fill_gradient(low = "white", high = masld_colors$up,
                        na.value = "#F0F0F0", limits = c(0, 1),
                        name = "Evidence\nstrength") +
    scale_color_manual(values = c(light = "white", dark = "gray30")) +
    labs(x = NULL, y = NULL,
         title = "Clinical drug targets \u00d7 evidence sources") +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
          axis.text.y = element_text(size = 5.5),
          legend.position = "right",
          legend.key.height = unit(0.6, "cm"),
          legend.key.width  = unit(0.25, "cm"),
          legend.title = element_text(size = 5.5),
          legend.text  = element_text(size = 5))

  # --- 7. Tier annotation strip -------------------------------------------
  tier_colors <- c(Strong = masld_colors$conserved,
                   Moderate = "#FFB300", Absent = masld_colors$ns)
  tier_dt <- unique(tile_long[, .(gene_label, atlas_support)])
  p_tier <- ggplot(tier_dt, aes(x = 1, y = gene_label, fill = atlas_support)) +
    geom_tile(color = "white", linewidth = 0.5, width = 0.8) +
    scale_fill_manual(values = tier_colors, name = "Atlas\nsupport") +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(legend.position = "right",
          legend.key.size = unit(0.25, "cm"),
          legend.title = element_text(size = 5.5),
          legend.text  = element_text(size = 5),
          axis.text = element_blank(),
          axis.ticks = element_blank())

  p_a <- p_heat + p_tier + plot_layout(widths = c(12, 1))
}

# ==========================================================================
# (b) LINCS L1000 top compounds + network proximity
# ==========================================================================
lincs <- load_lincs_compounds()
netprox <- load_network_proximity()

if (!is.null(lincs)) {
  # Aggregate top MOA classes by reversal score or count
  if ("moa" %in% names(lincs)) {
    moa_dt <- lincs[!is.na(moa) & nchar(moa) > 0]
    if (nrow(moa_dt) > 0) {
      # Count compounds per MOA
      moa_counts <- moa_dt[, .N, by = moa][order(-N)]
      top_moa <- head(moa_counts, 15)
      top_moa[, moa := factor(moa, levels = rev(top_moa$moa))]

      p_b_left <- ggplot(top_moa, aes(x = N, y = moa)) +
        geom_col(fill = masld_colors$up, width = 0.6) +
        geom_text(aes(label = N), hjust = -0.2, size = 1.8) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
        labs(x = "Compounds", y = NULL,
             title = "Top LINCS MOA classes") +
        theme_masld() +
        theme(axis.text.y = element_text(size = 4.5))

      # Combine with network proximity if available
      if (!is.null(netprox) && nrow(netprox) > 0) {
        # Show distribution of z-scores
        z_col <- intersect(c("z_combined", "z_closest", "z_shortest"),
                           names(netprox))[1]
        if (!is.na(z_col)) {
          netprox_plot <- copy(netprox)
          netprox_plot[, z_val := get(z_col)]
          netprox_plot[, sig := z_val < -2]

          p_b_right <- ggplot(netprox_plot, aes(x = z_val, fill = sig)) +
            geom_histogram(bins = 40, color = "white", linewidth = 0.1) +
            geom_vline(xintercept = -2, linetype = "dashed",
                       linewidth = 0.3, color = masld_colors$mr) +
            scale_fill_manual(values = c("TRUE" = masld_colors$up,
                                         "FALSE" = masld_colors$ns),
                              guide = "none") +
            annotate("text", x = -2, y = Inf,
                     label = "z < -2",
                     hjust = 1.1, vjust = 1.5, size = 2,
                     color = masld_colors$mr, fontface = "italic") +
            labs(x = paste0("Network proximity (", z_col, ")"),
                 y = "Drugs",
                 title = "Drug-disease network proximity") +
            theme_masld()

          p_b <- p_b_left | p_b_right
        } else {
          p_b <- p_b_left
        }
      } else {
        p_b <- p_b_left
      }
    }
  }
} else if (!is.null(netprox) && nrow(netprox) > 0) {
  # LINCS not available, show network proximity alone
  z_col <- intersect(c("z_combined", "z_closest", "z_shortest"),
                     names(netprox))[1]
  if (!is.na(z_col)) {
    netprox_plot <- copy(netprox)
    netprox_plot[, z_val := get(z_col)]
    netprox_plot[, sig := z_val < -2]

    p_b <- ggplot(netprox_plot, aes(x = z_val, fill = sig)) +
      geom_histogram(bins = 40, color = "white", linewidth = 0.1) +
      geom_vline(xintercept = -2, linetype = "dashed",
                 linewidth = 0.3, color = masld_colors$mr) +
      scale_fill_manual(values = c("TRUE" = masld_colors$up,
                                   "FALSE" = masld_colors$ns),
                        guide = "none") +
      labs(x = paste0("Network proximity (", z_col, ")"),
           y = "Drugs",
           title = "Drug-disease network proximity") +
      theme_masld()
  }
}

# ==========================================================================
# (c) Spatial SVG validation — Conserved enrichment in SVGs
# ==========================================================================
spatial_enr <- load_spatial_enrichment()
if (!is.null(spatial_enr)) {
  spatial_enr <- as.data.table(spatial_enr)

  # Remove row index column if present
  if ("V1" %in% names(spatial_enr) || "" %in% names(spatial_enr)) {
    idx_col <- intersect(c("V1", ""), names(spatial_enr))
    if (length(idx_col) > 0) spatial_enr[, (idx_col) := NULL]
  }

  if ("odds_ratio" %in% names(spatial_enr) && "gene_set" %in% names(spatial_enr)) {
    # Clean gene set names
    spatial_enr[, gene_set_label := gsub("_", " ", gene_set)]

    # Add significance annotation
    p_col <- intersect(c("fisher_padj_bh", "fisher_pval", "p_value"),
                       names(spatial_enr))[1]
    if (!is.na(p_col)) {
      spatial_enr[, sig := get(p_col) < 0.05]
      spatial_enr[, sig_label := ifelse(sig,
        paste0("p=", formatC(get(p_col), format = "e", digits = 1)),
        "NS")]
    } else {
      spatial_enr[, sig := TRUE]
      spatial_enr[, sig_label := ""]
    }

    spatial_enr[, gene_set_label := factor(gene_set_label,
      levels = rev(spatial_enr[order(odds_ratio)]$gene_set_label))]

    p_c <- ggplot(spatial_enr,
                  aes(x = odds_ratio, y = gene_set_label, fill = sig)) +
      geom_col(width = 0.6) +
      geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3,
                 color = "gray50") +
      geom_text(aes(label = paste0("OR=", round(odds_ratio, 2),
                                   " ", sig_label)),
                hjust = -0.05, size = 1.8, color = "gray30") +
      scale_fill_manual(values = c("TRUE" = masld_colors$conserved,
                                   "FALSE" = masld_colors$ns),
                        guide = "none") +
      scale_x_continuous(expand = expansion(mult = c(0, 0.4))) +
      labs(x = "Odds ratio (SVG enrichment)",
           y = NULL,
           title = "Spatial SVG enrichment in gene sets") +
      theme_masld()
  }
}

# ==========================================================================
# (d) Proteomics direction concordance — transcript vs protein logFC
# ==========================================================================
prot_conc <- load_protein_concordance_v2()
if (!is.null(prot_conc)) {
  # Need protein_logFC and dream_logFC (or equivalent transcript LFC)
  prot_lfc_col <- intersect(c("protein_logFC", "protein_lfc", "logFC_protein"),
                             names(prot_conc))[1]
  tx_lfc_col <- intersect(c("dream_logFC", "transcript_logFC", "logFC_transcript"),
                           names(prot_conc))[1]

  if (!is.na(prot_lfc_col) && !is.na(tx_lfc_col)) {
    prot_plot <- copy(prot_conc)
    prot_plot[, prot_lfc := get(prot_lfc_col)]
    prot_plot[, tx_lfc := get(tx_lfc_col)]
    prot_plot <- prot_plot[!is.na(prot_lfc) & !is.na(tx_lfc)]

    if (nrow(prot_plot) > 0) {
      # Direction concordance flag
      conc_col <- intersect(c("direction_concordant", "concordant"),
                             names(prot_plot))[1]
      if (!is.na(conc_col)) {
        prot_plot[, concordant := as.logical(get(conc_col))]
      } else {
        prot_plot[, concordant := sign(prot_lfc) == sign(tx_lfc)]
      }

      # Spearman rho
      rho_prot <- cor(prot_plot$tx_lfc, prot_plot$prot_lfc,
                      method = "spearman", use = "complete.obs")
      n_prot <- nrow(prot_plot)
      pct_conc <- round(100 * mean(prot_plot$concordant, na.rm = TRUE), 1)

      # Dataset faceting if available
      has_dataset <- "dataset" %in% names(prot_plot)

      p_d <- ggplot(prot_plot, aes(x = tx_lfc, y = prot_lfc,
                                    color = concordant)) +
        rasterize_layer(
          geom_point(alpha = 0.5, size = 0.6, shape = 16)
        ) +
        geom_hline(yintercept = 0, linewidth = 0.2, color = "gray60",
                   linetype = "dashed") +
        geom_vline(xintercept = 0, linewidth = 0.2, color = "gray60",
                   linetype = "dashed") +
        geom_smooth(method = "lm", se = FALSE, linewidth = 0.4,
                    color = "gray30", linetype = "solid") +
        scale_color_manual(
          values = c("TRUE" = proteomics_colors[["concordant"]],
                     "FALSE" = proteomics_colors[["discordant"]]),
          labels = c("TRUE" = "Concordant", "FALSE" = "Discordant"),
          name = "Direction") +
        annotate("text", x = Inf, y = -Inf,
                 label = paste0("rho = ", round(rho_prot, 3),
                                "\nn = ", comma(n_prot),
                                "\n", pct_conc, "% concordant"),
                 hjust = 1.1, vjust = -0.3, size = 2, fontface = "italic",
                 lineheight = 0.85) +
        labs(x = "Transcript logFC (integrated)",
             y = "Protein logFC",
             title = "Protein-transcript direction concordance") +
        theme_masld() +
        theme(legend.position = "bottom",
              legend.key.size = unit(0.2, "cm"))

      if (has_dataset) {
        p_d <- p_d + facet_wrap(~dataset, scales = "free")
      }
    }
  }
}

# ==========================================================================
# (e) Network convergence — 123 genes with 7/7 source convergence
# ==========================================================================
net_modules_path <- file.path(ME, "network_modules.csv")
net_prop_path <- file.path(ME, "network_propagation_scores.csv")

if (file.exists(net_prop_path)) {
  net_prop <- fread(net_prop_path)

  if ("network_convergence" %in% names(net_prop) &&
      "sources_active" %in% names(net_prop)) {
    # Compute distribution of network convergence scores
    conv_dist <- net_prop[, .N, by = network_convergence][order(network_convergence)]
    conv_dist[, pct := 100 * N / sum(N)]

    # Highlight 7/7 convergence
    conv_dist[, highlight := network_convergence == 7]

    p_e <- ggplot(conv_dist,
                  aes(x = factor(network_convergence), y = N,
                      fill = highlight)) +
      geom_col(width = 0.7) +
      geom_text(aes(label = N), vjust = -0.3, size = 1.8) +
      scale_fill_manual(values = c("TRUE" = masld_colors$conserved,
                                   "FALSE" = masld_colors$ns),
                        guide = "none") +
      labs(x = "Network convergence score (out of 7 sources)",
           y = "Number of genes",
           title = "Multi-source network convergence") +
      theme_masld()

    # Try to add enrichment annotation for 7/7 genes
    n_7of7 <- conv_dist[network_convergence == 7, N]
    if (length(n_7of7) > 0 && n_7of7 > 0) {
      # Check if Conserved enrichment data is available
      cc_genes <- load_concordance_atlas()
      if (!is.null(cc_genes)) {
        cc_symbols <- cc_genes[primary_category == "Conserved" |
                               dvc_category == "Conserved",
                               unique(human_symbol)]
        gene_col_net <- intersect(c("human_symbol", "symbol"),
                                   names(net_prop))[1]
        if (!is.na(gene_col_net)) {
          conv7 <- net_prop[network_convergence == 7]
          n_conv7_cc <- sum(conv7[[gene_col_net]] %in% cc_symbols)
          pct_cc <- round(100 * n_conv7_cc / nrow(conv7), 1)
          p_e <- p_e +
            annotate("text", x = Inf, y = Inf,
                     label = paste0(n_7of7, " genes at 7/7\n",
                                    pct_cc, "% Conserved\n",
                                    "(OR=4.85 for CC)"),
                     hjust = 1.1, vjust = 1.3, size = 2,
                     color = masld_colors$conserved,
                     lineheight = 0.85, fontface = "italic")
        }
      }
    }
  }
}

# ==========================================================================
# (f) Atlas-guided target prioritization — top convergence-score genes
# ==========================================================================
bayes_path <- file.path(ME, "bayesian_posterior.csv")
if (file.exists(bayes_path)) {
  bayes <- fread(bayes_path)

  if ("posterior_odds" %in% names(bayes) && "human_symbol" %in% names(bayes)) {
    # Top 30 genes by posterior odds
    top_genes <- head(bayes[order(-posterior_odds)], 30)

    # Identify source contributions (delta columns)
    # P0-E fix 2026-05-28: source numbering aligned to bayesian_posterior.csv
    # columns (S1 human bulk, S2 genetic, S3 essentiality, S4 epigenomic,
    # S5 spatial, S6 single-cell). File has only S1-S6 (no S7).
    delta_cols <- grep("^delta_S[1-6]", names(top_genes), value = TRUE)

    if (length(delta_cols) > 0) {
      # Melt for stacked bar chart showing evidence decomposition
      top_genes[, gene_rank := .I]
      top_genes[, gene_label := factor(human_symbol,
        levels = rev(top_genes$human_symbol))]

      bar_long <- melt(top_genes,
                       id.vars = c("human_symbol", "gene_label",
                                   "posterior_odds", "gene_rank"),
                       measure.vars = delta_cols,
                       variable.name = "source", value.name = "delta")

      # Clean source names
      # P0-E fix 2026-05-28: source numbering aligned to bayesian_posterior.csv
      # columns (S1 human bulk, S2 genetic, S3 essentiality, S4 epigenomic,
      # S5 spatial, S6 single-cell). Previous mapping mislabeled S2 onward
      # ("Mouse bulk"/shifted) and referenced a phantom S7. Labels kept
      # consistent with panel (a).
      bar_long[, source_label := fcase(
        grepl("S1", source), "S1: Human bulk",
        grepl("S2", source), "S2: Genetic",
        grepl("S3", source), "S3: Essentiality",
        grepl("S4", source), "S4: Epigenomic",
        grepl("S5", source), "S5: Spatial",
        grepl("S6", source), "S6: Single-cell",
        default = as.character(source)
      )]

      # Source colors
      source_colors <- c(
        "S1: Human bulk"   = masld_colors$deg,
        "S2: Genetic"      = masld_colors$mr,
        "S3: Essentiality" = masld_colors$gwas,
        "S4: Epigenomic"   = masld_colors$human_enriched,
        "S5: Spatial"      = masld_colors$conserved,
        "S6: Single-cell"  = masld_colors$sex
      )

      # Only show positive deltas (evidence supporting)
      bar_long[delta < 0, delta := 0]

      p_f <- ggplot(bar_long,
                    aes(x = delta, y = gene_label, fill = source_label)) +
        geom_col(width = 0.7, position = "stack") +
        scale_fill_manual(values = source_colors, name = "Source") +
        scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
        labs(x = "Evidence contribution (delta posterior odds)",
             y = NULL,
             title = "Top-ranked targets: evidence decomposition") +
        theme_masld() +
        theme(axis.text.y = element_text(size = 4.5),
              legend.position = "bottom",
              legend.key.size = unit(0.2, "cm"),
              legend.text = element_text(size = 5),
              legend.margin = margin(0, 0, 0, 0))

      # Highlight known drug targets
      known_targets <- c("THRB", "NR1H4", "PPARA", "PPARG", "PPARD",
                          "GLP1R", "FXR", "DGAT2")
      top_known <- top_genes[human_symbol %in% known_targets]
      if (nrow(top_known) > 0) {
        p_f <- p_f +
          annotate("text", x = Inf, y = -Inf,
                   label = paste0("Known targets in top 30:\n",
                                  paste(top_known$human_symbol,
                                        collapse = ", ")),
                   hjust = 1.1, vjust = -0.3, size = 1.8,
                   color = masld_colors$mr, fontface = "italic",
                   lineheight = 0.85)
      }
    } else {
      # Simpler waterfall without decomposition
      top_genes[, gene_label := factor(human_symbol,
        levels = rev(top_genes$human_symbol))]

      p_f <- ggplot(top_genes,
                    aes(x = posterior_odds, y = gene_label)) +
        geom_col(fill = masld_colors$up, width = 0.7) +
        labs(x = "Posterior odds",
             y = NULL,
             title = "Top-ranked targets (convergence score)") +
        theme_masld() +
        theme(axis.text.y = element_text(size = 4.5))
    }
  }
}

# ==========================================================================
# Compose
# Row 1: (a) drug validation matrix | (b) LINCS + proximity
# Row 2: (c) spatial SVG | (d) proteomics concordance
# Row 3: (e) network convergence | (f) target prioritization
# ==========================================================================
fig5 <- (p_a | p_b) /
        (p_c | p_d) /
        (p_e | p_f)

fig5 <- fig5 +
  plot_layout(heights = c(1.2, 0.9, 1.0)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

save_fig_tall(fig5, OUT, width = fig_full_width, height = 9.5)
cat("Saved:", OUT, "\n")
