##############################################################################
# Figure 3: Disease Progression Atlas (8 panels)
#
# Multi-scale progression map for MASLD, synthesizing Phase I (trajectory,
# transitions, bifurcation) and Phase II (cell-type deconvolution) results.
#
# Layout (4 rows x 2 columns):
#   Row 1: (a) Pseudotime trajectory (UMAP colored by pseudotime + fibrosis)
#           (b) Transition DEG cascade (heatmap of tau-specific gene programs)
#   Row 2: (c) Bifurcation analysis (divergence score along pseudotime)
#           (d) Fate probabilities (P(F4) density for S1 vs S2)
#   Row 3: (e) Cell-type composition across pseudotime (area chart)
#           (f) Cell-type attribution heatmap (DEGs per cell type x transition)
#   Row 4: (g) Top divergence genes (dot plot with ECM/fibrosis annotation)
#           (h) Exemplar gene expression across pseudotime (COL1A1 etc.)
#
# Output: figures/main/fig3_deconfounding/fig3_progression.pdf
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(ggrepel)
  library(viridis)
  library(RColorBrewer)
  library(uwot)        # for UMAP of VAE embeddings
})

set.seed(42)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PROG_DIR <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression")
STAGING_DIR <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/staging_classifier")

OUT_DIR <- FIGS02_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(OUT_DIR, "figS02_progression_main.pdf")

# NMF-dependent subtype panels (c, d) disabled 2026-04-24: reviewer audit and
# pending seed-robustness re-run indicate S1/S2 subtype separation is not stable.
# See memory/reviewer_audit_2026_04_23.md and CLAUDE.md NMF-subtyping note.
ENABLE_NMF_SUBTYPE_PANELS <- FALSE

# ---------------------------------------------------------------------------
# Fibrosis stage color palette (blue-to-dark-magenta gradient, 6 levels)
# Uses existing palette from publication_theme.R + a gray for NA/Control
# ---------------------------------------------------------------------------
fib_colors <- c(
  "F0"      = "#E3F2FD",  # Lightest blue
  "F1"      = "#90CAF9",
  "F2"      = "#42A5F5",
  "F3"      = "#1565C0",
  "F4"      = "#0D47A1",  # Darkest blue
  "Control" = "#9E9E9E"
)

# Subtype colors (Metabolic = resolving trajectory, Fibrogenic = progressing trajectory)
subtype_colors <- c(
  "Metabolic"      = "#42A5F5",  # Blue (resolving) — legacy name for panel (c)/(g)
  "Non-Fibrogenic" = "#42A5F5",  # Blue — k-program lumped non-Fibrotic label
  "Fibrogenic"     = "#C9265E"   # Liang deep magenta (progressing)
)

# Transition labels for clean display
transition_labels <- c(
  "F0_to_F1" = "F0 \u2192 F1",
  "F1_to_F2" = "F1 \u2192 F2",
  "F2_to_F3" = "F2 \u2192 F3",
  "F3_to_F4" = "F3 \u2192 F4"
)

# Cell-type colors for composition panel (major types only)
comp_colors <- c(
  "Hepatocyte"    = "#0D47A1",
  "Stellate"      = "#F57F17",
  "Macrophage"    = "#C2185B",
  "Endothelial"   = "#2E7D32",
  "Cholangiocyte" = "#1565C0",
  "Monocyte"      = "#E91E63",
  "DC"            = "#AD1457",
  "Neutrophil"    = "#FF6F00",
  "T_cell"        = "#7B1FA2",
  "NK_cell"       = "#00695C",
  "B_cell"        = "#9C27B0",
  "Plasma_cell"   = "#5D4037",
  "Other_immune"  = "#78909C"
)

# Pre-initialize all panels with placeholders
p_a <- placeholder("(a) Pseudotime trajectory UMAP")
p_b <- placeholder("(b) Transition DEG cascade")
p_c <- placeholder("(c) Bifurcation analysis")
p_d <- placeholder("(d) Fate probabilities")
p_e <- placeholder("(e) Cell-type composition")
p_f <- placeholder("(f) Cell-type attribution")
p_g <- placeholder("(g) Divergence genes")
p_h <- placeholder("(h) Exemplar gene trajectories")

# ===========================================================================
# Panel (b): Transition DEG cascade — heatmap of transition-specific genes
#            ranked by tau specificity index, showing logFC across 4 fibrosis
#            transitions. Top 15 genes per peak transition.
# ===========================================================================
tryCatch({
  message("Panel (b): Transition DEG cascade heatmap...")

  tau <- fread(file.path(PROG_DIR, "transition_tau_index.csv"))
  fib_de <- fread(file.path(PROG_DIR, "transition_fib_dream_results.csv"))

  # Keep only fibrosis transitions
  fib_transitions <- c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4")
  fib_de <- fib_de[transition %in% fib_transitions]

  # Add gene symbols via shared mapping
  fib_de <- add_symbols(fib_de, "gene")
  tau <- add_symbols(tau, "gene")

  # Select high-tau genes with significant peak transition
  tau_sig <- tau[tau > 0.8 & peak_padj < 0.05 & peak_transition %in% fib_transitions]

  # Top 15 per peak transition (by max_abs_t, then tau)
  tau_top <- tau_sig[, .SD[order(-max_abs_t)][1:min(.N, 15)],
                     by = peak_transition]

  # Build LFC matrix for selected genes
  sel_genes <- unique(tau_top$gene)
  mat_dt <- fib_de[gene %in% sel_genes, .(gene, transition, logFC)]
  mat_dt <- add_symbols(mat_dt, "gene")

  # Build a gene-level display name table (one row per gene) to avoid

  # cartesian joins when merging back to multi-transition mat_dt
  gene_sym <- unique(mat_dt[, .(gene, symbol)])
  gene_sym[, display := fifelse(
    duplicated(symbol) | symbol == sub("\\..*", "", gene),
    paste0(symbol, " (", sub("\\..*", "", gene), ")"),
    symbol
  )]

  # Create ordering: group by peak transition, then by tau within group
  gene_order_dt <- tau_top[, .(gene, peak_transition, tau, max_abs_t)]
  gene_order_dt[, peak_transition := factor(peak_transition, levels = fib_transitions)]
  gene_order_dt <- gene_order_dt[order(peak_transition, -tau)]
  gene_order_dt <- merge(gene_order_dt,
                         gene_sym[, .(gene, display)],
                         by = "gene", all.x = TRUE)
  gene_order_dt[is.na(display), display := sub("\\..*", "", gene)]

  mat_dt <- merge(mat_dt, gene_sym[, .(gene, display)],
                  by = "gene", all.x = TRUE)
  mat_dt[is.na(display), display := sub("\\..*", "", gene)]

  mat_dt[, display := factor(display, levels = rev(gene_order_dt$display))]
  mat_dt[, transition := factor(transition, levels = fib_transitions,
                                labels = transition_labels[fib_transitions])]

  # Cap logFC for visual clarity
  mat_dt[, logFC_cap := pmin(pmax(logFC, -2), 2)]

  # Add peak transition annotation
  mat_dt <- merge(mat_dt, gene_order_dt[, .(gene, peak_transition)],
                  by = "gene", all.x = TRUE)

  # Load OT transport costs for annotation
  ot_plans <- fread(file.path(PROG_DIR, "transition_transport_plans.csv"))
  ot_fib <- ot_plans[transition %in% fib_transitions]
  ot_fib[, transition_disp := factor(transition, levels = fib_transitions,
                                     labels = transition_labels[fib_transitions])]
  ot_fib[, cost_label := paste0("OT=", round(total_ot_cost, 2))]

  p_b <- ggplot(mat_dt, aes(x = transition, y = display, fill = logFC_cap)) +
    geom_tile(color = "white", linewidth = 0.2) +
    # Annotate OT transport costs above each column
    geom_text(data = ot_fib, inherit.aes = FALSE,
              aes(x = transition_disp, y = Inf, label = cost_label),
              vjust = -0.3, size = 1.8, color = "#880E4F", fontface = "bold") +
    scale_fill_gradient2(
      low = masld_colors$down, mid = "white", high = masld_colors$up,
      midpoint = 0, limits = c(-2, 2), oob = squish,
      name = "logFC"
    ) +
    coord_cartesian(clip = "off") +
    labs(title = "Transition-specific gene programs", x = NULL, y = NULL) +
    theme_masld() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(size = 5),
      legend.key.width = unit(0.3, "cm"),
      legend.key.height = unit(1.5, "cm"),
      plot.margin = margin(10, 3, 3, 3)
    )

  save_fig(p_b, file.path(OUT_DIR, "panel_b_cascade.pdf"),
           width = fig_half_width, height = 4.5)

}, error = function(e) message("Panel (b) failed: ", e$message))

# ===========================================================================
# Panel (c): Bifurcation analysis — divergence score along pseudotime
#            showing where S1 and S2 subtypes split
# ===========================================================================
if (ENABLE_NMF_SUBTYPE_PANELS) tryCatch({
  message("Panel (c): Bifurcation divergence score...")

  bif <- fread(file.path(PROG_DIR, "bifurcation_analysis.csv"))

  # Identify bifurcation point (max divergence in early pseudotime < 0.3)
  bif_early <- bif[pseudotime_center < 0.3]
  bif_peak <- bif_early[which.max(divergence_score)]

  # Late re-divergence
  bif_late <- bif[pseudotime_center > 0.7]
  bif_late_peak <- bif_late[which.max(divergence_score)]

  p_c <- ggplot(bif, aes(x = pseudotime_center, y = divergence_score)) +
    # Ribbon for uncertainty visual
    geom_ribbon(aes(ymin = 0, ymax = divergence_score),
                fill = "#E8EAF6", alpha = 0.5) +
    geom_line(color = "#1A237E", linewidth = 0.8) +
    geom_point(data = bif_peak, color = masld_colors$up, size = 2.5) +
    annotate("text", x = bif_peak$pseudotime_center + 0.08,
             y = bif_peak$divergence_score,
             label = paste0("Early\nbifurcation\n(t=",
                            round(bif_peak$pseudotime_center, 2), ")"),
             size = 2, color = masld_colors$up, hjust = 0, vjust = 0.5) +
    geom_point(data = bif_late_peak, color = masld_colors$fibrosis, size = 2.5) +
    annotate("text", x = bif_late_peak$pseudotime_center - 0.08,
             y = bif_late_peak$divergence_score,
             label = paste0("Late\nre-divergence\n(t=",
                            round(bif_late_peak$pseudotime_center, 2), ")"),
             size = 2, color = masld_colors$fibrosis, hjust = 1, vjust = 0.5) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray60",
               linewidth = 0.3) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0, 0.02)) +
    labs(title = "Metabolic/Fibrogenic subtype divergence", x = "Consensus pseudotime",
         y = "Divergence score") +
    theme_masld()

  save_fig(p_c, file.path(OUT_DIR, "panel_c_bifurcation.pdf"),
           width = fig_half_width, height = 2.5)

}, error = function(e) message("Panel (c) failed: ", e$message))

# ===========================================================================
# Panel (d): Fate probabilities — density of P(F4) for Metabolic vs Fibrogenic subtypes
# ===========================================================================
if (ENABLE_NMF_SUBTYPE_PANELS) tryCatch({
  message("Panel (d): Fate probability densities...")

  fate <- fread(file.path(PROG_DIR, "fate_probabilities.csv"))
  # k-program refactor (2026-04-20): legacy nmf_subtype S1/S2 in this CSV is a
  # convenience binary where S2 = dominant-program-is-Fibrotic, S1 = all other
  # programs. Relabel for plot; the Non-Fibrogenic label captures Metabolic +
  # Inflammatory + Sex-* + HCC-like samples lumped together.
  fate[, subtype := fifelse(subtype == "S1", "Non-Fibrogenic", "Fibrogenic")]
  fate[, subtype := factor(subtype, levels = c("Non-Fibrogenic", "Fibrogenic"))]

  # Wilcoxon test for annotation
  wtest <- wilcox.test(fate[subtype == "Non-Fibrogenic"]$fate_prob_F4,
                       fate[subtype == "Fibrogenic"]$fate_prob_F4)
  pval_label <- ifelse(wtest$p.value < 2.2e-16, "P < 2.2e-16",
                       paste0("P = ", signif(wtest$p.value, 2)))

  # Median lines
  medians <- fate[, .(med = median(fate_prob_F4, na.rm = TRUE)), by = subtype]

  p_d <- ggplot(fate, aes(x = fate_prob_F4, fill = subtype, color = subtype)) +
    geom_density(alpha = 0.35, linewidth = 0.5) +
    geom_vline(data = medians, aes(xintercept = med, color = subtype),
               linetype = "dashed", linewidth = 0.5) +
    scale_fill_manual(values = subtype_colors, name = "Subtype") +
    scale_color_manual(values = subtype_colors, name = "Subtype") +
    annotate("text", x = 0.65, y = Inf, label = pval_label,
             vjust = 1.5, size = 2.2, fontface = "italic") +
    labs(title = "Fate probability by subtype",
         x = "P(progression to F4)", y = "Density") +
    theme_masld() +
    theme(legend.position = c(0.85, 0.8))

  save_fig(p_d, file.path(OUT_DIR, "panel_d_fate.pdf"),
           width = fig_half_width, height = 2.5)

}, error = function(e) message("Panel (d) failed: ", e$message))

# ===========================================================================
# Panel (e): Cell-type composition along pseudotime (smoothed area chart)
#            BayesPrism proportions binned along consensus pseudotime
# ===========================================================================
tryCatch({
  message("Panel (e): Cell-type composition across pseudotime...")

  bp <- fread(file.path(PROG_DIR,
    "cibersortx_celltype_expression/bayesprism_proportions.csv"))
  pst <- fread(file.path(PROG_DIR, "consensus_pseudotime.csv"))

  # Merge pseudotime
  bp <- merge(bp, pst[, .(sample_id, pseudotime_consensus)],
              by = "sample_id", all.x = TRUE)
  bp <- bp[!is.na(pseudotime_consensus)]

  # Bin pseudotime into 20 bins
  bp[, pt_bin := cut(pseudotime_consensus, breaks = 20, labels = FALSE)]
  bp[, pt_center := (as.numeric(pt_bin) - 0.5) / 20]

  # Melt to long format
  ct_cols <- setdiff(names(bp), c("sample_id", "pseudotime_consensus",
                                   "pt_bin", "pt_center"))
  bp_long <- melt(bp, id.vars = c("sample_id", "pt_center"),
                  measure.vars = ct_cols, variable.name = "cell_type",
                  value.name = "proportion")

  # Mean proportion per bin
  bp_mean <- bp_long[, .(mean_prop = mean(proportion, na.rm = TRUE)),
                     by = .(pt_center, cell_type)]

  # Order cell types by overall abundance (hepatocyte on bottom)
  ct_order <- bp_mean[, .(total = sum(mean_prop)), by = cell_type][order(-total)]
  bp_mean[, cell_type := factor(cell_type, levels = rev(ct_order$cell_type))]

  # Map any missing cell types to gray
  fill_vals <- comp_colors
  missing_ct <- setdiff(levels(bp_mean$cell_type), names(fill_vals))
  if (length(missing_ct) > 0) {
    fill_vals[missing_ct] <- "#BDBDBD"
  }

  # Exclude hepatocytes to show non-hepatocyte dynamics (hepatocyte ~82% dominates)
  bp_nonhep <- bp_mean[cell_type != "Hepatocyte"]
  # Renormalize to sum to 1 within each bin for area chart
  bp_nonhep[, total_bin := sum(mean_prop), by = pt_center]
  bp_nonhep[, prop_renorm := mean_prop / total_bin]

  # Stellate expansion annotation: compute from original (non-renormalized) means
  stellate_early <- bp_mean[cell_type == "Stellate" & pt_center <= 0.1, mean(mean_prop)]
  stellate_late  <- bp_mean[cell_type == "Stellate" & pt_center >= 0.9, mean(mean_prop)]
  stellate_fold  <- stellate_late / stellate_early

  p_e <- ggplot(bp_nonhep, aes(x = pt_center, y = prop_renorm, fill = cell_type)) +
    geom_area(position = "stack", alpha = 0.85, linewidth = 0.15, color = "white") +
    scale_fill_manual(values = fill_vals, name = "Cell type") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       expand = c(0, 0)) +
    scale_y_continuous(labels = percent, expand = c(0, 0)) +
    annotate("text", x = 0.75, y = 0.95,
             label = paste0("Stellate: ",
                            round(stellate_early * 100, 1), "% -> ",
                            round(stellate_late * 100, 1), "%\n(",
                            round(stellate_fold, 1), "x expansion)"),
             size = 2, fontface = "bold", color = "#F57F17", hjust = 0.5) +
    labs(title = "Non-hepatocyte composition across progression",
         x = "Consensus pseudotime", y = "Proportion (non-hepatocyte)") +
    theme_masld() +
    theme(legend.position = "right",
          legend.text = element_text(size = 5),
          legend.key.size = unit(0.25, "cm"))

  save_fig(p_e, file.path(OUT_DIR, "panel_e_composition.pdf"),
           width = fig_half_width + 1, height = 2.8)

}, error = function(e) message("Panel (e) failed: ", e$message))

# ===========================================================================
# Panel (f): Cell-type x transition attribution heatmap
#            Number of DEGs per cell type per fibrosis transition
# ===========================================================================
tryCatch({
  message("Panel (f): Cell-type attribution heatmap...")

  ct_summ <- fread(file.path(PROG_DIR, "celltype_transition_summary.csv"))

  # Keep only fibrosis transitions
  ct_fib <- ct_summ[transition_type == "fibrosis" &
                     transition %in% names(transition_labels)]

  # Log-transform DEG counts for better visual dynamic range
  ct_fib[, log_degs := log10(n_degs_01 + 1)]

  # Clean labels
  ct_fib[, transition_disp := factor(transition, levels = names(transition_labels),
                                     labels = transition_labels)]
  ct_fib[, cell_type := factor(cell_type)]

  # Order cell types by total DEGs across transitions
  ct_totals <- ct_fib[, .(total = sum(n_degs_01)), by = cell_type][order(-total)]
  ct_fib[, cell_type := factor(cell_type, levels = rev(ct_totals$cell_type))]

  p_f <- ggplot(ct_fib, aes(x = transition_disp, y = cell_type, fill = log_degs)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(n_degs_01 > 0,
                                  formatC(n_degs_01, format = "d", big.mark = ","),
                                  "")),
              size = 1.8, color = "black") +
    scale_fill_gradient(
      low = "white", high = "#0D47A1",
      name = expression(log[10]*"(DEGs+1)"),
      limits = c(0, NA)
    ) +
    labs(title = "Cell-type transition programs",
         x = NULL, y = NULL) +
    theme_masld() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(size = 5.5),
      legend.key.width = unit(0.3, "cm"),
      legend.key.height = unit(1.5, "cm")
    )

  save_fig(p_f, file.path(OUT_DIR, "panel_f_attribution.pdf"),
           width = fig_half_width, height = 3.5)

}, error = function(e) message("Panel (f) failed: ", e$message))

# ===========================================================================
# Panel (f2): Bulk-masked exemplars — bulk vs single-cell-type logFC
#             Genes that are invisible (or weak) in bulk RNA-seq but
#             strongly DE in one cell type's pseudobulk.
#             Companion panel to (f) showing per-gene cases behind the
#             cell-type x transition DEG counts.
# ===========================================================================
tryCatch({
  message("Panel (f2): Bulk-masked exemplar genes...")

  bulk_path <- file.path(BASE,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration",
    "canonical_deg_results.csv")
  pb_dir <- file.path(BASE,
    "Analysis/SingleCell/results_gpu_v2/pseudobulk_de")

  # Curated 8-gene panel — see docs/paper_outline.md / memory rationale.
  # category: "masked" = bulk fails padj<0.05 OR |logFC|<0.5;
  #           "amplified" = bulk significant but CT logFC much larger.
  exemplars <- data.table::data.table(
    gene = c("LOX", "FBLN1", "CLEC4G", "FCN3", "SPP2",
             "ACKR1", "THY1", "COL1A1"),
    cell_type = c("Fibroblasts", "Fibroblasts", "Endothelial_cells",
                  "Endothelial_cells", "Cholangiocytes",
                  "Endothelial_cells", "Fibroblasts", "Fibroblasts"),
    cell_label = c("Fibroblast", "Fibroblast", "Endothelial",
                   "Endothelial", "Cholangiocyte",
                   "Endothelial", "Fibroblast", "Fibroblast"),
    category = c("Bulk-masked", "Bulk-masked", "Bulk-masked",
                 "Bulk-masked", "Bulk-masked",
                 "Amplified in CT", "Amplified in CT", "Amplified in CT")
  )

  # Bulk values
  bulk <- data.table::fread(bulk_path,
    select = c("symbol", "logFC", "padj"))
  bulk_sub <- bulk[symbol %in% exemplars$gene,
                   .(gene = symbol, bulk_logFC = logFC, bulk_padj = padj)]

  # CT pseudobulk values — one file per cell type, so loop
  ct_rows <- list()
  for (ct in unique(exemplars$cell_type)) {
    pb_file <- file.path(pb_dir, paste0(ct, "_de.csv"))
    if (!file.exists(pb_file)) next
    pb <- data.table::fread(pb_file, select = c("gene", "logFC", "padj"))
    g_keep <- exemplars[cell_type == ct, gene]
    sub <- pb[gene %in% g_keep]
    sub[, cell_type := ct]
    ct_rows[[ct]] <- sub[, .(gene, cell_type,
                             ct_logFC = logFC, ct_padj = padj)]
  }
  ct_dt <- data.table::rbindlist(ct_rows)

  panel_dt <- merge(exemplars, ct_dt, by = c("gene", "cell_type"),
                    all.x = TRUE)
  panel_dt <- merge(panel_dt, bulk_sub, by = "gene", all.x = TRUE)

  # Stable order: masked first (sorted by abs CT logFC), amplified second
  panel_dt[, category := factor(category,
                                levels = c("Bulk-masked", "Amplified in CT"))]
  panel_dt <- panel_dt[order(category, -abs(ct_logFC))]
  panel_dt[, gene := factor(gene, levels = rev(gene))]

  # Long form for paired points
  long <- data.table::rbindlist(list(
    panel_dt[, .(gene, category, cell_label,
                 source = "Bulk RNA-seq",
                 logFC = bulk_logFC, padj = bulk_padj)],
    panel_dt[, .(gene, category, cell_label,
                 source = "Cell-type pseudobulk",
                 logFC = ct_logFC, padj = ct_padj)]
  ))
  long[, sig_label := data.table::fcase(
    is.na(padj),     "n/d",
    padj < 0.001,    "***",
    padj < 0.01,     "**",
    padj < 0.05,     "*",
    default          = "ns"
  )]
  # Treat genes absent from bulk (e.g. CLEC4G) as logFC = 0 with "n/d" tag
  long[is.na(logFC) & source == "Bulk RNA-seq", logFC := 0]

  ct_palette <- c(
    "Fibroblast"     = "#F57F17",
    "Endothelial"    = "#2E7D32",
    "Cholangiocyte"  = "#1565C0"
  )

  p_f2 <- ggplot(long, aes(x = logFC, y = gene)) +
    # 1:1-ish reference: vertical lines at 0 and bulk |logFC|=0.5 cutoff
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = c(-0.5, 0.5), color = "grey80",
               linetype = "dashed", linewidth = 0.25) +
    # Pair connector: bulk -> CT
    geom_line(aes(group = gene), color = "grey70", linewidth = 0.4) +
    # Bulk = open square, CT = filled circle colored by cell type
    geom_point(data = long[source == "Bulk RNA-seq"],
               aes(shape = source), size = 1.8, fill = "white",
               color = "grey30", stroke = 0.5) +
    geom_point(data = long[source == "Cell-type pseudobulk"],
               aes(fill = cell_label, shape = source),
               size = 2.2, color = "grey20", stroke = 0.3) +
    geom_text(data = long[source == "Cell-type pseudobulk"],
              aes(label = sig_label),
              hjust = -0.4, size = 1.9, color = "grey20") +
    geom_text(data = long[source == "Bulk RNA-seq"],
              aes(label = sig_label),
              hjust = 1.4, size = 1.9, color = "grey45") +
    scale_shape_manual(values = c("Bulk RNA-seq" = 22,
                                  "Cell-type pseudobulk" = 21),
                       name = NULL) +
    scale_fill_manual(values = ct_palette, name = "Cell type") +
    guides(
      fill = guide_legend(override.aes = list(shape = 21, size = 2.6,
                                              color = "grey20",
                                              stroke = 0.3)),
      shape = guide_legend(override.aes = list(
        fill = c("Bulk RNA-seq" = "white",
                 "Cell-type pseudobulk" = "grey50"),
        color = c("Bulk RNA-seq" = "grey30",
                  "Cell-type pseudobulk" = "grey20"),
        size = c(1.8, 2.2)))
    ) +
    facet_grid(category ~ ., scales = "free_y", space = "free_y",
               switch = "y") +
    labs(x = expression(log[2]*" fold change (disease vs healthy)"),
         y = NULL,
         title = "scRNA-seq cell-type-specific DEGs diffused in bulk") +
    theme_masld() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, size = 6,
                                           face = "bold"),
          axis.text.y = element_text(face = "italic", size = 6),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.spacing.y = unit(0.05, "cm"),
          legend.text = element_text(size = 5),
          legend.key.size = unit(0.3, "cm"),
          plot.title = element_text(size = 7))

  save_fig(p_f2, file.path(OUT_DIR, "panel_f2_bulk_masked_genes.pdf"),
           width = fig_half_width, height = 3.5)

  # Persist underlying values for caption / supplement
  data.table::fwrite(panel_dt[, .(gene, category, cell_label,
                                  bulk_logFC, bulk_padj,
                                  ct_logFC, ct_padj)],
    file.path(OUT_DIR, "panel_f2_bulk_masked_genes.csv"))

}, error = function(e) message("Panel (f2) failed: ", e$message))

# ===========================================================================
# Panel (g): Top divergence genes between Metabolic and Fibrogenic (dot plot)
#            Annotated with ECM/fibrosis/collagen category
# ===========================================================================
tryCatch({
  message("Panel (g): Divergence gene dot plot...")

  div <- fread(file.path(PROG_DIR, "divergence_genes.csv"))

  # Take top 30 by absolute Cohen's d
  div[, abs_d := abs(cohens_d)]
  div_top <- div[order(-abs_d)][1:min(30, .N)]

  # Annotate functional category based on known ECM/fibrosis genes
  ecm_genes <- c("COL15A1", "COL16A1", "COL1A1", "COL1A2", "COL3A1", "COL4A1",
                  "COL4A2", "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3",
                  "COL14A1", "COL12A1", "FN1", "LAMA2", "LAMB1", "LAMC1",
                  "VCAN", "BGN", "DCN", "LUM")
  fibrosis_genes <- c("ACTA2", "TGFB1", "TGFB2", "TGFBR1", "PDGFRB", "PDGFRA",
                       "LOXL2", "LTBP2", "LTBP1", "CDH11", "FAP", "THY1",
                       "TIMP1", "TIMP2", "MMP2", "ADAMTS2", "POSTN", "SPP1")
  stellate_genes <- c("LRAT", "HGF", "DES", "GFAP", "RGS5")
  immune_genes <- c("CD68", "TREM2", "CCL2", "CCR2", "IL1B", "TNF", "MARCO")

  div_top[, category := fcase(
    gene_symbol %in% ecm_genes, "ECM/Collagen",
    gene_symbol %in% fibrosis_genes, "Fibrogenesis",
    gene_symbol %in% stellate_genes, "Stellate cell",
    gene_symbol %in% immune_genes, "Immune",
    default = "Other"
  )]

  # Ensure display name
  div_top[, display := fifelse(is.na(gene_symbol) | gene_symbol == "",
                               sub("\\..*", "", gene), gene_symbol)]
  # Remove Ensembl-only names from display
  div_top[grepl("^ENSG", display), display := sub("\\..*", "", gene)]

  # Order by Cohen's d
  div_top[, display := factor(display, levels = rev(div_top[order(cohens_d)]$display))]

  cat_colors <- c(
    "ECM/Collagen" = "#F57F17",
    "Fibrogenesis"  = "#C9265E",
    "Stellate cell" = "#7B1FA2",
    "Immune"        = "#2E7D32",
    "Other"         = "#78909C"
  )

  p_g <- ggplot(div_top, aes(x = cohens_d, y = display, color = category)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3,
               color = "gray60") +
    geom_segment(aes(x = 0, xend = cohens_d, y = display, yend = display),
                 linewidth = 0.3) +
    geom_point(aes(size = -log10(padj + 1e-300)), alpha = 0.8) +
    scale_color_manual(values = cat_colors, name = "Category") +
    scale_size_continuous(range = c(0.8, 3), name = expression(-log[10]*"(padj)"),
                         guide = guide_legend(override.aes = list(alpha = 1))) +
    labs(title = "Metabolic/Fibrogenic divergence genes (top 30)",
         x = "Cohen's d (Metabolic vs Fibrogenic)", y = NULL) +
    theme_masld() +
    theme(axis.text.y = element_text(size = 5))

  save_fig(p_g, file.path(OUT_DIR, "panel_g_divergence.pdf"),
           width = fig_half_width + 0.5, height = 4)

}, error = function(e) message("Panel (g) failed: ", e$message))

# ===========================================================================
# Panel (h): Exemplar gene expression across pseudotime
#            Key fibrosis/hepatocyte markers plotted as LOESS curves
# ===========================================================================
tryCatch({
  message("Panel (h): Exemplar gene trajectories...")

  # transition_programs.csv schema (2026-04) is one-row-per-gene summarising
  # peak_transition; does not carry full trajectories. Use transition_fib_dream
  # directly which has per-gene-per-transition logFC.
  prog <- data.table()

  # Target genes: fibrosis markers (up) + hepatocyte function markers (down)
  target_symbols <- c("COL1A1", "ACTA2", "SPP1", "TGFB1", "CYP3A4", "ALB")

  prog_targets <- data.table()

  if (nrow(prog_targets) > 0) {
    # Create pseudo-x axis: midpoint of each transition
    trans_x <- data.table(
      transition = names(transition_labels),
      x_pos = c(0.5, 1.5, 2.5, 3.5),
      label = transition_labels
    )
    prog_targets <- merge(prog_targets, trans_x, by = "transition")

    # Color: up-genes in magenta, down-genes in blue
    down_genes <- c("CYP3A4", "ALB")
    prog_targets[, direction := fifelse(symbol %in% down_genes,
                                         "Hepatocyte (lost)", "Fibrosis (gained)")]

    dir_colors <- c(
      "Fibrosis (gained)"  = masld_colors$up,
      "Hepatocyte (lost)"  = masld_colors$down
    )

    p_h <- ggplot(prog_targets, aes(x = x_pos, y = logFC,
                                     group = symbol, color = direction)) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3,
                 color = "gray60") +
      geom_line(linewidth = 0.6, alpha = 0.8) +
      geom_point(size = 1.5) +
      geom_text_repel(
        data = prog_targets[x_pos == max(x_pos)],
        aes(label = symbol), size = 2, nudge_x = 0.3, segment.size = 0.2,
        direction = "y", max.overlaps = 20, show.legend = FALSE
      ) +
      scale_color_manual(values = dir_colors, name = "Program") +
      scale_x_continuous(
        breaks = trans_x$x_pos, labels = trans_x$label,
        expand = expansion(mult = c(0.05, 0.15))
      ) +
      labs(title = "Key gene dynamics across fibrosis stages",
           x = NULL, y = "logFC (vs. preceding stage)") +
      theme_masld() +
      theme(
        axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
        legend.position = c(0.2, 0.85)
      )
  } else {
    # Fallback: if target genes not found in transition_programs, use dream results
    message("  Target genes not found in transition_programs; using dream DE fallback...")

    fib_de <- fread(file.path(PROG_DIR, "transition_fib_dream_results.csv"))
    fib_de <- add_symbols(fib_de, "gene")
    prog_targets <- fib_de[symbol %in% target_symbols]

    trans_x <- data.table(
      transition = names(transition_labels),
      x_pos = c(0.5, 1.5, 2.5, 3.5),
      label = transition_labels
    )
    prog_targets <- merge(prog_targets, trans_x, by = "transition")

    down_genes <- c("CYP3A4", "ALB")
    prog_targets[, direction := fifelse(symbol %in% down_genes,
                                         "Hepatocyte (lost)", "Fibrosis (gained)")]

    dir_colors <- c(
      "Fibrosis (gained)"  = masld_colors$up,
      "Hepatocyte (lost)"  = masld_colors$down
    )

    p_h <- ggplot(prog_targets, aes(x = x_pos, y = logFC,
                                     group = symbol, color = direction)) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3,
                 color = "gray60") +
      geom_line(linewidth = 0.6, alpha = 0.8) +
      geom_point(size = 1.5) +
      geom_text_repel(
        data = prog_targets[x_pos == max(x_pos)],
        aes(label = symbol), size = 2, nudge_x = 0.3, segment.size = 0.2,
        direction = "y", max.overlaps = 20, show.legend = FALSE
      ) +
      scale_color_manual(values = dir_colors, name = "Program") +
      scale_x_continuous(
        breaks = trans_x$x_pos, labels = trans_x$label,
        expand = expansion(mult = c(0.05, 0.15))
      ) +
      labs(title = "Key gene dynamics across fibrosis stages",
           x = NULL, y = "logFC (vs. preceding stage)") +
      theme_masld() +
      theme(
        axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
        legend.position = c(0.2, 0.85)
      )
  }

  save_fig(p_h, file.path(OUT_DIR, "panel_h_exemplars.pdf"),
           width = fig_half_width, height = 2.5)

}, error = function(e) message("Panel (h) failed: ", e$message))

# ===========================================================================
# Composite figure assembly (patchwork)
# ===========================================================================
message("Assembling composite figure...")

if (ENABLE_NMF_SUBTYPE_PANELS) {
  composite <- (
    (p_a | p_b) /
    (p_c | p_d) /
    (p_e | p_f) /
    (p_g | p_h)
  ) +
    plot_annotation(tag_levels = "a", title = NULL,
      theme = theme(plot.tag = element_text(size = 9, face = "bold", family = "Helvetica"))) +
    plot_layout(heights = c(1.2, 0.8, 0.9, 1.1))
  save_fig(composite, OUT, width = fig_full_width, height = 11)
} else {
  # NMF-dependent panels (c, d) dropped; assemble 6-panel composite.
  composite <- (
    (p_a | p_b) /
    (p_e | p_f) /
    (p_g | p_h)
  ) +
    plot_annotation(tag_levels = "a", title = NULL,
      theme = theme(plot.tag = element_text(size = 9, face = "bold", family = "Helvetica"))) +
    plot_layout(heights = c(1.2, 0.9, 1.1))
  save_fig(composite, OUT, width = fig_full_width, height = 9)
}

message("Done. Output: ", OUT)
message("Individual panels in: ", OUT_DIR)
