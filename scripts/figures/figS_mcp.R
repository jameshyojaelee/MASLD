#!/usr/bin/env Rscript
# figS_mcp.R — Supplementary figure panels for the cNMF multi-cellular program atlas.
#
# Panels (all under FIGS_MCP_DIR):
#   a: cNMF k-selection: cophenetic correlation and reconstruction error across k
#   b: cNMF consensus matrix at chosen k (3x3 panel)
#   c: Program recovery Jaccard across leave-one-dataset-out runs
#   d: Program stability across HVG counts (sensitivity to number of variable genes)
#   e: cNMF vs plain NMF cophenetic and cross-seed Jaccard comparison (matched configuration)
#   f: Multi-method factorization comparison (cNMF, plain NMF, scHPF, LIGER, MOFA+)
#   g: Program overlap heatmap: cNMF k=16 vs alternative methods (top-100 Jaccard)
#   h: Bulk k=6 vs scRNA cNMF program correspondence
#   i: Spatial coherence: Moran's I per program on Visium data (9,999 permutations)
#   j: Olink + DIA-MS proteomics concordance with cNMF program scores
#   k: COLOC enrichment per program with matched-random gene-set null
#   l: Mouse cross-species conservation scores with permutation null
#   m: Composition-adjustment controls: stage effects before and after cell-type fraction residualization
#   n: Exemplar program hero panel (top-scoring program across all modalities)
#   p: cNMF k=10 stage trajectories — z-score heatmap + raw mean-usage line plot
#      (k=10 primary manuscript k; direct analog of bulk figS08 panels c/d)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
})

# Config
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

MCP_DIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/mcp")
OUT <- FIGS_MCP_DIR

args <- commandArgs(trailingOnly = TRUE)
# Which panels to render; default: all
panels <- if (length(args) > 0) strsplit(args[1], ",")[[1]] else letters[1:14]

# --------------------------------------------------------------------------- #
# Panel A: k-selection stability + reconstruction error                       #
# --------------------------------------------------------------------------- #
panel_a <- function(name = "global") {
  kplot_f <- file.path(MCP_DIR, "cnmf_runs", name, name, sprintf("%s.k_selection.png", name))
  if (!file.exists(kplot_f)) kplot_f <- file.path(MCP_DIR, "cnmf_runs", name, sprintf("%s.k_selection.png", name))
  if (!file.exists(kplot_f)) {
    cat(sprintf("[figS_mcp] Panel A skipped: %s not found\n", kplot_f))
    return(NULL)
  }
  # Direct PDF conversion via magick
  if (requireNamespace("magick", quietly = TRUE)) {
    img <- magick::image_read(kplot_f)
    magick::image_write(img, path = file.path(OUT, sprintf("figSmcp_a_kselection_%s.pdf", name)),
                        format = "pdf")
  } else {
    file.copy(kplot_f, file.path(OUT, sprintf("figSmcp_a_kselection_%s.png", name)), overwrite = TRUE)
  }
  cat("[figS_mcp] panel a written\n")
}

# --------------------------------------------------------------------------- #
# Panel G: cNMF <-> DIALOGUE concordance heatmap                              #
# --------------------------------------------------------------------------- #
panel_g <- function() {
  f <- file.path(MCP_DIR, "integration/cnmf_dialogue_concordance.csv")
  if (!file.exists(f)) { cat("[figS_mcp] Panel G skipped (missing)\n"); return(NULL) }
  d <- fread(f)
  d[, neg_log10_q := -log10(pmax(hypergeom_q, 1e-10))]
  p <- ggplot(d, aes(factor(dialogue_meta_mcp), factor(cnmf_program))) +
    geom_tile(aes(fill = jaccard_top100)) +
    geom_text(data = d[hypergeom_q < 0.05], aes(label = "*"), size = GEOM_TEXT_6PT, color = "white") +
    scale_fill_gradient(low = "white", high = "#C2185B", limits = c(0, 0.5), oob = scales::squish) +
    facet_wrap(~ celltype, scales = "free", nrow = 2) +
    labs(x = "DIALOGUE meta-MCP", y = "cNMF program", fill = "Jaccard top-100") +
    theme_minimal(base_size = 6) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  ggsave(file.path(OUT, "figSmcp_g_cnmf_dialogue_concordance.pdf"),
         p, width = fig_full_width, height = 4.96, device = cairo_pdf)
  message("[caption] cNMF GEPs vs DIALOGUE MCPs (cross-tool concordance)")
  cat("[figS_mcp] panel g written\n")
}

# --------------------------------------------------------------------------- #
# Panel H: bulk k=6 <-> scRNA cNMF crossmap                                   #
# --------------------------------------------------------------------------- #
panel_h <- function(name = "global", k = 20) {
  f <- file.path(MCP_DIR, "integration", sprintf("bulk_sc_bridge_%s_k%d.tsv", name, k))
  if (!file.exists(f)) { cat("[figS_mcp] Panel H skipped (missing)\n"); return(NULL) }
  d <- fread(f)
  p <- ggplot(d, aes(bulk_program, factor(cnmf_program))) +
    geom_tile(aes(fill = pearson_r)) +
    scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C2185B", limits = c(-1, 1)) +
    labs(x = "Bulk k=6 program", y = "scRNA cNMF program", fill = "Pearson r") +
    theme_minimal(base_size = 6) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  ggsave(file.path(OUT, sprintf("figSmcp_h_bulk_sc_bridge_%s_k%d.pdf", name, k)),
         p, width = 6, height = 6, device = cairo_pdf)
  message("[caption] Bulk vs single-cell program correspondence")
  cat("[figS_mcp] panel h written\n")
}

# --------------------------------------------------------------------------- #
# Panel I: Spatial Visium Moran's I                                           #
# --------------------------------------------------------------------------- #
panel_i <- function(name = "global", k = 20) {
  f <- file.path(MCP_DIR, "validation/spatial", sprintf("spatial_moran_%s_k%d.tsv", name, k))
  if (!file.exists(f)) { cat("[figS_mcp] Panel I skipped (missing)\n"); return(NULL) }
  d <- fread(f)
  d[, program_num := as.integer(sub("program_", "", program))]
  p <- ggplot(d, aes(factor(program_num), morans_I)) +
    geom_boxplot(fill = "#F48FB1", alpha = 0.6, outlier.alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(x = "cNMF program", y = "Moran's I (spatial autocorrelation)") +
    theme_minimal(base_size = 6)
  ggsave(file.path(OUT, sprintf("figSmcp_i_spatial_moran_%s_k%d.pdf", name, k)),
         p, width = 7, height = 4, device = cairo_pdf)
  message("[caption] Spatial coherence of programs across Visium slides")
  cat("[figS_mcp] panel i written\n")
}

# --------------------------------------------------------------------------- #
# Panel K: COLOC enrichment                                                   #
# --------------------------------------------------------------------------- #
panel_k <- function(name = "global", k = 20) {
  f <- file.path(MCP_DIR, "validation/coloc", sprintf("coloc_enrichment_%s_k%d.tsv", name, k))
  if (!file.exists(f)) { cat("[figS_mcp] Panel K skipped (missing)\n"); return(NULL) }
  d <- fread(f)
  d[, neg_log10_q := -log10(pmax(qval, 1e-10))]
  p <- ggplot(d[set == "coloc_pp4_gt_0.5"],
              aes(factor(program), neg_log10_q)) +
    geom_col(fill = "#C2185B") +
    geom_hline(yintercept = -log10(0.05), linetype = 2, color = "gray40") +
    labs(x = "cNMF program", y = "-log10(q-value)") +
    theme_minimal(base_size = 6) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  ggsave(file.path(OUT, sprintf("figSmcp_k_coloc_enrichment_%s_k%d.pdf", name, k)),
         p, width = 8, height = 4, device = cairo_pdf)
  message("[caption] COLOC PP4 > 0.5 hypergeometric enrichment per program")
  cat("[figS_mcp] panel k written\n")
}

# --------------------------------------------------------------------------- #
# Panel O: Annotated cNMF k=16 program composition by disease stage           #
# (Promoted from a fig2d candidate; biological labels + lineage-grouped       #
#  colors so the cirrhosis lineage collapse reads without methods context.)   #
# --------------------------------------------------------------------------- #
panel_o <- function() {
  f <- file.path(MCP_DIR, "cnmf_annot/global/program_stage_usage.k16.tsv")
  if (!file.exists(f)) { cat("[figS_mcp] Panel O skipped (missing)\n"); return(NULL) }

  # Biological annotations from RESULTS_FINAL_k16.md.
  # `[b]` flags batch-confounded programs (P5, P9) -- shown but called out.
  bio <- data.table(
    program  = 1:16,
    label    = c("P1 Periportal hep.",          "P2 SH-peak hep.",
                 "P3 LSEC (sinusoidal)",         "P4 Healthy baseline hep.",
                 "P5 Hep secretory [b]",         "P6 Portal fibroblast",
                 "P7 Ductular reaction",         "P8 Macrophage (tissue-res.)",
                 "P9 APC / MHC-II [b]",          "P10 LSEC scavenger",
                 "P11 Late-emergent [a]",         "P12 Arterial endothelial",
                 "P13 Kupffer",                  "P14 SM / pericyte",
                 "P15 Cholangiocyte",            "P16 Activated HSC"),
    lineage  = c("Hepatocyte","Hepatocyte","Endothelial","Hepatocyte",
                 "Hepatocyte","Mesenchyme","Cholangiocyte","Immune",
                 "Immune","Endothelial","Endothelial","Endothelial",
                 "Immune","Mesenchyme","Cholangiocyte","Mesenchyme")
  )
  # Lineage-grouped palette (warm = hep, blues = endo, purple = immune,
  # browns = mesenchyme, greens = cholangiocyte). Order programs from
  # bottom (hepatocyte) -> top (cholangiocyte) so the cirrhosis collapse
  # is visible as warm bands shrinking and cool bands expanding.
  pal <- c(
    "P1 Periportal hep."         = "#B71C1C",
    "P2 SH-peak hep."            = "#E64A19",
    "P4 Healthy baseline hep."   = "#FF7043",
    "P5 Hep secretory [b]"       = "#FFAB91",
    "P3 LSEC (sinusoidal)"       = "#1565C0",
    "P10 LSEC scavenger"         = "#42A5F5",
    "P11 Late-emergent [a]"      = "#0277BD",
    "P12 Arterial endothelial"   = "#81D4FA",
    "P8 Macrophage (tissue-res.)"= "#6A1B9A",
    "P9 APC / MHC-II [b]"        = "#AB47BC",
    "P13 Kupffer"                = "#C2185B",
    "P6 Portal fibroblast"       = "#5D4037",
    "P14 SM / pericyte"          = "#8D6E63",
    "P16 Activated HSC"          = "#BF360C",
    "P7 Ductular reaction"       = "#2E7D32",
    "P15 Cholangiocyte"          = "#66BB6A"
  )
  stack_order <- names(pal)  # bottom -> top in stacked bar

  d <- fread(f)
  d <- merge(d, bio, by = "program")
  d_long <- melt(d,
    id.vars       = c("program","label","lineage"),
    measure.vars  = c("mean_stage_0","mean_stage_1","mean_stage_2","mean_stage_3"),
    variable.name = "stage", value.name = "mean_usage")
  # Cirrhosis (mean_stage_3) excluded: cNMF was run on contaminated atlas
  # including NPC-enriched donors; "Cirrhosis-spike" programs were confirmed
  # as artifacts of GSE136103 FACS-enriched donors in Phase 6.3 sensitivity
  # analysis (seg1 breakpoints lost in clean atlas). 3-stage axis is primary.
  d_long <- d_long[stage_raw != "mean_stage_3"]
  d_long[, stage := factor(stage,
    levels = c("mean_stage_0","mean_stage_1","mean_stage_2"),
    labels = c("Healthy","Steatosis","Steatohepatitis"))]
  tot <- d_long[, .(tot = sum(mean_usage)), by = stage]
  d_long <- merge(d_long, tot, by = "stage")
  d_long[, pct := 100 * mean_usage / tot]
  d_long[, label := factor(label, levels = stack_order)]

  p <- ggplot(d_long, aes(x = stage, y = pct, fill = label)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = pal, name = "cNMF k=16 program",
                      guide = guide_legend(ncol = 1, reverse = TRUE)) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Disease stage", y = "% mean program usage",
         caption = paste0("[a] Originally labeled 'Cirrhosis-spike'; reclassified as late-emergent after ",
                          "protocol-contamination remediation (2026-05-22) — breakpoint signal was driven by ",
                          "NPC-enriched GSE136103 donors (excluded). cNMF re-factorization on clean atlas pending.\n",
                          "[b] flags programs with high dataset-variance (P5: 76%; P9: 44%)")) +
    theme_minimal(base_size = 6) +
    theme(axis.text.x  = element_text(angle = 25, hjust = 1),
          legend.text  = element_text(size = 6),
          legend.title = element_text(size = 6),
          legend.key.size = unit(0.32, "cm"),
          plot.caption = element_text(size = 6, color = "gray35", hjust = 0))

  ggsave(file.path(OUT, "figSmcp_o_stage_composition_global_k16.pdf"),
         p, width = 7.4, height = 4.2, device = cairo_pdf)
  fwrite(d_long, file.path(OUT, "figSmcp_o_stage_composition_global_k16_data.csv"))
  message("[caption] cNMF k=16 program composition across disease stages (Healthy / Steatosis / Steatohepatitis)")
  cat("[figS_mcp] panel o written\n")
}

# --------------------------------------------------------------------------- #
# Panel P: cNMF k=10 stage trajectory heatmap + line plot                    #
# (k=10 = recommended primary manuscript k per RESULTS_FINAL_k16.md)         #
# Two sub-panels:                                                             #
#   p1: Z-scored mean-usage heatmap (programs × stages), programs sorted by  #
#       peak stage, row sidebar shows dominant cell-type lineage.             #
#   p2: Raw mean-usage line plot — one line per program, colored by lineage;  #
#       batch-confounded programs (P6) drawn dashed.                         #
# --------------------------------------------------------------------------- #
panel_p <- function(k = 10) {
  stage_f   <- file.path(MCP_DIR, sprintf("cnmf_annot/global/program_stage_usage.k%d.tsv", k))
  celltype_f <- file.path(MCP_DIR, sprintf("cnmf_annot/global/program_celltype_mean_usage.k%d.tsv", k))
  if (!file.exists(stage_f)) {
    cat("[figS_mcp] Panel P skipped: stage_usage file missing\n"); return(NULL)
  }

  # Biological annotations for k=10
  # Lineage derived from dominant cell-type in program_celltype_mean_usage.k10.tsv;
  # [b] = batch-confounded programs (k10-P6 maps to k16-P5 batch artifact per RESULTS_FINAL).
  # Stage direction (up/down/spike) from kruskal + spearman sign in program_phenotype_cor.k10.tsv.
  bio10 <- data.table(
    program   = 1:10,
    label     = c("P1 Periportal hep.",
                  "P2 SH-peak hep.",
                  "P3 LSEC (sinusoidal)",
                  "P4 Portal fibroblast",
                  "P5 Cholangiocyte",
                  "P6 Hep secretory [b]",
                  "P7 Macrophage",
                  "P8 Late-emergent [a]",
                  "P9 Stage-up (mixed)",
                  "P10 Activated HSC"),
    lineage   = c("Hepatocyte","Hepatocyte","Endothelial","Mesenchyme","Cholangiocyte",
                  "Hepatocyte","Immune","Immune","Endothelial","Mesenchyme"),
    batch_flag = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)
  )

  lineage_pal <- c(
    "Hepatocyte"    = "#E64A19",
    "Endothelial"   = "#1565C0",
    "Mesenchyme"    = "#5D4037",
    "Cholangiocyte" = "#2E7D32",
    "Immune"        = "#6A1B9A"
  )
  # Cirrhosis excluded — contamination artifact (see Panel O note); 3-stage primary
  stage_lvls <- c("Healthy", "Steatosis", "Steatohepatitis")

  # ---- load & merge ----
  su <- fread(stage_f)
  su <- merge(su, bio10, by = "program")
  dl <- melt(su,
    id.vars      = c("program","label","lineage","batch_flag","kruskal_p"),
    measure.vars = c("mean_stage_0","mean_stage_1","mean_stage_2","mean_stage_3"),
    variable.name = "stage_raw", value.name = "mean_usage")
  dl <- dl[stage_raw != "mean_stage_3"]   # remove Cirrhosis (contamination artifact)
  dl[, stage := factor(stage_raw,
    levels = c("mean_stage_0","mean_stage_1","mean_stage_2"),
    labels = stage_lvls)]

  # BH-correct kruskal p across programs for significance annotation
  pvals <- unique(su[, .(program, label, kruskal_p)])
  pvals[, kruskal_q := p.adjust(kruskal_p, "BH")]
  dl <- merge(dl, pvals[, .(program, kruskal_q)], by = "program")

  # ---- Sub-panel p1: z-score heatmap ----
  # Row z-score within each program (4 stage values -> z)
  hm <- dcast(dl, program + label + lineage + kruskal_q ~ stage, value.var = "mean_usage")
  for (s in stage_lvls) {
    vals <- hm[[s]]
    mn <- mean(vals); sd_v <- sd(vals)
    hm[[paste0("z_", s)]] <- if (sd_v > 0) (vals - mn) / sd_v else rep(0, length(vals))
  }
  # Sort programs: first by lineage (for grouping), then by stage of peak z-score
  hm[, peak_stage := stage_lvls[which.max(c(z_Healthy, z_Steatosis, z_Steatohepatitis))],
     by = program]
  hm[, peak_ord := match(peak_stage, stage_lvls)]
  hm[, lineage_ord := match(lineage, names(lineage_pal))]
  setorder(hm, lineage_ord, peak_ord)
  prog_order <- hm$label

  hz <- melt(hm, id.vars = c("program","label","lineage","kruskal_q","peak_stage"),
             measure.vars = paste0("z_", stage_lvls),
             variable.name = "stage_z", value.name = "z")
  hz[, stage := factor(gsub("^z_", "", stage_z), levels = stage_lvls)]
  hz[, label := factor(label, levels = prog_order)]
  hz[, sig_label := ifelse(kruskal_q < 0.05, "*", "")]

  # Row sidebar colors (lineage)
  sidebar_dt <- unique(hz[, .(label, lineage)])[, label := factor(label, levels = prog_order)]

  p1 <- ggplot(hz, aes(x = stage, y = label, fill = z)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sig_label), size = GEOM_TEXT_6PT, vjust = 0.8, color = "black") +
    scale_fill_gradient2(low = "#1565C0", mid = "grey95", high = "#C62828",
                         midpoint = 0, name = "Z-score\n(row-norm.)") +
    geom_tile(data = sidebar_dt,
              aes(x = -0.3, y = label, fill = NULL, color = lineage),
              width = 0.25, height = 0.85, inherit.aes = FALSE) +
    scale_color_manual(values = lineage_pal, name = "Lineage",
                       guide = guide_legend(override.aes = list(fill = lineage_pal,
                                                                 color = "white"))) +
    labs(x = NULL, y = NULL,
         caption = "* kruskal q < 0.05 (BH); [b] = batch-confounded program") +
    theme_masld(base_size = 6) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1),
          axis.text.y = element_text(size = 6),
          plot.caption = element_text(size = 6, color = "grey40", hjust = 0))

  message(sprintf("[caption] cNMF k=%d program usage across disease stages (z-scored)", k))
  save_fig(p1, file.path(OUT, sprintf("figSmcp_p1_stage_heatmap_k%d.pdf", k)),
           width = fig_half_width, height = 4.2)

  # ---- Sub-panel p2: raw mean-usage line plot ----
  dl[, label := factor(label, levels = rev(prog_order))]
  dl[, ltype := ifelse(batch_flag, "dashed", "solid")]

  p2 <- ggplot(dl, aes(x = stage, y = mean_usage, color = lineage,
                        group = label, linetype = ltype)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.2) +
    scale_color_manual(values = lineage_pal, name = "Lineage") +
    scale_linetype_identity(guide = guide_legend(title = NULL,
                                                  override.aes = list(color = "gray40")),
                            labels = c("dashed" = "batch-confounded", "solid" = "")) +
    facet_wrap(~ lineage, ncol = 3, scales = "free_y") +
    labs(x = NULL, y = "Mean cNMF program usage") +
    theme_masld(base_size = 6) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1),
          strip.background = element_rect(fill = "gray95", color = NA),
          panel.spacing = unit(4, "pt"),
          legend.position = "right")

  message(sprintf("[caption] cNMF k=%d stage trajectories by lineage", k))
  save_fig(p2, file.path(OUT, sprintf("figSmcp_p2_stage_lines_k%d.pdf", k)),
           width = fig_full_width, height = 4.0)

  cat(sprintf("[figS_mcp] panel p written (k=%d)\n", k))
}

# --------------------------------------------------------------------------- #
# Driver                                                                      #
# --------------------------------------------------------------------------- #
K_DEFAULT <- as.integer(Sys.getenv("MCP_K", "16"))
if ("a" %in% panels) panel_a()
if ("g" %in% panels) panel_g()
if ("h" %in% panels) panel_h(k = K_DEFAULT)
if ("i" %in% panels) panel_i(k = K_DEFAULT)
if ("k" %in% panels) panel_k(k = K_DEFAULT)
if ("o" %in% panels) panel_o()
if ("p" %in% panels) panel_p(k = 10)

cat(sprintf("[figS_mcp] DONE. Panels written to %s\n", OUT))
