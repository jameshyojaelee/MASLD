#!/usr/bin/env Rscript
# ============================================================================
# figS_composition_mediation.R
#
# Supplementary panel: formal causal-mediation decomposition of bulk
# F-transition DEGs into composition-mediated (via cell-type proportion shifts)
# vs cell-type-intrinsic (residual / direct effect).
#
# Sub-panels:
#   A. Stacked horizontal bar: top-30 F1->F2 bulk DEGs by |beta_TE|, segmented
#      composition vs intrinsic; bar coloured by top mediating CT.
#   B. Quadrant scatter (small multiples per major cell type): top-bulk DEGs
#      with x=Delta-pi_k(F1->F2) and y=median Delta-mu_g,k from BayesPrism Z.
#      Quadrants: expansion+activation (UR), shrinkage+activation (UL),
#      expansion+suppression (LR), shrinkage+suppression (LL).
#   C. Heatmap: top-30 DEGs (F1->F2 by intrinsic-fraction descending)
#      x 4 F-transitions, fill = MP (0=intrinsic, 1=composition).
#
# Inputs:
#   - mediation/F_transition_mediation_per_gene.csv
#   - mediation/mediation_per_celltype.csv
#   - mediation/quadrant_data_per_celltype.csv  (built here on-demand if absent)
#   - bayesprism_proportions.csv  (for Delta-pi)
#   - bayesprism per-CT expression (for Delta-mu)
#
# Outputs: FIGS_SENS_DIR/figS_composition_mediation/figS_composition_mediation.pdf
#
# Refs: Vellame 2021 PMID 34353365 (precedent), Imai & Keele 2010, Meng 2023
# PMID 36472568, BayesPrism PMID 35469013.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# Output directory: dedicated subdir under sensitivity
PANEL_BASE <- file.path(FIGS_SENS_DIR, "figS_composition_mediation")
dir.create(PANEL_BASE, recursive = TRUE, showWarnings = FALSE)

BASE_SIZE <- 6
LBL_PT    <- 6
LBL_SIZE  <- LBL_PT / ggplot2::.pt

theme_figS <- function() {
  theme_masld(base_size = BASE_SIZE) +
    theme(
      axis.title    = element_text(size = BASE_SIZE),
      axis.text     = element_text(size = BASE_SIZE),
      legend.title  = element_text(size = BASE_SIZE),
      legend.text   = element_text(size = BASE_SIZE),
      strip.text    = element_text(size = BASE_SIZE),
      plot.subtitle = element_blank(),
      plot.margin   = margin(3, 3, 3, 3)
    )
}

# Cell-type palette (matches fig2_panel_A0)
ct_palette <- c(
  "Hepatocyte"    = "#1B5E20",
  "Stellate"      = "#FF6F00",
  "Cholangiocyte" = "#00BCD4",
  "Macrophage"    = "#212121",
  "Endothelial"   = "#7CB342"
)
QUAD_CTS <- names(ct_palette)

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
INT     <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
MED_DIR <- file.path(INT, "results/mediation")
PROP_F  <- file.path(INT, "results/progression/cibersortx_celltype_expression",
                     "bayesprism_proportions.csv")
META_F  <- file.path(INT, "metadata/unified_metadata.csv")
QC_F    <- file.path(INT, "qc/sample_qc_report.csv")
CTEX_DIR <- file.path(INT, "results/progression/cibersortx_celltype_expression")

med_main <- fread(file.path(MED_DIR, "F_transition_mediation_per_gene.csv"))
med_ct   <- fread(file.path(MED_DIR, "mediation_per_celltype.csv"))

# ---------------------------------------------------------------------------
# Helper: build quadrant data on the fly (if not pre-computed)
# Computes Delta-pi_k between consecutive stages and median Delta-mu_g,k
# from BayesPrism per-CT expression, restricted to the union of top-200 DEGs.
# ---------------------------------------------------------------------------
QUAD_F <- file.path(MED_DIR, "quadrant_data_per_celltype.csv")
if (!file.exists(QUAD_F)) {
  message("[quad] Building quadrant_data_per_celltype.csv on the fly...")
  meta <- fread(META_F)
  qc   <- fread(QC_F)
  meta <- merge(meta, qc[, .(sample_id, pass_technical)], by = "sample_id",
                all.x = TRUE)
  meta <- meta[!is.na(pass_technical) & pass_technical == TRUE &
                 !is.na(fibrosis_stage) &
                 !dataset %in% c("GSE213621", "PRJNA512027")]  # coarse/dropped cohorts
  meta[, fibrosis_stage := as.integer(fibrosis_stage)]

  prop <- fread(PROP_F)
  prop_long <- melt(prop, id.vars = "sample_id",
                    variable.name = "cell_type", value.name = "pi")
  prop_long <- merge(prop_long, meta[, .(sample_id, fibrosis_stage)],
                     by = "sample_id")

  # Per-stage mean proportion per CT
  pi_stage <- prop_long[cell_type %in% QUAD_CTS,
                        .(pi_mean = mean(pi, na.rm = TRUE)),
                        by = .(cell_type, fibrosis_stage)]
  pi_stage[, fibrosis_stage := paste0("F", fibrosis_stage)]

  # Delta-pi at each consecutive transition
  pi_w <- dcast(pi_stage, cell_type ~ fibrosis_stage, value.var = "pi_mean")
  pi_w[, dpi_F1_vs_F0 := F1 - F0]
  pi_w[, dpi_F2_vs_F1 := F2 - F1]
  pi_w[, dpi_F3_vs_F2 := F3 - F2]
  pi_w[, dpi_F4_vs_F3 := F4 - F3]

  # Genes to tile: union of top-200 across transitions
  med_main[, ensembl_clean := sub("\\..*", "", gene)]
  gene_universe <- unique(med_main[, .(ensembl_clean, symbol)])

  # For each CT, load BayesPrism expression matrix (sample x gene_symbol)
  quad_rows <- list()
  for (ct in QUAD_CTS) {
    fct <- file.path(CTEX_DIR, sprintf("bayesprism_%s.csv.gz", ct))
    if (!file.exists(fct)) {
      message(sprintf("  [skip] %s not found", fct))
      next
    }
    message(sprintf("  [load] %s", fct))
    expr <- fread(fct)
    setnames(expr, names(expr)[1], "sample_id")
    # Long subset to gene universe (by symbol)
    sym_in <- intersect(names(expr), gene_universe$symbol)
    if (length(sym_in) == 0) next
    sub <- expr[, c("sample_id", sym_in), with = FALSE]
    # log2(x+1)
    for (j in sym_in) sub[, (j) := log2(get(j) + 1)]
    sub_long <- melt(sub, id.vars = "sample_id", variable.name = "symbol",
                     value.name = "mu")
    sub_long <- merge(sub_long, meta[, .(sample_id, fibrosis_stage)],
                      by = "sample_id")
    # Per-stage mean expression
    mu_stage <- sub_long[, .(mu_mean = mean(mu, na.rm = TRUE)),
                         by = .(symbol, fibrosis_stage)]
    mu_stage[, fibrosis_stage := paste0("F", fibrosis_stage)]
    mu_w <- dcast(mu_stage, symbol ~ fibrosis_stage, value.var = "mu_mean")
    # Delta-mu per transition
    if (all(c("F0","F1","F2","F3","F4") %in% names(mu_w))) {
      mu_w[, dmu_F1_vs_F0 := F1 - F0]
      mu_w[, dmu_F2_vs_F1 := F2 - F1]
      mu_w[, dmu_F3_vs_F2 := F3 - F2]
      mu_w[, dmu_F4_vs_F3 := F4 - F3]
    }
    # Reshape long
    dmu_cols <- grep("^dmu_", names(mu_w), value = TRUE)
    if (length(dmu_cols) == 0) next
    mu_l <- melt(mu_w[, c("symbol", dmu_cols), with = FALSE],
                 id.vars = "symbol", variable.name = "dmu_var",
                 value.name = "delta_mu")
    mu_l[, transition := sub("^dmu_", "", dmu_var)]
    mu_l[, dmu_var := NULL]
    mu_l[, cell_type := ct]
    # Attach delta_pi
    pi_l <- melt(pi_w[cell_type == ct], id.vars = "cell_type",
                 measure.vars = grep("^dpi_", names(pi_w), value = TRUE),
                 variable.name = "dpi_var", value.name = "delta_pi")
    pi_l[, transition := sub("^dpi_", "", dpi_var)]
    pi_l[, dpi_var := NULL]
    out <- merge(mu_l, pi_l[, .(transition, delta_pi)], by = "transition",
                 all.x = TRUE)
    quad_rows[[ct]] <- out
  }
  quad <- rbindlist(quad_rows, fill = TRUE)
  fwrite(quad, QUAD_F)
  message(sprintf("  Wrote: %s (%d rows)", QUAD_F, nrow(quad)))
} else {
  message("[quad] Reusing existing quadrant_data_per_celltype.csv")
}
quad <- fread(QUAD_F)

# ---------------------------------------------------------------------------
# Sub-panel A: Top-30 F1->F2 bulk DEGs, intrinsic vs composition stacked bar
# ---------------------------------------------------------------------------
PRIMARY_TRANS <- "F2_vs_F1"

# Top-30 by |beta_TE| at the F1->F2 boundary. Larger TEs yield more stable
# MP = IE/TE ratios (ratio is well-behaved when denominator is large).
a_top <- copy(med_main[transition == PRIMARY_TRANS & !is.na(MP)])
a_top[, abs_TE := abs(beta_TE)]
setorder(a_top, -abs_TE)
a_top <- head(a_top, 30)
message(sprintf("[A] Top-30 F1->F2: |beta_TE| range %.3f - %.3f",
                min(a_top$abs_TE), max(a_top$abs_TE)))

# Order bars by intrinsic fraction (1 - MP) descending
a_top[, intrinsic_frac := 1 - MP]
a_top[, label := ifelse(!is.na(symbol) & symbol != "", symbol, gene)]
a_top[, label := factor(label, levels = label[order(intrinsic_frac)])]

# Long form: composition + intrinsic segments
a_long <- rbind(
  a_top[, .(label, abs_TE, top_mediating_CT,
            component = "Composition", value = abs_TE * MP)],
  a_top[, .(label, abs_TE, top_mediating_CT,
            component = "Intrinsic",   value = abs_TE * (1 - MP))]
)
a_long[, component := factor(component, levels = c("Intrinsic", "Composition"))]

# Bar fill: composition segment uses CT colour for the top mediating CT;
# intrinsic segment is neutral grey. Implemented via interaction key.
a_long[, fill_key := ifelse(component == "Intrinsic", "Intrinsic",
                            paste0("Comp:", top_mediating_CT))]
fill_lookup <- c("Intrinsic" = "#9E9E9E",
                 setNames(ct_palette, paste0("Comp:", names(ct_palette))))
# Order in legend: Intrinsic first then CTs alphabetically
a_long[, fill_key := factor(fill_key,
                            levels = c("Intrinsic",
                                       paste0("Comp:", names(ct_palette))))]
# Stack order: Intrinsic first (left), Composition second (right)
a_long[, component := factor(component, levels = c("Intrinsic", "Composition"))]

p_A <- ggplot(a_long, aes(x = value, y = label, fill = fill_key,
                          group = component)) +
  geom_col(colour = "white", linewidth = 0.2, width = 0.78,
           position = position_stack()) +
  scale_fill_manual(values = fill_lookup, name = NULL,
                    labels = c("Intrinsic",
                               paste0("Composition (", names(ct_palette), ")"))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = "|beta_TE|  (composition + intrinsic)", y = NULL) +
  theme_figS() +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "right",
        legend.key.height = unit(7, "pt"),
        legend.key.width  = unit(7, "pt"))

# Save sub-panel data
fwrite(a_top,  file.path(PANEL_BASE, "figS_composition_mediation_A_top30.csv"))
fwrite(a_long, file.path(PANEL_BASE, "figS_composition_mediation_A_long.csv"))

# ---------------------------------------------------------------------------
# Sub-panel B: Quadrant scatter per cell type (F1->F2)
# ---------------------------------------------------------------------------
quad_b <- quad[transition == PRIMARY_TRANS & cell_type %in% QUAD_CTS &
                 is.finite(delta_mu) & is.finite(delta_pi)]
# Restrict to the top-bulk-DEG universe (union of all transitions = the genes
# we ran mediation on, but at F1->F2 we only need the F1->F2 200 list)
top_at_trans <- med_main[transition == PRIMARY_TRANS, unique(symbol)]
quad_b <- quad_b[symbol %in% top_at_trans]

# Highlight surviving-hep stress vs HSC expansion: top-10 by |dmu| within CT
quad_b[, abs_dmu := abs(delta_mu)]
quad_b[, hl := abs_dmu >= quantile(abs_dmu, 0.9, na.rm = TRUE),
       by = cell_type]

# CT-level Δπ (one number per CT); show as vertical reference line
ct_dpi <- unique(quad_b[, .(cell_type, delta_pi)])

# Fill order
quad_b[, cell_type := factor(cell_type, levels = QUAD_CTS)]

p_B <- ggplot(quad_b, aes(x = delta_mu, y = delta_pi, colour = cell_type)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3, linetype = 2) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3, linetype = 2) +
  geom_point(alpha = 0.55, size = 0.8, show.legend = FALSE) +
  ggrepel::geom_text_repel(
    data = quad_b[hl == TRUE],
    aes(label = symbol),
    size = LBL_SIZE * 0.75, colour = "grey15",
    max.overlaps = 20, segment.size = 0.2, segment.alpha = 0.5,
    min.segment.length = 0.1,
    show.legend = FALSE
  ) +
  facet_wrap(~ cell_type, scales = "free", ncol = 5) +
  scale_colour_manual(values = ct_palette) +
  labs(x = "Delta mu_g,k  (F1 -> F2 cell-type expression change, log2)",
       y = "Delta pi_k  (F1 -> F2 proportion change)") +
  theme_figS() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "plain", size = BASE_SIZE),
        panel.spacing = unit(6, "pt"))

fwrite(quad_b, file.path(PANEL_BASE, "figS_composition_mediation_B_quad.csv"))

# ---------------------------------------------------------------------------
# Sub-panel C: MP heatmap, top-30 F1->F2 DEGs x 4 transitions
# ---------------------------------------------------------------------------
sel_genes <- as.character(a_top$gene)  # ENSG-version IDs
sel_lbl   <- as.character(a_top$label)
names(sel_lbl) <- sel_genes

c_dat <- med_main[gene %in% sel_genes,
                   .(gene, transition, MP)]
c_dat[, label := factor(sel_lbl[gene], levels = levels(a_top$label))]
c_dat[, transition := factor(transition,
                             levels = c("F1_vs_F0", "F2_vs_F1",
                                        "F3_vs_F2", "F4_vs_F3"),
                             labels = c("F0->F1", "F1->F2", "F2->F3", "F3->F4"))]

p_C <- ggplot(c_dat, aes(x = transition, y = label, fill = MP)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(is.na(MP), "", sprintf("%.2f", MP))),
            size = LBL_SIZE * 0.8, colour = "grey20") +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C2185B",
                       midpoint = 0.5, limits = c(0, 1),
                       na.value = "grey90",
                       name = "MP\n(0=intrinsic\n1=composition)",
                       guide = guide_colorbar(barwidth = unit(0.25, "cm"),
                                              barheight = unit(2.5, "cm"))) +
  labs(x = NULL, y = NULL) +
  theme_figS() +
  theme(panel.grid = element_blank(),
        axis.line  = element_blank(),
        axis.ticks = element_blank())

fwrite(c_dat, file.path(PANEL_BASE, "figS_composition_mediation_C_heatmap.csv"))

# ---------------------------------------------------------------------------
# Compose
# ---------------------------------------------------------------------------
top_row    <- p_A | p_C
bottom_row <- p_B
p_full <- top_row / bottom_row +
  plot_layout(heights = c(1.0, 0.7))

message("[caption] A: Top-30 F1->F2 bulk DEGs: composition vs intrinsic split")
message("[caption] B: Quadrant: cell-type proportion change vs cell-type expression change")
message("[caption] C: Mediation proportion across F-transitions")

out_pdf <- file.path(PANEL_BASE, "figS_composition_mediation.pdf")
save_fig(p_full, out_pdf, width = fig_full_width, height = fig_full_width * 10 / 13)
message("Wrote: ", out_pdf)

# ---------------------------------------------------------------------------
# Stand-alone summary report (printed to log)
# ---------------------------------------------------------------------------
mp_bins <- med_main[!is.na(MP), .(
  pct_high  = mean(MP > 0.7) * 100,
  pct_mid   = mean(MP >= 0.3 & MP <= 0.7) * 100,
  pct_low   = mean(MP < 0.3) * 100,
  n         = .N
), by = transition]
print(mp_bins)
fwrite(mp_bins, file.path(PANEL_BASE, "figS_composition_mediation_mp_summary.csv"))

message("Done.")
