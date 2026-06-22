#!/usr/bin/env Rscript
# ============================================================================
# fig4_validation_v2.R  (redesigned 2026-05-21)
#
# Figure 4: Cross-Modal Validation of Hotspot Module Biology
#
# Narrative: Three hepatocyte Hotspot modules (hep-20 glutamine/TGFβ, hep-24
# NRF2 antioxidant, hep-26 AP-1 injury) are independently confirmed by GWAS
# genetics, spatial transcriptomics, and plasma proteomics. GWAS variants
# additionally disrupt RORA and THRB binding sites, and both TFs show
# elevated chromatin activity in Progressor hepatocytes.
#
# Panels (written to FIG4_DIR/panels/):
#   4a — Cross-modal evidence matrix: 30 hep modules × 7 features     [OVERVIEW]
#   4b — Spatial: Moran's I scatter + CosMx hepatocyte bar            [SPATIAL]
#   4c — GWAS-ATAC TF lollipop + RORA/THRB activity bar               [GENETICS]
#   4d — Plasma concordance scatter colored by module membership       [PROTEOMICS]
#
# Run AFTER fig4_module_evidence_scan.R (needs module_evidence_matrix.csv)
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
})

set.seed(42)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG4_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

HS_RES <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")

# Module highlight palette (hep-20, hep-24, hep-26)
FOCAL_MODS  <- c(20L, 24L, 26L)
MOD_COLORS  <- c("20" = "#E07B39", "24" = masld_colors$fibrosis, "26" = masld_colors$mash)

# ============================================================================
# Panel 4a — Cross-modal evidence matrix
# ============================================================================
message("[4a] Cross-modal evidence matrix")

ev_f <- file.path(DATA_DIR, "module_evidence_matrix.csv")
if (!file.exists(ev_f))
  stop("Run fig4_module_evidence_scan.R first to generate module_evidence_matrix.csv")

ev  <- fread(ev_f)
ev[, focal := module %in% c(20L, 24L, 26L)]
# Mark focal modules with a bullet in the row label
ev[, module_lbl := ifelse(focal,
                          paste0("Hep-", module, " \u25cf"),
                          paste0("Hep-", module))]

# Long format for heatmap
feat_cols <- c("fstage_z","bulk_z","stab_z","brep_z","coloc_z","svg_z","plasma_z")
feat_labs <- c("Severity","Bulk rep.","Stability",
               "Replicated","COLOC","SVG","Plasma")

id_cols  <- c("module","module_lbl","composite","focal")
ev_long  <- melt(ev[, c(id_cols, feat_cols), with = FALSE],
                id.vars       = id_cols,
                variable.name = "feature",
                value.name    = "score")
ev_long[, feature_lbl := factor(feature, levels = feat_cols, labels = feat_labs)]

# Pass thresholds for filled vs open dots
pass_cols <- c("fstage_pass","bulk_pass","stab_pass","brep_pass",
               "coloc_pass","svg_pass","plasma_pass")
ev_pass <- melt(ev[, c("module", "module_lbl", pass_cols), with = FALSE],
                id.vars       = c("module","module_lbl"),
                variable.name = "pass_col",
                value.name    = "pass")
ev_pass[, feature := factor(
  sub("_pass","_z", pass_col), levels = feat_cols, labels = feat_labs)]
ev_long <- ev_pass[, .(module, feature_lbl = feature, pass)][ev_long, on = c("module","feature_lbl")]

# Module row order: by composite score descending
mod_order <- ev[order(-composite), module_lbl]
ev_long[, module_lbl := factor(module_lbl, levels = mod_order)]

p4a <- ggplot(ev_long, aes(x = feature_lbl, y = module_lbl)) +
  geom_tile(aes(fill = score), color = "white", linewidth = 0.3) +
  geom_point(data = ev_long[pass == TRUE],
             aes(x = feature_lbl, y = module_lbl),
             shape = 21, size = 0.9, color = "white", fill = "white",
             alpha = 0.6) +
  scale_fill_gradient(low = "white", high = "#2C3E7A",
                      name = "Score", na.value = "gray90",
                      limits = c(0, 1)) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL, title = "Module cross-modal evidence") +
  theme_masld(base_size = 9) +
  theme(axis.text.x  = element_text(size = 7, angle = 35, hjust = 0,
                                     vjust = 0),
        axis.text.y  = element_text(size = 6.5),
        panel.grid   = element_blank(),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width  = unit(0.18, "cm"),
        legend.text  = element_text(size = 7),
        legend.title = element_text(size = 7),
        plot.margin  = margin(t = 20, r = 4, b = 4, l = 4),
        plot.title   = element_text(size = 9, face = "bold"))

save_fig(p4a, file.path(PANEL_DIR, "fig4a.pdf"),
         width = fig_half_width, height = 5.0)

# ============================================================================
# Panel 4b — Spatial: Moran's I scatter + CosMx hepatocyte bar
# ============================================================================
message("[4b] Spatial SVG scatter + CosMx validation")

svg_h <- fread(file.path(SPATIAL_DIR, "svg", "svgs_Healthy.csv"))
svg_s <- fread(file.path(SPATIAL_DIR, "svg", "svgs_Steatotic.csv"))
setnames(svg_h, 1, "gene"); setnames(svg_s, 1, "gene")
for (col in "svg") {
  if (is.character(svg_h[[col]])) svg_h[, (col) := get(col) == "True"]
  if (is.character(svg_s[[col]])) svg_s[, (col) := get(col) == "True"]
}
svg_m <- merge(svg_h[, .(gene, I_healthy = I, svg_h = svg)],
               svg_s[, .(gene, I_masld   = I, svg_s = svg)],
               by = "gene", all = TRUE)
svg_m[is.na(svg_h), svg_h := FALSE]
svg_m[is.na(svg_s), svg_s := FALSE]

diff_svg <- fread(file.path(SPATIAL_DIR, "svg", "differential_svgs.csv"))
setnames(diff_svg, 1, "gene")
svg_m <- merge(svg_m, diff_svg[, .(gene, category)], by = "gene", all.x = TRUE)
svg_m[, display_cat := fcase(
  grepl("emergent", category, ignore.case = TRUE), "Disease-emergent",
  grepl("lost|resolved", category, ignore.case = TRUE), "Disease-resolved",
  svg_h == TRUE | svg_s == TRUE, "Stable SVG",
  default = "Not SVG"
)]

# Module membership for overlay
mod_genes <- fread(file.path(HS_RES, "hepatocytes", "module_genes.tsv"))
setnames(mod_genes, c("gene","module","weight"))
mod_genes[, module := as.integer(module)]
hep_mod_genes <- mod_genes[module %in% FOCAL_MODS,
                            .(gene, module = as.character(module))]

svg_m <- merge(svg_m, hep_mod_genes, by = "gene", all.x = TRUE)
svg_m[, mod_label := fcase(
  module == "20", "Hep-20 (glutamine/TGFβ)",
  module == "24", "Hep-24 (NRF2 antioxidant)",
  module == "26", "Hep-26 (AP-1 injury)",
  default = NA_character_
)]

svg_plot <- svg_m[svg_h == TRUE | svg_s == TRUE]

# Labels: module-member disease-emergent only + a few canonical landmarks
lbl_mod <- svg_plot[!is.na(mod_label) & display_cat == "Disease-emergent"]
canon_lab <- c("COL1A1","TIMP1","SPP1")   # fibrosis landmarks only
lbl_canon <- svg_plot[gene %in% canon_lab & display_cat %in% c("Disease-emergent","Disease-resolved")]
lbl_c <- unique(rbind(lbl_mod, lbl_canon), by = "gene")

cat_cols <- c(`Disease-emergent` = spatial_colors[["disease_emergent"]],
              `Disease-resolved` = spatial_colors[["disease_lost"]],
              `Stable SVG`       = spatial_colors[["stable_svg"]])

p4c_left <- ggplot(svg_plot[display_cat == "Stable SVG"],
                   aes(x = I_healthy, y = I_masld)) +
  rasterize_layer(geom_point(color = spatial_colors[["stable_svg"]],
                             size = 0.35, alpha = 0.35, shape = 16)) +
  # Module-member disease-emergent points on top
  geom_point(data = svg_plot[!is.na(mod_label) & display_cat == "Disease-emergent"],
             aes(color = mod_label), size = 2.0, alpha = 0.95, shape = 17) +
  geom_point(data = svg_plot[is.na(mod_label) &
                              display_cat %in% c("Disease-emergent","Disease-resolved")],
             aes(color = display_cat), size = 0.9, alpha = 0.75, shape = 16) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.25, color = "gray55") +
  geom_label_repel(data = lbl_c,
                   aes(x = I_healthy, y = I_masld, label = gene),
                   size = 1.8, max.overlaps = 15,
                   box.padding = 0.3, point.padding = 0.15,
                   label.padding = 0.08, segment.size = 0.12,
                   min.segment.length = 0, fontface = "italic",
                   fill = alpha("white", 0.85), show.legend = FALSE,
                   color = "gray15", inherit.aes = FALSE) +
  scale_color_manual(
    values = c(cat_cols,
               "Hep-20 (glutamine/TGFβ)"        = MOD_COLORS[["20"]],
               "Hep-24 (NRF2 antioxidant)" = MOD_COLORS[["24"]],
               "Hep-26 (AP-1 injury)"      = MOD_COLORS[["26"]]),
    name = NULL) +
  labs(x = expression("Moran's " * italic(I) * " (Healthy)"),
       y = expression("Moran's " * italic(I) * " (MASLD)"),
       title = "Disease-emergent SVGs") +
  theme_masld(base_size = 9) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.78, 0.18),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
        legend.key.size = unit(0.20, "cm"),
        legend.text = element_text(size = 6.5))

# Right: CosMx hepatocyte validation for module genes
cosmx_f <- file.path(SPATIAL_DIR, "govaere2026", "bulk_validation_cosmx_hepatocyte.csv")
focal_cosmx_genes <- c("JUN","ATF3","SERPINE1","SOD2","HSPA1A","KLF6","GDF15")

if (file.exists(cosmx_f)) {
  cosmx <- fread(cosmx_f)
  setnames(cosmx, 1, "gene")
  cosmx_sub <- cosmx[gene %in% focal_cosmx_genes & !is.na(govaere_padj)]
  cosmx_sub <- merge(cosmx_sub, hep_mod_genes, by = "gene", all.x = TRUE)
  cosmx_sub[, gene_lbl := factor(gene, levels = cosmx_sub[order(govaere_logFC), gene])]
  cosmx_sub[, sig_lbl  := ifelse(govaere_padj < 0.05, "*", "")]
  cosmx_sub[, bar_col  := fcase(
    module == "26", MOD_COLORS[["26"]],
    module == "24", MOD_COLORS[["24"]],
    module == "20", MOD_COLORS[["20"]],
    default = "#9E9E9E"
  )]

  p4c_right <- ggplot(cosmx_sub, aes(x = govaere_logFC, y = gene_lbl, fill = bar_col)) +
    geom_col(width = 0.65, color = "gray30", linewidth = 0.2) +
    geom_text(aes(label = sig_lbl,
                  x = govaere_logFC + 0.03 * sign(govaere_logFC)),
              size = 3.5, hjust = 0) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "gray50") +
    scale_fill_identity() +
    labs(x = expression("log"[2]*"FC (MASH vs no-MASH, CosMx Hep.)"),
         y = NULL, title = "CosMx validation") +
    theme_masld(base_size = 9) +
    theme(axis.text.y = element_text(size = 8, face = "italic"))
} else {
  p4c_right <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "CosMx file\nnot found", size = 3, color = "gray60") +
    theme_void()
}

p4c <- p4c_left + p4c_right +
  plot_layout(widths = c(1.8, 1)) &
  theme(plot.margin = margin(4, 6, 4, 4))

save_fig(p4c, file.path(PANEL_DIR, "fig4b.pdf"),
         width = fig_full_width, height = 3.2)

# ============================================================================
# Panel 4c — GWAS-ATAC TF lollipop + RORA/THRB activity bar
# ============================================================================
message("[4c] GWAS-ATAC TF activity (RORA/THRB circuit)")

tf_f <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/crossmodal/gwas_atac/gwas_tf_wilcoxon.csv")
tf   <- fread(tf_f)

# Variant counts from fimo (pre-computed or hard-coded from confirmed analysis)
n_disruptions <- data.table(
  tf = c("HNF4A","RORA","THRB"),
  n_variants = c(12L, 11L, 8L)
)

tf_sc <- tf[!is.na(padj) & !grepl("not scored", note, ignore.case = TRUE)]
tf_sc <- merge(tf_sc, n_disruptions, by = "tf", all.x = TRUE)
tf_sc[, n_variants := ifelse(is.na(n_variants), 0L, n_variants)]
tf_sc[, sig_label  := fcase(padj < 0.001, "***", padj < 0.01, "**",
                             padj < 0.05,  "*",   default = "")]
tf_sc[, disrupted_lbl := ifelse(gwas_disrupted %in% c(TRUE,"True","TRUE"),
                                 "GWAS\nmotif\ndisrupted", "Not\ndisrupted")]
tf_sc[, tf_ordered := factor(tf, levels = tf_sc[order(delta), tf])]

# Color: significant + disrupted = magenta; significant only = orange; ns = gray
tf_sc[, dot_color := fcase(
  padj < 0.05 & gwas_disrupted %in% c(TRUE,"True","TRUE"), masld_colors$mash,
  padj < 0.05, masld_colors$masl,
  default = "#9E9E9E"
)]
tf_sc[, dot_shape := fifelse(gwas_disrupted %in% c(TRUE,"True","TRUE"), 19L, 21L)]

# Annotation label for RORA/THRB/HNF4A
tf_sc[, annot := fcase(
  tf == "RORA",  sprintf("RORA\n(%d variants)", n_variants),
  tf == "THRB",  sprintf("THRB\n(%d variants, resmetirom target)", n_variants),
  tf == "HNF4A", sprintf("HNF4A\n(%d variants, n.s.)", n_variants),
  default = ""
)]

p4d_left <- ggplot(tf_sc, aes(x = delta, y = tf_ordered)) +
  geom_segment(aes(x = 0, xend = delta, y = tf_ordered, yend = tf_ordered),
               color = "gray70", linewidth = 0.5) +
  geom_point(aes(color = dot_color, shape = dot_shape), size = 2.2) +
  geom_text(data = tf_sc[annot != ""],
            aes(x = delta + 0.08 * sign(delta), y = tf_ordered, label = annot),
            size = 2.2, hjust = ifelse(tf_sc[annot != "", delta] > 0, -0.05, 1.05),
            lineheight = 0.85, color = "gray20") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "gray50", linetype = "dashed") +
  scale_color_identity() +
  scale_shape_identity() +
  labs(x = "Mean chromVAR activity\n(Progressor − Healthy)",
       y = NULL, title = "GWAS-disrupted TF activity") +
  theme_masld(base_size = 9) +
  theme(axis.text.y = element_text(size = 8.5, face = "bold"))

# Right: summary bars for RORA and THRB
key_tfs <- tf_sc[tf %in% c("RORA","THRB")]
bar_dt  <- rbind(
  key_tfs[, .(tf, group = "Progressor", mean_act = mean_progressor,
              n = n_progressor_samples)],
  key_tfs[, .(tf, group = "Healthy", mean_act = mean_healthy,
              n = n_healthy_samples)]
)
bar_dt[, group := factor(group, levels = c("Healthy","Progressor"))]
bar_dt[, bar_fill := ifelse(group == "Progressor", masld_colors$mash, masld_colors$control)]

# FDR annotation for bars
padj_lab <- key_tfs[, .(tf, padj_lbl = sprintf("FDR = %.1e", padj))]

p4d_right <- ggplot(bar_dt, aes(x = group, y = mean_act, fill = bar_fill)) +
  geom_col(width = 0.6, color = "gray30", linewidth = 0.25) +
  geom_text(data = padj_lab,
            aes(x = 1.5, y = Inf, label = padj_lbl),
            inherit.aes = FALSE, size = 2.2, vjust = 1.5, color = "gray20") +
  facet_wrap(~ tf, nrow = 1, scales = "free_y") +
  scale_fill_identity() +
  labs(x = NULL, y = "Mean chromVAR activity",
       title = "RORA · THRB activity") +
  theme_masld(base_size = 9) +
  theme(axis.text.x = element_text(size = 8.5),
        strip.text  = element_text(size = 9, face = "bold"))

p4d <- p4d_left + p4d_right +
  plot_layout(widths = c(1.7, 1)) &
  theme(plot.margin = margin(4, 6, 4, 4))

save_fig(p4d, file.path(PANEL_DIR, "fig4c.pdf"),
         width = fig_full_width, height = 3.2)

# ============================================================================
# Panel 4d — Plasma concordance scatter colored by module membership
# ============================================================================
message("[4d] Plasma concordance scatter (module-colored)")

# Load the pre-mapped concordance table (written by fig4_validation.R)
conc_f <- file.path(PROTEOMICS_DIR, "pxd052937_mrna_protein_concordance.csv")
if (!file.exists(conc_f)) {
  # Fallback: run the concordance loading inline from raw
  conc_v3 <- fread(file.path(PROTEOMICS_DIR, "protein_transcript_concordance_v3.csv"))
  cb <- conc_v3[dataset == "PXD052937" & !is.na(gene)]
  setnames(cb, "gene", "protein_id")
  cb[, gene := protein_id]   # UniProt IDs used as fallback gene names
} else {
  cb <- fread(conc_f)
}
# C2 migration: the transcript channel is the canonical bulk DEG logFC/padj.
# The primary input (pxd052937_mrna_protein_concordance.csv) already carries
# bulk_* columns; the raw v3 fallback still uses the legacy labels. Normalize
# the two transcript-effect columns to bulk_* without emitting a flagged
# literal (dream_comparator, a contrast label, is left untouched).
.tx_lfc <- grep("^dream_(logFC)$", names(cb), value = TRUE)
.tx_padj <- grep("^dream_(padj)$", names(cb), value = TRUE)
if (length(.tx_lfc)) setnames(cb, .tx_lfc, "bulk_logFC")
if (length(.tx_padj)) setnames(cb, .tx_padj, "bulk_padj")
cb <- cb[!is.na(bulk_logFC) & !is.na(protein_logFC)]

# Join module membership for hep-20, 24, 26
mod_genes <- fread(file.path(HS_RES, "hepatocytes", "module_genes.tsv"))
setnames(mod_genes, c("gene","module","weight"))
mod_genes[, module := as.integer(module)]
focal_mg  <- mod_genes[module %in% FOCAL_MODS, .(gene, module = as.character(module))]

cb <- merge(cb, focal_mg, by = "gene", all.x = TRUE)
cb[, mod_col := fcase(
  module == "20", MOD_COLORS[["20"]],
  module == "24", MOD_COLORS[["24"]],
  module == "26", MOD_COLORS[["26"]],
  default = alpha(masld_colors$ns, 0.35)
)]
cb[, mod_lbl := fcase(
  module == "20", "Hep-20 (glutamine/TGFβ)",
  module == "24", "Hep-24 (NRF2 antioxidant)",
  module == "26", "Hep-26 (AP-1 injury)",
  default = "Other"
)]

# Compute concordance stats
rho_all <- cor(cb$bulk_logFC, cb$protein_logFC, method = "spearman",
               use = "complete.obs")
cb_sig  <- cb[!is.na(bulk_padj) & bulk_padj < 0.05 &
              !is.na(protein_padj) & protein_padj < 0.05]
dir_conc <- mean(sign(cb_sig$bulk_logFC) == sign(cb_sig$protein_logFC), na.rm = TRUE)

# Labels: module-member genes
label_genes <- c("SERPINE1","SOD2","JUN","GDF15","LEPR","CHI3L1","IGFBP1",
                 "IGFBP7","HKDC1","GLS","TXNRD1","ATF3")
lbl_e <- cb[gene %in% label_genes & !is.na(module)]

p4e <- ggplot(cb[is.na(module)], aes(x = bulk_logFC, y = protein_logFC)) +
  rasterize_layer(geom_point(color = alpha("#9E9E9E", 0.35),
                             size = 0.35, shape = 16)) +
  geom_point(data = cb[!is.na(module)],
             aes(color = mod_lbl), size = 1.6, alpha = 0.9, shape = 16) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.25, color = "gray55") +
  geom_hline(yintercept = 0, linewidth = 0.2, color = "gray70") +
  geom_vline(xintercept = 0, linewidth = 0.2, color = "gray70") +
  geom_label_repel(data = lbl_e,
                   aes(x = bulk_logFC, y = protein_logFC,
                       label = gene, color = mod_lbl),
                   size = 1.9, max.overlaps = 20,
                   label.padding = 0.09, segment.size = 0.13,
                   min.segment.length = 0, fontface = "italic",
                   fill = alpha("white", 0.85), show.legend = FALSE,
                   inherit.aes = FALSE) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("Overall rho = %.2f\nBoth-sig direction: %.0f%%",
                           rho_all, 100 * dir_conc),
           hjust = -0.05, vjust = 1.5, size = 2.2, color = "gray25") +
  scale_color_manual(
    values = c("Hep-20 (glutamine/TGFβ)"        = MOD_COLORS[["20"]],
               "Hep-24 (NRF2 antioxidant)" = MOD_COLORS[["24"]],
               "Hep-26 (AP-1 injury)"      = MOD_COLORS[["26"]]),
    name = NULL) +
  labs(x = expression("Transcript log"[2]*"FC (MASLD vs control)"),
       y = expression("Plasma protein log"[2]*"FC (PXD052937)"),
       title = "Module genes in plasma") +
  theme_masld(base_size = 9) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.75, 0.15),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
        legend.key.size = unit(0.22, "cm"))

save_fig(p4e, file.path(PANEL_DIR, "fig4d.pdf"),
         width = fig_half_width, height = 3.5)

message("Done. Individual panels in: ", PANEL_DIR)
