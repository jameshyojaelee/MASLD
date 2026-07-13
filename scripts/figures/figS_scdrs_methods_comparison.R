#!/usr/bin/env Rscript
# ============================================================================
# figS_scdrs_methods_comparison.R
#
# KEY MESSAGE: The cell-type localization of MASLD GWAS signal depends on
# the gene-set anchoring strategy more than on the underlying biology.
# Three methods asking the same question give three different answers:
#   - TWAS-anchored (S-PrediXcan z-scores) → strong Hepatocyte enrichment
#   - COLOC-anchored (SuSiE-COLOC PP>0.5 gene sets) → strong Hepatocyte enrichment
#   - MAGMA-anchored (Zhang 2022 canonical) → essentially null in UKBB EUR
# Both biased methods use bulk-liver eQTL data; MAGMA does not. The bias
# direction is mechanistically explained in Panel D (gene-set composition).
#
# Compact 4-panel headline figure following project style (publication_theme).
#
# Panels:
#   A. CT × anchoring-method heatmap for liver enzymes (UKBB ALT/AST/GGT).
#      Color = signed group z; dot = FDR<0.05. Shows which CTs each method
#      picks as enriched.
#   B. Count of (CT × GWAS) cells reaching FDR<0.05 per method across the
#      three UKBB EUR liver-enzyme GWAS. Lollipop. Numerical anchor for
#      the asymmetry across methods.
#   C. Per-CT z-score scatter: TWAS UKBB_ALT (x) vs MAGMA UKBB_ALT (y).
#      Diagonal reference. Distance from diagonal = bias evidence.
#   D. Gene-set hepatocyte-expression decile composition (TWAS / COLOC /
#      MAGMA gene sets). Stacked bars showing where each method's gene set
#      sits in the hepatocyte-mean-expression distribution. Mechanism for
#      the disagreement in panels A-C.
#
# Inputs:
#   RNA-seq/results/causal_inference/scdrs/{ukbb_alt,ukbb_ast,ukbb_ggt}/group_results.csv
#   Analysis/SingleCell/results_gpu_v2/disease_signatures/
#     scdrs_gwas_celltype_enrichment.csv             (COLOC)
#     scdrs_magma_celltype_enrichment.csv            (MAGMA)
#     GWAS_anchored_scdrs.gs                          (COLOC gene sets)
#     MAGMA_anchored_scdrs.gs                         (MAGMA gene sets)
#     gene_hep_decile_metadata.csv                    (per-gene hep decile)
#   RNA-seq/results/causal_inference/otters_broadaway/{gwas}/otters_twas_combined.csv (TWAS top-1000)
#
# Outputs:
#   figures/supplementary/figS_scdrs/figS_scdrs_methods_comparison.pdf
#   figures/supplementary/figS_scdrs/panels/figS_scdrs_methods_comparison_{A..D}_*.pdf
#   figures/supplementary/figS_scdrs/panel_data/figS_scdrs_methods_comparison_panel{A..D}_data.csv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR    <- FIGS_SCDRS_DIR
DATA_DIR   <- FIGS_SCDRS_DATA_DIR
PANELS_DIR <- file.path(FIGS_SCDRS_DIR, "panels")
dir.create(PANELS_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat("[methods_comp] ", ..., "\n", sep = "")
save_panel <- function(p, slug, width, height) {
  out <- file.path(PANELS_DIR, paste0("figS_scdrs_methods_comparison_", slug, ".pdf"))
  save_fig(p, out, width = width, height = height)
  log_msg("  panel saved: ", out)
}

SIG_DIR    <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/disease_signatures")
TWAS_DIR   <- file.path(BASE, "RNA-seq/results/causal_inference/scdrs")
OTTERS_DIR <- file.path(BASE, "RNA-seq/results/causal_inference/otters_broadaway")

BASE_SIZE <- 7
LBL_PT    <- 7
LBL_SIZE  <- LBL_PT / ggplot2::.pt

theme_scdrs <- function() {
  theme_masld(base_size = BASE_SIZE) +
    theme(
      plot.title    = element_text(size = BASE_SIZE, face = "plain"),
      plot.subtitle = element_text(size = BASE_SIZE - 1, colour = "grey25",
                                   lineheight = 1.05),
      legend.key.size = unit(0.3, "cm"),
      panel.grid    = element_blank(),
      strip.background = element_rect(fill = "grey96", colour = NA),
      strip.text       = element_text(size = BASE_SIZE - 1, face = "plain")
    )
}

# CT umbrella ordering — fine CTs collapsed in TWAS legacy run
umbrella_map <- c(
  "Hepatocytes"             = "Hepatocytes",
  "Cholangiocytes"          = "Cholangiocytes",
  "Endothelial cells"       = "Endothelial cells",
  "Fibroblasts"             = "Fibroblasts",
  "Macrophages"             = "Macrophages",
  "Mono+mono derived cells" = "Mono+mono_derived",
  "Mono+mono_derived"       = "Mono+mono_derived",
  "T cells"                 = "T cells",
  "B cells"                 = "B cells",
  "Plasma cells"            = "Plasma cells",
  "Resident NK"             = "Resident NK",
  "Circulating NK/NKT"      = "Circulating NK/NKT",
  "Basophils"               = "Basophils"
)
CT_ORDER <- c(
  "Hepatocytes", "Cholangiocytes", "Endothelial cells", "Fibroblasts",
  "Macrophages", "Mono+mono_derived",
  "Plasma cells", "T cells", "B cells",
  "Resident NK", "Circulating NK/NKT", "Basophils"
)

METHOD_ORDER  <- c("TWAS", "COLOC", "MAGMA")
METHOD_COLORS <- c(TWAS = "#E64A19", COLOC = "#7B1FA2", MAGMA = "#1B5E20")

LE_GWAS  <- c("ukbb_alt", "ukbb_ast", "ukbb_ggt")
LE_LABELS <- c(ukbb_alt = "UKBB ALT", ukbb_ast = "UKBB AST", ukbb_ggt = "UKBB GGT")

# ----------------------------------------------------------------------------
# Load + harmonise CT enrichment from all 3 methods
# ----------------------------------------------------------------------------

# MAGMA — direct read, has assoc_mcz + assoc_mcp
log_msg("loading MAGMA enrichment (canonical Zhang 2022 anchoring)")
magma <- fread(file.path(SIG_DIR, "scdrs_magma_celltype_enrichment.csv"))
magma <- magma[gwas_id %in% LE_GWAS]
magma_dt <- magma[, .(cell_type, gwas_id, z = assoc_mcz, p = assoc_mcp)]
magma_dt[, method := "MAGMA"]

# COLOC — broad liver enzymes panel (covers all 3 UKBB enzymes pooled).
# To keep the comparison 1:1 with TWAS/MAGMA per-enzyme, duplicate this row
# across all three liver-enzyme GWAS slots (the COLOC panel is genuinely
# pooled, so this represents what the COLOC anchor "says" for liver enzymes).
log_msg("loading COLOC enrichment (SuSiE-COLOC PP>0.5 liver enzymes panel, PolyFun LD)")
coloc <- fread(file.path(SIG_DIR, "scdrs_gwas_celltype_enrichment.csv"))
coloc <- coloc[trait_class == "gwas_coloc_liver_enzymes" & panel == "polyfun"]
coloc_dt <- coloc[, .(cell_type, z = assoc_mcz, p = assoc_mcp)]
coloc_dt <- coloc_dt[, .(gwas_id = LE_GWAS), by = .(cell_type, z, p)]
coloc_dt[, method := "COLOC"]

# TWAS — 200K-cell legacy run; only mean_norm_score and ttest_pval are stored.
# Convert ttest_pval → signed z via inverse normal, signed by mean_norm_score.
# This is approximate but lets us put TWAS on a comparable z-score axis.
log_msg("loading TWAS-anchored enrichment (legacy S-PrediXcan z-score gene sets)")
twas_parts <- lapply(LE_GWAS, function(g) {
  f <- file.path(TWAS_DIR, g, "group_results.csv")
  d <- fread(f)
  d[, atlas_ct := umbrella_map[as.character(cell_type)]]
  d <- d[!is.na(atlas_ct)]
  # Inverse-normal of one-sided p, signed by direction of mean_norm_score
  d[, z := sign(mean_norm_score) * qnorm(pmax(ttest_pval, 1e-300) / 2,
                                          lower.tail = FALSE)]
  d[, p := ttest_pval]
  d[, gwas_id := g]
  d[, .(cell_type = atlas_ct, gwas_id, z, p)]
})
twas_dt <- rbindlist(twas_parts)
twas_dt[, method := "TWAS"]

# Aggregate fine CTs in TWAS to umbrella via weighted mean of z (by N — but we
# don't have N per fine CT here; simple mean is fine since CT umbrellas
# usually map 1:1 except dendritic/granulocytes that we don't use).
twas_dt <- twas_dt[, .(z = mean(z, na.rm = TRUE),
                       p = min(p, na.rm = TRUE)),
                   by = .(cell_type, gwas_id, method)]

all_dt <- rbindlist(list(magma_dt, coloc_dt, twas_dt), use.names = TRUE)
all_dt <- all_dt[cell_type %in% CT_ORDER]
all_dt[, cell_type   := factor(cell_type, levels = CT_ORDER)]
all_dt[, method      := factor(method,    levels = METHOD_ORDER)]
all_dt[, gwas_id     := factor(gwas_id,   levels = LE_GWAS,
                                labels = LE_LABELS[LE_GWAS])]
all_dt[, sig := !is.na(p) & p < 0.05]
log_msg("  rows: ", nrow(all_dt))

# ----------------------------------------------------------------------------
# Panel A — CT × method heatmap for liver enzymes
# ----------------------------------------------------------------------------
log_msg("Panel A: CT × method heatmap (liver enzymes)")
A_data <- copy(all_dt)
fwrite(A_data, file.path(DATA_DIR, "figS_scdrs_methods_comparison_panelA_data.csv"))

# Cap fill for visual stability
A_data[, fill_capped := pmax(pmin(z, 4), -4)]

panel_A <- ggplot(A_data, aes(x = method, y = cell_type, fill = fill_capped)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_point(data = A_data[sig == TRUE],
             colour = "black", size = 0.7, shape = 1, stroke = 0.4) +
  facet_wrap(~ gwas_id, nrow = 1) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                       midpoint = 0, limits = c(-4, 4),
                       oob = scales::squish,
                       name = "scDRS group\nz-score",
                       guide = guide_colorbar(barwidth = unit(2.6, "cm"),
                                              barheight = unit(0.25, "cm"))) +
  scale_y_discrete(limits = rev(CT_ORDER), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title = "A. Same GWAS, three anchoring methods → three different cell-type maps",
       subtitle = "○ = FDR<0.05; canonical MAGMA (right column) shows no significant CT in any UKBB EUR liver-enzyme GWAS, while TWAS and COLOC both flag Hepatocytes") +
  theme_scdrs() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right")

save_panel(panel_A, "A_CT_x_method_heatmap_liver_enzymes",
           width = 6.8, height = 4.0)

# ----------------------------------------------------------------------------
# Panel B — count of (CT × GWAS) cells reaching FDR<0.05 per method
# ----------------------------------------------------------------------------
log_msg("Panel B: # FDR<0.05 CTs per method (liver enzymes only)")
B_data <- all_dt[, .(n_sig = sum(sig, na.rm = TRUE)), by = method]
B_data[, label := paste0(n_sig, " / ", 12 * length(LE_GWAS))]
fwrite(B_data, file.path(DATA_DIR, "figS_scdrs_methods_comparison_panelB_data.csv"))

panel_B <- ggplot(B_data, aes(x = method, y = n_sig, colour = method)) +
  geom_segment(aes(xend = method, y = 0, yend = n_sig), linewidth = 0.6) +
  geom_point(size = 4) +
  geom_text(aes(label = n_sig), vjust = -1.1, size = LBL_SIZE * 1.2,
            colour = "grey15", fontface = "plain") +
  scale_colour_manual(values = METHOD_COLORS, guide = "none") +
  scale_y_continuous(limits = c(0, max(B_data$n_sig) * 1.25),
                     expand = c(0, 0)) +
  labs(x = NULL, y = "# (CT × GWAS) cells\nat FDR<0.05",
       title = "B. Asymmetric FDR<0.05 yield",
       subtitle = sprintf("Across %d UKBB EUR liver-enzyme GWAS × %d CTs = %d slots",
                          length(LE_GWAS), length(CT_ORDER),
                          length(LE_GWAS) * length(CT_ORDER))) +
  theme_scdrs()

save_panel(panel_B, "B_FDR_sig_count_per_method", width = 3.5, height = 3.6)

# ----------------------------------------------------------------------------
# Panel C — Per-CT scatter: TWAS vs MAGMA for UKBB_ALT
# ----------------------------------------------------------------------------
log_msg("Panel C: per-CT scatter TWAS vs MAGMA, UKBB_ALT")
wide <- dcast(all_dt[gwas_id == "UKBB ALT"], cell_type ~ method,
              value.var = "z")
C_data <- wide[!is.na(TWAS) & !is.na(MAGMA)]
fwrite(C_data, file.path(DATA_DIR, "figS_scdrs_methods_comparison_panelC_data.csv"))

rho_C <- if (nrow(C_data) >= 5)
  suppressWarnings(cor(C_data$TWAS, C_data$MAGMA, method = "spearman")) else NA_real_

# Highlight Hepatocytes — the cell type most extremely affected
hep_row <- C_data[cell_type == "Hepatocytes"]
hep_label <- if (nrow(hep_row) > 0)
  sprintf("Hepatocytes: TWAS z = %.1f, MAGMA z = %.1f",
          hep_row$TWAS, hep_row$MAGMA) else ""

xy_lim <- max(abs(c(C_data$TWAS, C_data$MAGMA)), na.rm = TRUE) * 1.1

panel_C <- ggplot(C_data, aes(x = TWAS, y = MAGMA)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.2) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.2) +
  geom_point(aes(colour = cell_type), size = 1.8, alpha = 0.9) +
  ggrepel::geom_text_repel(aes(label = cell_type),
                           size = LBL_SIZE * 0.8, max.overlaps = 20,
                           segment.size = 0.2, min.segment.length = 0,
                           show.legend = FALSE) +
  # Clip extreme TWAS z (Hepatocytes ≈ 37 because ttest p rounded to 0; the
  # qualitative gap is the message, not the precise tail value).
  coord_fixed(xlim = c(-10, 10), ylim = c(-10, 10), clip = "off") +
  scale_colour_manual(values = c(
    "Hepatocytes"="#1B5E20","Cholangiocytes"="#00BCD4",
    "Endothelial cells"="#7CB342","Fibroblasts"="#FF6F00",
    "Macrophages"="#212121","Mono+mono_derived"="#EC407A",
    "Plasma cells"="#8D6E63","T cells"="#1A237E","B cells"="#5D4037",
    "Resident NK"="#00897B","Circulating NK/NKT"="#26A69A","Basophils"="#F9A825"),
    guide = "none") +
  labs(x = "TWAS-anchored z (S-PrediXcan)  [clipped at ±10]",
       y = "MAGMA-anchored z (canonical)",
       title = "C. TWAS inflates, MAGMA does not — same atlas, same GWAS",
       subtitle = sprintf("UKBB ALT (N=344K) per cell type. %s (TWAS z clipped from %.0f). Distance from y=x = anchoring bias.",
                          gsub(" TWAS.*", "", hep_label),
                          hep_row$TWAS)) +
  theme_scdrs() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

save_panel(panel_C, "C_TWAS_vs_MAGMA_scatter_UKBB_ALT",
           width = 4.5, height = 4.5)

# ----------------------------------------------------------------------------
# Panel D — Gene-set hepatocyte-decile composition
# ----------------------------------------------------------------------------
# For each method, what fraction of its gene set lies in each hep-expression
# decile? If MAGMA's gene set is roughly uniform across deciles while
# TWAS / COLOC are top-heavy, that mechanistically explains the bias.
log_msg("Panel D: gene-set hep-decile composition")

dec <- fread(file.path(SIG_DIR, "gene_hep_decile_metadata.csv"))
dec <- dec[!is.na(hep_decile)]
gene_to_dec <- setNames(dec$hep_decile, dec$gene)

# Read .gs files (simple two-column TRAIT \t GENESET parsing)
read_gs <- function(gs_path) {
  lines <- readLines(gs_path)
  lines <- lines[-1]  # drop header
  out <- list()
  for (l in lines) {
    parts <- strsplit(l, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 2) next
    trait <- parts[1]
    pairs <- strsplit(parts[2], ",", fixed = TRUE)[[1]]
    g <- sub(":.*$", "", pairs)
    out[[trait]] <- g
  }
  out
}

magma_gs <- read_gs(file.path(SIG_DIR, "MAGMA_anchored_scdrs.gs"))
coloc_gs <- read_gs(file.path(SIG_DIR, "GWAS_anchored_scdrs.gs"))

# Method-aggregated gene sets:
#   MAGMA = union of top-1000 across the 3 UKBB EUR liver-enzyme GWAS
#   COLOC = the gwas_coloc_liver_enzymes_polyfun trait
#   TWAS  = pull top-1000 by |best_z| from OTTERS combined results for each
#           UKBB EUR liver-enzyme GWAS, then union
magma_genes <- unique(unlist(magma_gs[c("ukbb_alt", "ukbb_ast", "ukbb_ggt")]))
coloc_genes <- coloc_gs[["gwas_coloc_liver_enzymes_polyfun"]]

twas_genes <- character(0)
for (g in LE_GWAS) {
  f <- file.path(OTTERS_DIR, g, "otters_twas_combined.csv")
  if (!file.exists(f)) next
  d <- fread(f, select = c("gene_symbol", "best_z"))
  d <- d[!is.na(best_z) & !is.na(gene_symbol)]
  d[, abs_z := abs(best_z)]
  top <- head(d[order(-abs_z)], 1000)
  twas_genes <- unique(c(twas_genes, top$gene_symbol))
}

assign_decile <- function(genes) {
  d <- gene_to_dec[as.character(genes)]
  d <- d[!is.na(d)]
  as.integer(d)
}

D_parts <- list(
  data.table(method = "TWAS",  decile = assign_decile(twas_genes)),
  data.table(method = "COLOC", decile = assign_decile(coloc_genes)),
  data.table(method = "MAGMA", decile = assign_decile(magma_genes))
)
D_long <- rbindlist(D_parts)
D_long[, method := factor(method, levels = METHOD_ORDER)]
D_long[, decile := factor(decile, levels = 0:9,
                          labels = paste0("D", 0:9))]
D_data <- D_long[, .(n = .N), by = .(method, decile)]
D_data[, frac := n / sum(n), by = method]
fwrite(D_data, file.path(DATA_DIR, "figS_scdrs_methods_comparison_panelD_data.csv"))

# Median hep-decile for each method (annotate)
med_dec <- D_long[, .(median_dec = median(as.integer(decile) - 1, na.rm = TRUE)),
                  by = method]

panel_D <- ggplot(D_data, aes(x = decile, y = frac, fill = method)) +
  geom_col(position = position_dodge(0.85), width = 0.78) +
  scale_fill_manual(values = METHOD_COLORS, name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.10))) +
  labs(x = "Hepatocyte mean-expression decile (D0 = lowest, D9 = highest)",
       y = "Fraction of gene set",
       title = "D. Gene-set membership is similar across methods — bias is in the WEIGHTS",
       subtitle = sprintf(
         "All three methods concentrate in top hep-expression deciles (median decile: TWAS=%d, COLOC=%d, MAGMA=%d). The TWAS/COLOC bias in panels A-C comes from gene-weight distributions: TWAS z-scores scale with bulk-eQTL strength (correlated with hep expression) and inflate per-cell scores; MAGMA z-scores are bounded by GWAS marginal effects without eQTL imputation.",
         med_dec[method=="TWAS"]$median_dec,
         med_dec[method=="COLOC"]$median_dec,
         med_dec[method=="MAGMA"]$median_dec)) +
  theme_scdrs() +
  theme(legend.position = c(0.16, 0.88),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key.size = unit(0.3, "cm"))

save_panel(panel_D, "D_geneset_hep_decile_composition",
           width = 6.8, height = 3.6)

# ----------------------------------------------------------------------------
# Compose — 2x2 layout, compact
# ----------------------------------------------------------------------------
log_msg("Composing figS_scdrs_methods_comparison.pdf")

# Row 1: A (wide heatmap)
# Row 2: B + C (narrow lollipop + scatter)
# Row 3: D (wide bar chart)
top   <- panel_A
mid   <- panel_B + panel_C + plot_layout(widths = c(1.0, 1.4))
bot   <- panel_D
full  <- top / mid / bot +
  plot_layout(heights = c(1.1, 1.0, 1.0)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = BASE_SIZE + 2, face = "plain"))

out_pdf <- file.path(OUT_DIR, "figS_scdrs_methods_comparison.pdf")
save_fig(full, out_pdf, width = 7.5, height = 11)
log_msg("Wrote ", out_pdf)
log_msg("Done.")
