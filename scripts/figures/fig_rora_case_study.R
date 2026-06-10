#!/usr/bin/env Rscript
# fig_rora_case_study.R — RORA multi-modal case study (Panels B-F)
# Panel A is the locus zoom from figS09_locus_zoom.R (fig3b)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- FIGS_RORA_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── shared aesthetics ─────────────────────────────────────────────────────────
BLUE   <- "#1565C0"
RED    <- "#C2185B"
COL_F  <- masld_colors$female
COL_M  <- masld_colors$male
COL_NS <- "#CFD8DC"

PANEL_THEME <- theme_masld(base_size = 8) +
  theme(
    plot.title       = element_text(size = 8.5, face = "bold", hjust = 0,
                                    margin = margin(b = 3)),
    plot.subtitle    = element_text(size = 6.8, colour = "#555555",
                                    margin = margin(b = 5)),
    axis.title       = element_text(size = 7.5),
    axis.text        = element_text(size = 7),
    legend.text      = element_text(size = 6.5),
    legend.title     = element_text(size = 7),
    legend.key.size  = unit(0.32, "cm"),
    panel.grid.minor = element_blank()
  )

# =============================================================================
# Panel B: Cross-trait × cross-ancestry COLOC tile
# =============================================================================
cat("-- Panel B: Comprehensive COLOC --\n")

coloc_raw <- data.frame(
  trait    = c("GGT","GGT","GGT","GGT",
               "AST","AST","AST","AST",
               "ALT","ALT","ALT","ALT",
               "NAFLD","NAFLD","NAFLD","NAFLD"),
  ancestry = c("EUR","EAS","CSA","AFR",
               "EUR","EAS","CSA","AFR",
               "EUR","EAS","CSA","AFR",
               "EUR","EAS","CSA","AFR"),
  pp4      = c(0.995, 0.997, 0.995, 0.012,  # GGT
               0.989, 0.353, 0.079, 0.006,  # AST
               0.997, 0.166, 0.161, 0.006,  # ALT
               0.777,    NA,    NA,    NA)   # NAFLD (UKBB EUR only)
)
coloc_raw$trait    <- factor(coloc_raw$trait,
                             levels = rev(c("GGT","AST","ALT","NAFLD")))
coloc_raw$ancestry <- factor(coloc_raw$ancestry,
                             levels = c("EUR","EAS","CSA","AFR"))
coloc_raw$tier <- with(coloc_raw,
                       ifelse(is.na(pp4),           "no data",
                       ifelse(pp4 >= 0.8,           "strong",
                       ifelse(pp4 >= 0.5,           "suggestive", "ns"))))
coloc_raw$label <- with(coloc_raw,
                        ifelse(is.na(pp4), "–",
                               sprintf("%.2f", pp4)))

pB <- ggplot(coloc_raw, aes(x = ancestry, y = trait, fill = pp4)) +
  geom_tile(colour = "white", linewidth = 1.4) +
  geom_text(aes(label = label,
                colour = pp4 > 0.5,
                fontface = ifelse(!is.na(pp4) & pp4 >= 0.8, "bold", "plain")),
            size = 2.7) +
  scale_fill_gradientn(
    colours = c("#FFFFFF", "#E3F2FD", "#90CAF9", BLUE, "#0D47A1"),
    values  = rescale(c(0, 0.3, 0.6, 0.8, 1)),
    limits  = c(0, 1), na.value = "#EEEEEE",
    name    = "COLOC\nPP4", breaks = c(0, 0.5, 0.8, 1.0)
  ) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#555555"),
                      guide = "none") +
  scale_x_discrete(labels = c(EUR = "EUR\n(UKBB)", EAS = "EAS\n(BBJ)",
                              CSA = "CSA\n(PanUKBB)", AFR = "AFR\n(PanUKBB)")) +
  labs(x = NULL, y = "GWAS trait",
       title = "Cross-trait, cross-ancestry colocalization",
       subtitle = "Liver enzymes replicate in EUR; GGT replicates in 3 of 4 ancestries") +
  coord_equal() +
  PANEL_THEME +
  theme(
    panel.grid.major = element_blank(),
    legend.position  = "right",
    axis.text.x      = element_text(lineheight = 0.9)
  )

ggsave(file.path(OUT_DIR, "rora_crossancestry_coloc_pp4.pdf"), pB,
       width = 4.5, height = 3.2, useDingbats = FALSE)
cat("  Saved rora_crossancestry_coloc_pp4\n")

# =============================================================================
# Panel C: GWAS-ATAC motif disruption — clean horizontal bar
# =============================================================================
cat("-- Panel C: GWAS-ATAC motif disruption --\n")

motif_raw <- fread(file.path(BASE,
  "GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv"))
rora_m <- motif_raw[tf_name == "RORA",
                    .(alleleDiff = alleleDiff[1],
                      max_pip    = max(max_pip, na.rm = TRUE),
                      priority_score = max(priority_score, na.rm = TRUE)),
                    by = SNP_id]
# Only keep variants with credible GWAS causal support (PIP >= 0.2)
rora_m <- rora_m[max_pip >= 0.2]
rora_m <- rora_m[order(max_pip)]
cat("  Variants retained (PIP >= 0.2):", nrow(rora_m), "\n")
rora_m[, snp_label := gsub("^(\\d+:\\d+).+", "chr\\1", SNP_id)]
rora_m[, snp_lab_f := factor(SNP_id, levels = SNP_id)]
# Disruption direction
rora_m[, direction := ifelse(alleleDiff < 0,
                              "weaker RORA binding",
                              "stronger RORA binding")]

pC <- ggplot(rora_m, aes(x = alleleDiff, y = snp_lab_f, fill = max_pip)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "black") +
  geom_col(width = 0.65, colour = NA) +
  # label the GWAS max_pip on the right (most causal variants get emphasis)
  geom_text(aes(x = ifelse(alleleDiff > 0, alleleDiff + 0.08, alleleDiff - 0.08),
                label = sprintf("PIP %.3f", max_pip),
                hjust = ifelse(alleleDiff > 0, 0, 1)),
            size = 2.2, colour = "#333333") +
  scale_fill_gradientn(
    colours = c("#CFD8DC", "#90CAF9", BLUE),
    limits = c(0, 1), name = "GWAS\nmax PIP",
    breaks = c(0, 0.5, 1)
  ) +
  scale_y_discrete(labels = function(x) gsub("^(\\d+:\\d+).+", "chr\\1", x)) +
  scale_x_continuous(
    limits = c(-2.6, 2.6),
    breaks = c(-2, -1, 0, 1, 2),
    labels = c("-2", "-1", "0", "+1", "+2")
  ) +
  labs(x = "Motif score change (alt − ref allele)", y = NULL,
       title = "MASLD variants disrupt RORA binding at distal loci",
       subtitle = "3 credible MASLD variants (PIP >= 0.2, non-RORA loci) · regulon padj = 0.033") +
  PANEL_THEME +
  theme(
    axis.text.y    = element_text(size = 6.5, family = "mono"),
    legend.position = c(0.95, 0.3),
    legend.background = element_blank()
  )

ggsave(file.path(OUT_DIR, "rora_gwas_atac_motif.pdf"), pC,
       width = 4.5, height = 3.0, useDingbats = FALSE)
cat("  Saved rora_gwas_atac_motif\n")

# =============================================================================
# Panel D: scRNA cell-type expression — where is RORA expressed?
# =============================================================================
cat("-- Panel D: Cell-type expression --\n")

ct_files <- list.files(
  file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudobulk_de"),
  pattern = "_de\\.csv$", full.names = TRUE)
ct_df <- rbindlist(lapply(ct_files, function(f) {
  ct <- gsub("_de\\.csv", "", basename(f))
  if (ct == "allcell_pseudobulk") return(NULL)
  x <- fread(f)
  if (!("AveExpr" %in% names(x) && "gene" %in% names(x))) return(NULL)
  rora <- x[gene == "RORA"]
  if (nrow(rora) == 0) return(NULL)
  data.table(celltype = ct, ave_expr = rora$AveExpr[1])
}))
# Clean cell-type labels
ct_df[, ct_clean := gsub("_", " ", celltype)]
ct_df[, ct_clean := gsub("Mono\\+mono derived cells",
                          "Mono/mono-derived", ct_clean)]
ct_df[, is_hep := celltype == "Hepatocytes"]
ct_df <- ct_df[order(ave_expr)]
ct_df[, ct_clean := factor(ct_clean, levels = ct_clean)]

pD <- ggplot(ct_df, aes(x = ave_expr, y = ct_clean,
                         fill = is_hep)) +
  geom_col(width = 0.7, colour = NA) +
  geom_text(aes(label = sprintf("%.1f", ave_expr),
                x = ave_expr + 0.15),
            hjust = 0, size = 2.3, colour = "#333333") +
  scale_fill_manual(values = c(`TRUE` = BLUE, `FALSE` = COL_NS),
                    guide = "none") +
  scale_x_continuous(limits = c(0, 13),
                     breaks = c(0, 3, 6, 9, 12),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "RORA expression (scRNA pseudobulk, log-CPM)", y = NULL,
       title = "Hepatocyte-dominant expression",
       subtitle = "11 cell types from integrated liver scRNA atlas") +
  PANEL_THEME +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 6.8,
                                face = ifelse(levels(ct_df$ct_clean)
                                              == "Hepatocytes",
                                              "bold", "plain"),
                                colour = ifelse(levels(ct_df$ct_clean)
                                                == "Hepatocytes",
                                                BLUE, "black"))
  )

ggsave(file.path(OUT_DIR, "rora_celltype_expression.pdf"), pD,
       width = 4.2, height = 3.2, useDingbats = FALSE)
cat("  Saved rora_celltype_expression\n")

# =============================================================================
# Panel E: Disease progression + sex — simple vertical bars
# =============================================================================
cat("-- Panel E: Progression + sex --\n")

prog_df <- data.frame(
  group = factor(c("NAFL","NASH","F2\ninflection",
                   "Adv.\nFibrosis","Cirrhosis"),
                 levels = c("NAFL","NASH","F2\ninflection",
                            "Adv.\nFibrosis","Cirrhosis")),
  logFC = c(-0.402, -0.407, -0.236, -0.216, -0.240),
  padj  = c(0.014, 3.97e-08, 1.28e-09, 1.02e-06, 0.032),
  panel = "Disease progression"
)
prog_df$stars <- as.character(cut(prog_df$padj,
                                   breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                   labels = c("***","**","*","")))

sex_df <- data.frame(
  group = factor(c("Female","Male"),
                 levels = c("Female","Male")),
  logFC = c(-0.462, -0.242),
  padj  = c(NA, NA),
  panel = "Sex stratification"
)
sex_df$stars <- ""

pE_prog <- ggplot(prog_df, aes(x = group, y = logFC)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "black") +
  geom_col(fill = BLUE, width = 0.6, colour = NA) +
  geom_text(aes(label = stars, y = logFC - 0.015),
            vjust = 1, size = 3.3, colour = "#333333") +
  geom_text(aes(label = sprintf("%.2f", logFC),
                y = 0.01),
            vjust = 0, size = 2.4, colour = "#333333") +
  scale_y_continuous(limits = c(-0.55, 0.08),
                     breaks = c(-0.5, -0.3, -0.1, 0),
                     expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = NULL, y = "log2 fold change",
       title = "Disease progression",
       subtitle = "NAFL/NASH vs Ctrl (4 cohorts)\nF2/Adv Fib (7 cohorts) · Cirrhosis (6 cohorts)") +
  PANEL_THEME +
  theme(axis.text.x = element_text(size = 6.8, lineheight = 0.9),
        panel.grid.major.x = element_blank())

pE_sex <- ggplot(sex_df, aes(x = group, y = logFC,
                              fill = group)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "black") +
  geom_col(width = 0.6, colour = NA) +
  geom_text(aes(label = sprintf("%.2f", logFC), y = 0.01),
            vjust = 0, size = 2.4, colour = "#333333") +
  scale_fill_manual(values = c(Female = COL_F, Male = COL_M),
                    guide = "none") +
  scale_y_continuous(limits = c(-0.55, 0.08),
                     breaks = c(-0.5, -0.3, -0.1, 0),
                     expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = NULL, y = NULL,
       title = "Sex stratification",
       subtitle = "delta = 0.22 (interaction p = 0.28)") +
  PANEL_THEME +
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_blank())

pE <- pE_prog + pE_sex + plot_layout(widths = c(5, 2))

ggsave(file.path(OUT_DIR, "rora_progression_sex.pdf"), pE,
       width = 5.5, height = 3.0, useDingbats = FALSE)
cat("  Saved rora_progression_sex\n")

# =============================================================================
# Panel F: Visium — spatial LOCALIZATION of RORA (NOT a disease contrast)
# =============================================================================
# CRITICAL: the disease axis in this Visium data is fully confounded with study.
# All Healthy + Steatotic sections come from GSE192741 (JBO*); all MASLD-spectrum
# sections come from Vu 2025 (VLP*). There is zero within-study disease-vs-control
# overlap, and per-spot log1p normalization does NOT remove cross-study batch.
# A disease-vs-control comparison here is therefore uninterpretable. This panel is
# used for in-situ LOCALIZATION only; the disease direction for RORA (DOWN) is
# carried by the bulk RNA-seq (Panel E) — a within-design, cohort-adjusted contrast.
cat("-- Panel F: Visium spatial localization (no disease claim) --\n")

visium <- fread(file.path(BASE,
  "Analysis/Spatial/results/rora_case_study/rora_visium.csv"))

# Cohort of origin is the only honest grouping (JBO = GSE192741, VLP = Vu 2025).
visium[, cohort := ifelse(grepl("^JBO", sample_id), "GSE192741", "Vu 2025")]

# Representative, density-matched sections — one per cohort. NOT chosen to
# maximize a disease contrast (the prior "highest-mean MASLD" cherry-pick is gone).
rep_ids <- c("JBO018", "VLP119_A")
vis_rep <- visium[sample_id %in% rep_ids]
vis_rep[, facet_lab := factor(
  sprintf("%s · %s", cohort, sample_id),
  levels = c("GSE192741 · JBO018", "Vu 2025 · VLP119_A"))]
vis_rep[, x_norm := (x - min(x)) / (max(x) - min(x)), by = sample_id]
vis_rep[, y_norm := 1 - (y - min(y)) / (max(y) - min(y)), by = sample_id]

global_cap <- quantile(vis_rep$RORA, 0.99)
vis_rep[, expr_c := pmin(RORA, global_cap)]

viridis_pal <- c("#F7F7F7", "#440154", "#31688E", "#35B779", "#FDE725")

pF_spatial <- ggplot(vis_rep, aes(x = x_norm, y = y_norm, colour = expr_c)) +
  geom_point(size = 1.1, alpha = 0.9, shape = 16, stroke = 0) +
  scale_colour_gradientn(
    colours = viridis_pal,
    limits  = c(0, global_cap),
    name    = "RORA",
    breaks  = c(0, global_cap),
    labels  = c("low", "high")
  ) +
  facet_wrap(~ facet_lab, nrow = 1) +
  coord_equal() +
  labs(title = "RORA expression in human liver Visium",
       subtitle = "Patchy lobular hepatocyte expression · two independent cohorts") +
  theme_void(base_size = 8) +
  theme(
    strip.text        = element_text(size = 7.5, face = "bold",
                                     margin = margin(b = 2)),
    legend.position   = "right",
    legend.title      = element_text(size = 7, face = "bold"),
    legend.text       = element_text(size = 6.5),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width  = unit(0.22, "cm"),
    plot.title        = element_text(size = 8, face = "bold", hjust = 0.5,
                                     margin = margin(b = 2)),
    plot.subtitle     = element_text(size = 6.5, colour = "#555555", hjust = 0.5,
                                     margin = margin(b = 4)),
    plot.background   = element_rect(fill = "white", colour = NA)
  )

# Right: per-section mean RORA, grouped by COHORT (not disease). Showing the
# cross-cohort technical offset explicitly is the honest framing — RORA is
# robustly expressed in both datasets, and the between-cohort difference is a
# batch/platform effect, NOT a disease signal.
sample_means <- visium[, .(mean_rora = mean(RORA),
                           sd_rora   = sd(RORA),
                           n_spots   = .N),
                       by = .(sample_id, cohort)]
sample_means[, cohort := factor(cohort, levels = c("GSE192741", "Vu 2025"))]

cohort_cols <- c("GSE192741" = "#3949AB", "Vu 2025" = "#00897B")

pF_quant <- ggplot(sample_means, aes(x = cohort, y = mean_rora,
                                      fill = cohort)) +
  geom_boxplot(width = 0.55, alpha = 0.35, outlier.shape = NA,
               colour = "#444444", linewidth = 0.35) +
  geom_jitter(aes(colour = cohort, size = n_spots),
              width = 0.15, alpha = 0.85, shape = 16) +
  scale_fill_manual(values = cohort_cols, guide = "none") +
  scale_colour_manual(values = cohort_cols, guide = "none") +
  scale_size_continuous(range = c(1.5, 3.5), name = "n spots",
                        breaks = c(500, 1000, 2000)) +
  labs(x = NULL, y = "Mean RORA per section",
       title = "Per-section expression by cohort",
       subtitle = sprintf(paste0("Localization only (n = %d sections) · ",
                                  "cross-cohort difference is technical,\n",
                                  "not a disease contrast; see Panel E for ",
                                  "disease direction"),
                          length(unique(visium$sample_id)))) +
  PANEL_THEME +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position  = "none"
  )

pF <- pF_spatial + pF_quant + plot_layout(widths = c(2, 1.3))

ggsave(file.path(OUT_DIR, "rora_visium_spatial.pdf"), pF,
       width = 7.0, height = 3.2, useDingbats = FALSE)
cat("  Saved rora_visium_spatial\n")

# =============================================================================
# Individual panels only — no composite. The PI assembles Fig 3 in Illustrator.
# =============================================================================

# Clean up old orphan files
for (f in c("panelC_progression_logFC.pdf", "panelD_sex_logFC.pdf",
            "panelE_visium_spatial.pdf", "panelF_twas.pdf",
            "panelF_twas_hyprcoloc.pdf",
            "panelD_meta_subtype_activity.pdf",
            "panelF_visium_spatial.pdf")) {
  old <- file.path(OUT_DIR, f)
  if (file.exists(old)) file.remove(old)
}

cat("\nDone. Outputs in:", OUT_DIR, "\n")
cat("  rora_crossancestry_coloc_pp4.pdf — 4 traits x 4 ancestries COLOC\n")
cat("  rora_gwas_atac_motif.pdf         — Motif disruption, signed bars\n")
cat("  rora_celltype_expression.pdf     — scRNA cell-type expression\n")
cat("  rora_progression_sex.pdf         — Vertical bars: progression + sex\n")
cat("  rora_visium_spatial.pdf          — Spatial localization (by cohort; no disease claim)\n")
