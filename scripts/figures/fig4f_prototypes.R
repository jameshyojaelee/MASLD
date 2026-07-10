#!/usr/bin/env Rscript
##############################################################################
# fig4f_prototypes.R
# Three alternative designs for Fig 4 drug-target genetic validation panel.
# Gene universe: 8 clinical anchors + 11 novel druggable genes with
# SuSiE PP4 >= 0.5 from the 50-GWAS SuSiE-COLOC portfolio (MVP 5-ancestry expansion; was 23-GWAS pre-2026-07-05).
#
# Outputs (figures/main/fig4_validation/panels/):
#   fig4f_opt1_lollipop.pdf  -- lollipop by best SuSiE PP4
#   fig4f_opt2_scatter.pdf   -- dream logFC vs PP4 scatter
#   fig4f_opt3_forest.pdf    -- per-GWAS PP4 forest/heatmap
#   fig4f_compare.pdf        -- side-by-side comparison sheet
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

PANEL_DIR <- file.path(FIG4_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Gene universe + clinical status tiers
# ============================================================================
# Clinical anchors: include even if COLOC weak (shows null vs genetic support)
clinical_tiers <- data.table(
  gene = c("THRB",    "NR1H4",            "PPARA",        "PPARG",
           "GLP1R",   "SCD",              "HSD17B13",     "DGAT2"),
  drug = c("Resmetirom", "Obeticholic acid", "Elafibranor", "Lanifibranor",
           "Semaglutide", "Aramchol",       "Rapirosiran",  "ION224"),
  status = c("FDA approved", "Phase 3 (rejected)", "Phase 3 (failed)",
             "Phase 3", "FDA approved", "Phase 3",
             "Phase 3 planned", "Phase 2b")
)
clinical_tiers[, tier := fcase(
  grepl("FDA", status),         "FDA approved",
  grepl("Phase 3", status),     "Phase 3",
  grepl("Phase 2", status),     "Phase 2/earlier",
  default = "Other"
)]
clinical_tiers[, tier := factor(tier,
  levels = c("FDA approved", "Phase 3", "Phase 2/earlier", "Novel (druggable)"))]

# Novel druggable genes with SuSiE PP4 >= 0.5 (excluding GGT1 = cis artifact)
# HKDC1 = progression driver headline (F2->F3); ADH4 = pericentral-causal
novel_genes <- c("HKDC1", "ADH4", "CHI3L1", "GAS6", "PKN3",
                 "CTSD", "F2RL1", "ST14", "ALDH1B1", "PKM", "VDR")

novel_tiers <- data.table(
  gene   = novel_genes,
  drug   = NA_character_,
  status = "Novel (druggable)",
  tier   = factor("Novel (druggable)", levels = levels(clinical_tiers$tier))
)

gene_meta <- rbind(clinical_tiers, novel_tiers)

tier_colors <- c(
  "FDA approved"      = "#880E4F",  # deep magenta
  "Phase 3"           = "#E91E63",  # magenta
  "Phase 2/earlier"   = "#F4A674",  # Liang warm peach
  "Novel (druggable)" = "#1565C0"   # deep blue
)

# ============================================================================
# Pull COLOC + DE evidence from multi-evidence atlas
# ============================================================================
atlas_f <- file.path(BASE,
  "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
atlas <- fread(atlas_f,
  select = c("human_symbol", "bulk_logFC", "bulk_padj", "is_conserved",
             "coloc_susie_best_pp4", "coloc_susie_best_gwas",
             "coloc_susie_n_gwas_h4_05", "coloc_susie_n_gwas_h4_08",
             "coloc_abf_best_pp4", "coloc_abf_best_gwas",
             "best_liver_enzyme_pp4"))
setnames(atlas, "human_symbol", "gene")

gene_tbl <- merge(gene_meta, atlas, by = "gene", all.x = TRUE)

# Prefer SuSiE PP4; fall back to ABF for anchors where SuSiE failed
gene_tbl[, coloc_pp4 := coloc_susie_best_pp4]
gene_tbl[is.na(coloc_pp4) | coloc_pp4 < coloc_abf_best_pp4,
         coloc_pp4 := coloc_abf_best_pp4]
gene_tbl[, coloc_gwas := coloc_susie_best_gwas]
gene_tbl[is.na(coloc_gwas) | coloc_pp4 == coloc_abf_best_pp4,
         coloc_gwas := coloc_abf_best_gwas]

# Human-friendly GWAS labels
shorten_gwas <- function(x) {
  x <- gsub("^UKBB_", "UKBB ", x)
  x <- gsub("^BBJ_",  "BBJ ",  x)
  x <- gsub("^FinnGen_", "FinnGen ", x)
  x <- gsub("^PanUKBB_CSA_", "PanUKBB-CSA ", x)
  x <- gsub("^2023_36280732_NAFLD_",  "NAFLD:",   x)
  x <- gsub("^2019_31311600_NAFLD_",  "NAFLD:",   x)
  x <- gsub("^2020_32514122_",         "",        x)
  x <- gsub("_EUR$|_EAS$|_AFR$|_SAS$", "", x)
  x
}
gene_tbl[, coloc_gwas_lbl := shorten_gwas(coloc_gwas)]

cat("\n--- Gene panel (n=", nrow(gene_tbl), ") ---\n", sep = "")
print(gene_tbl[, .(gene, tier, drug, coloc_pp4, coloc_gwas_lbl,
                   bulk_logFC, is_conserved)])

# ============================================================================
# OPTION 1: Lollipop by best SuSiE PP4
# ============================================================================
cat("\n[Option 1] Lollipop\n")
d1 <- copy(gene_tbl)
setorder(d1, -coloc_pp4, na.last = TRUE)
d1[is.na(coloc_pp4), coloc_pp4 := 0]
d1[, gene := factor(gene, levels = rev(d1$gene))]
d1[, abs_lfc := abs(bulk_logFC)]
d1[is.na(abs_lfc), abs_lfc := 0]

p1 <- ggplot(d1, aes(x = coloc_pp4, y = gene, color = tier)) +
  geom_segment(aes(x = 0, xend = coloc_pp4, yend = gene),
               linewidth = 0.4) +
  geom_point(aes(size = abs_lfc)) +
  geom_vline(xintercept = c(0.5, 0.9), linetype = "dashed",
             linewidth = 0.25, color = "gray55") +
  geom_text(aes(label = coloc_gwas_lbl),
            hjust = -0.15, size = GEOM_TEXT_6PT, color = "black") +
  scale_color_manual(values = tier_colors, name = NULL, drop = FALSE) +
  scale_size_continuous(name = expression("|log"[2]*"FC|"),
                        range = c(0.8, 3.2), limits = c(0, 1.6)) +
  scale_x_continuous(limits = c(0, 1.35), breaks = c(0, 0.5, 0.9, 1.0),
                     labels = c("0", "0.5", "0.9", "1")) +
  annotate("text", x = 0.5,  y = 0.5, label = "PP4 = 0.5",
           hjust = -0.05, vjust = 0, size = GEOM_TEXT_6PT, color = "black") +
  annotate("text", x = 0.9,  y = 0.5, label = "PP4 = 0.9",
           hjust = -0.05, vjust = 0, size = GEOM_TEXT_6PT, color = "black") +
  labs(x = "Best SuSiE colocalization PP4 (across 50 GWAS)", y = NULL) +
  theme_masld() +
  theme(axis.text.y = element_text(size = 6, face = "italic"),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.key.size = unit(0.25, "cm")) +
  guides(color = guide_legend(override.aes = list(size = 2), nrow = 2),
         size  = guide_legend(nrow = 1))

message("[caption] Drug-target genetic validation, option 1: lollipop by best SuSiE PP4. Point size = transcriptional effect magnitude; label = best GWAS trait.")

save_fig(p1, file.path(PANEL_DIR, "fig4f_opt1_lollipop.pdf"),
         width = fig_half_width * 1.1, height = 4.2)

# ============================================================================
# OPTION 2: Scatter of bulk_logFC vs best SuSiE PP4
# ============================================================================
cat("[Option 2] Scatter\n")
d2 <- copy(gene_tbl)
d2[is.na(bulk_logFC), bulk_logFC := 0]
d2[is.na(coloc_pp4), coloc_pp4 := 0]

# Label all anchors + all novels with PP4 > 0.5
d2[, do_label := tier != "Novel (druggable)" | coloc_pp4 >= 0.5]

p2 <- ggplot(d2, aes(x = bulk_logFC, y = coloc_pp4, color = tier)) +
  geom_hline(yintercept = c(0.5, 0.9), linetype = "dashed",
             linewidth = 0.25, color = "gray55") +
  geom_vline(xintercept = 0, linewidth = 0.2, color = "gray70") +
  geom_point(size = 2.2, alpha = 0.9) +
  geom_label_repel(data = d2[do_label == TRUE],
                   aes(label = gene),
                   size = GEOM_TEXT_6PT, max.overlaps = 30,
                   label.padding = 0.1, segment.size = 0.2,
                   min.segment.length = 0, fontface = "italic",
                   fill = alpha("white", 0.85), show.legend = FALSE) +
  scale_color_manual(values = tier_colors, name = NULL, drop = FALSE) +
  scale_y_continuous(limits = c(-0.02, 1.05),
                     breaks = c(0, 0.5, 0.9, 1.0),
                     labels = c("0", "0.5", "0.9", "1")) +
  annotate("text", x = Inf, y = 0.9, label = "PP4 = 0.9 (strong)",
           hjust = 1.05, vjust = -0.3, size = GEOM_TEXT_6PT, color = "black") +
  annotate("text", x = Inf, y = 0.5, label = "PP4 = 0.5 (canonical)",
           hjust = 1.05, vjust = -0.3, size = GEOM_TEXT_6PT, color = "black") +
  labs(x = expression("Transcript log"[2]*"FC (dream; MASLD vs control)"),
       y = "Best SuSiE colocalization PP4") +
  theme_masld() +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.25, "cm")) +
  guides(color = guide_legend(override.aes = list(size = 2.5), nrow = 2))

message("[caption] Drug-target genetic validation, option 2: scatter of transcript logFC vs best SuSiE PP4. Top-right / top-left = genetically validated + transcriptionally regulated.")

save_fig(p2, file.path(PANEL_DIR, "fig4f_opt2_scatter.pdf"),
         width = fig_half_width * 1.1, height = 4.2)

# ============================================================================
# OPTION 3: Per-GWAS PP4 heatmap (forest-plot feel)
# ============================================================================
cat("[Option 3] Per-GWAS heatmap\n")

susie_f <- file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv")
susie <- fread(susie_f, select = c("gene", "gwas_name", "PP.H4.susie", "PP.H4.abf"))

# Collapse per (gene, gwas): take max of SuSiE (fall back to ABF)
susie[, pp4 := pmax(PP.H4.susie, PP.H4.abf, na.rm = TRUE)]
susie_long <- susie[gene %in% gene_tbl$gene,
                    .(pp4 = max(pp4, na.rm = TRUE)),
                    by = .(gene, gwas_name)]
susie_long[is.infinite(pp4), pp4 := NA_real_]

# Select an informative GWAS subset: 12 traits covering enzymes, NAFLD, HCC,
# cross-ancestry. Preserve ordering so enzymes come first.
gwas_keep <- c(
  # UKBB liver enzymes
  "UKBB_ALT", "UKBB_AST", "UKBB_GGT",
  # BBJ liver enzymes (East Asian)
  "BBJ_ALT",  "BBJ_AST",  "BBJ_GGT",
  # Steatosis / PDFF
  "UKBB_PDFF",
  # NAFLD meta-analyses
  "2023_36280732_NAFLD_UKBB_EUR",
  "2023_36280732_NAFLD_Intermountain_EUR",
  "2019_31311600_NAFLD_EUR",
  # FinnGen
  "FinnGen_NAFLD", "FinnGen_NASH",
  # Cross-ancestry
  "PanUKBB_CSA_ALT", "PanUKBB_CSA_GGT", "PanUKBB_AFR_ALT"
)

d3 <- susie_long[gwas_name %in% gwas_keep]
# Fill missing (gene, gwas) pairs with PP4 = 0 for full-matrix visual
all_pairs <- CJ(gene = gene_tbl$gene, gwas_name = gwas_keep)
d3 <- merge(all_pairs, d3, by = c("gene", "gwas_name"), all.x = TRUE)
d3[is.na(pp4), pp4 := 0]

# Order genes by best PP4 across these GWAS
gene_order <- d3[, .(best = max(pp4, na.rm = TRUE)), by = gene][order(best)]
d3[, gene := factor(gene, levels = gene_order$gene)]

# Friendly GWAS labels + ordering
gwas_lbl <- c(
  "UKBB_ALT"   = "UKBB ALT",   "UKBB_AST" = "UKBB AST", "UKBB_GGT" = "UKBB GGT",
  "BBJ_ALT"    = "BBJ ALT",    "BBJ_AST"  = "BBJ AST",  "BBJ_GGT"  = "BBJ GGT",
  "UKBB_PDFF"  = "UKBB PDFF",
  "2023_36280732_NAFLD_UKBB_EUR"          = "NAFLD:UKBB",
  "2023_36280732_NAFLD_Intermountain_EUR" = "NAFLD:Intermountain",
  "2019_31311600_NAFLD_EUR"               = "NAFLD:Anstee",
  "FinnGen_NAFLD" = "FinnGen NAFLD",
  "FinnGen_NASH"  = "FinnGen NASH",
  "PanUKBB_CSA_ALT" = "PanUKBB-CSA ALT",
  "PanUKBB_CSA_GGT" = "PanUKBB-CSA GGT",
  "PanUKBB_AFR_ALT" = "PanUKBB-AFR ALT"
)
d3[, gwas_lbl := factor(gwas_lbl[gwas_name], levels = gwas_lbl)]

# Tier annotation on y-axis
tier_tbl <- gene_meta[, .(gene, tier)]
d3 <- merge(d3, tier_tbl, by = "gene", all.x = TRUE)

p3 <- ggplot(d3, aes(x = gwas_lbl, y = gene, fill = pp4)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(data = d3[pp4 >= 0.5],
            aes(label = sprintf("%.2f", pp4)),
            size = GEOM_TEXT_6PT, color = "white", fontface = "plain") +
  scale_fill_gradientn(
    colours = c("#F5F5F5", "#CFD8DC", "#42A5F5", "#C9265E", "#A01753"),
    values  = scales::rescale(c(0, 0.25, 0.5, 0.8, 1.0)),
    limits  = c(0, 1),
    breaks  = c(0, 0.5, 0.9, 1.0),
    name    = "PP4"
  ) +
  labs(x = NULL, y = NULL) +
  theme_masld() +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6, face = "italic"),
        legend.key.size = unit(0.3, "cm"))

message("[caption] Drug-target genetic validation, option 3: per-GWAS SuSiE PP4 heatmap; tile annotated when PP4 >= 0.5.")

save_fig(p3, file.path(PANEL_DIR, "fig4f_opt3_forest.pdf"),
         width = fig_full_width * 0.7, height = 4.5)

# ============================================================================
# Side-by-side comparison sheet
# ============================================================================
cat("[compare] Assembling comparison sheet\n")

message("[caption] Fig 4f alternatives: drug-target genetic validation. Same gene universe; compare visual clarity and information density.")

compare <- (p1 | p2) / p3 +
  plot_layout(heights = c(1, 1.2)) &
  theme(plot.tag = element_text(size = 6, face = "plain"))

save_fig(compare, file.path(PANEL_DIR, "fig4f_compare.pdf"),
         width = fig_full_width, height = 9)

cat("\nDone.\n")
cat("  Option 1 (lollipop):  ", file.path(PANEL_DIR, "fig4f_opt1_lollipop.pdf"), "\n")
cat("  Option 2 (scatter):   ", file.path(PANEL_DIR, "fig4f_opt2_scatter.pdf"), "\n")
cat("  Option 3 (heatmap):   ", file.path(PANEL_DIR, "fig4f_opt3_forest.pdf"), "\n")
cat("  Comparison sheet:     ", file.path(PANEL_DIR, "fig4f_compare.pdf"), "\n")
