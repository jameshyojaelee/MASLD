#!/usr/bin/env Rscript
# KEY MESSAGE: HKDC1 anchors a previously GLYCEMIC locus to MASLD via liver-enzyme
# colocalization (SuSiE 0.99, UKBB ALT/AST) and sits in a hepatocyte antioxidant
# co-expression module. The new contribution is MASLD-anchoring a glycemic locus,
# NOT a fibrosis-trajectory claim (the prior trajectory sub-panel was a scVI
# F-stage leakage artifact: on honest documented Andrews staging, n = 39, F0-F3,
# the module score is FLAT and non-monotonic — DROPPED, see ledger HKDC1_trajectory=DROP).
#
# Panels (revised, two-panel):
#   i.  COLOC PP.H4 lollipop: HKDC1 vs context genes. Axis is "COLOC PP.H4"
#       (NOT "SuSiE-") because some plotted genes (THRB, SERPINE1, NR1H4) are
#       abf-only — they have no convergent SuSiE posterior.
#   ii. hep antioxidant co-expression module weight lollipop (top genes).
#
# Backbone (confirmatory, in legend): mouse Hkdc1 KO rescues MASH (Xu/Khan 2025);
# the ALT/AST GWAS locus is female-biased (Pazoki 2021).
# Caveat (in legend): the co-expression module is HYPOTHESIS-GENERATING, NOT
# NRF2-target proof — HKDC1's published mechanism is metabolic/mitochondrial.
#
# COLOC numbers are read from the canonical per-GWAS coloc file
# (susie_coloc_all_gwas.csv), NOT the method-inconsistent atlas *_coloc_pp4
# convenience columns (HARD GATE).
#
# Output: figures/main/fig5_molecular_context/fig4g_hkdc1_circuit.pdf
# (relettered e->g 2026-07-08: now sits after both cyp3a4 panels, per the Fig4 lineup)
# Env:    rnaseq

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── Data ────────────────────────────────────────────────────────────────────
# (1) Canonical per-GWAS COLOC (NEVER use atlas *_coloc_pp4 convenience cols).
coloc <- read.csv(
  file.path(BASE, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"),
  stringsAsFactors = FALSE)

# (2) Atlas — for bulk_logFC only (direction colouring of module genes).
atlas <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/integration/multi_evidence_atlas_with_spatial.csv"),
  stringsAsFactors = FALSE)

# (3) Hepatocyte antioxidant co-expression module — cirrhosis-excluded run
#     (nocirrhosis module 17 == canonical module 24; HKDC1 is a member).
mod_genes <- read.table(
  file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules/hepatocytes_nocirrhosis/module_genes.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ── Per-gene best COLOC PP.H4 (SuSiE where convergent, else coloc.abf) ────────
context_genes <- c("HKDC1", "THRB", "RORA", "FADS2", "SERPINE1", "NR1H4")

best_coloc <- coloc %>%
  filter(gene %in% context_genes) %>%
  group_by(gene) %>%
  summarise(
    susie_pp4 = suppressWarnings(max(PP.H4.susie, na.rm = TRUE)),
    susie_gwas = ifelse(all(is.na(PP.H4.susie)), NA_character_,
                        gwas_name[which.max(replace(PP.H4.susie, is.na(PP.H4.susie), -Inf))]),
    abf_pp4   = suppressWarnings(max(PP.H4.abf, na.rm = TRUE)),
    abf_gwas  = gwas_name[which.max(PP.H4.abf)],
    .groups = "drop") %>%
  mutate(
    susie_pp4 = ifelse(is.infinite(susie_pp4), NA_real_, susie_pp4),
    has_susie = !is.na(susie_pp4),
    # Plot SuSiE posterior where it converged, otherwise coloc.abf.
    pp4       = ifelse(has_susie, susie_pp4, abf_pp4),
    method    = ifelse(has_susie, "SuSiE", "coloc.abf"),
    best_gwas = ifelse(has_susie, susie_gwas, abf_gwas))

# Order: focal HKDC1 in the middle of the strong band; weak (abf-only) at bottom.
gene_order <- c("NR1H4", "SERPINE1", "FADS2", "RORA", "HKDC1", "THRB")
best_coloc <- best_coloc %>%
  mutate(
    is_focal = gene == "HKDC1",
    gene_f   = factor(gene, levels = gene_order))

# ── hep antioxidant module lollipop (top genes by weight) ─────────────────────
nrf2_canonical <- c("HKDC1", "SOD2", "TXNRD1", "AKR1B10", "SQSTM1",
                    "AKR1C1", "AKR1C2", "ME1", "BACH1", "PTGR1")

hepmod <- mod_genes %>%
  filter(module == 17) %>%
  left_join(atlas %>% select(human_symbol, bulk_logFC),
            by = c("gene" = "human_symbol")) %>%
  mutate(
    is_nrf2 = gene %in% nrf2_canonical,
    direction = case_when(
      !is.na(bulk_logFC) & bulk_logFC > 0 ~ "up",
      !is.na(bulk_logFC) & bulk_logFC < 0 ~ "down",
      TRUE                                 ~ "ns"))

# Protein-coding focus for the panel: drop the lncRNA/antisense module members
# (LINC*, *HG, MALAT1, PVT1, LUCAT1, ENSG ids) so the antioxidant program reads
# clearly; HKDC1 + the canonical antioxidant genes are protein-coding.
hepmod_pc <- hepmod %>%
  filter(!grepl("^(LINC|ENSG|GS1-|MIR)", gene),
         !gene %in% c("MALAT1", "PVT1", "LUCAT1"))

top_mod <- hepmod_pc %>%
  slice_max(weight, n = 14) %>%
  mutate(gene = factor(gene, levels = gene[order(weight)]))

# ── Panel i: COLOC PP.H4 bar chart ───────────────────────────────────────────
# abf-only genes (no convergent SuSiE posterior) shown at reduced alpha
p_coloc <- ggplot(best_coloc, aes(x = pp4, y = gene_f)) +
  geom_vline(xintercept = 0.5, linewidth = 0.3, linetype = "dashed",
             color = "gray70") +
  geom_col(aes(fill = is_focal, alpha = has_susie), width = 0.65) +
  scale_fill_manual(values = c(`FALSE` = "#9E9E9E", `TRUE` = masld_colors[["up"]]),
                    guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.45),
                     labels = c(`TRUE` = "SuSiE", `FALSE` = "coloc.abf"),
                     name = "Method",
                     guide = guide_legend(override.aes = list(fill = "#9E9E9E"))) +
  scale_x_continuous(limits = c(0, 1.08), breaks = c(0, 0.5, 1.0),
                     expand = expansion(mult = c(0, 0))) +
  coord_cartesian(clip = "off") +
  labs(x = "COLOC PP.H4", y = NULL) +
  theme_masld_compact() +
  theme(
    axis.text.y = element_text(size = 6, face = "italic", color = "black"),
    legend.position = "bottom",
    legend.key.size = unit(0.28, "lines"))

# ── Panel ii: hepatocyte antioxidant co-expression module ─────────────────────
p_module <- ggplot(top_mod, aes(x = weight, y = gene)) +
  geom_col(aes(fill = direction), width = 0.65) +
  scale_fill_manual(
    values = c(up = masld_colors[["up"]], down = masld_colors[["down"]],
               ns = "#9E9E9E"),
    breaks = c("up", "down", "ns"),
    labels = c("Up in MASLD", "Down", "n.s."),
    name   = "Bulk direction") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.06))) +
  coord_cartesian(clip = "off") +
  labs(x = "Co-expression weight", y = NULL) +
  theme_masld_compact() +
  theme(
    axis.text.y = element_text(size = 6, face = "italic", color = "black"),
    legend.position = "bottom",
    legend.key.size = unit(0.28, "lines"))

# ── Legend stats → stdout (NOT on the panel; PI directive) ────────────────────
hk <- best_coloc %>% filter(gene == "HKDC1")
hk_bulk <- atlas$bulk_logFC[match("HKDC1", atlas$human_symbol)]
hk_prot <- atlas$best_protein_logFC[match("HKDC1", atlas$human_symbol)]
hk_prot_padj <- atlas$best_protein_padj[match("HKDC1", atlas$human_symbol)]

message("==================== HKDC1 circuit — figure legend (stats) ====================")
message(sprintf("Panel i (COLOC PP.H4): HKDC1 SuSiE PP.H4 = 0.992 (UKBB_AST) / coloc.abf = 0.993 (UKBB_ALT); plotted = %.3f (%s, %s).",
                hk$pp4, hk$method, hk$best_gwas))
message("  Axis labelled 'COLOC PP.H4' (NOT SuSiE-) because THRB / SERPINE1 / NR1H4 are abf-only")
message("  (no convergent SuSiE posterior; shown as open rings). Context: THRB abf 1.000 (UKBB_GGT),")
message("  RORA SuSiE 0.998 (BBJ_GGT), FADS2 SuSiE 0.948 (BBJ_ALT), SERPINE1 abf 0.138, NR1H4 abf 0.208.")
message(sprintf("Panel ii: HKDC1 in the hepatocyte antioxidant co-expression module (cirrhosis-excluded Hotspot run; nocirr module 17 == canonical module 24); %d module members, %d protein-coding shown.",
                nrow(hepmod), nrow(hepmod_pc)))
message(sprintf("HKDC1 bulk mRNA logFC = %.2f (limma-voom-qw C2). Protein logFC = %.2f (DIA-MS, padj = %.2f, n.s. — not significant).",
                hk_bulk, hk_prot, hk_prot_padj))
message("CAVEAT (bake into legend): the co-expression module is HYPOTHESIS-GENERATING (antioxidant CO-EXPRESSION),")
message("  NOT NRF2-target proof — HKDC1's published mechanism is metabolic/mitochondrial. The novel layer here is")
message("  MASLD-anchoring of a previously GLYCEMIC locus via liver-enzyme colocalization.")
message("BACKBONE (confirmatory): mouse Hkdc1 KO rescues MASH (Xu/Khan 2025); ALT/AST GWAS locus female-biased (Pazoki 2021).")
message("DROPPED: fibrosis-trajectory sub-panel — scVI F-stage leakage artifact (ledger HKDC1_trajectory=DROP);")
message("  on honest documented Andrews SAF staging (n = 39, F0-F3) the module score is FLAT/non-monotonic, n.s.")
message("===============================================================================")

# ── Assemble (two panels) ─────────────────────────────────────────────────────
p_out <- (p_coloc | p_module) +
  plot_layout(widths = c(1, 1.25))

out <- file.path(FIG4_DIR, "panels", "fig4g_hkdc1_circuit.pdf")
if (!dir.exists(FIG4_DIR)) dir.create(FIG4_DIR, recursive = TRUE)
cairo_pdf(out, width = fig_full_width * 0.66, height = 2.05, family = "Helvetica")
print(p_out)
invisible(dev.off())
message("Saved: ", out)
