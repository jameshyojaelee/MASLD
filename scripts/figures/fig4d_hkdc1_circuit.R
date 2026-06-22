#!/usr/bin/env Rscript
# KEY MESSAGE: HKDC1 is connected from a fine-mapped GWAS variant through a
# hepatocyte NRF2 co-expression module to protein-level validation, demonstrating
# a complete genetic-to-functional circuit at the F2→F3 fibrosis transition.
#
# Panels:
#   i.   SuSiE-COLOC PP.H4 bar for HKDC1 vs context genes
#   ii.  hep-24 module lollipop (top genes by weight)
#   iii. Stage-resolved bulk expression heatmap-style (available as logFC vs control)
#   iv.  Protein logBF comparison: HKDC1 vs atlas mean
#
# Note: Full GWAS locus zoom available via:
#   Rscript scripts/figures/figS09_locus_zoom.R <locus_id_for_HKDC1>
#
# Module source = CIRRHOSIS-EXCLUDED hepatocyte Hotspot run
# (hepatocytes_nocirrhosis, module 17 -> canonical module 24; n = 192 donors).
# Trajectory betas for the caption: F_stage +0.47 (q 3.5e-10),
# disease_stage +0.75 (q 2.7e-12).
#
# Output: figures/main/fig4_validation/hkdc1_nrf2_circuit.pdf
# Env:    rnaseq

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── Data ──────────────────────────────────────────────────────────────────────
atlas <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/integration/multi_evidence_atlas_with_spatial.csv"),
  stringsAsFactors = FALSE)

# hep-24 module genes — sourced from the cirrhosis-excluded run
# (nocirrhosis module 17 == canonical module 24; n = 192 donors).
mod_genes <- read.table(
  file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules/hepatocytes_nocirrhosis/module_genes.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)
hep24 <- mod_genes %>%
  filter(module == 17) %>%
  arrange(desc(weight)) %>%
  left_join(atlas %>% select(human_symbol, bulk_logFC, bulk_padj),
            by = c("gene" = "human_symbol")) %>%
  mutate(
    is_nrf2_canonical = gene %in% c("HKDC1", "SOD2", "TXNRD1", "AKR1B10",
                                     "SQSTM1", "AKR1C1", "AKR1C2", "ME1"),
    direction = case_when(
      !is.na(bulk_logFC) & bulk_logFC > 0 ~ "up",
      !is.na(bulk_logFC) & bulk_logFC < 0 ~ "down",
      TRUE                                 ~ "ns"
    )
  )

# Module metadata from all_modules.tsv
mod_meta <- read.table(
  file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules/all_modules.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
) %>% filter(cell_type == "hepatocytes", module == 24)

# HKDC1 from atlas
hkdc1_row <- atlas %>% filter(human_symbol == "HKDC1")

# ── hep-24 module score vs fibrosis trajectory (cirrhosis-excluded run) ───────
# Per-donor module score (nocirr module 17 == canonical 24) joined to the
# augmented per-donor F-stage. Restrict to the stage-analysis donors
# (exclude_stage_analysis == FALSE) so the plotted set matches the n = 192
# the trajectory betas were fit on (F_stage beta +0.47, q 3.5e-10).
donor_scores <- read.table(
  file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules/hepatocytes_nocirrhosis/donor_scores.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
) %>% filter(module == 17) %>% select(sample, mod17_score = score)

donor_meta <- read.table(
  file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, quote = ""
) %>% select(sample, F_stage_augmented, exclude_stage_analysis)

traj <- donor_scores %>%
  inner_join(donor_meta, by = "sample") %>%
  filter(exclude_stage_analysis == "FALSE", !is.na(F_stage_augmented)) %>%
  mutate(F_stage = factor(paste0("F", F_stage_augmented),
                          levels = paste0("F", 0:4)))
n_traj <- nrow(traj)

# ── Panel i: COLOC PP.H4 context (HKDC1 vs NR triad) ────────────────────────
coloc_context <- atlas %>%
  filter(human_symbol %in% c("HKDC1", "THRB", "RORA", "NR1H4", "SERPINE1", "FADS2")) %>%
  select(human_symbol, coloc_susie_best_pp4, bulk_logFC) %>%
  mutate(
    is_focal = human_symbol == "HKDC1",
    gene     = factor(human_symbol,
                      levels = c("THRB", "RORA", "HKDC1", "FADS2", "SERPINE1", "NR1H4"))
  )

p_coloc <- ggplot(coloc_context,
                  aes(x = coloc_susie_best_pp4, y = gene, color = is_focal)) +
  geom_vline(xintercept = 0.5, linewidth = 0.3, linetype = "dashed", color = "gray60") +
  geom_segment(aes(x = 0, xend = coloc_susie_best_pp4, y = gene, yend = gene),
               linewidth = 0.5) +
  geom_point(aes(size = is_focal)) +
  scale_color_manual(values = c(`FALSE` = "gray60", `TRUE` = masld_colors[["up"]]),
                     guide = "none") +
  scale_size_manual(values = c(`FALSE` = 1.5, `TRUE` = 2.8), guide = "none") +
  scale_x_continuous(limits = c(0, 1.05), breaks = c(0, 0.5, 1.0)) +
  coord_cartesian(clip = "off") +
  labs(x = "COLOC PP.H4", y = NULL,
       title = "Genetic evidence") +   # SuSiE where available, else coloc.abf (SERPINE1/NR1H4)
  theme_masld() + theme_pub() +
  theme(axis.text.y = element_text(
    face = ifelse(levels(coloc_context$gene) == "HKDC1", "bold.italic", "italic"),
    size = PUB_AXIS_TEXT
  ))

# ── Panel ii: hep-24 module lollipop (top 15 genes) ─────────────────────────
top_hep24 <- hep24 %>% slice_max(weight, n = 15) %>%
  mutate(gene = factor(gene, levels = gene[order(weight)]))

# Canonical NRF2 members are highlighted by larger magenta dots; the y-axis
# already names every gene, so no duplicate right-side labels (keeps it clean).
p_module <- ggplot(top_hep24, aes(x = weight, y = gene, color = direction)) +
  geom_segment(aes(x = 0, xend = weight, y = gene, yend = gene), linewidth = 0.5) +
  geom_point(aes(size = is_nrf2_canonical)) +
  scale_color_manual(
    values = c(up = masld_colors[["up"]], down = masld_colors[["down"]], ns = "gray60"),
    labels = c(up = "Up in MASLD", down = "Down", ns = "n.s."),
    name   = "Bulk direction"
  ) +
  scale_size_manual(values = c(`TRUE` = 2.8, `FALSE` = 1.2), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(clip = "off") +
  labs(x = "Module weight", y = NULL,
       title = "hep-24 module") +
  theme_masld() + theme_pub() +
  theme(axis.text.y   = element_text(size = PUB_AXIS_TEXT, face = "italic"),
        legend.position = "bottom")

# ── Panel iii: hep-24 module score vs fibrosis trajectory ───────────────────
# Per-donor points (jittered) + per-stage boxplot; the rising trend across
# F0..F4 is the scientific point. The trajectory betas live in the figure
# legend (emitted to stdout below), NOT on the panel (PI directive).
p_traj <- ggplot(traj, aes(x = F_stage, y = mod17_score)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.3,
               fill = NA, color = "gray55") +
  geom_jitter(width = 0.16, height = 0, size = 0.7, alpha = 0.55,
              color = masld_colors[["up"]]) +
  stat_summary(fun = mean, geom = "line", aes(group = 1),
               color = masld_colors[["up"]], linewidth = 0.7) +
  stat_summary(fun = mean, geom = "point",
               color = masld_colors[["up"]], size = 1.6) +
  coord_cartesian(clip = "off") +
  labs(x = "Fibrosis stage", y = "hep-24 module score",
       title = "Fibrosis trajectory") +
  theme_masld() + theme_pub() +
  theme(axis.title.y = element_text(size = PUB_AXIS_TITLE))

# ── Assemble ──────────────────────────────────────────────────────────────────
# Stats (mRNA logFC, COLOC PP.H4, protein logFC) belong in the legend, not the
# panel subtitle (PI directive). Emit them to stdout for the caption.
message(sprintf("[hkdc1 legend] HKDC1: mRNA logFC = %.2f; COLOC PP.H4 = %.3f; protein logFC = %.2f",
                hkdc1_row$bulk_logFC, hkdc1_row$coloc_susie_best_pp4,
                hkdc1_row$best_protein_logFC))
message(sprintf("[hkdc1 legend] hep-24 module trajectory (cirrhosis-excluded, n = %d): ", n_traj),
        "F_stage beta +0.47 (q 3.5e-10); disease_stage beta +0.75 (q 2.7e-12)")

p_out <- (p_coloc | p_module | p_traj) +
  plot_layout(widths = c(1, 1.5, 1.2))

out <- file.path(FIG4_DIR, "hkdc1_nrf2_circuit.pdf")
cairo_pdf(out, width = fig_full_width * 1.18, height = 2.8)
print(p_out)
dev.off()
message("Saved: ", out)
