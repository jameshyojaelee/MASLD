#!/usr/bin/env Rscript
# ============================================================================
# figS_hotspot_landscape.R
#
# KEY MESSAGE: The Hotspot autocorrelation module catalog (202 modules across
# 7 cell types) is structurally honest: large-effect progression modules
# concentrate at high LOO stability, and the OR-criterion saturation against
# cNMF k=16 / bulk NMF k=6 / Hallmark / SCENIC+ / curated panels is shown
# transparently. Four panels (A-D).
#
# A — Catalog: 202 modules per cell type, stacked by class
# B — Phenotype: 50 progression modules ordered by disease β
# C — Stability vs significance: per-CT scatter; 0.5 stability floor
# D — Method honesty: where each module's top match landed
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- FIGS_HOTSPOT_PANELS_DIR
DATA_DIR  <- FIGS_HOTSPOT_DATA_DIR
PREFIX    <- "landscape_"

HS_RES <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")
mods <- fread(file.path(HS_RES, "all_modules.tsv"))

# CT display names
CT_DISPLAY <- c(global="Global", hepatocytes="Hepatocytes",
                macrophages="Macrophages", fibroblasts="Fibroblasts",
                endothelial_cells="Endothelial cells",
                cholangiocytes="Cholangiocytes", tcells="T cells")
mods[, ct := CT_DISPLAY[cell_type]]

# Order cell types by total module count (decreasing)
CT_ORDER <- mods[, .N, by = ct][order(-N), ct]
mods[, ct := factor(ct, levels = CT_ORDER)]

# Map ct to ct_palette where present, otherwise gray
CT_FILL <- sapply(CT_ORDER, function(x) {
  if (x %in% names(ct_palette)) ct_palette[[x]] else "#7E57C2"
})
names(CT_FILL) <- CT_ORDER

# ----------------------------------------------------------------------------
# Panel A — Catalog by class
# ----------------------------------------------------------------------------
mods[, class := dplyr::case_when(
  progression_module & bulk_replicated & stability_score >= 0.5 ~ "Gold (prog+bulk+stab)",
  progression_module                                            ~ "Progression",
  activity_module                                               ~ "Activity (NAS-only)",
  bulk_replicated                                               ~ "Bulk-replicated",
  TRUE                                                          ~ "Other"
)]
class_levels <- c("Gold (prog+bulk+stab)", "Progression",
                  "Activity (NAS-only)", "Bulk-replicated", "Other")
mods[, class := factor(class, levels = class_levels)]
CLASS_FILL <- c("Gold (prog+bulk+stab)" = "#A01753",
                "Progression"           = "#C9265E",
                "Activity (NAS-only)"   = "#F4A674",
                "Bulk-replicated"       = "#90A4AE",
                "Other"                 = "#E0E0E0")

count_by_class <- mods[, .N, by = .(ct, class)]
pA <- ggplot(count_by_class, aes(ct, N, fill = class)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = CLASS_FILL, name = "Module class") +
  labs(x = NULL, y = "Modules") +
  theme_masld() + theme_pub() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "right")

# ----------------------------------------------------------------------------
# Panel B — 50 progression modules ordered by disease β
# ----------------------------------------------------------------------------
prog <- mods[progression_module == TRUE]
prog[, mod_label := sprintf("%s__%d", cell_type, module)]
prog <- prog[order(disease_stage_beta)]
prog[, mod_label := factor(mod_label, levels = mod_label)]
prog[, signed_logq := sign(disease_stage_beta) *
                      pmin(-log10(pmax(disease_stage_q, 1e-20)), 15)]

pB <- ggplot(prog, aes(mod_label, ct, fill = signed_logq)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C9265E",
                       midpoint = 0,
                       name = expression(sign(beta) %*% -log[10](q[disease]))) +
  labs(x = NULL, y = NULL) +
  theme_masld() + theme_pub() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.key.width = unit(0.4, "cm"))

# ----------------------------------------------------------------------------
# Panel C — Stability vs significance scatter
# ----------------------------------------------------------------------------
sc <- mods[, .(cell_type = ct, stability_score, disease_stage_q,
               bulk_replicated, progression_module)]
sc[, neg_log_q := -log10(pmax(disease_stage_q, 1e-20))]
sc[, role := dplyr::case_when(
  progression_module & bulk_replicated ~ "Progression + bulk",
  progression_module                    ~ "Progression only",
  bulk_replicated                       ~ "Bulk-replicated only",
  TRUE                                  ~ "Neither"
)]
ROLE_SHAPE <- c("Progression + bulk"  = 21,
                "Progression only"    = 22,
                "Bulk-replicated only"= 24,
                "Neither"             = 25)

pC <- ggplot(sc, aes(stability_score, neg_log_q,
                     fill = cell_type, shape = role)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted",
             linewidth = 0.3, color = "gray40") +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             linewidth = 0.3, color = "gray40") +
  geom_point(size = 1.4, alpha = 0.85, stroke = 0.2, color = "black") +
  scale_fill_manual(values = CT_FILL, name = "Cell type") +
  scale_shape_manual(values = ROLE_SHAPE, name = "Role") +
  labs(x = "LOO stability (49/49)",
       y = expression(-log[10](q[disease]))) +
  theme_masld() + theme_pub() +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "gray60")))

# ----------------------------------------------------------------------------
# Panel D — Best-match panel breakdown
# ----------------------------------------------------------------------------
PANEL_FILL <- c("cnmf_k16"      = "#1565C0",
                "hallmark"      = "#0288D1",
                "scenic_hep"    = "#7B1FA2",
                "curated_liver" = "#C9265E",
                "bulk_nmf_k6"   = "#F57F17")
bm <- mods[, .N, by = .(ct, best_match_panel)]
pD <- ggplot(bm, aes(ct, N, fill = best_match_panel)) +
  geom_col(position = "stack", width = 0.75,
           color = "white", linewidth = 0.2) +
  scale_fill_manual(values = PANEL_FILL, name = "Best-match panel",
                    breaks = names(PANEL_FILL)) +
  labs(x = NULL, y = "Modules") +
  theme_masld() + theme_pub() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "right")

# ----------------------------------------------------------------------------
# Save each panel + composite (2x2)
# ----------------------------------------------------------------------------
save_fig(pA, file.path(PANEL_DIR, paste0(PREFIX, "A_catalog_by_class.pdf")),
         width = fig_half_width, height = 2.2)
save_fig(pB, file.path(PANEL_DIR, paste0(PREFIX, "B_progression_phenotype.pdf")),
         width = fig_full_width, height = 1.8)
save_fig(pC, file.path(PANEL_DIR, paste0(PREFIX, "C_stability_vs_significance.pdf")),
         width = fig_half_width, height = 2.2)
save_fig(pD, file.path(PANEL_DIR, paste0(PREFIX, "D_best_match_panel.pdf")),
         width = fig_half_width, height = 2.2)

composite <- (pA + pD) / (pB) / (pC) +
  patchwork::plot_layout(heights = c(1, 1, 1.3)) +
  patchwork::plot_annotation(tag_levels = "a")
save_fig(composite,
         file.path(FIGS_HOTSPOT_DIR, "figS_hotspot_landscape.pdf"),
         width = fig_full_width, height = 6.5)

# Data tables for transparency
fwrite(count_by_class, file.path(DATA_DIR, paste0(PREFIX, "A_catalog_data.csv")))
fwrite(prog[, .(mod_label, cell_type = ct, module, disease_stage_beta,
                disease_stage_q, signed_logq)],
       file.path(DATA_DIR, paste0(PREFIX, "B_progression_data.csv")))
fwrite(sc, file.path(DATA_DIR, paste0(PREFIX, "C_stability_data.csv")))
fwrite(bm, file.path(DATA_DIR, paste0(PREFIX, "D_bestmatch_data.csv")))

cat("Wrote figS_hotspot_landscape: 4 panels + composite\n")
