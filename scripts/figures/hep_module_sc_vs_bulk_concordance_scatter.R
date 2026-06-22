#!/usr/bin/env Rscript
# ============================================================================
# hep_module_sc_vs_bulk_concordance_scatter.R
# Fig 3 supplement (figS_hotspot) — hepatocyte Hotspot module sc-vs-bulk
# concordance.
#
# Claim: of 30 hepatocyte Hotspot modules, those tracking the leak-safe
# disease stage (disease_stage_q<0.05) and replicating by direction in bulk
# align — sc stage-beta (gaining/declining) matches bulk mean log2FC
# (up/down). Quadrant I (gain+bulk-up) and III (decline+bulk-down) carry the
# coordinated programs.
#
# Scatter: x = disease_stage_beta (leak-safe coarse 4-bin axis; NEVER
# F_stage_beta = scVI leakage). y = mean_dream_logFC (bulk). Solid = bulk
# replicated; open circle when disease_stage_q >= 0.05 (non-progression).
# hep-19 is cirrhosis-dependent -> dashed ring + footnote.
#
# Output: figures/supplementary/figS_hotspot/panels/
#         hep_module_sc_vs_bulk_concordance_scatter.pdf
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- FIGS_HOTSPOT_PANELS_DIR
OUT_PDF <- file.path(OUT_DIR, "hep_module_sc_vs_bulk_concordance_scatter.pdf")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load + filter to hepatocytes
# ---------------------------------------------------------------------------
mods <- fread(file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/hotspot_modules/all_modules.tsv"))
d <- mods[cell_type == "hepatocytes"]

d[, mod_lab := paste0("hep-", module)]
d[, progresses := disease_stage_q < 0.05]
d[, bulk_rep := bulk_replicated == TRUE]

# hep-19 is cirrhosis-dependent (see project notes) -> flag for dashed ring.
CIRR_DEP <- "hep-19"
d[, cirr_dep := mod_lab == CIRR_DEP]

# HKDC1 lives in hepatocyte module 24 (module_genes.tsv). Label hero modules.
hero_labels <- c("hep-19", "hep-20", "hep-24", "hep-26", "hep-27")
d[, do_label := mod_lab %in% hero_labels]
d[mod_lab == "hep-24", mod_lab_show := "hep-24 (HKDC1)"]
d[mod_lab != "hep-24", mod_lab_show := mod_lab]

# ---------------------------------------------------------------------------
# Hero numbers
# ---------------------------------------------------------------------------
n_total <- nrow(d)
n_prog  <- sum(d$progresses)
n_rep   <- sum(d$bulk_rep)
cat(sprintf("[hero] Hepatocyte Hotspot modules: total=%d | disease_stage_q<0.05=%d | bulk_replicated=%d\n",
            n_total, n_prog, n_rep))
# Concordance of the progressing modules: do sc-beta and bulk-LFC share sign?
prog <- d[progresses == TRUE]
conc <- sum(sign(prog$disease_stage_beta) == sign(prog$mean_dream_logFC))
cat(sprintf("[hero] Among %d stage-progressing modules, %d share sc-beta / bulk-logFC sign (%.0f%% directional concordance)\n",
            nrow(prog), conc, 100 * conc / nrow(prog)))

dir.create(FIGS_HOTSPOT_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
fwrite(d[, .(cell_type, module, mod_lab, disease_stage_beta, disease_stage_q,
             mean_dream_logFC, pct_concordant_up, bulk_replicated,
             progresses, cirr_dep)],
       file.path(FIGS_HOTSPOT_DATA_DIR, "hep_module_sc_vs_bulk_concordance_scatter.csv"))

# ---------------------------------------------------------------------------
# Aesthetics:
#   fill color   = bulk direction sign (magenta up / blue down)
#   open vs solid= disease_stage_q < 0.05 (solid) vs >= 0.05 (open/hollow)
#   ring linetype= solid normally, dashed for cirrhosis-dependent hep-19
# ---------------------------------------------------------------------------
d[, dir := ifelse(mean_dream_logFC >= 0, "Bulk up", "Bulk down")]
dir_colors <- c("Bulk up" = "#C9265E", "Bulk down" = "#1565C0")
d[, fillcol := ifelse(progresses, dir_colors[dir], "white")]

rng_x <- max(abs(d$disease_stage_beta)) * 1.12
rng_y <- max(abs(d$mean_dream_logFC)) * 1.18

p <- ggplot(d, aes(x = disease_stage_beta, y = mean_dream_logFC)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "#9E9E9E") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "#9E9E9E") +
  # progressing modules: solid fill by direction
  geom_point(data = d[progresses == TRUE & cirr_dep == FALSE],
             aes(fill = dir), shape = 21, color = "grey15",
             size = 2.6, stroke = 0.5) +
  # non-progressing: hollow / open circle
  geom_point(data = d[progresses == FALSE],
             shape = 21, fill = "white", color = "grey45",
             size = 2.2, stroke = 0.5) +
  # cirrhosis-dependent hep-19: heavier ring + outer dashed open ring marker
  geom_point(data = d[cirr_dep == TRUE],
             aes(fill = dir), shape = 21, color = "grey15",
             size = 2.6, stroke = 0.5) +
  geom_point(data = d[cirr_dep == TRUE],
             shape = 1, color = "grey15", size = 4.2, stroke = 0.5) +
  geom_text_repel(data = d[do_label == TRUE],
                  aes(label = mod_lab_show),
                  size = 2.0, color = "grey15", segment.size = 0.2,
                  min.segment.length = 0, max.overlaps = Inf,
                  box.padding = 0.3, seed = 42) +
  scale_fill_manual(values = dir_colors, name = "Bulk direction") +
  scale_x_continuous(limits = c(-rng_x, rng_x),
                     name = "sc disease-stage slope (coarse 4-bin axis)") +
  scale_y_continuous(limits = c(-rng_y, rng_y),
                     name = expression("Bulk mean log"[2]*" fold change")) +
  labs(title = "Hepatocyte modules: single-cell stage track vs bulk direction",
       caption = paste0("Solid = bulk-replicated stage modules; open = non-progressing (q≥0.05). ",
                        "hep-19 dashed = cirrhosis-dependent.")) +
  theme_masld(base_size = 7) +
  theme(
    plot.title    = element_text(size = 7.3, face = "bold", margin = margin(b = 5)),
    plot.caption  = element_text(size = 4.8, color = "grey35", hjust = 0),
    legend.position = "right",
    legend.text   = element_text(size = 5.5),
    legend.title  = element_text(size = 6),
    legend.key.size = unit(0.16, "cm")
  )

ggsave(OUT_PDF, p,
       width  = 90 / 25.4,
       height = 78 / 25.4,
       units  = "in",
       device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))
