#!/usr/bin/env Rscript
# ============================================================================
# nmf_programs.R  — RETIRED 2026-06-12 (superseded by nmf_programs_lines.R).
#   BUG: the hardcoded program-label map below (P1=Stable, ... P6=Quiescent-
#   Parenchyma_3) is stale vs the current program_labels.csv (P1=Progression-
#   Inflammatory, P5=Stable_1, P6=Stable_2). The relabel looked for "Stable"
#   (no longer a label), so Stable_1 + Stable_2 (291 samples) fell to NA and
#   were dropped, leaving only 4 of 6 programs with renormalised (wrong) fracs.
#   The canonical NMF panel is now the dynamic, label-driven nmf_programs_lines.R.
# Fig 2 panel f — k=6 NMF dominant-program composition across fibrosis stages
#
# Stacked bar showing fraction of samples assigned to each of 6 NMF programs
# at each fibrosis stage F0..F4.
#
# Sources:
#   - RNA-seq/results/subtypes/nmf_assignments.csv (sample_id, dominant_program)
#   - RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv (fibrosis_stage)
#
# Output: figures/main/fig3_RNAseq/panels/nmf_programs.pdf
#   sized 60 x 55 mm.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG2_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
OUT_PDF   <- file.path(PANEL_DIR, "nmf_programs.pdf")

# ---------------------------------------------------------------------------
# Stacked bar: dominant_program x fibrosis_stage
# ---------------------------------------------------------------------------
nmf <- fread(file.path(BASE, "RNA-seq/results/subtypes/nmf_assignments.csv"))
meta <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"))
meta_use <- meta[, .(sample_id, fibrosis_stage, condition)]
d <- merge(nmf[, .(sample_id, dominant_program, dominant_program_code)],
           meta_use, by = "sample_id", all.y = TRUE)

d[, stage := as.character(fibrosis_stage)]
d[condition == "Control", stage := "Control"]
d[condition != "Control" & !is.na(fibrosis_stage), stage := paste0("F", fibrosis_stage)]

d <- d[stage %in% c("F0", "F1", "F2", "F3", "F4")]
d[, fibrosis_stage := factor(stage, levels = c("F0", "F1", "F2", "F3", "F4"))]

# Program labels — v2 clean gene pool (pseudogene-stripped)
# P1 Stable (GABRG2/RBFOX1/MUC17, neural-like, stage_rho=-0.4)
# P2 Progression-Inflammatory (AKR1B10/CCL20/TREM2, stage_rho=+0.9)
# P3 Fibrogenic (DES/CNN1/MYOCD, stage_rho=+0.9)
# P4 Quiescent-Parenchyma_1 (CCDC144A/NEUROD6, stage_rho=-0.8)
# P5 Quiescent-Parenchyma_2 (CALML6/SCX/CIMIP2A, stage_rho=-1.0)
# P6 Quiescent-Parenchyma_3 (TIMD4/CXCR1-2/FCGR3B, innate immune, stage_rho=-1.0)
prog_levels <- c("Stromal",
                  "Inflammatory",
                  "Fibrogenic",
                  "Quiescent-1",
                  "Quiescent-2",
                  "Kupffer-cell")
prog_codes  <- c("P1", "P2", "P3", "P4", "P5", "P6")
prog_pal <- c(
  "Stromal"       = "#F4A674",  # warm peach (masld_colors$nafl — mild/background)
  "Inflammatory"  = "#C9265E",  # Liang magenta (masld_colors$nash — primary disease)
  "Fibrogenic"    = "#1565C0",  # deep blue (fibrosis stage F3/F4 gradient)
  "Quiescent-1"   = "#9E9E9E",  # control gray
  "Quiescent-2"   = "#BDBDBD",  # light gray
  "Kupffer-cell"  = "#00695C"   # teal (masld_colors$conserved — protective/resident)
)
d[dominant_program == "Stable",                   dominant_program := "Stromal"]
d[dominant_program == "Progression-Inflammatory", dominant_program := "Inflammatory"]
d[dominant_program == "Quiescent-Parenchyma_1",   dominant_program := "Quiescent-1"]
d[dominant_program == "Quiescent-Parenchyma_2",   dominant_program := "Quiescent-2"]
d[dominant_program == "Quiescent-Parenchyma_3",   dominant_program := "Kupffer-cell"]
d[, dominant_program := factor(dominant_program, levels = prog_levels)]
d <- d[!is.na(dominant_program)]

# Compute dominant-program fraction per stage
comp <- d[, .N, by = .(fibrosis_stage, dominant_program)]
comp[, frac := N / sum(N), by = fibrosis_stage]
n_per_stage <- d[, .N, by = fibrosis_stage]
setorder(comp, fibrosis_stage, dominant_program)

# Hero numbers
p2_f1 <- comp[fibrosis_stage == "F1" & dominant_program == "Inflammatory", frac]
p2_f2 <- comp[fibrosis_stage == "F2" & dominant_program == "Inflammatory", frac]
p3_f2 <- comp[fibrosis_stage == "F2" & dominant_program == "Fibrogenic",   frac]
p3_f3 <- comp[fibrosis_stage == "F3" & dominant_program == "Fibrogenic",   frac]
cat(sprintf("[hero] Inflammatory F1=%.1f%% F2=%.1f%% (%.1fx); Fibrogenic F2=%.1f%% F3=%.1f%% (%.1fx)\n",
            100*p2_f1, 100*p2_f2, p2_f2/p2_f1,
            100*p3_f2, 100*p3_f3, p3_f3/p3_f2))

p_stack <- ggplot(comp, aes(x = fibrosis_stage, y = frac,
                              fill = dominant_program)) +
  geom_col(width = 0.78, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = prog_pal, name = NULL) +
  scale_y_continuous(name = "Fraction of samples",
                     labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1), expand = c(0, 0),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(x = NULL,
       title = "NMF programs across fibrosis stage") +
  theme_masld(base_size = 7) +
  theme(
    plot.title    = element_text(size = 7.3, face = "bold", margin = margin(b = 6)),
    plot.subtitle = element_text(size = 5.5, color = "grey35"),
    axis.text.x   = element_text(size = 6.5, face = "bold"),
    legend.position = "right",
    legend.text   = element_text(size = 5),
    legend.key.size = unit(0.17, "cm")
  )

ggsave(OUT_PDF, p_stack,
       width  = 80 / 25.4,
       height = 70 / 25.4,
       units  = "in",
       device = cairo_pdf)
fwrite(comp, file.path(DATA_DIR, "nmf_composition.csv"))
cat(sprintf("[saved] %s\n", OUT_PDF))
