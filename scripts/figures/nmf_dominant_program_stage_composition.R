#!/usr/bin/env Rscript
# ============================================================================
# nmf_dominant_program_stage_composition.R
# Fig 3 (RNA-seq) — winner-take-all NMF dominant-program composition across
# fibrosis stage.
#
# 100% stacked bar: x = documented fibrosis stage F0..F4, fill = dominant
# NMF program (biological label). Quiescent programs dominate ~2/3 of F0;
# the Inflammatory-EMT program rises with stage; the Fibrotic-ECM  # was: Progression-Inflammatory / Fibrogenic
# program expands sharply in advanced fibrosis.
#
# Controls excluded (NMF fitted on disease samples only); NA-stage excluded.
# Labels read from program_labels.csv (NEVER hardcode P-code -> label).
#
# Output: figures/main/fig3_RNAseq/panels/fig3k_nmf_dominant_program_stage_composition.pdf
#   sized 80 x 70 mm
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
OUT_PDF   <- file.path(PANEL_DIR, "fig3k_nmf_dominant_program_stage_composition.pdf")

# ---------------------------------------------------------------------------
# Load and merge
# ---------------------------------------------------------------------------
nmf  <- fread(file.path(BASE, "RNA-seq/results/subtypes/nmf_assignments.csv"))
meta <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"))
prog_lab <- fread(file.path(BASE, "RNA-seq/results/subtypes/program_labels.csv"))
setorder(prog_lab, program_code)

# Resolve each sample's dominant program label from its code via the
# authoritative label table (do not trust the free-text dominant_program col).
nmf <- merge(
  nmf[, .(sample_id, dominant_program_code)],
  prog_lab[, .(program_code, biological_label)],
  by.x = "dominant_program_code", by.y = "program_code"
)

# Clean staged disease samples only (excludes coarse-staged GSE213621 [F0F1/F3F4
# grouped, not true Kleiner] + dropped GSE/PRJNA512027; provides the F0-F4 `stage`
# factor). See load_figure_data.R::staged_disease_meta(). `meta` read above is no
# longer used for stage binning.
stm <- staged_disease_meta()
d <- merge(nmf, stm[, .(sample_id, stage)], by = "sample_id")

# ---------------------------------------------------------------------------
# Program color / order (label-driven). Magenta = Inflammatory-EMT,  # was: Progression-Inflammatory
# blue = Fibrotic-ECM, neutrals for quiescent/stable programs.  # was: Fibrogenic
# ---------------------------------------------------------------------------
prog_meta <- data.table(program_code = prog_lab$program_code, label = prog_lab$biological_label)
neutral_pool <- c("#9E9E9E", "#BDBDBD", "#6D6D6D", "#424242", "#00695C")
prog_meta[, color := NA_character_]
prog_meta[grepl("Inflammatory", label, ignore.case = TRUE), color := "#C9265E"]
prog_meta[program_code == "P3", color := "#1565C0"]  # was: grepl("Fibrogenic")
rest <- which(is.na(prog_meta$color))
prog_meta$color[rest] <- neutral_pool[seq_along(rest)]
prog_colors <- setNames(prog_meta$color, prog_meta$label)

# ---------------------------------------------------------------------------
# Composition: % of samples per stage assigned to each dominant program
# ---------------------------------------------------------------------------
comp <- d[, .N, by = .(stage, biological_label)]
comp[, frac := N / sum(N), by = stage]
comp[, label := factor(biological_label, levels = prog_meta$label)]
setorder(comp, stage, label)

fwrite(comp, file.path(DATA_DIR, "nmf_dominant_program_stage_composition.csv"))

# ---------------------------------------------------------------------------
# Hero numbers
# ---------------------------------------------------------------------------
inflam_lab <- prog_meta$label[grepl("Inflammatory", prog_meta$label, ignore.case = TRUE)][1]
fibro_lab  <- prog_meta$label[prog_meta$program_code == "P3"][1]  # was: grepl("Fibrogenic")
getpct <- function(lab, st) {
  v <- comp[label == lab & stage == st, frac]
  if (length(v) == 0) 0 else 100 * v
}
cat("[hero] Inflammatory-EMT dominant %% per stage:\n")  # was: "Progression-Inflammatory"
for (st in paste0("F", 0:4)) cat(sprintf("    %s = %.1f%%\n", st, getpct(inflam_lab, st)))
cat("[hero] Fibrotic-ECM dominant %% per stage:\n")  # was: "Fibrogenic"
for (st in paste0("F", 0:4)) cat(sprintf("    %s = %.1f%%\n", st, getpct(fibro_lab, st)))
cat(sprintf("[hero] Inflammatory-EMT F0=%.0f%% -> F3=%.0f%%; Fibrotic-ECM F2=%.0f%% -> F4=%.0f%%\n",  # was: "Progression-Inflammatory"/"Fibrogenic"
            getpct(inflam_lab, "F0"), getpct(inflam_lab, "F3"),
            getpct(fibro_lab, "F2"),  getpct(fibro_lab, "F4")))

# ---------------------------------------------------------------------------
# Transition annotations
# ---------------------------------------------------------------------------
ann <- data.table(
  x     = c("F3", "F4"),
  y     = c(0.5, 0.5),
  txt   = c(sprintf("Inflammatory\n%.0f%%→%.0f%%", getpct(inflam_lab, "F0"), getpct(inflam_lab, "F3")),
            sprintf("Fibrotic-ECM\n%.0f%%→%.0f%%",   getpct(fibro_lab, "F2"),  getpct(fibro_lab, "F4")))  # was: "Fibrogenic"
)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
p <- ggplot(comp, aes(x = stage, y = frac, fill = label)) +
  geom_col(width = 0.82, color = "white", linewidth = 0.15) +
  geom_text(data = ann, aes(x = x, y = y, label = txt),
            inherit.aes = FALSE, size = 1.9, fontface = "plain",
            color = "white", lineheight = 0.9) +
  scale_fill_manual(values = prog_colors, name = NULL) +
  scale_y_continuous(
    name   = "Samples (dominant program)",
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = "Fibrosis stage",
       title = "NMF dominant-program composition reorganizes with fibrosis") +
  theme_masld(base_size = 7) +
  theme(
    plot.title      = element_text(size = 7.0, face = "plain", margin = margin(b = 6)),
    axis.text.x     = element_text(size = 6.5, face = "plain"),
    legend.position = "right",
    legend.text     = element_text(size = 5.2),
    legend.key.size = unit(0.16, "cm")
  )

ggsave(OUT_PDF, p,
       width  = 88 / 25.4,
       height = 70 / 25.4,
       units  = "in",
       device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))
message("[fig2d caption] Dominant-program composition is a HARD argmax on the k=6 NMF (seed 42). ",
        "This partition is seed-fragile: cross-seed dominant-program ARI ≈ 0.50 (only ~half of donors ",
        "retain their program label across seeds); k=4 is more reproducible (cophenetic 0.875 vs k=6 0.757). ",
        "Read the cascade percentages qualitatively — one-seed illustrations of the trajectory, not exact composition.")
