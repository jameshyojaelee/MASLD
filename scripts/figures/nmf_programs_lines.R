#!/usr/bin/env Rscript
# ============================================================================
# nmf_programs_lines.R
# Fig 2 panel f (alt) — k=6 NMF program activity across fibrosis stages
#
# Line plot of mean normalised program score ± SE across F0..F4.
# Each sample's 6 raw scores are normalised to sum to 1 (relative composition)
# before averaging, preserving the continuous signal that winner-take-all
# dominant-program composition discards.
#
# Controls excluded: NMF was fitted on disease samples only; controls carry
# no valid program scores.
#
# Output: figures/main/fig3_RNAseq/panels/nmf_programs_lines.pdf
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
OUT_PDF   <- file.path(PANEL_DIR, "nmf_programs_lines.pdf")

# ---------------------------------------------------------------------------
# Load and merge
# ---------------------------------------------------------------------------
nmf  <- fread(file.path(BASE, "RNA-seq/results/subtypes/nmf_assignments.csv"))
meta <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"))

# ---------------------------------------------------------------------------
# Program labels read from the authoritative program_labels.csv (Script 44) so
# labels/order always match the current canonical k-program fit. NEVER hardcode
# P-code -> label (the fit can be re-run at a different k / with relabelling).
# ---------------------------------------------------------------------------
prog_lab <- fread(file.path(BASE, "RNA-seq/results/subtypes/program_labels.csv"))
setorder(prog_lab, program_code)          # P1, P2, ... order
score_cols <- prog_lab$program_code       # bare per-program score cols in nmf_assignments

d <- merge(
  nmf[, c("sample_id", score_cols), with = FALSE],
  meta[, .(sample_id, fibrosis_stage, condition)],
  by = "sample_id"
)

# Exclude controls (no valid NMF scores) and samples without fibrosis stage
d <- d[condition != "Control" & !is.na(fibrosis_stage)]
d[, stage := factor(paste0("F", fibrosis_stage), levels = paste0("F", 0:4))]

# ---------------------------------------------------------------------------
# Normalise each sample's scores to relative composition (sum to 1)
# ---------------------------------------------------------------------------
d[, row_sum := rowSums(.SD), .SDcols = score_cols]
for (p in score_cols) d[, (paste0(p, "_norm")) := get(p) / row_sum]

# ---------------------------------------------------------------------------
# Program metadata (label-driven). Magenta = Progression-Inflammatory, blue =
# Fibrogenic, distinct neutrals for the quiescent/stable programs.
# ---------------------------------------------------------------------------
prog_meta <- data.table(
  code  = paste0(prog_lab$program_code, "_norm"),
  label = prog_lab$biological_label
)
neutral_pool <- c("#9E9E9E", "#BDBDBD", "#6D6D6D", "#424242", "#00695C")
prog_meta[, color := NA_character_]
prog_meta[grepl("Inflammatory", label, ignore.case = TRUE), color := "#C9265E"]
prog_meta[grepl("Fibrogenic",  label, ignore.case = TRUE), color := "#1565C0"]
rest <- which(is.na(prog_meta$color))
prog_meta$color[rest] <- neutral_pool[seq_along(rest)]

# ---------------------------------------------------------------------------
# Mean ± SE per stage per program
# ---------------------------------------------------------------------------
norm_cols <- paste0(score_cols, "_norm")
long <- melt(d[, c("stage", norm_cols), with = FALSE],
             id.vars     = "stage",
             variable.name = "code",
             value.name    = "score")
long <- merge(long, prog_meta, by = "code")
long[, label := factor(label, levels = prog_meta$label)]

summary_dt <- long[, .(
  mean = mean(score),
  se   = sd(score) / sqrt(.N),
  n    = .N
), by = .(stage, label, color)]

# Emit hero numbers (reference the actual canonical labels)
inflam_lab <- prog_meta$label[grepl("Inflammatory", prog_meta$label, ignore.case = TRUE)][1]
fibro_lab  <- prog_meta$label[grepl("Fibrogenic",  prog_meta$label, ignore.case = TRUE)][1]
p_inflam <- summary_dt[label == inflam_lab]
p_fibro  <- summary_dt[label == fibro_lab]
cat(sprintf("[hero] %s F0=%.1f%% -> F4=%.1f%% (+%.1f pp)\n", inflam_lab,
            100 * p_inflam[stage == "F0", mean],
            100 * p_inflam[stage == "F4", mean],
            100 * (p_inflam[stage == "F4", mean] - p_inflam[stage == "F0", mean])))
cat(sprintf("[hero] %s F0=%.1f%% -> F4=%.1f%% (+%.1f pp)\n", fibro_lab,
            100 * p_fibro[stage == "F0", mean],
            100 * p_fibro[stage == "F4", mean],
            100 * (p_fibro[stage == "F4", mean] - p_fibro[stage == "F0", mean])))

# Save summary table
fwrite(summary_dt, file.path(DATA_DIR, "nmf_programs_lines.csv"))

# ---------------------------------------------------------------------------
# Color map (stable order)
# ---------------------------------------------------------------------------
prog_colors <- setNames(prog_meta$color, prog_meta$label)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
p_lines <- ggplot(summary_dt,
                  aes(x = stage, y = mean, color = label, group = label)) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se, fill = label),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 1.6) +
  scale_color_manual(values = prog_colors, name = NULL) +
  scale_fill_manual( values = prog_colors, name = NULL) +
  scale_y_continuous(
    name   = "Mean relative program score",
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(x     = "Fibrosis stage",
       title = "NMF program activity across fibrosis stage") +
  theme_masld(base_size = 7) +
  theme(
    plot.title      = element_text(size = 7.3, face = "bold", margin = margin(b = 6)),
    axis.text.x     = element_text(size = 6.5, face = "bold"),
    legend.position = "right",
    legend.text     = element_text(size = 5.5),
    legend.key.size = unit(0.18, "cm")
  )

ggsave(OUT_PDF, p_lines,
       width  = 80 / 25.4,
       height = 70 / 25.4,
       units  = "in",
       device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))
