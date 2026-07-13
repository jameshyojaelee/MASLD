#!/usr/bin/env Rscript
# fig4_variant_motif_disruption.R
# ============================================================================
# Fig 4 Panel B — Sequence-level variant -> TF-motif disruption at MASLD loci.
#
# WHY THIS PANEL IS VALID (and what it does NOT claim):
#   Fig 4 is the VALIDATION figure. The n=18 multiome cohort (5 control / 13
#   disease) cannot make a per-cell or per-condition disease-DE claim — that is
#   pseudoreplicated (Squair 2021) and the donor-level regulon/chromVAR DE is
#   honestly 0-significant (sc_underpowered = TRUE). This panel sidesteps all of
#   that: it is a SEQUENCE-LEVEL, sample-size-INDEPENDENT measurement. For each
#   GWAS fine-mapped credible-set variant we compute the allelic difference in
#   the TF position-weight-matrix match score (motifbreakR `alleleDiff`) between
#   the reference and alternate alleles. This is a property of the genome
#   sequence + the PWM, not of any cohort. It CORROBORATES the well-powered
#   Fig 2 (GWAS/COLOC) + Fig 3 (bulk RNA-seq, n=846) disease master regulators
#   by showing that the same disease-regulon TF motifs are physically disrupted
#   by the fine-mapped risk alleles.
#
#   `alleleDiff` is the SIGNED motif-score difference (alt vs ref): positive =
#   the alternate (risk) allele STRENGTHENS the TF motif (predicted gain);
#   negative = the alternate allele WEAKENS it (predicted loss). |alleleDiff| is
#   the disruption MAGNITUDE. This is DESCRIPTIVE of disruption magnitude — there
#   is NO significance test on this axis and NO p-value is shown (Refpvalue /
#   Altpvalue are not populated in the source table; max_pip is the variant
#   fine-mapping posterior, used only to weight point size).
#
# Source (read directly; do not trust column comments):
#   GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv
#     cols used: SNP_id, tf_name, alleleDiff, max_pip, effect,
#                motif_in_disease_regulon, ref_genome, alt_genome
#
# Filter: motif_in_disease_regulon == TRUE, restricted to the disease master-
#   regulator TFs HNF4A, RORA, THRB, TCF7L2 (Fig 2/Fig 3 hits).
#
# Output: figures/main/fig4_validation/fig4b_variant_motif_disruption.pdf (PDF only)
# ============================================================================

suppressMessages({
  library(data.table)
  library(ggplot2)
})

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

MOTIF_FILE <- file.path(BASE_DIR,
  "GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv")
OUT_DIR <- file.path(BASE_DIR, "figures/main/fig4_validation/_supp")  # demoted to supp (2026-06-22): redundant with atac_rora_motif
OUT_PDF <- file.path(OUT_DIR, "fig4b_variant_motif_disruption.pdf")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Disease master regulators (Fig 2 / Fig 3 hits). HNF4A/RORA/THRB/TCF7L2 required.
DISEASE_TFS <- c("HNF4A", "RORA", "THRB", "TCF7L2")

# Publication theme (theme_masld, save_fig, fig_half_width, etc.)
theme_file <- file.path(BASE_DIR, "scripts/figures/publication_theme.R")
if (file.exists(theme_file)) {
  source(theme_file)
} else {
  theme_masld <- function(base_size = 7) theme_classic(base_size = base_size)
  save_fig <- function(plot, filename, width = 7, height = 5, dpi = 300)
    ggsave(filename, plot, width = width, height = height, device = cairo_pdf)
  fig_half_width <- 88 / 25.4
}

# ── 1. Load + filter ────────────────────────────────────────────────────────
cat("Reading", MOTIF_FILE, "\n")
m <- fread(MOTIF_FILE)
cat("  total variant-motif rows:", nrow(m), "\n")

reg <- m[motif_in_disease_regulon == TRUE & tf_name %in% DISEASE_TFS]
# Source table carries exact duplicate rows per (SNP, TF); collapse to unique.
reg <- unique(reg, by = c("SNP_id", "tf_name"))
cat("  disease-regulon rows for", paste(DISEASE_TFS, collapse = "/"), ":",
    nrow(reg), "unique (SNP, TF) pairs\n")

stopifnot(nrow(reg) > 0)
stopifnot(all(DISEASE_TFS %in% reg$tf_name))

# Allelic direction (descriptive only — sign of the motifbreakR score diff).
reg[, direction := ifelse(alleleDiff >= 0, "Alt strengthens motif (gain)",
                                           "Alt weakens motif (loss)")]
reg[, direction := factor(direction,
      levels = c("Alt strengthens motif (gain)", "Alt weakens motif (loss)"))]

# Order TF facets by the strongest single disruption magnitude they carry,
# then order variants within a TF by signed alleleDiff for a clean strip.
tf_order <- reg[, .(m = max(abs(alleleDiff))), by = tf_name][order(-m), tf_name]
reg[, tf_name := factor(tf_name, levels = tf_order)]
setorder(reg, tf_name, alleleDiff)
reg[, snp_y := factor(SNP_id, levels = unique(SNP_id))]

# ── 2. Console: paste the exact rows the panel shows ────────────────────────
cat("\n--- Variant -> TF disruption rows shown in panel ---\n")
show <- reg[, .(SNP_id, tf_name = as.character(tf_name),
                alleleDiff = round(alleleDiff, 3),
                max_pip = round(max_pip, 4), effect,
                ref = ref_genome, alt = alt_genome)]
setorder(show, tf_name, -max_pip)
print(show, nrows = 100)

# ── 3. Panel ────────────────────────────────────────────────────────────────
# Strip / dot panel (NOT a lollipop — house rule): per disease-regulon TF, each
# fine-mapped variant placed at its SIGNED allelic motif-score difference,
# point size = fine-mapping posterior (max_pip), colour = predicted direction.
# Zero line = no allelic motif change. Distance from zero = disruption magnitude.
GAIN_COL <- "#C2185B"   # magenta — alt allele strengthens the motif
LOSS_COL <- "#1565C0"   # blue    — alt allele weakens the motif

xr <- range(reg$alleleDiff)
xpad <- 0.55
x_lim <- c(min(xr[1], 0) - xpad, max(xr[2], 0) + xpad)

# NOTE: dots only (NO geom_segment stem) — a stem-to-zero + ball is a lollipop,
# which is a hard house-rule ban. Position relative to the zero reference line
# already encodes both direction (left/right) and disruption magnitude (distance).
p <- ggplot(reg, aes(x = alleleDiff, y = snp_y)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey55", linewidth = 0.3) +
  geom_point(aes(colour = direction, size = max_pip)) +
  facet_grid(tf_name ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = c("Alt strengthens motif (gain)" = GAIN_COL,
                                 "Alt weakens motif (loss)"     = LOSS_COL),
                      name = NULL, drop = FALSE) +
  scale_size_continuous(range = c(1.0, 3.6), limits = c(0, 1),
                        breaks = c(0.1, 0.5, 0.9),
                        name = "Variant\nfine-map PIP") +
  scale_x_continuous(limits = x_lim, breaks = seq(-2, 2, 1)) +
  labs(
    x = expression(paste("Allelic motif-score difference (alt - ref), motifbreakR ", Delta)),
    y = NULL,
    caption = "Sequence-level, sample-size-independent; magnitude only, not a test"
  ) +
  theme_masld(base_size = 7) +
  theme(
    axis.text.y      = element_text(size = 6, colour = "black"),
    strip.text.y.left = element_text(angle = 0, face = "italic", size = 6,
                                     colour = "black"),
    strip.placement  = "outside",
    panel.spacing.y  = unit(2, "pt"),
    plot.caption     = element_text(size = 6, colour = "black", hjust = 0),
    legend.position  = "right",
    legend.box       = "vertical",
    legend.key.size  = unit(0.28, "cm")
  ) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 2.4)),
         size   = guide_legend(order = 2))

save_fig(p, OUT_PDF, width = fig_half_width + 1.2, height = 4.4)
cat("\nSaved:", OUT_PDF, "\n")

# Panel-level table for the manuscript ledger (PDF is the deliverable; this is
# a side artifact under the same fig4 dir, not a separate figure).
fwrite(show, file.path(OUT_DIR, "fig4b_variant_motif_disruption_rows.csv"))
cat("Wrote row table:", file.path(OUT_DIR, "fig4b_variant_motif_disruption_rows.csv"), "\n")
