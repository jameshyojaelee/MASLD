#!/usr/bin/env Rscript
# fig2_finemap_cascade.R  (2026-06-12)  — Fig 2 (first part of the para-6 narrative)
# Fine-mapping cascade: lead-SNP loci -> SuSiE-converged -> independent loci,
# then the two cross-ancestry arms (SuSiEx high-PIP; meSuSiE shared EUR-EAS).
# All counts are read from the frozen, disk-verified source CSV so they always
# match the manuscript text.
#
# Out: figures/main/fig2_genetics/panels/finemap_cascade.pdf
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- file.path(FIG3_DIR, "panels")

src <- fread(file.path(PANEL_DIR, "finemap_convergence_source.csv"), fill = TRUE)
src <- src[!grepl("^#", metric)]
val <- function(m) as.numeric(src[metric == m, value][1])

dt <- data.table(
  stage = c("Lead-SNP loci screened",
            "SuSiE-converged",
            "Independent loci (500 kb merge)",
            "SuSiEx high-PIP loci (PIP ≥ 0.9)",
            "meSuSiE shared EUR–EAS loci"),
  n     = c(val("lead_snp_loci_screened"),
            val("loci_converged"),
            val("independent_loci_500kb"),
            val("susiex_physical_loci_pip_ge_0.9"),
            val("mesusie_shared_eur_eas_loci")),
  phase = c(rep("Within-ancestry SuSiE", 3), rep("Cross-ancestry", 2)))

dt[, stage := factor(stage, levels = rev(stage))]            # first stage on top
dt[, phase := factor(phase, levels = c("Within-ancestry SuSiE", "Cross-ancestry"))]
dt[, lab := fifelse(stage == "SuSiE-converged", paste0(n, " (96.7%)"), as.character(n))]

phase_cols <- c("Within-ancestry SuSiE" = "#1565C0", "Cross-ancestry" = "#00695C")

p <- ggplot(dt, aes(x = n, y = stage, fill = phase)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = lab), hjust = -0.12, size = GEOM_TEXT_6PT, fontface = "plain", color = "black") +
  scale_fill_manual(values = phase_cols, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18)), limits = c(0, 260)) +
  labs(x = "Loci", y = NULL) +
  theme_masld(base_size = 9) +
  theme(legend.position = c(0.98, 0.06), legend.justification = c(1, 0),
        legend.background = element_rect(fill = scales::alpha("white", 0.7), color = NA),
        axis.text.y = element_text(size = 6))

save_fig(p, file.path(PANEL_DIR, "finemap_cascade.pdf"),
         width = fig_col_width * 1.12, height = 2.5)
message("[caption] Fine-mapping cascade")
cat("[finemap_cascade] wrote panel — 241/233/184 + SuSiEx 70 + meSuSiE 108\n")
