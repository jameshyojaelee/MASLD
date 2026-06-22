#!/usr/bin/env Rscript
# figS_cas13_library_core_vs_mash_overlap.R
# ---------------------------------------------------------------------------
# v8 Cas13 library composition. The human spine = integrated CORE (disease-vs-
# control, canonical limma-voom QW C2) UNION COHORT-REPLICATED (DE in >=2 of 5
# cohorts); the library adds mouse-confirmed, COLOC and positive-control tiers.
# miRNA and the MASH tier are retired (the legacy "core_vs_mash" filename is kept
# for runner compatibility; the figure now contrasts core vs cohort-replicated).
# v8 PC gate = mouse-hepatocyte scRNA substrate (Panel D shows the substrate split).
#
# Individual panels (PDF only):
#   A  cas13_v8_core_vs_cohort_euler.pdf  area-proportional 2-set Euler of the
#        human-spine structure: integrated Core vs Cohort-replicated (>=2 cohorts).
#   B  cas13_v8_evidence_upset.pdf        UpSet over library mouse genes by evidence
#        axis: Core / Cohort-replicated / Mouse cross-diet(>=3) / COLOC / Pos. control.
#   C  cas13_v8_tier_composition.pdf      stacked bar: tier x biotype gene counts.
#   D  cas13_v8_mouse_hep_substrate.pdf   PC genes by mouse-hepatocyte substrate
#        (high / ambient_suspect / absent); the v8 gate axis + exemption book-keeping.
#
# Output dir: Cas13_Library_Design/figures/ (FIGS_CAS13LIB_DIR).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggforce); library(UpSetR)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
OUT_DIR <- FIGS_CAS13LIB_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

LIB <- file.path(BASE, "Cas13_Library_Design/data/cas13_library_v3.0.csv")

COL_CORE   <- "#e35070"; COL_COHORT <- "#f0a050"; COL_MOUSE <- "#4baeef"
COL_COLOC  <- "#8e6fc7"; COL_CTRL   <- "#9E9E9E"
BIOTYPE_COLS <- c(protein_coding = "#5480C2", lncRNA = "#D4834A")

lib <- fread(LIB)

# ── Panel A: Core vs Cohort-replicated (mouse-gene level, from the library) ────
nC  <- sum(lib$is_core); nCo <- sum(lib$is_cohort_replicated)
nO  <- sum(lib$is_core & lib$is_cohort_replicated)
nCo_new <- sum(lib$is_cohort_replicated & !lib$is_core)
message(sprintf("Core=%d  Cohort-rep=%d  overlap=%d  cohort-only(net-new)=%d", nC, nCo, nO, nCo_new))

lens_area <- function(r1, r2, d) {
  if (d >= r1 + r2) return(0); if (d <= abs(r1 - r2)) return(pi * min(r1, r2)^2)
  a1 <- acos(pmin(1, pmax(-1, (d^2 + r1^2 - r2^2)/(2*d*r1))))
  a2 <- acos(pmin(1, pmax(-1, (d^2 + r2^2 - r1^2)/(2*d*r2))))
  r1^2*(a1 - sin(2*a1)/2) + r2^2*(a2 - sin(2*a2)/2)
}
solve_d <- function(r1, r2, target) {
  if (target <= 0) return((r1 + r2)*1.04)
  if (target >= pi*min(r1,r2)^2) return(abs(r1 - r2))
  uniroot(function(d) lens_area(r1, r2, d) - target, lower = abs(r1-r2)+1e-6, upper = r1+r2-1e-6)$root
}
rC <- sqrt(nC/pi); rCo <- sqrt(nCo/pi); d <- solve_d(rC, rCo, nO)
cxC <- -d/2; cxCo <- d/2
circles <- data.frame(x0 = c(cxC, cxCo), y0 = c(0,0), r = c(rC, rCo),
                      set = c("Core (integrated)", "Cohort-replicated (>=2)"))
labs <- data.frame(x = c(-(rC+rCo)/2, (rC-rCo)/2, (rC+rCo)/2), y = 0,
                   lab = c(nC - nO, nO, nCo - nO))
pA <- ggplot() +
  geom_circle(data = circles, aes(x0 = x0, y0 = y0, r = r, fill = set), colour = NA, alpha = 0.55) +
  geom_text(data = labs, aes(x, y, label = lab), size = 3.0, fontface = "bold") +
  annotate("text", x = cxC, y = rC + 1.5, label = sprintf("Core\n%d", nC),
           colour = COL_CORE, size = 2.6, fontface = "bold", lineheight = 0.9) +
  annotate("text", x = cxCo, y = rCo + 1.5, label = sprintf("Cohort >=2\n%d", nCo),
           colour = "#b06a1f", size = 2.6, fontface = "bold", lineheight = 0.9) +
  scale_fill_manual(values = c("Core (integrated)" = COL_CORE,
                               "Cohort-replicated (>=2)" = COL_COHORT), guide = "none") +
  coord_fixed(clip = "off") +
  labs(title = "Human spine: integrated Core vs cohort-replicated",
       subtitle = sprintf("%d cohort-replicated genes are net-new vs the integrated core", nCo_new)) +
  theme_void(base_size = 7) +
  theme(plot.title = element_text(size = 7, hjust = 0.5),
        plot.subtitle = element_text(size = 6.5, hjust = 0.5, colour = "grey30"),
        plot.margin = margin(14, 6, 6, 6))
save_fig(pA, file.path(OUT_DIR, "cas13_v8_core_vs_cohort_euler.pdf"),
         width = fig_half_width, height = fig_half_width * 0.95)
message("wrote cas13_v8_core_vs_cohort_euler.pdf")

# ── Panel B: UpSet over 5 evidence axes ───────────────────────────────────────
sets <- lapply(list(
  `Core`              = lib[is_core == TRUE, gene_id_mouse],
  `Cohort-replicated` = lib[is_cohort_replicated == TRUE, gene_id_mouse],
  `Mouse cross-diet`  = lib[n_diets_up >= 3, gene_id_mouse],
  `COLOC`             = lib[has_coloc == TRUE, gene_id_mouse],
  `Positive control`  = lib[is_positive_control == TRUE, gene_id_mouse]), unique)
pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
pdf_device(file.path(OUT_DIR, "cas13_v8_evidence_upset.pdf"), width = fig_full_width, height = 3.6)
print(upset(fromList(sets), order.by = "freq", nsets = 5, nintersects = 24,
            sets.bar.color = c(COL_CORE, COL_COHORT, COL_MOUSE, COL_COLOC, COL_CTRL),
            main.bar.color = "#37474F", matrix.color = "#37474F",
            mainbar.y.label = "library genes", sets.x.label = "evidence-axis size", text.scale = 1.1))
dev.off()
message("wrote cas13_v8_evidence_upset.pdf")

# ── Panel C: tier x biotype composition ───────────────────────────────────────
tier_lvls <- c("core", "cohort_replicated", "mouse_confirmed", "coloc", "positive_control")
tier_labs <- c("Core\n(integrated)", "Cohort\nreplicated", "Mouse\nconfirmed", "COLOC", "Positive\ncontrol")
comp <- lib[, .N, by = .(tier, biotype)]
comp[, tier := factor(tier, levels = tier_lvls, labels = tier_labs)]
comp[, biotype := factor(biotype, levels = c("protein_coding", "lncRNA"))]
totals <- comp[, .(N = sum(N)), by = tier]
pC <- ggplot(comp, aes(tier, N, fill = biotype)) +
  geom_col(width = 0.7) +
  geom_text(data = totals, aes(tier, N, label = N), vjust = -0.4, inherit.aes = FALSE,
            size = 2.6, fontface = "bold") +
  scale_fill_manual(values = BIOTYPE_COLS, name = "biotype") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "library genes",
       title = sprintf("Cas13 library v8 composition (n = %d)", nrow(lib))) +
  theme_masld() + theme_pub() +
  theme(legend.position = "right", axis.text.x = element_text(size = 6.5, lineheight = 0.85))
save_fig(pC, file.path(OUT_DIR, "cas13_v8_tier_composition.pdf"),
         width = fig_half_width * 1.3, height = fig_half_width * 0.9)
message("wrote cas13_v8_tier_composition.pdf")

# ── Panel D: mouse-hepatocyte substrate composition (the v8 gate axis) ─────────
sub_lvls <- c("high", "ambient_suspect", "absent")
sub_cols <- c(high = "#2E9A86", ambient_suspect = "#9FBFB5", absent = "#B0413E")
pcd  <- lib[biotype == "protein_coding"]
dsub <- pcd[, .N, by = mouse_hep_substrate]
dsub[, mouse_hep_substrate := factor(mouse_hep_substrate, levels = sub_lvls)]
pD <- ggplot(dsub, aes(mouse_hep_substrate, N, fill = mouse_hep_substrate)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = N), vjust = -0.4, size = 2.8, fontface = "bold") +
  scale_fill_manual(values = sub_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(x = "mouse-hepatocyte substrate", y = "protein-coding targets",
       title = sprintf("v8 gate: mouse-hepatocyte substrate (PC, n = %d)", nrow(pcd)),
       subtitle = sprintf("%d hep-failing PC genes kept via positive-control / COLOC>=0.9 exemption (flagged mouse_untestable)",
                          sum(lib$mouse_untestable))) +
  theme_masld() + theme_pub() +
  theme(plot.subtitle = element_text(size = 6, colour = "grey35"),
        axis.text.x = element_text(size = 6.5))
save_fig(pD, file.path(OUT_DIR, "cas13_v8_mouse_hep_substrate.pdf"),
         width = fig_half_width * 1.25, height = fig_half_width * 0.9)
message("wrote cas13_v8_mouse_hep_substrate.pdf")

message("figS_cas13_library_core_vs_mash_overlap.R complete.")
