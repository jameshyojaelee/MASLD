#!/usr/bin/env Rscript
# figS_cas13_library_core_vs_mash_overlap.R
# ---------------------------------------------------------------------------
# v6 Cas13 library: how much of the MASH-vs-MASL (NASH>NAFL) signature is
# additive over the core Disease-vs-Control spine, and what the 4-tier library
# is made of. Both DEG sources are limma-voom + metafor (REML), ashr-shrunk (06b).
#
# Individual panels (PDF only, assembled in Illustrator):
#   A  cas13_v6_core_vs_mash_euler.pdf   area-proportional 2-set Euler of the
#        human-DEG overlap: Core DvC UP vs MASH-vs-MASL UP (lfsr<0.05 & shrunk>0.2).
#   B  cas13_v6_evidence_upset.pdf       UpSet over library mouse genes by evidence
#        axis: core spine / MASH-vs-MASL / mouse cross-diet(>=3) / miRNA-conserved.
#   C  cas13_v6_tier_composition.pdf     stacked bar: tier x biotype gene counts.
#
# Output dir: Cas13_Library_Design/figures/ (FIGS_CAS13LIB_DIR).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggforce)
  library(UpSetR)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- FIGS_CAS13LIB_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results")
ASHR <- file.path(INT, "integration/meta_results_ashr.csv")
MASH <- file.path(INT, "disease_signatures/nafl_vs_nash_meta_ashr.csv")
LIB  <- file.path(BASE, "Cas13_Library_Design/data/cas13_library_v3.0.csv")

LFSR_THR  <- 0.05
SHRUNK_THR <- 0.20    # matches HUMAN_SHRUNK_THR in rebuild_cas13_library.R
strip_v <- function(x) sub("[.][0-9]+$", "", x)

# tier colors (palette1; distinct categorical, control-free)
COL_CORE  <- "#e35070"   # magenta  -- disease-vs-control core
COL_MASH  <- "#b771e6"   # purple   -- MASH-vs-MASL progression
COL_MOUSE <- "#4baeef"   # blue     -- mouse cross-diet
COL_MIR   <- "#30d796"   # green    -- miRNA conserved
BIOTYPE_COLS <- c(protein_coding = "#5480C2", lncRNA = "#D4834A", miRNA = "#30d796")

# ── Human-DEG UP sets (the actual DEG overlap, pre-ortholog) ───────────────────
core <- fread(ASHR, select = c("gene", "shrunk_logFC", "lfsr"))
mash <- fread(MASH, select = c("gene", "shrunk_logFC", "lfsr"))
core_up <- unique(strip_v(core[!is.na(lfsr) & lfsr < LFSR_THR & shrunk_logFC > SHRUNK_THR, gene]))
mash_up <- unique(strip_v(mash[!is.na(lfsr) & lfsr < LFSR_THR & shrunk_logFC > SHRUNK_THR, gene]))
nC <- length(core_up); nM <- length(mash_up)
nO <- length(intersect(core_up, mash_up))
nM_new <- length(setdiff(mash_up, core_up))
message(sprintf("Core UP=%d  MASH UP=%d  overlap=%d  MASH-only(net-new)=%d (Jaccard=%.3f)",
                nC, nM, nO, nM_new, nO / length(union(core_up, mash_up))))

# ── Panel A: area-proportional 2-set Euler (Core left, MASH right) ─────────────
lens_area <- function(r1, r2, d) {
  if (d >= r1 + r2) return(0)
  if (d <= abs(r1 - r2)) return(pi * min(r1, r2)^2)
  a1 <- acos(pmin(1, pmax(-1, (d^2 + r1^2 - r2^2) / (2 * d * r1))))
  a2 <- acos(pmin(1, pmax(-1, (d^2 + r2^2 - r1^2) / (2 * d * r2))))
  r1^2 * (a1 - sin(2 * a1) / 2) + r2^2 * (a2 - sin(2 * a2) / 2)
}
solve_d <- function(r1, r2, target) {
  if (target <= 0)               return((r1 + r2) * 1.04)
  if (target >= pi * min(r1, r2)^2) return(abs(r1 - r2))
  uniroot(function(d) lens_area(r1, r2, d) - target,
          lower = abs(r1 - r2) + 1e-6, upper = r1 + r2 - 1e-6)$root
}
rC <- sqrt(nC / pi); rM <- sqrt(nM / pi)
d  <- solve_d(rC, rM, nO)
cxC <- -d / 2; cxM <- d / 2
circles <- data.frame(
  x0 = c(cxC, cxM), y0 = c(0, 0), r = c(rC, rM),
  set = c("Core (Disease vs Control)", "MASH vs MASL"))
labs <- data.frame(
  x = c(-(rC + rM) / 2, (rC - rM) / 2, (rC + rM) / 2), y = 0,
  lab = c(nC - nO, nO, nM - nO))
pA <- ggplot() +
  geom_circle(data = circles, aes(x0 = x0, y0 = y0, r = r, fill = set),
              colour = NA, alpha = 0.55) +
  geom_text(data = labs, aes(x, y, label = lab), size = 3.0, fontface = "bold") +
  annotate("text", x = cxC, y = rC + 1.5, label = sprintf("Core\n%d", nC),
           colour = COL_CORE, size = 2.6, fontface = "bold", lineheight = 0.9) +
  annotate("text", x = cxM, y = rM + 1.5, label = sprintf("MASH vs MASL\n%d", nM),
           colour = "#7a3fb0", size = 2.6, fontface = "bold", lineheight = 0.9) +
  scale_fill_manual(values = c("Core (Disease vs Control)" = COL_CORE,
                               "MASH vs MASL" = COL_MASH), guide = "none") +
  coord_fixed(clip = "off") +
  labs(title = "Human DEG overlap (limma-voom + metafor, UP, lfsr<0.05 & shrunk>0.2)",
       subtitle = sprintf("%d MASH-vs-MASL genes are net-new vs the core spine", nM_new)) +
  theme_void(base_size = 7) +
  theme(plot.title = element_text(size = 7, hjust = 0.5),
        plot.subtitle = element_text(size = 6.5, hjust = 0.5, colour = "grey30"),
        plot.margin = margin(14, 6, 6, 6))
save_fig(pA, file.path(OUT_DIR, "cas13_v6_core_vs_mash_euler.pdf"),
         width = fig_half_width, height = fig_half_width * 0.95)
message("wrote cas13_v6_core_vs_mash_euler.pdf")

# ── Library evidence axes (mouse-gene level) ──────────────────────────────────
lib <- fread(LIB)
sets <- list(
  `Core spine`          = lib[has_human_de == TRUE, gene_id_mouse],
  `MASH vs MASL`        = lib[is_mash_deg == TRUE, gene_id_mouse],
  `Mouse cross-diet>=3` = lib[n_diets_up >= 3, gene_id_mouse],
  `miRNA conserved`     = lib[biotype == "miRNA", gene_id_mouse])
sets <- lapply(sets, unique)

# ── Panel B: UpSet over the four evidence axes ────────────────────────────────
pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
pdf_device(file.path(OUT_DIR, "cas13_v6_evidence_upset.pdf"),
           width = fig_full_width, height = 3.6)
print(upset(fromList(sets), order.by = "freq", nsets = 4, nintersects = 20,
            sets.bar.color = c(COL_CORE, COL_MASH, COL_MOUSE, COL_MIR),
            main.bar.color = "#37474F", matrix.color = "#37474F",
            mainbar.y.label = "library genes", sets.x.label = "tier size",
            text.scale = 1.1))
dev.off()
message("wrote cas13_v6_evidence_upset.pdf")

# ── Panel C: tier x biotype composition ───────────────────────────────────────
tier_lvls <- c("human", "mash_progression", "mouse_confirmed", "mirna_conserved")
comp <- lib[, .N, by = .(tier, biotype)]
comp[, tier := factor(tier, levels = tier_lvls,
                      labels = c("Core\n(Disease vs Control)", "MASH\nprogression",
                                 "Mouse\nconfirmed", "miRNA\nconserved"))]
comp[, biotype := factor(biotype, levels = c("protein_coding", "lncRNA", "miRNA"))]
totals <- comp[, .(N = sum(N)), by = tier]
pC <- ggplot(comp, aes(tier, N, fill = biotype)) +
  geom_col(width = 0.7) +
  geom_text(data = totals, aes(tier, N, label = N), vjust = -0.4,
            inherit.aes = FALSE, size = 2.6, fontface = "bold") +
  scale_fill_manual(values = BIOTYPE_COLS, name = "biotype") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "library genes",
       title = sprintf("Cas13 library v6 composition (n = %d)", nrow(lib))) +
  theme_masld() + theme_pub() +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 6.5, lineheight = 0.85))
save_fig(pC, file.path(OUT_DIR, "cas13_v6_tier_composition.pdf"),
         width = fig_half_width * 1.25, height = fig_half_width * 0.9)
message("wrote cas13_v6_tier_composition.pdf")

message("figS_cas13_library_core_vs_mash_overlap.R complete.")
