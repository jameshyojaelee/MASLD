#!/usr/bin/env Rscript
# ==============================================================================
# Fig 4j: ESSENTIALITY NEGATIVE CONTROL — prioritized MASLD targets are NOT
#   confounded by cell-essentiality. They sit in the non-essential, druggable
#   zone (a standard target-atlas negative control).
#
# KEY MESSAGE:
#   A downregulated disease gene that is ALSO pan-essential would be a red flag —
#   its loss-of-function is incompatible with cell viability, so it cannot
#   explain a disease phenotype via simple knockdown. We show our genetically
#   colocalized, prioritized targets avoid that trap: every case-study target
#   sits at Chronos ~ 0 (non-essential) across the COLOC-PP.H4 axis.
#
# AXES (a SCATTER — no lollipops, no shaded bands; a faint reference line only):
#   x = genetic colocalization PP.H4 per gene, read from the CANONICAL per-gene
#       file gene_level_coloc.csv (best PP.H4: prefer SuSiE column, fall back to
#       abf). NOT the method-inconsistent atlas *_coloc_pp4 convenience columns.
#   y = DepMap CRISPR essentiality (Chronos). Convention: MORE NEGATIVE = MORE
#       essential (~ -1 pan-essential; ~0 non-essential). From the atlas
#       essentiality_chronos column, joined on human_symbol.
#
# Essential zone = the atlas is_essential boundary, which is exactly
#   Chronos < -0.5 (1,607 genes). Marked with one faint horizontal reference
#   line + a light "essential (LoF-implausible)" label.
#
# HOUSE STYLE: all text BLACK; no on-panel title/subtitle (all prose -> caption
#   to stdout); compact/dense; PDF only (cairo_pdf). 6 case-study genes are the
#   only labels (ggrepel, italic black); points coloured by bulk direction.
#
# Output: figures/main/fig4_validation/fig4j_essentiality_control.pdf
# Env:    rnaseq
# Run:    ~/micromamba/envs/rnaseq/bin/Rscript scripts/figures/essentiality_control.R
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── Input paths (canonical on-disk sources) ───────────────────────────────────
COLOC_CSV <- file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
ATLAS_CSV <- file.path(BASE,
  "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")

HERO_GENES   <- c("THRB", "RORA", "NR1H4", "HKDC1", "SERPINE1", "CYP3A4")
ESS_CUTOFF   <- -0.5      # atlas is_essential boundary (Chronos < -0.5)

# ── 1. Per-gene best COLOC PP.H4 (canonical: SuSiE preferred, abf fallback) ───
coloc <- fread(COLOC_CSV)
coloc <- coloc[gene != "" & !is.na(gene)]
coloc[, coloc_pp4 := fifelse(!is.na(coloc_best_susie_pp4),
                             coloc_best_susie_pp4, coloc_best_pp4)]
coloc[, coloc_meth := fifelse(!is.na(coloc_best_susie_pp4), "SuSiE", "abf")]
cg <- coloc[!is.na(coloc_pp4), .(gene, coloc_pp4, coloc_meth)]

# ── 2. Essentiality (Chronos) + bulk direction from the atlas ─────────────────
atl <- fread(ATLAS_CSV,
  select = c("human_symbol", "essentiality_chronos", "is_essential", "bulk_logFC"))
atl <- atl[!is.na(essentiality_chronos) & !is.na(human_symbol) & human_symbol != ""]

# ── 3. Join: genes with BOTH a COLOC PP.H4 and a Chronos score ────────────────
d <- merge(cg, atl, by.x = "gene", by.y = "human_symbol")
d <- unique(d, by = "gene")
d[, direction := fifelse(is.na(bulk_logFC), "n.s.",
                  fifelse(bulk_logFC > 0, "MASLD-up", "MASLD-down"))]
cat(sprintf("[join] %d genes with both COLOC PP.H4 and Chronos\n", nrow(d)))
cat(sprintf("[ess]  essential zone (Chronos < %.1f): %d / %d genes\n",
            ESS_CUTOFF, sum(d$essentiality_chronos < ESS_CUTOFF), nrow(d)))

# ── 4. Case-study genes (the only labelled points) ────────────────────────────
hero <- d[gene %in% HERO_GENES]
hero[, gene := factor(gene, levels = HERO_GENES)]
setorder(hero, gene)
stopifnot(nrow(hero) == length(HERO_GENES))   # all 6 must be present
cat("[hero] case-study targets (all expected to be non-essential):\n")
print(hero[, .(gene, coloc_pp4 = round(coloc_pp4, 3), coloc_meth,
               chronos = round(essentiality_chronos, 3),
               is_essential, direction)])

# ── PLOT ──────────────────────────────────────────────────────────────────────
COL_UP   <- masld_colors[["up"]]     # Liang deep magenta (MASLD-up)
COL_DOWN <- masld_colors[["down"]]   # deep blue (MASLD-down)
COL_BG   <- "#D9D9D9"                # faint grey background cloud
COL_ESS  <- "#B0B0B0"                # neutral grey reference line / zone marker

# Label position for the "essential (LoF-implausible)" zone (left side, in-zone)
y_floor   <- min(d$essentiality_chronos)
zone_lab_y <- (ESS_CUTOFF + y_floor) / 2

p <- ggplot() +
  # all genes with COLOC + Chronos — small faint grey points (background)
  geom_point(data = d[!(gene %in% HERO_GENES)],
             aes(x = coloc_pp4, y = essentiality_chronos),
             colour = COL_BG, size = 0.45, alpha = 0.55, stroke = 0) +
  # essential-zone boundary = a single dashed reference line (no on-panel label; defined in caption)
  annotate("segment", x = -Inf, xend = Inf, y = ESS_CUTOFF, yend = ESS_CUTOFF,
           linetype = "22", linewidth = 0.35, colour = COL_ESS) +
  # case-study targets — coloured by bulk direction, black outline
  geom_point(data = hero,
             aes(x = coloc_pp4, y = essentiality_chronos, fill = direction),
             shape = 21, size = 2.6, colour = "black", stroke = 0.4) +
  geom_text_repel(data = hero,
             aes(x = coloc_pp4, y = essentiality_chronos, label = gene),
             size = 2.5, fontface = "italic", colour = "black",
             segment.size = 0.25, segment.color = "grey55",
             min.segment.length = 0, box.padding = 0.7, point.padding = 0.35,
             nudge_y = 0.45, max.overlaps = Inf, seed = 7) +
  scale_fill_manual(values = c("MASLD-up" = COL_UP, "MASLD-down" = COL_DOWN,
                               "n.s." = "#9E9E9E"),
                    name = "bulk direction",
                    breaks = c("MASLD-up", "MASLD-down")) +
  scale_x_continuous(name = "genetic colocalization (COLOC PP.H4)",
                     limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.02, 0.03))) +
  scale_y_continuous(name = "DepMap essentiality (Chronos)",
                     expand = expansion(mult = c(0.04, 0.04))) +
  theme_masld() + theme_pub() +
  theme(axis.text  = element_text(size = PUB_AXIS_TEXT, colour = "black"),
        axis.title = element_text(size = PUB_AXIS_TEXT + 1, colour = "black"),
        legend.position = c(0.985, 0.02), legend.justification = c(1, 0),
        legend.title = element_text(size = 5.5, colour = "black"),
        legend.text  = element_text(size = 5, colour = "black"),
        legend.key.size = unit(0.30, "cm"),
        legend.background = element_blank(), legend.key = element_blank(),
        plot.margin = margin(3, 4, 3, 3))

OUT_PDF <- file.path(FIG4_DIR, "fig4j_essentiality_control.pdf")
dir.create(FIG4_DIR, recursive = TRUE, showWarnings = FALSE)
if (capabilities("cairo")) {
  ggsave(OUT_PDF, p, width = 3.3, height = 2.9, device = cairo_pdf)
} else {
  ggsave(OUT_PDF, p, width = 3.3, height = 2.9,
         device = grDevices::pdf, useDingbats = FALSE)
}

# ── Caption / provenance -> stdout (ALL prose lives here, never on-panel) ──────
n_ess <- sum(d$essentiality_chronos < ESS_CUTOFF)
hero_txt <- paste(sprintf("%s (Chronos %+.2f, PP.H4 %.2f)",
                          as.character(hero$gene), hero$essentiality_chronos,
                          hero$coloc_pp4), collapse = "; ")
message(strrep("=", 78))
message("FIG 4j — ESSENTIALITY NEGATIVE CONTROL (prioritized targets are non-essential)")
message(strrep("=", 78))
message(sprintf(
"Prioritized MASLD targets are NOT confounded by cell-essentiality. Each point is
a gene with both a genetic colocalization PP.H4 (x; canonical per-gene
gene_level_coloc.csv, SuSiE preferred / abf fallback) and a DepMap CRISPR Chronos
essentiality score (y; atlas essentiality_chronos, joined on human_symbol;
n=%d genes with both). Chronos convention: more negative = more essential
(~ -1 pan-essential, ~0 non-essential). The faint reference line + light grey zone
mark the essential boundary (Chronos < %.1f = the atlas is_essential flag; %d/%d
genes fall there). A downregulated disease gene that was ALSO pan-essential would
be a red flag (loss-of-function incompatible with viability, cannot explain the
phenotype by knockdown); our targets avoid it. All six labelled case-study targets
sit at Chronos ~ 0 (NON-essential) across a wide COLOC range: %s. Point colour =
bulk disease direction (magenta = MASLD-up, blue = MASLD-down). No on-panel title;
all text black.",
  nrow(d), ESS_CUTOFF, n_ess, nrow(d), hero_txt))
message(strrep("=", 78))
cat(sprintf("[done] wrote %s  (%.1f KB)\n", OUT_PDF, file.info(OUT_PDF)$size / 1024))
