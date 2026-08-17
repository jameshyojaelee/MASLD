#!/usr/bin/env Rscript
# ==============================================================================
# Fig 4j: ESSENTIALITY NEGATIVE CONTROL — prioritized MASLD targets avoid a
#   generic liver-cell viability liability in DepMap.
#
# KEY MESSAGE:
#   A target that is also essential in liver cancer cell lines would raise a
#   generic cell-viability concern for inhibition. All six case-study targets
#   sit at Chronos ~ 0 (non-essential), across a wide COLOC-PP.H4 range.
#
# AXES (compact quadrant scatter with a binned reference distribution):
#   x = genetic colocalization PP.H4 per gene, read from the CANONICAL per-gene
#       file gene_level_coloc.csv (best PP.H4: prefer SuSiE column, fall back to
#       abf). NOT the method-inconsistent atlas *_coloc_pp4 convenience columns.
#   y = DepMap CRISPR essentiality (Chronos). Convention: MORE NEGATIVE = MORE
#       essential (~ -1 pan-essential; ~0 non-essential). From the atlas
#       essentiality_chronos column, joined on human_symbol.
#
# Essential zone = the atlas is_essential boundary, which is exactly
#   Chronos < -0.5. The display is focused on -1.5 to 0.25; out-of-range counts
#   are disclosed on-panel and the full values remain in the source data.
#
# HOUSE STYLE: all text BLACK; no on-panel title/subtitle (all prose -> caption
#   to stdout); compact/dense; PDF only (cairo_pdf). 6 case-study genes are the
#   only labels (fixed 6 pt, italic black); points show canonical DEG status.
#
# Output: figures/main/fig5_molecular_context/fig4j_essentiality_control.pdf
# Env:    rnaseq
# Run:    ~/micromamba/envs/rnaseq/bin/Rscript scripts/figures/essentiality_control.R
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
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

# ── 2. Essentiality (Chronos) + canonical bulk DEG status from the atlas ──────
atl <- fread(ATLAS_CSV,
  # bulk_padj added 2026-08-12 for the canonical gate.
  select = c("human_symbol", "essentiality_chronos", "is_essential",
             "bulk_logFC", "bulk_padj", "bulk_treat_fdr"))
atl <- atl[!is.na(essentiality_chronos) & !is.na(human_symbol) & human_symbol != ""]

# ── 3. Join: genes with BOTH a COLOC PP.H4 and a Chronos score ────────────────
d <- merge(cg, atl, by.x = "gene", by.y = "human_symbol")
d <- unique(d, by = "gene")
d[, bulk_status := fifelse(is_canonical_deg(d) & bulk_logFC > 0, "DEG up",
                    fifelse(is_canonical_deg(d) & bulk_logFC < 0, "DEG down",
                            "Not DEG"))]
d[, bulk_status := factor(bulk_status,
                           levels = c("DEG up", "DEG down", "Not DEG"))]
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
               is_essential, bulk_status)])

# ── PLOT ──────────────────────────────────────────────────────────────────────
COL_UP   <- masld_colors[["up"]]     # Liang deep magenta (MASLD-up)
COL_DOWN <- masld_colors[["down"]]   # deep blue (MASLD-down)
COL_NS   <- "#E3E3E3"                # non-DEG case-study targets
COL_BG   <- "#8F8F8F"                # binned reference distribution
COL_REF  <- "#8A8A8A"                # cutoff guides / annotations

Y_MIN <- -1.50
Y_MAX <-  0.25
PP4_CUTOFF <- 0.50

# Bin only the visible reference values. Out-of-range values are counted below
# rather than silently squished into the boundary bins.
d_bg <- d[!(gene %in% HERO_GENES) &
            essentiality_chronos >= Y_MIN & essentiality_chronos <= Y_MAX]
n_clip_low  <- sum(d$essentiality_chronos < Y_MIN)
n_clip_high <- sum(d$essentiality_chronos > Y_MAX)

# Fixed label positions prevent the nearly coincident RORA/THRB/HKDC1 labels
# from changing between renders. Segments fan out from the true coordinates.
label_pos <- data.table(
  gene = HERO_GENES,
  label_x = c(0.94, 0.73, 0.29, 0.84, 0.11, 0.53),
  label_y = c(0.18, 0.18, -0.20, -0.20, 0.18, 0.18),
  line_x  = c(0.96, 0.82, 0.25, 0.90, 0.13, 0.55),
  line_y  = c(0.13, 0.14, -0.15, -0.15, 0.13, 0.13)
)
hero <- merge(hero, label_pos, by = "gene", all.x = TRUE, sort = FALSE)
hero[, gene := factor(gene, levels = HERO_GENES)]
setorder(hero, gene)
# RORA and THRB are effectively coincident (both PP.H4 ~1.00, Chronos ~0.05).
# Shift RORA slightly left for display and connect it to its true coordinate.
hero[, plot_x := coloc_pp4]
hero[gene == "RORA", plot_x := 0.975]

p <- ggplot() +
  # Light zones make the two independent decision boundaries explicit.
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = Y_MIN, ymax = ESS_CUTOFF,
           fill = "#F4F4F4", colour = NA) +
  annotate("rect", xmin = PP4_CUTOFF, xmax = Inf, ymin = ESS_CUTOFF, ymax = Y_MAX,
           fill = "#F2F7F4", colour = NA) +
  # All reference genes are summarized as bins to avoid a 13k-point cloud.
  geom_bin_2d(data = d_bg,
              aes(x = coloc_pp4, y = essentiality_chronos,
                  fill = after_stat(count)),
              binwidth = c(0.02, 0.04), colour = NA) +
  scale_fill_gradient(low = "#F0F0F0", high = COL_BG, guide = "none") +
  geom_hline(yintercept = ESS_CUTOFF, linetype = "22",
             linewidth = 0.35, colour = COL_REF) +
  geom_vline(xintercept = PP4_CUTOFF, linetype = "22",
             linewidth = 0.35, colour = COL_REF) +
  # Case-study targets — canonical DEG class, black outline.
  geom_segment(data = hero[gene == "RORA"],
               aes(x = coloc_pp4, y = essentiality_chronos,
                   xend = plot_x, yend = essentiality_chronos),
               linewidth = 0.25, linetype = "22", colour = "grey50") +
  geom_point(data = hero, aes(x = plot_x, y = essentiality_chronos),
             shape = 21, size = 2.8, fill = "black", colour = "black", stroke = 0) +
  geom_point(data = hero,
             aes(x = plot_x, y = essentiality_chronos, colour = bulk_status),
             shape = 16, size = 2.15) +
  geom_segment(data = hero,
               aes(x = plot_x, y = essentiality_chronos,
                   xend = line_x, yend = line_y),
               linewidth = 0.25, colour = "grey50") +
  geom_text(data = hero,
            aes(x = label_x, y = label_y, label = gene),
            size = GEOM_TEXT_6PT, fontface = "italic", colour = "black") +
  annotate("text", x = 0.015, y = ESS_CUTOFF - 0.04,
           label = "Essential", hjust = 0, vjust = 1,
           size = GEOM_TEXT_6PT, colour = COL_REF) +
  annotate("text", x = PP4_CUTOFF + 0.015, y = Y_MIN + 0.04,
           label = "Coloc-supported", hjust = 0, vjust = 0,
           size = GEOM_TEXT_6PT, colour = COL_REF) +
  annotate("text", x = 0.015, y = Y_MIN + 0.04,
           label = sprintf("%s genes < %.1f; %s > %.2f",
                           format(n_clip_low, big.mark = ","), Y_MIN,
                           format(n_clip_high, big.mark = ","), Y_MAX),
           hjust = 0, vjust = 0, size = GEOM_TEXT_6PT, colour = COL_REF) +
  scale_colour_manual(values = c("DEG up" = COL_UP, "DEG down" = COL_DOWN,
                                 "Not DEG" = COL_NS),
                      name = "Bulk status",
                      breaks = c("DEG up", "DEG down", "Not DEG")) +
  scale_x_continuous(name = "Coloc probability (PP.H4)",
                     breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.02, 0.03))) +
  scale_y_continuous(name = "Liver-cell essentiality (Chronos)",
                     breaks = c(-1.5, -1.0, -0.5, 0),
                     expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(-0.02, 1.03), ylim = c(Y_MIN, Y_MAX),
                  expand = FALSE) +
  theme_masld_compact() +
  theme(legend.position = c(0.985, 0.025), legend.justification = c(1, 0),
        legend.key.size = unit(0.30, "cm"),
        legend.background = element_blank(), legend.key = element_blank(),
        plot.margin = margin(3, 4, 3, 3))

# fig4j_essentiality_control.pdf RETIRED 2026-07-07 (do NOT re-create):
# the essentiality negative-control panel was dropped from Fig 4. The PDF is no
# longer written. Caption/provenance below is kept for reference.
# OUT_PDF <- file.path(FIG4_DIR, "panels", "fig4j_essentiality_control.pdf")
# dir.create(FIG4_DIR, recursive = TRUE, showWarnings = FALSE)
# if (capabilities("cairo")) {
#   ggsave(OUT_PDF, p, width = 2.8, height = 2.05, device = cairo_pdf)
# } else {
#   ggsave(OUT_PDF, p, width = 2.8, height = 2.05,
#          device = grDevices::pdf, useDingbats = FALSE)
# }

# ── Caption / provenance -> stdout (ALL prose lives here, never on-panel) ──────
n_ess <- sum(d$essentiality_chronos < ESS_CUTOFF)
hero_txt <- paste(sprintf("%s (Chronos %+.2f, PP.H4 %.2f)",
                          as.character(hero$gene), hero$essentiality_chronos,
                          hero$coloc_pp4), collapse = "; ")
message(strrep("=", 78))
message("FIG 4j — ESSENTIALITY NEGATIVE CONTROL (prioritized targets are non-essential)")
message(strrep("=", 78))
message(sprintf(
"Prioritized MASLD targets do not show a generic liver-cell viability liability.
The background distribution summarizes genes with both a genetic colocalization
PP.H4 (x; canonical per-gene
gene_level_coloc.csv, SuSiE preferred / abf fallback) and a DepMap CRISPR Chronos
essentiality score averaged across liver cancer cell lines (y; atlas
essentiality_chronos, joined on human_symbol;
n=%d genes with both). Chronos convention: more negative = more essential
(~ -1 pan-essential, ~0 non-essential). Dashed guides mark the canonical COLOC
threshold (PP.H4 >= %.1f) and essentiality boundary (Chronos < %.1f; %d/%d genes
fall below the latter). The displayed y-range is %.2f to %.2f; %d lower and %d
upper values are omitted and counted on-panel. All six labelled case-study targets
sit at Chronos ~ 0 (non-essential) across a wide COLOC range: %s. Point colour =
canonical bulk DEG status (FDR < 0.05 at lfc=0.25; magenta = DEG up, blue = DEG
down, grey = not DEG). RORA is displaced slightly left, with a dashed connector to
its true coordinate, to separate it from the nearly coincident THRB point. Lack of
liver-cell essentiality reduces concern that target
inhibition would cause nonspecific loss of viability; it does not establish
organism-level safety. No on-panel title; all text black.",
  nrow(d), PP4_CUTOFF, ESS_CUTOFF, n_ess, nrow(d), Y_MIN, Y_MAX,
  n_clip_low, n_clip_high, hero_txt))
message(strrep("=", 78))
# fig4j PDF write RETIRED 2026-07-07 — OUT_PDF no longer defined; nothing written.
cat("[done] fig4j_essentiality_control panel RETIRED — no PDF written.\n")
