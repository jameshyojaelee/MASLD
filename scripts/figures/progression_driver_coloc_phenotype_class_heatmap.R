#!/usr/bin/env Rscript
# ============================================================================
# progression_driver_coloc_phenotype_class_heatmap.R
# Fig 2 (genetics) supplement — Fig S2N (relocated from Fig 3 supp 2026-07-08;
# S2M is the meSuSiE shared/ancestry-specific panel, so this appends as S2N).
# Output dir figS04_coloc. Progression-driver colocalization by GWAS phenotype class.
#
# Among stage-progression drivers with strong hepatic colocalization
# (coloc_best_pp4>0.5), the genetic anchor is almost exclusively a liver-enzyme
# endophenotype (ALT / AST / GGT), not diagnosed-patient onset or progression
# GWAS. Only HKDC1 colocalizes across all phenotype classes.
#
# Rows  = drivers with coloc_best_pp4>0.5 (deduped per gene, max).
# Cols  = GWAS phenotype classes {enzyme, onset, progression}.
#         steatosis dropped (max PP.H4 only ~0.31 — not a driver axis).
# Fill  = colocalization PP.H4.
#
# INTEGRITY: per-class PP.H4 columns come from progression_driver_genetics.csv
# (derived from gene_level_coloc.csv). Atlas coloc_susie_best_pp4 NEVER used.
# Assert all PP.H4 <= 1.
#
# Output: figures/main/fig2_genetics/panels/
#         FigS2N_progression_driver_coloc_phenotype_class_heatmap.pdf
#         (co-located + FigS2N_-prefixed to match the other FigS2 panels, 2026-07-07)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- file.path(FIG3_DIR, "panels")   # fig2_genetics/panels (co-located with FigS2A-M)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_PDF  <- file.path(OUT_DIR, "FigS2N_progression_driver_coloc_phenotype_class_heatmap.pdf")
DATA_CSV <- file.path(OUT_DIR, "FigS2N_progression_driver_coloc_phenotype_class_heatmap_source.csv")

# ---------------------------------------------------------------------------
# Load + dedup per gene (file is per-transition; per-gene PP.H4 cols are
# gene-level aggregates -> take max across rows)
# ---------------------------------------------------------------------------
d <- fread(file.path(BASE,
  "RNA-seq/results/stratified_causal/progression_driver_genetics.csv"))

pcols <- c("coloc_best_pp4", "pp4_best_enzyme", "pp4_best_onset",
           "pp4_best_progression", "pp4_best_steatosis")
g <- d[, lapply(.SD, function(x) suppressWarnings(max(x, na.rm = TRUE))),
       by = gene_symbol, .SDcols = pcols]
for (cc in pcols) g[is.infinite(get(cc)), (cc) := NA_real_]

stopifnot(max(unlist(g[, ..pcols]), na.rm = TRUE) <= 1)

# Drivers with strong colocalization
core <- g[coloc_best_pp4 > 0.5]
setorder(core, -coloc_best_pp4)

# ---------------------------------------------------------------------------
# Long format: enzyme / onset / progression (drop steatosis)
# ---------------------------------------------------------------------------
class_map <- c(pp4_best_enzyme = "Enzyme\n(ALT/AST/GGT)",
               pp4_best_onset = "Onset\n(diagnosis)",
               pp4_best_progression = "Progression")
long <- melt(core, id.vars = "gene_symbol",
             measure.vars = names(class_map),
             variable.name = "pheno", value.name = "pp4")
long[, pheno := factor(class_map[as.character(pheno)], levels = class_map)]
long[, gene_symbol := factor(gene_symbol, levels = rev(core$gene_symbol))]
long[is.na(pp4), pp4 := 0]

# ---------------------------------------------------------------------------
# Hero numbers
# ---------------------------------------------------------------------------
core[, n_classes_lit := (pp4_best_enzyme > 0.5 & !is.na(pp4_best_enzyme)) +
                        (pp4_best_onset > 0.5 & !is.na(pp4_best_onset)) +
                        (pp4_best_progression > 0.5 & !is.na(pp4_best_progression))]
multi <- core[n_classes_lit > 1, gene_symbol]
cat(sprintf("[hero] %d progression drivers with coloc_best_pp4>0.5\n", nrow(core)))
cat("[hero] genes:", paste(core$gene_symbol, collapse = ", "), "\n")
cat("[hero] light >1 phenotype class (PP.H4>0.5):",
    ifelse(length(multi), paste(multi, collapse = ", "), "none"), "\n")

fwrite(core, DATA_CSV)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
p <- ggplot(long, aes(x = pheno, y = gene_symbol, fill = pp4)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(pp4 > 0.5, sprintf("%.2f", pp4), "")),
            size = 1.8, color = "white", fontface = "plain") +
  scale_fill_gradient(low = "#F3E1EA", high = "#C9265E",
                      limits = c(0, 1), name = "PP.H4",
                      breaks = c(0, 0.5, 1)) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title = "Progression drivers anchor on liver-enzyme GWAS") +
  theme_masld(base_size = 7) +
  theme(
    plot.title      = element_text(size = 7.1, face = "plain", margin = margin(b = 5)),
    axis.text.x.top = element_text(size = 6, face = "plain"),
    axis.text.y     = element_text(size = 6, face = "italic"),
    axis.line       = element_blank(),
    axis.ticks      = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.3, "cm")
  )

ggsave(OUT_PDF, p,
       width  = 70 / 25.4,
       height = max(60, 6 * nrow(core) + 24) / 25.4,
       units  = "in",
       device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))
