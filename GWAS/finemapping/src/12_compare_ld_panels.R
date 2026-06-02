#!/usr/bin/env Rscript
# 12_compare_ld_panels.R
# Compare SuSiE-COLOC results across LD panels (UKBB v1 snapshot / 1KG EUR / TOP-LD EUR).
# Output: concordance metrics + supplementary figure.
#
# Inputs:
#   - results/susie_coloc_ukbb_sghatan_v1_2026-04-21/<gwas>/susie_coloc_chr<N>.csv (frozen v1)
#   - results/susie_coloc_1kg/<gwas>/susie_coloc_chr<N>.csv (Phase 4 production output)
#   - results/susie_coloc_topld/<gwas>/susie_coloc_chr<N>.csv (Phase 3 partial output, after TOP-LD build)
# Outputs:
#   - results/ld_panel_comparison/ld_panel_concordance.csv
#   - figures/supplementary/figS04_coloc/figS_ld_panel_comparison.pdf
#
# Usage:
#   Rscript 12_compare_ld_panels.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

source(file.path(PROJ, "scripts/figures/load_figure_data.R"))
source(file.path(PROJ, "scripts/figures/publication_theme.R"))

FM_DIR <- file.path(PROJ, "GWAS/finemapping")
OUT_CSV <- file.path(FM_DIR, "results/ld_panel_comparison/ld_panel_concordance.csv")
OUT_PDF <- file.path(FIGS04_DIR, "figS_ld_panel_comparison.pdf")
dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
load_panel <- function(suffix_dir) {
  root <- file.path(FM_DIR, "results", suffix_dir)
  if (!dir.exists(root)) return(NULL)
  files <- list.files(root, pattern = "susie_coloc_chr[0-9]+\\.csv$",
                      recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) return(NULL)
  rows <- lapply(files, function(f) {
    dt <- try(fread(f), silent = TRUE)
    if (inherits(dt, "try-error") || nrow(dt) == 0) return(NULL)
    dt
  })
  dt <- rbindlist(rows, fill = TRUE, use.names = TRUE)
  dt[]
}

join_gwas_gene <- function(dt, panel_label) {
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  dt[, panel := panel_label]
  dt[, max_pp4 := pmax(PP.H4.abf, PP.H4.susie, na.rm = TRUE)]
  dt[, .(panel, gwas_name, gene, ensembl, chr, n_snps, method,
         PP.H4.abf, PP.H4.susie, max_pp4, top_snp)]
}

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
cat("Loading panels...\n")
p_ukbb  <- load_panel("susie_coloc_ukbb_sghatan_v1_2026-04-21")
p_1kg   <- load_panel("susie_coloc_1kg")
p_topld <- load_panel("susie_coloc_topld")

cat("  UKBB v1:  ", if (is.null(p_ukbb))  "MISSING" else paste(nrow(p_ukbb),  "rows"), "\n")
cat("  1KG:      ", if (is.null(p_1kg))   "MISSING" else paste(nrow(p_1kg),   "rows"), "\n")
cat("  TOP-LD:   ", if (is.null(p_topld)) "MISSING" else paste(nrow(p_topld), "rows"), "\n")

tidy <- rbindlist(list(
  join_gwas_gene(p_ukbb,  "UKBB_v1"),
  join_gwas_gene(p_1kg,   "1KG"),
  join_gwas_gene(p_topld, "TOP-LD")
), fill = TRUE)

if (nrow(tidy) == 0) {
  cat("\nNo panel data found — nothing to compare. Exiting.\n")
  quit(status = 0)
}

# Wide: one row per (gwas, gene), columns per panel
wide <- dcast(tidy, gwas_name + gene + ensembl + chr ~ panel,
              value.var = c("max_pp4", "PP.H4.susie", "PP.H4.abf", "n_snps"),
              fun.aggregate = function(x) x[1])

fwrite(wide, OUT_CSV)
cat("\nWrote concordance table: ", OUT_CSV, " (", nrow(wide), " rows)\n", sep = "")

# ---------------------------------------------------------------------------
# Summary stats
# ---------------------------------------------------------------------------
panels_present <- unique(tidy$panel)
cat("\nPanels with data:", paste(panels_present, collapse = ", "), "\n")

if ("UKBB_v1" %in% panels_present && "1KG" %in% panels_present) {
  common <- wide[!is.na(max_pp4_UKBB_v1) & !is.na(max_pp4_1KG)]
  if (nrow(common) > 0) {
    r2 <- cor(common$max_pp4_UKBB_v1, common$max_pp4_1KG, use = "pairwise.complete.obs")^2
    rho <- cor(common$max_pp4_UKBB_v1, common$max_pp4_1KG, method = "spearman", use = "pairwise.complete.obs")
    cat(sprintf("\nUKBB↔1KG (n=%d): r² = %.3f, Spearman ρ = %.3f\n", nrow(common), r2, rho))
    cat(sprintf("  Both PP4>0.5: %d  | UKBB-only: %d  | 1KG-only: %d\n",
                sum(common$max_pp4_UKBB_v1 > 0.5 & common$max_pp4_1KG > 0.5),
                sum(common$max_pp4_UKBB_v1 > 0.5 & !(common$max_pp4_1KG > 0.5)),
                sum(!(common$max_pp4_UKBB_v1 > 0.5) & common$max_pp4_1KG > 0.5)))
  }
}

if ("UKBB_v1" %in% panels_present && "TOP-LD" %in% panels_present) {
  common <- wide[!is.na(max_pp4_UKBB_v1) & !is.na(`max_pp4_TOP-LD`)]
  if (nrow(common) > 0) {
    r2 <- cor(common$max_pp4_UKBB_v1, common$`max_pp4_TOP-LD`, use = "pairwise.complete.obs")^2
    rho <- cor(common$max_pp4_UKBB_v1, common$`max_pp4_TOP-LD`, method = "spearman", use = "pairwise.complete.obs")
    cat(sprintf("\nUKBB↔TOP-LD (n=%d): r² = %.3f, Spearman ρ = %.3f\n", nrow(common), r2, rho))
  }
}

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
plots <- list()

# Panel (a) UKBB vs 1KG scatter
if ("UKBB_v1" %in% panels_present && "1KG" %in% panels_present) {
  d <- wide[!is.na(max_pp4_UKBB_v1) & !is.na(max_pp4_1KG)]
  plots$a <- ggplot(d, aes(x = max_pp4_UKBB_v1, y = max_pp4_1KG)) +
    geom_point(alpha = 0.4, size = 0.5, color = "#2171B5") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0.5, linetype = "dotted", color = "grey60") +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey60") +
    labs(x = "Max PP.H4 (UKBB v1 snapshot)", y = "Max PP.H4 (1KG EUR)",
         title = "(a) UKBB v1 vs 1KG EUR",
         subtitle = sprintf("n = %d gene×GWAS comparisons", nrow(d))) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_masld(base_size = 8)
}

# Panel (b) UKBB vs TOP-LD scatter
if ("UKBB_v1" %in% panels_present && "TOP-LD" %in% panels_present) {
  d <- wide[!is.na(max_pp4_UKBB_v1) & !is.na(`max_pp4_TOP-LD`)]
  plots$b <- ggplot(d, aes(x = max_pp4_UKBB_v1, y = `max_pp4_TOP-LD`)) +
    geom_point(alpha = 0.4, size = 0.5, color = "#D62728") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0.5, linetype = "dotted", color = "grey60") +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey60") +
    labs(x = "Max PP.H4 (UKBB v1)", y = "Max PP.H4 (TOP-LD EUR)",
         title = "(b) UKBB v1 vs TOP-LD EUR") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_masld(base_size = 8)
}

# Panel (c) n_snps per panel (boxplot)
plots$c <- ggplot(tidy, aes(x = panel, y = n_snps, fill = panel)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.15, outlier.size = 0.4) +
  scale_y_log10() +
  labs(x = NULL, y = "n overlapping SNPs per test",
       title = "(c) LD-ref overlap with GWAS sumstats") +
  theme_masld(base_size = 8) +
  theme(legend.position = "none")

# Panel (d) UpSet-style bar chart of PP4>0.5 genes by panel
hi_by_panel <- tidy[max_pp4 > 0.5, .N, by = panel]
plots$d <- ggplot(hi_by_panel, aes(x = panel, y = N, fill = panel)) +
  geom_col() +
  geom_text(aes(label = N), vjust = -0.3, size = 3) +
  labs(x = NULL, y = "# gene×GWAS with PP.H4 > 0.5",
       title = "(d) PP.H4 > 0.5 hits per panel") +
  theme_masld(base_size = 8) +
  theme(legend.position = "none")

composite <- wrap_plots(plots, ncol = 2) +
  plot_annotation(
    title = "Figure S | LD panel sensitivity for SuSiE-COLOC (UKBB v1 / 1KG EUR / TOP-LD EUR)",
    theme = theme(plot.title = element_text(size = 11, face = "bold"))
  )

dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)
ggsave(OUT_PDF, composite, width = 10, height = 8, device = cairo_pdf)
cat("\nWrote figure: ", OUT_PDF, "\n", sep = "")
cat("\nDone.\n")
