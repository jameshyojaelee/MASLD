#!/usr/bin/env Rscript
# NOTE (Phase 10 swap, 2026-05-06): This script intentionally loads the FROZEN
# UKBB v1 sghatan snapshot as a sensitivity-comparison panel. v1 is no longer
# the production EUR LD reference — PolyFun replaced it (concordance r=0.986).
# v1 archive: archive/ld_panel_v1_sghatan_2026-05-06/
#
# presentation_ld_panel_comparison.R
# Presentation-grade "1KG ≈ UKBB" concordance figure.
# Larger fonts, fewer panels, clearer message than the supplementary version.
#
# Output: figures/presentation/ld_panel_concordance_presentation.pdf (+ per-panel PDFs)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(ggrastr)
})

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(PROJ, "scripts/figures/publication_theme.R"))

FM_DIR <- file.path(PROJ, "GWAS/finemapping")
OUT_DIR <- file.path(PROJ, "figures/supplementary/figS04_coloc/ld_panel_2way_ukbb_vs_1kg/coloc")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "../../ld_panel_3way_eur/coloc"),
           recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load panel outputs
# ---------------------------------------------------------------------------
load_panel <- function(subdir) {
  root <- file.path(FM_DIR, "results", subdir)
  files <- list.files(root, pattern = "susie_coloc_chr[0-9]+\\.csv$",
                      recursive = TRUE, full.names = TRUE)
  rows <- lapply(files, function(f) {
    dt <- try(fread(f), silent = TRUE)
    if (inherits(dt, "try-error") || nrow(dt) == 0) return(NULL)
    dt
  })
  dt <- rbindlist(rows, fill = TRUE)
  dt[, max_pp4 := pmax(PP.H4.abf, PP.H4.susie, na.rm = TRUE)]
  dt[is.infinite(max_pp4), max_pp4 := NA]
  dt[]
}

cat("Loading panels...\n")
v1    <- load_panel("../../archive/ld_panel_v1_sghatan_2026-05-06/susie_coloc_v1_2026-04-21")
kg    <- load_panel("susie_coloc_1kg")
topld <- load_panel("susie_coloc_topld")
cat("  UKBB v1:  ", nrow(v1),    " rows\n", sep = "")
cat("  1KG EUR:  ", nrow(kg),    " rows\n", sep = "")
cat("  TOP-LD :  ", if (is.null(topld)) "MISSING" else paste(nrow(topld), "rows"), "\n", sep = "")

# Only 17 EUR GWAS in EUR-only comparison (1kg runs are 17 EUR only)
eur_gwas <- unique(kg$gwas_name)
v1 <- v1[gwas_name %in% eur_gwas]

# Merge by (gwas_name, ensembl) — gene-level comparison
m <- merge(
  v1[, .(gwas_name, ensembl, gene, pp4_ukbb = max_pp4,
         susie_ukbb = PP.H4.susie, abf_ukbb = PP.H4.abf,
         method_ukbb = method, n_snps_ukbb = n_snps, top_snp_ukbb = top_snp)],
  kg[, .(gwas_name, ensembl, pp4_1kg = max_pp4,
         susie_1kg = PP.H4.susie, abf_1kg = PP.H4.abf,
         method_1kg = method, n_snps_1kg = n_snps, top_snp_1kg = top_snp)],
  by = c("gwas_name", "ensembl")
)
m <- m[!is.na(pp4_ukbb) & !is.na(pp4_1kg)]
cat("  paired gene×GWAS rows: ", nrow(m), "\n", sep = "")

pearson <- cor(m$pp4_ukbb, m$pp4_1kg)
spearman <- cor(m$pp4_ukbb, m$pp4_1kg, method = "spearman")
both_hi <- nrow(m[pp4_ukbb > 0.5 & pp4_1kg > 0.5])
ukbb_only <- nrow(m[pp4_ukbb > 0.5 & pp4_1kg <= 0.5])
kg_only   <- nrow(m[pp4_ukbb <= 0.5 & pp4_1kg > 0.5])
cat(sprintf("\n  r = %.3f (r² = %.3f), Spearman ρ = %.3f\n", pearson, pearson^2, spearman))
cat(sprintf("  PP4>0.5:  both = %d, UKBB-only = %d, 1KG-only = %d\n\n",
            both_hi, ukbb_only, kg_only))

# ---------------------------------------------------------------------------
# Panel A — main scatter (UKBB v1 on x, 1KG on y)
# ---------------------------------------------------------------------------
# Label only a short list of iconic MASLD genes
hallmark <- c("RORA", "THRB", "HKDC1", "MTTP", "GCKR", "CFLAR")
lbl <- m[gene %in% hallmark & pp4_ukbb > 0.3,
         .(gwas_name, gene, pp4_ukbb, pp4_1kg)][order(-pp4_ukbb)]
lbl <- unique(lbl, by = "gene")

p_a <- ggplot(m, aes(x = pp4_ukbb, y = pp4_1kg)) +
  geom_abline(slope = 1, intercept = 0, color = "#B0B0B0", linetype = "dashed", linewidth = 0.8) +
  rasterise(geom_point(alpha = 0.25, size = 0.9, color = "#1565C0"), dpi = 200) +
  geom_hline(yintercept = 0.5, color = "grey60", linetype = "dotted", linewidth = 0.5) +
  geom_vline(xintercept = 0.5, color = "grey60", linetype = "dotted", linewidth = 0.5) +
  geom_point(data = lbl, aes(x = pp4_ukbb, y = pp4_1kg),
             color = "black", size = 4) +
  geom_text_repel(data = lbl,
                  aes(x = pp4_ukbb, y = pp4_1kg, label = gene),
                  size = 8, color = "black", fontface = "bold",
                  box.padding = 0.8, max.overlaps = 20,
                  segment.color = "grey30", segment.size = 0.4) +
  annotate("text", x = 0.05, y = 0.92,
           label = sprintf("r = %.3f\nρ = %.3f\nn = %s",
                           pearson, spearman, format(nrow(m), big.mark = ",")),
           size = 8, hjust = 0, lineheight = 1.15, fontface = "plain") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    x = "PP.H4 — UKBB EUR",
    y = "PP.H4 — 1KG EUR",
    title = "Gene-level SuSiE-COLOC concordance"
  ) +
  theme_masld(base_size = 20) +
  theme(plot.title = element_text(size = 26, face = "bold"),
        axis.title = element_text(size = 22),
        axis.text  = element_text(size = 18))

# ---------------------------------------------------------------------------
# Panel B — PP4>0.5 set agreement (bar)
# ---------------------------------------------------------------------------
bar_dat <- data.frame(
  category = factor(c("Both panels agree", "UKBB only", "1KG only"),
                    levels = c("Both panels agree", "UKBB only", "1KG only")),
  n = c(both_hi, ukbb_only, kg_only)
)

p_b <- ggplot(bar_dat, aes(x = category, y = n, fill = category)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = n), vjust = -0.3, size = 9, fontface = "bold") +
  scale_fill_manual(values = c("#7B1FA2", "#F4511E", "#42A5F5"), guide = "none") +
  labs(
    x = NULL, y = "# gene × GWAS at PP.H4 > 0.5",
    title = sprintf("Top-hit reproducibility (%.1f%%)",
                    100 * both_hi / (both_hi + ukbb_only))
  ) +
  theme_masld(base_size = 20) +
  theme(plot.title = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        panel.grid.major.x = element_blank()) +
  expand_limits(y = max(bar_dat$n) * 1.15)

# ---------------------------------------------------------------------------
# Panel C — Hallmark-gene head-to-head table
# ---------------------------------------------------------------------------
# Select the top hallmark hits in both panels, show side-by-side
hl_tbl <- m[gene %in% hallmark & pmax(pp4_ukbb, pp4_1kg) > 0.3][order(-pp4_ukbb)]
hl_tbl <- unique(hl_tbl, by = "gene")
hl_tbl <- hl_tbl[1:min(10, nrow(hl_tbl)),
                 .(gene,
                   pp4_ukbb = sprintf("%.3f", pp4_ukbb),
                   pp4_1kg  = sprintf("%.3f", pp4_1kg),
                   Δ = sprintf("%+.4f", pp4_1kg - pp4_ukbb))]
hl_tbl <- hl_tbl[!is.na(gene)]

p_c <- ggplot(hl_tbl) +
  geom_text(aes(x = 0, y = rev(seq_len(nrow(hl_tbl))), label = gene),
            hjust = 0, size = 8, fontface = "bold") +
  geom_text(aes(x = 1.2, y = rev(seq_len(nrow(hl_tbl))), label = pp4_ukbb),
            hjust = 1, size = 7.5) +
  geom_text(aes(x = 2.4, y = rev(seq_len(nrow(hl_tbl))), label = pp4_1kg),
            hjust = 1, size = 7.5) +
  geom_text(aes(x = 3.3, y = rev(seq_len(nrow(hl_tbl))), label = Δ),
            hjust = 1, size = 7, color = "grey40") +
  annotate("text", x = c(0, 1.2, 2.4, 3.3), y = nrow(hl_tbl) + 0.8,
           label = c("Gene", "UKBB v1", "1KG EUR", "Δ"),
           hjust = c(0, 1, 1, 1), size = 8, fontface = "bold") +
  scale_x_continuous(limits = c(-0.1, 3.5)) +
  scale_y_continuous(limits = c(0.5, nrow(hl_tbl) + 1.3)) +
  labs(title = "Hallmark genes — identical PP.H4") +
  theme_void(base_size = 20) +
  theme(plot.title = element_text(size = 24, face = "bold"),
        plot.margin = margin(10, 10, 10, 10))

# ---------------------------------------------------------------------------
# Compose + save
# ---------------------------------------------------------------------------
composite <- (p_a | (p_b / p_c)) + plot_layout(widths = c(1.35, 1)) +
  plot_annotation(
    title = "1KG EUR reproduces sghatan UKBB fine-mapping",
    theme = theme(plot.title = element_text(size = 30, face = "bold"))
  )

# ---------------------------------------------------------------------------
# 3-way panel (UKBB / 1KG / TOP-LD) when TOP-LD has data
# ---------------------------------------------------------------------------
if (!is.null(topld) && nrow(topld) > 0) {
  cat("\n=== 3-way comparison (UKBB / 1KG / TOP-LD) ===\n")
  topld <- topld[gwas_name %in% eur_gwas]
  m3 <- merge(
    m,
    topld[, .(gwas_name, ensembl, pp4_topld = max_pp4)],
    by = c("gwas_name", "ensembl"), all.x = TRUE
  )
  m3 <- m3[!is.na(pp4_topld)]
  cat("  3-way paired rows: ", nrow(m3), "\n", sep = "")
  if (nrow(m3) > 0) {
    p_topld_ukbb <- cor(m3$pp4_ukbb, m3$pp4_topld)
    p_topld_1kg  <- cor(m3$pp4_1kg,  m3$pp4_topld)
    cat(sprintf("  UKBB↔TOP-LD r = %.3f\n", p_topld_ukbb))
    cat(sprintf("  1KG↔TOP-LD  r = %.3f\n", p_topld_1kg))

    # Long-form for facetted scatter
    long <- rbind(
      m3[, .(comparison = sprintf("UKBB vs 1KG    (r = %.3f)", pearson),
             x = pp4_ukbb,  y = pp4_1kg,   gene)],
      m3[, .(comparison = sprintf("UKBB vs TOP-LD (r = %.3f)", p_topld_ukbb),
             x = pp4_ukbb,  y = pp4_topld, gene)],
      m3[, .(comparison = sprintf("1KG  vs TOP-LD (r = %.3f)", p_topld_1kg),
             x = pp4_1kg,   y = pp4_topld, gene)]
    )
    lbl3 <- m3[gene %in% hallmark, .(gene, pp4_ukbb, pp4_1kg, pp4_topld)]
    lbl3 <- unique(lbl3, by = "gene")

    p_3way <- ggplot(long, aes(x = x, y = y)) +
      geom_abline(slope = 1, intercept = 0,
                  color = "#B0B0B0", linetype = "dashed", linewidth = 0.6) +
      rasterise(geom_point(alpha = 0.25, size = 0.7, color = "#1565C0"), dpi = 200) +
      geom_hline(yintercept = 0.5, color = "grey60", linetype = "dotted", linewidth = 0.4) +
      geom_vline(xintercept = 0.5, color = "grey60", linetype = "dotted", linewidth = 0.4) +
      facet_wrap(~ comparison, ncol = 3) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(x = "PP.H4 (panel on x)", y = "PP.H4 (panel on y)",
           title = "3-way LD-panel concordance (gene-level SuSiE-COLOC)") +
      theme_masld(base_size = 18) +
      theme(plot.title  = element_text(size = 24, face = "bold"),
            strip.text  = element_text(size = 16, face = "bold"),
            axis.title  = element_text(size = 18),
            axis.text   = element_text(size = 14))

    ggsave(file.path(OUT_DIR, "../../ld_panel_3way_eur/coloc/3way_facet_scatter_simple.pdf"), p_3way,
           width = 18, height = 7, device = cairo_pdf)
    cat("  Wrote: ld_panel_3way_concordance.pdf\n")
  }
}

out_composite <- file.path(OUT_DIR, "concordance_presentation.pdf")
ggsave(out_composite, composite, width = 17, height = 9, device = cairo_pdf)
cat("Wrote:", out_composite, "\n")

# Save each panel separately for slide-deck flexibility
ggsave(file.path(OUT_DIR, "concordance_scatter.pdf"), p_a,
       width = 9, height = 9, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "concordance_bar.pdf"), p_b,
       width = 6.5, height = 5, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "concordance_hallmarks.pdf"), p_c,
       width = 6.5, height = 5, device = cairo_pdf)

# Emit a compact summary text
sink(file.path(OUT_DIR, "concordance_summary.txt"))
cat("LD Panel Concordance — UKBB v1 (sghatan, N≈337K) vs 1KG EUR (ours, N=379)\n")
cat("Generated: ", as.character(Sys.time()), "\n\n", sep = "")
cat("Total paired gene × GWAS comparisons: ", format(nrow(m), big.mark = ","), "\n", sep = "")
cat("Pearson r:      ", sprintf("%.4f", pearson), "\n", sep = "")
cat("Pearson r²:     ", sprintf("%.4f", pearson^2), "\n", sep = "")
cat("Spearman ρ:     ", sprintf("%.4f", spearman), "\n\n", sep = "")
cat("At PP.H4 > 0.5:\n")
cat("  Both panels agree:  ", both_hi, "\n", sep = "")
cat("  UKBB-only (lost in 1KG):  ", ukbb_only, "\n", sep = "")
cat("  1KG-only (new in 1KG):    ", kg_only, "\n", sep = "")
cat("  Reproducibility: ", sprintf("%.1f%%", 100*both_hi/(both_hi+ukbb_only)), "\n\n", sep = "")
cat("Top hallmark genes (UKBB v1 vs 1KG):\n")
print(hl_tbl)
sink()

cat("\nAll presentation artifacts written to: ", OUT_DIR, "\n", sep = "")
