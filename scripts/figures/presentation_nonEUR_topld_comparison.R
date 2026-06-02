#!/usr/bin/env Rscript
# presentation_nonEUR_topld_comparison.R
# Non-EUR (EAS/AFR/SAS) 2-way comparison: 1KG-per-ancestry vs TOP-LD-per-ancestry.
# No UKBB v1 column for non-EUR (sghatan is EUR-only).
#
# Inputs:
#   1KG baseline:  GWAS/finemapping/results/susie_coloc/<gwas>/susie_coloc_chr*.csv
#                  (uses 1kg_eas / 1kg_afr / 1kg_sas via LD_PANEL=1kg dispatch)
#   TOP-LD:        GWAS/finemapping/results/susie_coloc_topld/<gwas>/susie_coloc_chr*.csv
#                  (Phase 7b output)
#
# GWAS lists per ancestry:
#   EAS: 5 GWAS (BBJ_ALT, BBJ_AST, BBJ_GGT, 2020_32514122_Cirrhosis_EAS, 2020_32514122_HCC_EAS)
#   AFR: 3 GWAS (PanUKBB_AFR_ALT/AST/GGT)
#   SAS: 3 GWAS (PanUKBB_CSA_ALT/AST/GGT)
#
# Outputs in figures/supplementary/figS09_multi_ancestry/ld_panel_nonEUR/:
#   eas/coloc/{scatter,hits_bar}.pdf
#   afr/coloc/{scatter,hits_bar}.pdf
#   sas/coloc/{scatter,hits_bar}.pdf
#   all_ancestries/coloc_concordance_per_ancestry.pdf
#   summary.txt

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
OUT_DIR <- file.path(PROJ, "figures/supplementary/figS09_multi_ancestry/ld_panel_nonEUR")
for (anc in c("eas", "afr", "sas", "all_ancestries"))
  dir.create(file.path(OUT_DIR, anc, "coloc"), recursive = TRUE, showWarnings = FALSE)

ANCESTRY_GWAS <- list(
  EAS = c("BBJ_ALT", "BBJ_AST", "BBJ_GGT",
          "2020_32514122_Cirrhosis_EAS", "2020_32514122_HCC_EAS"),
  AFR = c("PanUKBB_AFR_ALT", "PanUKBB_AFR_AST", "PanUKBB_AFR_GGT"),
  SAS = c("PanUKBB_CSA_ALT", "PanUKBB_CSA_AST", "PanUKBB_CSA_GGT")
)
ANC_COLORS <- c(EAS = "#F4511E", AFR = "#00695C", SAS = "#7B1FA2")
PANEL_COLORS <- c("1KG" = "#7B1FA2", "TOP-LD" = "#C2185B")

hallmark <- c("RORA","THRB","HKDC1","MTTP","GCKR","CFLAR","EPHA2",
              "TNFSF10","SLC39A8","ADH4","DGAT2","LPL","PNPLA3","TM6SF2",
              "APOE","HSD17B13","MARC1","TRIB1","HFE")

# ----------------------------------------------------------------------------
# Load per-GWAS COLOC for one panel
# ----------------------------------------------------------------------------
load_panel <- function(subdir, gwas_filter) {
  root <- file.path(FM_DIR, "results", subdir)
  if (!dir.exists(root)) return(NULL)
  gwas_dirs <- intersect(gwas_filter, list.files(root))
  if (length(gwas_dirs) == 0) return(NULL)
  rows <- lapply(gwas_dirs, function(g) {
    files <- list.files(file.path(root, g),
                        pattern = "susie_coloc_chr[0-9]+\\.csv$",
                        full.names = TRUE)
    if (length(files) == 0) return(NULL)
    rbindlist(lapply(files, function(f) {
      dt <- try(fread(f), silent = TRUE)
      if (inherits(dt, "try-error") || nrow(dt) == 0) return(NULL)
      dt
    }), fill = TRUE)
  })
  dt <- rbindlist(rows, fill = TRUE)
  dt[, max_pp4 := pmax(PP.H4.abf, PP.H4.susie, na.rm = TRUE)]
  dt[is.infinite(max_pp4), max_pp4 := NA]
  dt[]
}

per_anc_results <- list()
all_summary_lines <- c(
  paste0("Non-EUR 1KG vs TOP-LD comparison (", Sys.time(), ")"),
  paste0("EAS GWAS: ", paste(ANCESTRY_GWAS$EAS, collapse = ", ")),
  paste0("AFR GWAS: ", paste(ANCESTRY_GWAS$AFR, collapse = ", ")),
  paste0("SAS GWAS: ", paste(ANCESTRY_GWAS$SAS, collapse = ", ")),
  ""
)

for (anc in names(ANCESTRY_GWAS)) {
  cat(sprintf("\n=== %s comparison ===\n", anc))

  kg <- load_panel("susie_coloc",       ANCESTRY_GWAS[[anc]])
  tl <- load_panel("susie_coloc_topld", ANCESTRY_GWAS[[anc]])

  if (is.null(kg) || is.null(tl)) {
    cat(sprintf("  [WARN] %s: missing one panel (1KG=%s, TOP-LD=%s); skipping\n",
                anc, !is.null(kg), !is.null(tl)))
    next
  }
  cat(sprintf("  1KG rows: %d   TOP-LD rows: %d\n", nrow(kg), nrow(tl)))

  m <- merge(
    kg[, .(gwas_name, ensembl, gene, pp4_kg = max_pp4)],
    tl[, .(gwas_name, ensembl, pp4_tl = max_pp4)],
    by = c("gwas_name", "ensembl"), all = FALSE
  )
  m <- m[!is.na(pp4_kg) & !is.na(pp4_tl)]
  cat(sprintf("  Paired rows: %d\n", nrow(m)))

  r   <- cor(m$pp4_kg, m$pp4_tl)
  rho <- cor(m$pp4_kg, m$pp4_tl, method = "spearman")
  hi_kg <- sum(m$pp4_kg > 0.5)
  hi_tl <- sum(m$pp4_tl > 0.5)
  hi_both <- sum(m$pp4_kg > 0.5 & m$pp4_tl > 0.5)

  cat(sprintf("  r = %.3f   rho = %.3f\n", r, rho))
  cat(sprintf("  PP.H4 > 0.5: 1KG=%d  TOP-LD=%d  both=%d\n", hi_kg, hi_tl, hi_both))

  # Per-ancestry scatter
  lbl_pts <- m[gene %in% hallmark & pmax(pp4_kg, pp4_tl) > 0.5]
  lbl_pts <- unique(lbl_pts, by = "gene")

  p_scatter <- ggplot(m, aes(x = pp4_kg, y = pp4_tl)) +
    geom_abline(slope = 1, intercept = 0, color = "#B0B0B0",
                linetype = "dashed", linewidth = 0.8) +
    rasterise(geom_point(alpha = 0.30, size = 0.9, color = ANC_COLORS[anc]), dpi = 220) +
    geom_hline(yintercept = 0.5, color = "grey60", linetype = "dotted") +
    geom_vline(xintercept = 0.5, color = "grey60", linetype = "dotted") +
    geom_point(data = lbl_pts, color = "black", size = 2.6) +
    geom_text_repel(data = lbl_pts, aes(label = gene), size = 4.5,
                    color = "black", fontface = "bold", max.overlaps = 25) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "PP.H4 (1KG-EUR per-ancestry)",
         y = "PP.H4 (TOP-LD per-ancestry)",
         title = sprintf("%s: 1KG ↔ TOP-LD COLOC concordance\n%d genes paired, r = %.3f, ρ = %.3f",
                         anc, nrow(m), r, rho)) +
    theme_masld(base_size = 16) +
    theme(plot.title = element_text(size = 20, face = "bold"),
          axis.title = element_text(size = 18))
  ggsave(file.path(OUT_DIR, tolower(anc), "coloc/scatter.pdf"),
         p_scatter, width = 9, height = 9, device = cairo_pdf)
  cat(sprintf("  Wrote: %s/coloc/scatter.pdf\n", tolower(anc)))

  # Per-ancestry hits bar
  hits <- data.table(
    panel = factor(rep(c("1KG", "TOP-LD"), each = 3), levels = c("1KG", "TOP-LD")),
    threshold = factor(rep(c("PP.H4 > 0.5", "PP.H4 > 0.8", "PP.H4 > 0.9"), 2),
                       levels = c("PP.H4 > 0.5", "PP.H4 > 0.8", "PP.H4 > 0.9")),
    n = c(sum(m$pp4_kg > 0.5), sum(m$pp4_kg > 0.8), sum(m$pp4_kg > 0.9),
          sum(m$pp4_tl > 0.5), sum(m$pp4_tl > 0.8), sum(m$pp4_tl > 0.9))
  )
  p_hits <- ggplot(hits, aes(x = panel, y = n, fill = panel)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = n), vjust = -0.35, size = 5, fontface = "bold") +
    facet_wrap(~ threshold, ncol = 3, scales = "free_y") +
    scale_fill_manual(values = PANEL_COLORS, guide = "none") +
    labs(x = NULL, y = "# gene × GWAS",
         title = sprintf("%s: COLOC hit counts per panel", anc)) +
    theme_masld(base_size = 16) +
    theme(plot.title = element_text(size = 20, face = "bold"),
          strip.text = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size = 13)) +
    expand_limits(y = max(hits$n) * 1.20)
  ggsave(file.path(OUT_DIR, tolower(anc), "coloc/hits_bar.pdf"),
         p_hits, width = 12, height = 5.5, device = cairo_pdf)

  per_anc_results[[anc]] <- m
  per_anc_results[[anc]][, ancestry := anc]

  all_summary_lines <- c(all_summary_lines,
    sprintf("=== %s ===", anc),
    sprintf("  Paired rows: %d", nrow(m)),
    sprintf("  Pearson r:   %.3f", r),
    sprintf("  Spearman rho:%.3f", rho),
    sprintf("  PP.H4 > 0.5: 1KG=%d  TOP-LD=%d  both=%d", hi_kg, hi_tl, hi_both),
    sprintf("  PP.H4 > 0.8: 1KG=%d  TOP-LD=%d", sum(m$pp4_kg > 0.8), sum(m$pp4_tl > 0.8)),
    ""
  )
}

# ----------------------------------------------------------------------------
# Combined per-ancestry facet figure
# ----------------------------------------------------------------------------
if (length(per_anc_results) > 0) {
  combined <- rbindlist(per_anc_results, fill = TRUE)
  combined[, ancestry := factor(ancestry, levels = c("EAS", "AFR", "SAS"))]

  ann <- combined[, .(r = cor(pp4_kg, pp4_tl),
                       rho = cor(pp4_kg, pp4_tl, method = "spearman"),
                       n = .N), by = ancestry]
  ann[, label := sprintf("r = %.3f\nρ = %.3f\nn = %d", r, rho, n)]

  lbl_combined <- combined[gene %in% hallmark & pmax(pp4_kg, pp4_tl) > 0.5]
  lbl_combined <- unique(lbl_combined, by = c("gene", "ancestry"))

  p_facet <- ggplot(combined, aes(x = pp4_kg, y = pp4_tl)) +
    geom_abline(slope = 1, intercept = 0, color = "#B0B0B0",
                linetype = "dashed", linewidth = 0.8) +
    rasterise(geom_point(alpha = 0.30, size = 0.7, aes(color = ancestry)), dpi = 220) +
    geom_hline(yintercept = 0.5, color = "grey60", linetype = "dotted", linewidth = 0.4) +
    geom_vline(xintercept = 0.5, color = "grey60", linetype = "dotted", linewidth = 0.4) +
    geom_point(data = lbl_combined, color = "black", size = 2.0) +
    geom_text_repel(data = lbl_combined, aes(label = gene), size = 3.6,
                    color = "black", fontface = "bold", max.overlaps = 20) +
    geom_text(data = ann, aes(label = label), x = 0.05, y = 0.92,
              hjust = 0, vjust = 1, size = 5, fontface = "bold", color = "grey20") +
    facet_wrap(~ ancestry, ncol = 3) +
    scale_color_manual(values = ANC_COLORS, guide = "none") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "PP.H4 (1KG-per-ancestry)",
         y = "PP.H4 (TOP-LD-per-ancestry)",
         title = "Non-EUR COLOC concordance: 1KG ↔ TOP-LD by ancestry") +
    theme_masld(base_size = 16) +
    theme(plot.title = element_text(size = 22, face = "bold"),
          strip.text = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 18))
  ggsave(file.path(OUT_DIR, "all_ancestries/coloc_concordance_per_ancestry.pdf"),
         p_facet, width = 18, height = 7, device = cairo_pdf)
  cat("\nWrote: all_ancestries/coloc_concordance_per_ancestry.pdf\n")

  writeLines(all_summary_lines, file.path(OUT_DIR, "all_ancestries/summary.txt"))
  cat("Wrote: all_ancestries/summary.txt\n")
}

cat("\nALL DONE. Outputs in:", OUT_DIR, "\n")
