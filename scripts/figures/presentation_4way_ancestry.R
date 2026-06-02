#!/usr/bin/env Rscript
# presentation_4way_ancestry.R
# Compare 2-ancestry (EUR+EAS) vs 4-ancestry (EUR+EAS+AFR+SAS) cross-ancestry
# fine-mapping. Both runs use 1KG LD across all arms.
#
# Inputs:
#   results/susiex_1kg/susiex_gene_summary_1kg.csv          — 2-way (EUR+EAS)
#   results/susiex_4way/susiex_gene_summary_4way.csv        — 4-way (EUR+EAS+AFR+SAS)
#   results/mesusie_1kg/mesusie_gene_summary.csv            — 2-way MESuSiE
#   results/mesusie_4way/mesusie_gene_summary_4way.csv      — 4-way MESuSiE
#
# Outputs (figures/presentation/4way_ancestry_expansion/):
#   susiex/{pip_scatter.pdf,hits_bar.pdf,top_snp_agreement.pdf,summary.txt}
#   mesusie/{pip_scatter.pdf,hits_bar.pdf,summary.txt}

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel)
})

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(PROJ, "scripts/figures/publication_theme.R"))

FM_DIR <- file.path(PROJ, "GWAS/finemapping")
OUT    <- file.path(PROJ, "figures/presentation/4way_ancestry_expansion")
for (sub in c("susiex", "mesusie")) dir.create(file.path(OUT, sub), recursive = TRUE, showWarnings = FALSE)

PANEL_COLS <- c("2-way (EUR+EAS)" = "#377EB8", "4-way (EUR+EAS+AFR+SAS)" = "#E41A1C")

hallmark <- c("RORA","THRB","HKDC1","MTTP","GCKR","CFLAR","EPHA2","TNFSF10",
              "SLC39A8","ADH4","DGAT2","LPL","PNPLA3","TM6SF2","APOE","HSD17B13",
              "MARC1","TRIB1","HFE","SERPINA1")

write_lines_safe <- function(lines, path) {
  con <- file(path, open = "w"); writeLines(lines, con); close(con)
}

# ============================================================================
# SuSiEX 2-way vs 4-way
# ============================================================================
cat("=== SuSiEX 2-way vs 4-way ===\n")

sx2_path <- file.path(FM_DIR, "results/susiex_1kg/susiex_gene_summary_1kg.csv")
sx4_path <- file.path(FM_DIR, "results/susiex_4way/susiex_gene_summary_4way.csv")

if (!file.exists(sx2_path) || !file.exists(sx4_path)) {
  cat("MISSING:\n  ", sx2_path, " exists=", file.exists(sx2_path),
      "\n  ", sx4_path, " exists=", file.exists(sx4_path), "\n")
} else {
  sx2 <- fread(sx2_path)[, .(ENSG, GeneSymbol, trait_pairs,
                              pip_2way = susiex_max_pip,
                              cs_2way  = susiex_cs_size_joint,
                              n_cs_2way = susiex_n_cs)]
  sx4 <- fread(sx4_path)[, .(ENSG, GeneSymbol, trait_pairs,
                              pip_4way = susiex_max_pip,
                              cs_4way  = susiex_cs_size_joint,
                              n_cs_4way = susiex_n_cs)]
  m_sx <- merge(sx2, sx4, by = c("ENSG","GeneSymbol","trait_pairs"), all = FALSE)
  cat(sprintf("  Paired genes: %d\n", nrow(m_sx)))

  # PIP scatter (gene-level)
  pearson <- if (nrow(m_sx) >= 3) cor(m_sx$pip_2way, m_sx$pip_4way, use = "complete.obs") else NA
  spearman <- if (nrow(m_sx) >= 3) cor(m_sx$pip_2way, m_sx$pip_4way, method = "spearman", use = "complete.obs") else NA
  m_sx[, hl := GeneSymbol %in% hallmark]
  p_sx <- ggplot(m_sx, aes(pip_2way, pip_4way)) +
    geom_abline(slope = 1, linetype = "dashed", color = "grey60") +
    geom_point(alpha = 0.5, size = 1.4) +
    geom_point(data = m_sx[hl == TRUE], color = "#E41A1C", size = 2.2) +
    geom_text_repel(data = m_sx[hl == TRUE & (pip_2way > 0.3 | pip_4way > 0.3)],
                    aes(label = GeneSymbol), size = 3, max.overlaps = 20) +
    labs(x = "Max PIP — 2-way (EUR+EAS)",
         y = "Max PIP — 4-way (EUR+EAS+AFR+SAS)",
         title = sprintf("SuSiEX gene-level max PIP (n=%d genes)", nrow(m_sx)),
         subtitle = sprintf("Pearson r = %.3f  ·  Spearman ρ = %.3f", pearson, spearman)) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_pub()
  ggsave(file.path(OUT, "susiex/pip_scatter.pdf"), p_sx, width = 6.5, height = 6.5)

  # Hits bar at PIP > 0.5 / 0.8 / 0.9
  hits <- data.table(
    panel = c(rep("2-way (EUR+EAS)", 3), rep("4-way (EUR+EAS+AFR+SAS)", 3)),
    threshold = rep(c("> 0.5", "> 0.8", "> 0.9"), 2),
    n_genes = c(sum(m_sx$pip_2way > 0.5, na.rm = TRUE),
                sum(m_sx$pip_2way > 0.8, na.rm = TRUE),
                sum(m_sx$pip_2way > 0.9, na.rm = TRUE),
                sum(m_sx$pip_4way > 0.5, na.rm = TRUE),
                sum(m_sx$pip_4way > 0.8, na.rm = TRUE),
                sum(m_sx$pip_4way > 0.9, na.rm = TRUE))
  )
  p_hits <- ggplot(hits, aes(threshold, n_genes, fill = panel)) +
    geom_col(position = position_dodge(0.7), width = 0.6) +
    geom_text(aes(label = n_genes), position = position_dodge(0.7), vjust = -0.4, size = 3.2) +
    scale_fill_manual(values = PANEL_COLS, name = NULL) +
    labs(x = "Max PIP threshold", y = "Genes",
         title = "SuSiEX hit count by PIP threshold") +
    theme_pub() + theme(legend.position = "top")
  ggsave(file.path(OUT, "susiex/hits_bar.pdf"), p_hits, width = 6.5, height = 4.5)

  # CS-size comparison (smaller = tighter, harder fine-mapping is more confident)
  cs_dt <- m_sx[!is.na(cs_2way) & !is.na(cs_4way)]
  if (nrow(cs_dt) > 0) {
    p_cs <- ggplot(cs_dt, aes(cs_2way, cs_4way)) +
      geom_abline(slope = 1, linetype = "dashed", color = "grey60") +
      geom_point(alpha = 0.5, size = 1.4) +
      scale_x_log10() + scale_y_log10() +
      labs(x = "Joint CS size — 2-way", y = "Joint CS size — 4-way",
           title = sprintf("SuSiEX joint credible-set size (n=%d genes)", nrow(cs_dt))) +
      coord_equal() + theme_pub()
    ggsave(file.path(OUT, "susiex/cs_size_scatter.pdf"), p_cs, width = 6.5, height = 6.5)
  }

  # Summary text
  summary_lines <- c(
    "SuSiEX 2-way (EUR+EAS) vs 4-way (EUR+EAS+AFR+SAS) — gene-level max PIP",
    sprintf("Loci anchor:                 140 EUR×EAS shared loci × 3 GWAS pairs (ALT/AST/GGT)"),
    sprintf("Paired genes (in both runs): %d", nrow(m_sx)),
    sprintf("Pearson r:                   %.3f", pearson),
    sprintf("Spearman rho:                %.3f", spearman),
    "",
    "Hit counts (max PIP):",
    sprintf("  > 0.5  : 2-way=%d  4-way=%d  delta=%+d",
            hits$n_genes[hits$panel=="2-way (EUR+EAS)" & hits$threshold=="> 0.5"],
            hits$n_genes[hits$panel=="4-way (EUR+EAS+AFR+SAS)" & hits$threshold=="> 0.5"],
            hits$n_genes[hits$panel=="4-way (EUR+EAS+AFR+SAS)" & hits$threshold=="> 0.5"] -
            hits$n_genes[hits$panel=="2-way (EUR+EAS)" & hits$threshold=="> 0.5"]),
    sprintf("  > 0.8  : 2-way=%d  4-way=%d  delta=%+d",
            hits$n_genes[hits$panel=="2-way (EUR+EAS)" & hits$threshold=="> 0.8"],
            hits$n_genes[hits$panel=="4-way (EUR+EAS+AFR+SAS)" & hits$threshold=="> 0.8"],
            hits$n_genes[hits$panel=="4-way (EUR+EAS+AFR+SAS)" & hits$threshold=="> 0.8"] -
            hits$n_genes[hits$panel=="2-way (EUR+EAS)" & hits$threshold=="> 0.8"]),
    sprintf("  > 0.9  : 2-way=%d  4-way=%d  delta=%+d",
            hits$n_genes[hits$panel=="2-way (EUR+EAS)" & hits$threshold=="> 0.9"],
            hits$n_genes[hits$panel=="4-way (EUR+EAS+AFR+SAS)" & hits$threshold=="> 0.9"],
            hits$n_genes[hits$panel=="4-way (EUR+EAS+AFR+SAS)" & hits$threshold=="> 0.9"] -
            hits$n_genes[hits$panel=="2-way (EUR+EAS)" & hits$threshold=="> 0.9"])
  )
  write_lines_safe(summary_lines, file.path(OUT, "susiex/summary.txt"))
  cat("  Wrote SuSiEX figures + summary to ", file.path(OUT, "susiex"), "\n", sep = "")
}

# ============================================================================
# MESuSiE 2-way vs 4-way
# ============================================================================
cat("\n=== MESuSiE 2-way vs 4-way ===\n")

ms2_path <- file.path(FM_DIR, "results/mesusie_1kg/mesusie_gene_summary.csv")
ms4_path <- file.path(FM_DIR, "results/mesusie_4way/mesusie_gene_summary_4way.csv")

if (!file.exists(ms2_path) || !file.exists(ms4_path)) {
  cat("MISSING:\n  ", ms2_path, " exists=", file.exists(ms2_path),
      "\n  ", ms4_path, " exists=", file.exists(ms4_path), "\n")
} else {
  ms2 <- fread(ms2_path)[, .(ENSG, GeneSymbol, trait_pairs, pip_2way = mesusie_max_pip)]
  ms4 <- fread(ms4_path)[, .(ENSG, GeneSymbol, trait_pairs, pip_4way = mesusie_max_pip)]
  m_ms <- merge(ms2, ms4, by = c("ENSG","GeneSymbol","trait_pairs"), all = FALSE)
  cat(sprintf("  Paired genes: %d\n", nrow(m_ms)))

  pearson  <- if (nrow(m_ms) >= 3) cor(m_ms$pip_2way, m_ms$pip_4way, use = "complete.obs") else NA
  spearman <- if (nrow(m_ms) >= 3) cor(m_ms$pip_2way, m_ms$pip_4way, method = "spearman", use = "complete.obs") else NA

  m_ms[, hl := GeneSymbol %in% hallmark]
  p_ms <- ggplot(m_ms, aes(pip_2way, pip_4way)) +
    geom_abline(slope = 1, linetype = "dashed", color = "grey60") +
    geom_point(alpha = 0.5, size = 1.4) +
    geom_point(data = m_ms[hl == TRUE], color = "#E41A1C", size = 2.2) +
    geom_text_repel(data = m_ms[hl == TRUE & (pip_2way > 0.3 | pip_4way > 0.3)],
                    aes(label = GeneSymbol), size = 3, max.overlaps = 20) +
    labs(x = "Max PIP — 2-way (EUR+EAS)",
         y = "Max PIP — 4-way (EUR+EAS+AFR+SAS)",
         title = sprintf("MESuSiE gene-level max PIP (n=%d genes)", nrow(m_ms)),
         subtitle = sprintf("Pearson r = %.3f  ·  Spearman ρ = %.3f", pearson, spearman)) +
    coord_equal(xlim = c(0,1), ylim = c(0,1)) + theme_pub()
  ggsave(file.path(OUT, "mesusie/pip_scatter.pdf"), p_ms, width = 6.5, height = 6.5)

  hits <- data.table(
    panel = c(rep("2-way (EUR+EAS)", 3), rep("4-way (EUR+EAS+AFR+SAS)", 3)),
    threshold = rep(c("> 0.5", "> 0.8", "> 0.9"), 2),
    n_genes = c(sum(m_ms$pip_2way > 0.5, na.rm = TRUE),
                sum(m_ms$pip_2way > 0.8, na.rm = TRUE),
                sum(m_ms$pip_2way > 0.9, na.rm = TRUE),
                sum(m_ms$pip_4way > 0.5, na.rm = TRUE),
                sum(m_ms$pip_4way > 0.8, na.rm = TRUE),
                sum(m_ms$pip_4way > 0.9, na.rm = TRUE))
  )
  p_hits <- ggplot(hits, aes(threshold, n_genes, fill = panel)) +
    geom_col(position = position_dodge(0.7), width = 0.6) +
    geom_text(aes(label = n_genes), position = position_dodge(0.7), vjust = -0.4, size = 3.2) +
    scale_fill_manual(values = PANEL_COLS, name = NULL) +
    labs(x = "Max PIP threshold", y = "Genes",
         title = "MESuSiE hit count by PIP threshold") +
    theme_pub() + theme(legend.position = "top")
  ggsave(file.path(OUT, "mesusie/hits_bar.pdf"), p_hits, width = 6.5, height = 4.5)

  summary_lines <- c(
    "MESuSiE 2-way (EUR+EAS) vs 4-way (EUR+EAS+AFR+SAS) — gene-level max PIP",
    sprintf("Loci anchor:                 140 EUR×EAS shared loci × 3 GWAS pairs (ALT/AST/GGT)"),
    sprintf("Paired genes (in both runs): %d", nrow(m_ms)),
    sprintf("Pearson r:                   %.3f", pearson),
    sprintf("Spearman rho:                %.3f", spearman),
    "",
    "Hit counts (max PIP):",
    sprintf("  > 0.5  : 2-way=%d  4-way=%d  delta=%+d",
            hits$n_genes[hits$panel=="2-way (EUR+EAS)" & hits$threshold=="> 0.5"],
            hits$n_genes[hits$panel=="4-way (EUR+EAS+AFR+SAS)" & hits$threshold=="> 0.5"],
            hits$n_genes[hits$panel=="4-way (EUR+EAS+AFR+SAS)" & hits$threshold=="> 0.5"] -
            hits$n_genes[hits$panel=="2-way (EUR+EAS)" & hits$threshold=="> 0.5"]),
    sprintf("  > 0.8  : 2-way=%d  4-way=%d  delta=%+d",
            hits$n_genes[hits$panel=="2-way (EUR+EAS)" & hits$threshold=="> 0.8"],
            hits$n_genes[hits$panel=="4-way (EUR+EAS+AFR+SAS)" & hits$threshold=="> 0.8"],
            hits$n_genes[hits$panel=="4-way (EUR+EAS+AFR+SAS)" & hits$threshold=="> 0.8"] -
            hits$n_genes[hits$panel=="2-way (EUR+EAS)" & hits$threshold=="> 0.8"]),
    sprintf("  > 0.9  : 2-way=%d  4-way=%d  delta=%+d",
            hits$n_genes[hits$panel=="2-way (EUR+EAS)" & hits$threshold=="> 0.9"],
            hits$n_genes[hits$panel=="4-way (EUR+EAS+AFR+SAS)" & hits$threshold=="> 0.9"],
            hits$n_genes[hits$panel=="4-way (EUR+EAS+AFR+SAS)" & hits$threshold=="> 0.9"] -
            hits$n_genes[hits$panel=="2-way (EUR+EAS)" & hits$threshold=="> 0.9"])
  )
  write_lines_safe(summary_lines, file.path(OUT, "mesusie/summary.txt"))
  cat("  Wrote MESuSiE figures + summary to ", file.path(OUT, "mesusie"), "\n", sep = "")
}
cat("\nDone. Figures in: ", OUT, "\n", sep = "")
