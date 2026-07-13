#!/usr/bin/env Rscript
# presentation_susiex_mesusie_panels.R
# Cross-LD-panel comparison for SuSiEX (4 variants) and MESuSiE (6 variants).
# Both methods are at the same 140 EUR×EAS shared loci × 3 GWAS pairs.
#
# SuSiEX panels:
#   v1   = sghatan UKBB v1   (results/susiex/susiex_gene_summary.csv)
#   1kg  = 1KG EUR + 1KG EAS  (results/susiex_1kg/susiex_gene_summary_1kg.csv)
#   topld     = TOP-LD EUR + 1KG EAS (results/susiex_topld/susiex_gene_summary_topld.csv)
#   topld_full = TOP-LD EUR + TOP-LD EAS (results/susiex_topld_full/susiex_gene_summary_topld_full.csv)
#
# MESuSiE panels (same 4 + 2 PolyFun variants):
#   v1, 1kg, topld, topld_full, polyfun, polyfun_full
#
# Outputs in figures/supplementary/figS04_coloc/ld_panel_4way_eur/{susiex,mesusie}/:
#   pip_scatter_pairwise.pdf — pairwise PIP scatter (gene-level max PIP)
#   hits_bar.pdf             — gene count at PIP>0.5 / 0.8 / 0.9
#   summary.txt              — pairwise r/rho

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
OUT_DIR <- file.path(PROJ, "figures/supplementary/figS04_coloc/ld_panel_4way_eur")
for (sub in c("susiex", "mesusie")) {
  dir.create(file.path(OUT_DIR, sub), recursive = TRUE, showWarnings = FALSE)
}

PANEL_COLORS <- c("UKBB v1" = "#0D47A1", "1KG" = "#7B1FA2",
                  "TOP-LD" = "#C2185B", "TOP-LD-full" = "#F4511E",
                  "PolyFun" = "#42A5F5", "PolyFun-full" = "#E377C2")

hallmark <- c("RORA","THRB","HKDC1","MTTP","GCKR","CFLAR","EPHA2",
              "TNFSF10","SLC39A8","ADH4","DGAT2","LPL","PNPLA3","TM6SF2",
              "APOE","HSD17B13","MARC1","TRIB1","HFE","SERPINA1")

# ============================================================================
# SuSiEX 4-way (gene-level)
# ============================================================================
cat("\n=== SuSiEX 4-way comparison ===\n")

load_susiex <- function(path, label) {
  if (!file.exists(path)) { cat("  MISSING:", label, "\n"); return(NULL) }
  dt <- fread(path)
  setnames(dt, "susiex_max_pip", paste0("pip_", label))
  dt[, .(ENSG, GeneSymbol, trait_pairs, pip = get(paste0("pip_", label)))]
}

sx_v1   <- load_susiex(file.path(FM_DIR, "results/susiex/susiex_gene_summary.csv"), "v1")
sx_1kg  <- load_susiex(file.path(FM_DIR, "results/susiex_1kg/susiex_gene_summary_1kg.csv"), "1kg")
sx_tl   <- load_susiex(file.path(FM_DIR, "results/susiex_topld/susiex_gene_summary_topld.csv"), "topld")
sx_tlf  <- load_susiex(file.path(FM_DIR, "results/susiex_topld_full/susiex_gene_summary_topld_full.csv"), "topld_full")

if (!is.null(sx_v1) && !is.null(sx_1kg) && !is.null(sx_tl) && !is.null(sx_tlf)) {
  setnames(sx_v1,  "pip", "v1")
  setnames(sx_1kg, "pip", "kg")
  setnames(sx_tl,  "pip", "tl")
  setnames(sx_tlf, "pip", "tlf")

  m_sx <- Reduce(function(a,b) merge(a, b, by = c("ENSG","GeneSymbol","trait_pairs"), all = FALSE),
                  list(sx_v1, sx_1kg[, .(ENSG, GeneSymbol, trait_pairs, kg)],
                       sx_tl[,  .(ENSG, GeneSymbol, trait_pairs, tl)],
                       sx_tlf[, .(ENSG, GeneSymbol, trait_pairs, tlf)]))
  cat(sprintf("  4-way paired rows: %d\n", nrow(m_sx)))

  pair_specs <- list(
    list("UKBB v1 vs 1KG",        "v1",  "kg"),
    list("UKBB v1 vs TOP-LD",     "v1",  "tl"),
    list("UKBB v1 vs TOP-LD-full","v1",  "tlf"),
    list("1KG vs TOP-LD",         "kg",  "tl"),
    list("1KG vs TOP-LD-full",    "kg",  "tlf"),
    list("TOP-LD vs TOP-LD-full", "tl",  "tlf")
  )
  for (p in pair_specs) {
    r  <- cor(m_sx[[p[[2]]]], m_sx[[p[[3]]]])
    rho <- cor(m_sx[[p[[2]]]], m_sx[[p[[3]]]], method = "spearman")
    cat(sprintf("  %-30s r = %.3f   rho = %.3f\n", p[[1]], r, rho))
  }

  # Pairwise scatter (6 panels)
  panel_long <- function(name, x_col, y_col) {
    d <- m_sx[, .(GeneSymbol, x = get(x_col), y = get(y_col))]
    r <- cor(d$x, d$y)
    d[, comparison := sprintf("%s   (r = %.3f)", name, r)]
    d
  }
  long_sx <- rbindlist(lapply(pair_specs, function(p) panel_long(p[[1]], p[[2]], p[[3]])))
  long_sx[, comparison := factor(comparison, levels = unique(comparison))]
  lbl_sx <- long_sx[GeneSymbol %in% hallmark & pmax(x, y, na.rm = TRUE) > 0.5]
  lbl_sx <- unique(lbl_sx, by = c("GeneSymbol", "comparison"))

  p_sx_scatter <- ggplot(long_sx, aes(x = x, y = y)) +
    geom_abline(slope = 1, intercept = 0, color = "#B0B0B0",
                linetype = "dashed", linewidth = 0.8) +
    rasterise(geom_point(alpha = 0.25, size = 0.7, color = "#1565C0"), dpi = 220) +
    geom_hline(yintercept = 0.5, color = "grey60", linetype = "dotted", linewidth = 0.4) +
    geom_vline(xintercept = 0.5, color = "grey60", linetype = "dotted", linewidth = 0.4) +
    geom_point(data = lbl_sx, color = "black", size = 2.2) +
    geom_text_repel(data = lbl_sx, aes(label = GeneSymbol),
                    size = 4.2, color = "black", fontface = "plain",
                    box.padding = 0.4, max.overlaps = 25,
                    segment.color = "grey30", segment.size = 0.3) +
    facet_wrap(~ comparison, ncol = 3) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "max SuSiEX PIP (panel on x)",
         y = "max SuSiEX PIP (panel on y)",
         title = "Gene-level SuSiEX max-PIP concordance — 4 panels (6 pairwise)") +
    theme_masld(base_size = 16) +
    theme(plot.title = element_text(size = 22, face = "plain"),
          strip.text = element_text(size = 14, face = "plain"),
          axis.title = element_text(size = 18))
  ggsave(file.path(OUT_DIR, "susiex/pip_scatter_pairwise.pdf"), p_sx_scatter,
         width = 18, height = 12, device = cairo_pdf)
  cat("  Wrote: susiex/pip_scatter_pairwise.pdf\n")

  # Hits bar
  hits_sx <- data.table(
    panel = factor(rep(c("UKBB v1","1KG","TOP-LD","TOP-LD-full"), each = 3),
                   levels = c("UKBB v1","1KG","TOP-LD","TOP-LD-full")),
    threshold = factor(rep(c("PIP > 0.5","PIP > 0.8","PIP > 0.9"), 4),
                       levels = c("PIP > 0.5","PIP > 0.8","PIP > 0.9")),
    n = c(
      sum(m_sx$v1  > 0.5), sum(m_sx$v1  > 0.8), sum(m_sx$v1  > 0.9),
      sum(m_sx$kg  > 0.5), sum(m_sx$kg  > 0.8), sum(m_sx$kg  > 0.9),
      sum(m_sx$tl  > 0.5), sum(m_sx$tl  > 0.8), sum(m_sx$tl  > 0.9),
      sum(m_sx$tlf > 0.5), sum(m_sx$tlf > 0.8), sum(m_sx$tlf > 0.9)
    )
  )
  p_sx_hits <- ggplot(hits_sx, aes(x = panel, y = n, fill = panel)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = n), vjust = -0.35, size = 5.5, fontface = "plain") +
    facet_wrap(~ threshold, ncol = 3, scales = "free_y") +
    scale_fill_manual(values = PANEL_COLORS, guide = "none") +
    labs(x = NULL, y = "# gene × trait", title = "SuSiEX hits at increasing PIP thresholds (4-way)") +
    theme_masld(base_size = 16) +
    theme(plot.title = element_text(size = 22, face = "plain"),
          strip.text = element_text(size = 16, face = "plain"),
          axis.text.x = element_text(size = 13, angle = 25, hjust = 1)) +
    expand_limits(y = max(hits_sx$n) * 1.18)
  ggsave(file.path(OUT_DIR, "susiex/hits_bar.pdf"), p_sx_hits,
         width = 16, height = 7, device = cairo_pdf)
  cat("  Wrote: susiex/hits_bar.pdf\n")

  # Summary text
  fileConn <- file(file.path(OUT_DIR, "susiex/summary.txt"))
  writeLines(c(
    paste0("SuSiEX 4-panel comparison (", Sys.time(), ")"),
    "",
    paste0("4-way paired rows: ", nrow(m_sx)),
    "",
    "Pairwise concordance:",
    sapply(pair_specs, function(p) {
      r  <- cor(m_sx[[p[[2]]]], m_sx[[p[[3]]]])
      rho <- cor(m_sx[[p[[2]]]], m_sx[[p[[3]]]], method = "spearman")
      sprintf("  %-30s r = %.3f   rho = %.3f", p[[1]], r, rho)
    }),
    "",
    "Hit counts (PIP > 0.5 / 0.8 / 0.9):",
    sprintf("  UKBB v1:      %d / %d / %d", sum(m_sx$v1>0.5),  sum(m_sx$v1>0.8),  sum(m_sx$v1>0.9)),
    sprintf("  1KG:          %d / %d / %d", sum(m_sx$kg>0.5),  sum(m_sx$kg>0.8),  sum(m_sx$kg>0.9)),
    sprintf("  TOP-LD:       %d / %d / %d", sum(m_sx$tl>0.5),  sum(m_sx$tl>0.8),  sum(m_sx$tl>0.9)),
    sprintf("  TOP-LD-full:  %d / %d / %d", sum(m_sx$tlf>0.5), sum(m_sx$tlf>0.8), sum(m_sx$tlf>0.9))
  ), fileConn)
  close(fileConn)
  cat("  Wrote: susiex/summary.txt\n")
} else {
  cat("  Skipping SuSiEX panels — at least one input missing.\n")
}

# ============================================================================
# MESuSiE 6-way (gene-level max PIP across shared/EUR/EAS)
# ============================================================================
cat("\n=== MESuSiE 6-way comparison ===\n")

load_mesusie <- function(path, label) {
  if (!file.exists(path)) { cat("  MISSING:", label, "\n"); return(NULL) }
  dt <- fread(path)
  cols <- c("ENSG", "GeneSymbol", "trait_pairs", "mesusie_max_pip")
  if (!all(cols %in% names(dt))) { cat("  WARN:", label, "missing cols\n"); return(NULL) }
  dt[, .(ENSG, GeneSymbol, trait_pairs, pip = mesusie_max_pip)]
}

ms_v1  <- load_mesusie(file.path(FM_DIR, "results/mesusie/mesusie_gene_summary.csv"), "v1")
ms_kg  <- load_mesusie(file.path(FM_DIR, "results/mesusie_1kg/mesusie_gene_summary.csv"), "1kg")
ms_tl  <- load_mesusie(file.path(FM_DIR, "results/mesusie_topld/mesusie_gene_summary.csv"), "topld")
ms_tlf <- load_mesusie(file.path(FM_DIR, "results/mesusie_topld_full/mesusie_gene_summary.csv"), "topld_full")
ms_pf  <- load_mesusie(file.path(FM_DIR, "results/mesusie_polyfun/mesusie_gene_summary.csv"), "polyfun")
ms_pff <- load_mesusie(file.path(FM_DIR, "results/mesusie_polyfun_full/mesusie_gene_summary.csv"), "polyfun_full")

panels_avail <- list(v1 = ms_v1, kg = ms_kg, tl = ms_tl, tlf = ms_tlf, pf = ms_pf, pff = ms_pff)
panels_avail <- panels_avail[!sapply(panels_avail, is.null)]
cat(sprintf("  Panels available: %d / 6 (%s)\n",
            length(panels_avail), paste(names(panels_avail), collapse = ", ")))

if (length(panels_avail) >= 2) {
  for (k in names(panels_avail)) setnames(panels_avail[[k]], "pip", k)

  m_ms <- Reduce(function(a, b) merge(a, b, by = c("ENSG","GeneSymbol","trait_pairs"), all = FALSE),
                  lapply(names(panels_avail), function(k) {
                    panels_avail[[k]][, c("ENSG","GeneSymbol","trait_pairs", k), with = FALSE]
                  }))
  cat(sprintf("  paired rows across %d panels: %d\n", length(panels_avail), nrow(m_ms)))

  pretty_label <- c(v1 = "UKBB v1", kg = "1KG", tl = "TOP-LD",
                    tlf = "TOP-LD-full", pf = "PolyFun", pff = "PolyFun-full")
  pretty_color <- c(v1 = "#0D47A1", kg = "#7B1FA2", tl = "#C2185B",
                    tlf = "#F4511E", pf = "#42A5F5", pff = "#E377C2")
  panel_keys <- names(panels_avail)
  pair_specs <- combn(panel_keys, 2, simplify = FALSE)

  pair_long <- function(p) {
    d <- m_ms[, .(GeneSymbol, x = get(p[1]), y = get(p[2]))]
    r <- cor(d$x, d$y)
    d[, comparison := sprintf("%s vs %s   (r = %.3f)",
                              pretty_label[p[1]], pretty_label[p[2]], r)]
    d
  }
  long_ms <- rbindlist(lapply(pair_specs, pair_long))
  long_ms[, comparison := factor(comparison, levels = unique(comparison))]
  lbl_ms <- long_ms[GeneSymbol %in% hallmark & pmax(x, y, na.rm = TRUE) > 0.5]
  lbl_ms <- unique(lbl_ms, by = c("GeneSymbol", "comparison"))

  p_ms_scatter <- ggplot(long_ms, aes(x = x, y = y)) +
    geom_abline(slope = 1, intercept = 0, color = "#B0B0B0",
                linetype = "dashed", linewidth = 0.8) +
    rasterise(geom_point(alpha = 0.20, size = 0.6, color = "#1565C0"), dpi = 200) +
    geom_hline(yintercept = 0.5, color = "grey60", linetype = "dotted", linewidth = 0.4) +
    geom_vline(xintercept = 0.5, color = "grey60", linetype = "dotted", linewidth = 0.4) +
    geom_point(data = lbl_ms, color = "black", size = 1.8) +
    geom_text_repel(data = lbl_ms, aes(label = GeneSymbol),
                    size = 3.2, color = "black", fontface = "plain",
                    box.padding = 0.3, max.overlaps = 20,
                    segment.color = "grey30", segment.size = 0.3) +
    facet_wrap(~ comparison, ncol = 5, scales = "fixed") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "max MESuSiE PIP (x)",
         y = "max MESuSiE PIP (y)",
         title = sprintf("Gene-level MESuSiE max-PIP concordance — %d panels (%d pairwise)",
                         length(panels_avail), length(pair_specs))) +
    theme_masld(base_size = 12) +
    theme(plot.title = element_text(size = 18, face = "plain"),
          strip.text = element_text(size = 9, face = "plain"),
          axis.title = element_text(size = 14))
  n_facet_rows <- ceiling(length(pair_specs) / 5)
  ggsave(file.path(OUT_DIR, "mesusie/pip_scatter_pairwise.pdf"), p_ms_scatter,
         width = 22, height = 4 * n_facet_rows + 2, device = cairo_pdf, limitsize = FALSE)
  cat("  Wrote: mesusie/pip_scatter_pairwise.pdf\n")

  # Hits bar (6-panel × 3 thresholds)
  hits_ms <- rbindlist(lapply(panel_keys, function(k) {
    pip <- m_ms[[k]]
    data.table(
      panel = factor(pretty_label[k], levels = pretty_label[panel_keys]),
      threshold = c("PIP > 0.5","PIP > 0.8","PIP > 0.9"),
      n = c(sum(pip > 0.5), sum(pip > 0.8), sum(pip > 0.9))
    )
  }))
  hits_ms[, threshold := factor(threshold, levels = c("PIP > 0.5","PIP > 0.8","PIP > 0.9"))]

  p_ms_hits <- ggplot(hits_ms, aes(x = panel, y = n, fill = panel)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = n), vjust = -0.35, size = 5, fontface = "plain") +
    facet_wrap(~ threshold, ncol = 3, scales = "free_y") +
    scale_fill_manual(values = pretty_color[panel_keys], guide = "none",
                      labels = pretty_label[panel_keys]) +
    labs(x = NULL, y = "# gene × trait",
         title = sprintf("MESuSiE hits at increasing PIP thresholds (%d-way)", length(panels_avail))) +
    theme_masld(base_size = 16) +
    theme(plot.title = element_text(size = 22, face = "plain"),
          strip.text = element_text(size = 16, face = "plain"),
          axis.text.x = element_text(size = 12, angle = 28, hjust = 1)) +
    expand_limits(y = max(hits_ms$n) * 1.18)
  ggsave(file.path(OUT_DIR, "mesusie/hits_bar.pdf"), p_ms_hits,
         width = 16, height = 7, device = cairo_pdf)
  cat("  Wrote: mesusie/hits_bar.pdf\n")

  fileConn <- file(file.path(OUT_DIR, "mesusie/summary.txt"))
  writeLines(c(
    paste0("MESuSiE ", length(panels_avail), "-panel comparison (", Sys.time(), ")"),
    "",
    paste0("Panels: ", paste(pretty_label[panel_keys], collapse = ", ")),
    paste0("Paired rows: ", nrow(m_ms)),
    "",
    "Pairwise concordance:",
    sapply(pair_specs, function(p) {
      r <- cor(m_ms[[p[1]]], m_ms[[p[2]]])
      rho <- cor(m_ms[[p[1]]], m_ms[[p[2]]], method = "spearman")
      sprintf("  %-30s r = %.3f   rho = %.3f",
              paste(pretty_label[p[1]], "vs", pretty_label[p[2]]), r, rho)
    }),
    "",
    "Hit counts (PIP > 0.5 / 0.8 / 0.9):",
    sapply(panel_keys, function(k) {
      sprintf("  %-13s %d / %d / %d", pretty_label[k],
              sum(m_ms[[k]] > 0.5), sum(m_ms[[k]] > 0.8), sum(m_ms[[k]] > 0.9))
    })
  ), fileConn)
  close(fileConn)
  cat("  Wrote: mesusie/summary.txt\n")
}

cat("\nALL DONE. Outputs in:", OUT_DIR, "\n")
