#!/usr/bin/env Rscript
# presentation_3way_full.R
# Comprehensive 3-way LD-panel comparison for both COLOC (gene-level) and
# standalone fine-mapping (per-variant PIP). Three panels: UKBB v1 (sghatan,
# N≈337K) / 1KG EUR (N=379) / TOP-LD EUR (N=13,160).
#
# Outputs (figures/presentation/3way/):
#   coloc_3way_scatter.pdf            — 3 facetted scatters (UKBB×1KG, UKBB×TOPLD, 1KG×TOPLD)
#   coloc_3way_hits_bar.pdf           — PP4>0.5/0.8/0.9 hit counts per panel
#   coloc_3way_upset.pdf              — Set agreement (PP4>0.5 sets)
#   coloc_3way_hallmarks.pdf          — Hallmark genes head-to-head
#   fm_3way_scatter.pdf               — Per-variant SuSiE PIP, 3 facets
#   fm_3way_pip_hits_bar.pdf          — PIP>0.5 hit counts per panel
#   3way_summary.txt                  — text summary of all numbers

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
OUT_DIR <- file.path(PROJ, "figures/supplementary/figS04_coloc/ld_panel_3way_eur")
for (sub in c("", "coloc", "fm", "susiex"))
  dir.create(file.path(OUT_DIR, sub), recursive = TRUE, showWarnings = FALSE)

# ==========================================================================
# COLOC (gene-level) 3-way
# ==========================================================================
load_coloc <- function(subdir) {
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

cat("\n=== Loading COLOC outputs ===\n")
v1    <- load_coloc("susie_coloc_ukbb_sghatan_v1_2026-04-21")
kg    <- load_coloc("susie_coloc_1kg")
topld <- load_coloc("susie_coloc_topld")
cat(sprintf("  UKBB v1: %d rows\n  1KG EUR: %d rows\n  TOP-LD : %d rows\n",
            nrow(v1), nrow(kg), nrow(topld)))

# Restrict to 17 EUR GWAS (1kg/topld panels are EUR-only)
eur_gwas <- unique(kg$gwas_name)
v1    <- v1[gwas_name %in% eur_gwas]

m_coloc <- Reduce(function(a,b) merge(a, b, by = c("gwas_name","ensembl"), all = FALSE),
  list(
    v1[, .(gwas_name, ensembl, gene, pp4_ukbb = max_pp4)],
    kg[, .(gwas_name, ensembl, pp4_1kg = max_pp4)],
    topld[, .(gwas_name, ensembl, pp4_topld = max_pp4)]
  ))
m_coloc <- m_coloc[!is.na(pp4_ukbb) & !is.na(pp4_1kg) & !is.na(pp4_topld)]
cat(sprintf("  3-way paired rows: %d\n", nrow(m_coloc)))

# r and ρ
panels <- list(
  list("UKBB vs 1KG",   "pp4_ukbb", "pp4_1kg"),
  list("UKBB vs TOP-LD","pp4_ukbb", "pp4_topld"),
  list("1KG  vs TOP-LD","pp4_1kg",  "pp4_topld")
)
for (p in panels) {
  r <- cor(m_coloc[[p[[2]]]], m_coloc[[p[[3]]]])
  cat(sprintf("  %-20s r = %.3f\n", p[[1]], r))
}

# Hallmark genes
hallmark <- c("RORA", "THRB", "HKDC1", "MTTP", "GCKR", "CFLAR", "EPHA2",
              "TNFSF10", "SLC39A8", "ADH4", "DGAT2", "LPL")

# ----- Panel: 3-way scatter, faceted -----
panel_long <- function(name, x_col, y_col) {
  d <- m_coloc[, .(gene, x = get(x_col), y = get(y_col))]
  r <- cor(d$x, d$y)
  d[, comparison := sprintf("%s   (r = %.3f)", name, r)]
  d
}
long_coloc <- rbind(
  panel_long("UKBB vs 1KG",    "pp4_ukbb", "pp4_1kg"),
  panel_long("UKBB vs TOP-LD", "pp4_ukbb", "pp4_topld"),
  panel_long("1KG  vs TOP-LD", "pp4_1kg",  "pp4_topld")
)
lbl_long <- long_coloc[gene %in% hallmark & pmax(x, y, na.rm = TRUE) > 0.5]
lbl_long <- unique(lbl_long, by = c("gene", "comparison"))

p_coloc_scatter <- ggplot(long_coloc, aes(x = x, y = y)) +
  geom_abline(slope = 1, intercept = 0,
              color = "#B0B0B0", linetype = "dashed", linewidth = 0.8) +
  rasterise(geom_point(alpha = 0.25, size = 0.7, color = "#1565C0"), dpi = 220) +
  geom_hline(yintercept = 0.5, color = "grey60", linetype = "dotted", linewidth = 0.4) +
  geom_vline(xintercept = 0.5, color = "grey60", linetype = "dotted", linewidth = 0.4) +
  geom_point(data = lbl_long, color = "black", size = 2.5) +
  geom_text_repel(data = lbl_long, aes(label = gene),
                  size = 5, color = "black", fontface = "plain",
                  box.padding = 0.5, max.overlaps = 25,
                  segment.color = "grey30", segment.size = 0.3) +
  facet_wrap(~ comparison, ncol = 3) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "PP.H4 (panel on x)", y = "PP.H4 (panel on y)",
       title = "Gene-level SuSiE-COLOC concordance across 3 LD panels") +
  theme_masld(base_size = 18) +
  theme(plot.title = element_text(size = 26, face = "plain"),
        strip.text = element_text(size = 18, face = "plain"),
        axis.title = element_text(size = 20),
        axis.text  = element_text(size = 16))

ggsave(file.path(OUT_DIR, "coloc/scatter.pdf"), p_coloc_scatter,
       width = 20, height = 8, device = cairo_pdf)
cat("  Wrote: coloc_3way_scatter.pdf\n")

# ----- Panel: hit counts per panel -----
hits_coloc <- data.table(
  panel = factor(rep(c("UKBB v1", "1KG EUR", "TOP-LD EUR"), each = 3),
                 levels = c("UKBB v1", "1KG EUR", "TOP-LD EUR")),
  threshold = factor(rep(c("PP.H4 > 0.5", "PP.H4 > 0.8", "PP.H4 > 0.9"), 3),
                     levels = c("PP.H4 > 0.5", "PP.H4 > 0.8", "PP.H4 > 0.9")),
  n = c(
    sum(m_coloc$pp4_ukbb  > 0.5), sum(m_coloc$pp4_ukbb  > 0.8), sum(m_coloc$pp4_ukbb  > 0.9),
    sum(m_coloc$pp4_1kg   > 0.5), sum(m_coloc$pp4_1kg   > 0.8), sum(m_coloc$pp4_1kg   > 0.9),
    sum(m_coloc$pp4_topld > 0.5), sum(m_coloc$pp4_topld > 0.8), sum(m_coloc$pp4_topld > 0.9)
  )
)
p_coloc_hits <- ggplot(hits_coloc, aes(x = panel, y = n, fill = panel)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), vjust = -0.35, size = 6, fontface = "plain") +
  facet_wrap(~ threshold, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = c("UKBB v1" = "#0D47A1", "1KG EUR" = "#7B1FA2",
                                "TOP-LD EUR" = "#C2185B"), guide = "none") +
  labs(x = NULL, y = "# gene × GWAS",
       title = "Hit counts at increasing PP.H4 thresholds (3-way)") +
  theme_masld(base_size = 18) +
  theme(plot.title = element_text(size = 26, face = "plain"),
        strip.text = element_text(size = 18, face = "plain"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        panel.grid.major.x = element_blank()) +
  expand_limits(y = max(hits_coloc$n) * 1.18)
ggsave(file.path(OUT_DIR, "coloc/hits_bar.pdf"), p_coloc_hits,
       width = 16, height = 7, device = cairo_pdf)
cat("  Wrote: coloc_3way_hits_bar.pdf\n")

# ----- Panel: 3-way set agreement (PP4>0.5) -----
hi <- m_coloc[, .(uk = pp4_ukbb > 0.5, kg = pp4_1kg > 0.5, tl = pp4_topld > 0.5)]
sets <- data.table(
  category = factor(c(
    "All 3 agree", "UKBB+1KG only", "UKBB+TOP-LD only", "1KG+TOP-LD only",
    "UKBB only", "1KG only", "TOP-LD only"
  ), levels = c("All 3 agree", "UKBB+1KG only", "UKBB+TOP-LD only", "1KG+TOP-LD only",
                "UKBB only", "1KG only", "TOP-LD only")),
  n = c(
    sum(hi$uk & hi$kg & hi$tl),
    sum(hi$uk & hi$kg & !hi$tl),
    sum(hi$uk & !hi$kg & hi$tl),
    sum(!hi$uk & hi$kg & hi$tl),
    sum(hi$uk & !hi$kg & !hi$tl),
    sum(!hi$uk & hi$kg & !hi$tl),
    sum(!hi$uk & !hi$kg & hi$tl)
  )
)
p_coloc_set <- ggplot(sets, aes(x = category, y = n, fill = category)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = n), vjust = -0.3, size = 6, fontface = "plain") +
  scale_fill_manual(values = c("All 3 agree" = "#0D47A1",
                                "UKBB+1KG only" = "#1F77B4",
                                "UKBB+TOP-LD only" = "#9467BD",
                                "1KG+TOP-LD only" = "#FF7F0E",
                                "UKBB only" = "#7F7F7F",
                                "1KG only" = "#7F7F7F",
                                "TOP-LD only" = "#7F7F7F"),
                    guide = "none") +
  labs(x = NULL, y = "# gene × GWAS at PP.H4 > 0.5",
       title = "3-way set agreement at PP.H4 > 0.5") +
  theme_masld(base_size = 18) +
  theme(plot.title = element_text(size = 26, face = "plain"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        panel.grid.major.x = element_blank()) +
  expand_limits(y = max(sets$n) * 1.18)
ggsave(file.path(OUT_DIR, "coloc/set_agreement.pdf"), p_coloc_set,
       width = 12, height = 7, device = cairo_pdf)
cat("  Wrote: coloc_3way_set_agreement.pdf\n")

# ----- Panel: hallmark gene table -----
hl_tbl <- m_coloc[gene %in% hallmark][order(-pmax(pp4_ukbb, pp4_1kg, pp4_topld))]
hl_tbl <- unique(hl_tbl, by = "gene")[1:min(12, .N)]
hl_tbl <- hl_tbl[!is.na(gene),
                 .(gene,
                   UKBB    = sprintf("%.3f", pp4_ukbb),
                   `1KG`   = sprintf("%.3f", pp4_1kg),
                   `TOPLD` = sprintf("%.3f", pp4_topld))]

N_HL <- nrow(hl_tbl)
hl_tbl[, row_y := rev(seq_len(N_HL))]
p_hl <- ggplot(hl_tbl) +
  geom_text(aes(x = 0,   y = row_y, label = gene),
            hjust = 0, size = 7, fontface = "plain") +
  geom_text(aes(x = 1.2, y = row_y, label = UKBB),
            hjust = 1, size = 7) +
  geom_text(aes(x = 2.4, y = row_y, label = `1KG`),
            hjust = 1, size = 7) +
  geom_text(aes(x = 3.6, y = row_y, label = TOPLD),
            hjust = 1, size = 7) +
  annotate("text", x = c(0, 1.2, 2.4, 3.6), y = N_HL + 0.8,
           label = c("Gene", "UKBB v1", "1KG EUR", "TOP-LD EUR"),
           hjust = c(0, 1, 1, 1), size = 8, fontface = "plain") +
  scale_x_continuous(limits = c(-0.1, 3.8)) +
  scale_y_continuous(limits = c(0.5, N_HL + 1.3)) +
  labs(title = "Hallmark genes — 3-way head-to-head") +
  theme_void(base_size = 20) +
  theme(plot.title = element_text(size = 24, face = "plain"),
        plot.margin = margin(15, 15, 15, 15))
ggsave(file.path(OUT_DIR, "coloc/hallmarks.pdf"), p_hl,
       width = 9, height = 8, device = cairo_pdf)
cat("  Wrote: coloc_3way_hallmarks.pdf\n")

# ==========================================================================
# Standalone fine-mapping (per-variant PIP) 3-way
# ==========================================================================
cat("\n=== Loading standalone FM outputs ===\n")
load_fm <- function(panel_dir) {
  registry <- fread(file.path(FM_DIR, "config/gwas_registry.tsv"))
  studies  <- registry[ancestry == "EUR"]$study_name
  files <- unlist(lapply(studies, function(s) {
    list.files(file.path(FM_DIR, "output", s, panel_dir, "susie"),
               pattern = "\\.tsv$", full.names = TRUE)
  }))
  if (length(files) == 0) return(NULL)
  rows <- lapply(files, function(f) {
    dt <- try(fread(f), silent = TRUE)
    if (inherits(dt, "try-error") || nrow(dt) == 0) return(NULL)
    # path is .../output/<study>/<panel_dir>/susie/<file>.tsv → study = 3 dirs up
    dt[, study := basename(dirname(dirname(dirname(f))))]
    dt[, locus_id := gsub("_cov0\\.95_.*", "", basename(f))]
    dt[, converged_flag := !grepl("_notconverged", basename(f))]
    dt
  })
  dt <- rbindlist(rows, fill = TRUE)
  dt[]
}
fm_v1    <- load_fm("EUR_0.5Mb")
fm_kg    <- load_fm("1kg_eur_0.5Mb")
fm_topld <- load_fm("topld_eur_0.5Mb")

cat(sprintf("  UKBB v1: %s rows\n", if (is.null(fm_v1)) "MISSING" else format(nrow(fm_v1))))
cat(sprintf("  1KG EUR: %s rows\n", if (is.null(fm_kg)) "MISSING" else format(nrow(fm_kg))))
cat(sprintf("  TOP-LD : %s rows\n", if (is.null(fm_topld)) "MISSING" else format(nrow(fm_topld))))

if (!is.null(fm_topld) && nrow(fm_topld) > 0) {
  # Pair-specific 2-way joins (don't require the third panel to also have the
  # variant — that strict 3-way intersection drops most credible-set candidates,
  # esp. when TopLD has lower variant density). Each pair gets its own legitimate
  # sample of converged-in-both variants.
  v1_dt    <- fm_v1[converged_flag == TRUE,   .(study, locus_id, SNP, PIP_ukbb  = PIP)]
  kg_dt    <- fm_kg[converged_flag == TRUE,   .(study, locus_id, SNP, PIP_1kg   = PIP)]
  topld_dt <- fm_topld[converged_flag == TRUE,.(study, locus_id, SNP, PIP_topld = PIP)]

  m_uv_kg_full    <- merge(v1_dt, kg_dt,    by = c("study", "locus_id", "SNP"), all = FALSE)
  m_uv_topld_full <- merge(v1_dt, topld_dt, by = c("study", "locus_id", "SNP"), all = FALSE)
  m_kg_topld_full <- merge(kg_dt, topld_dt, by = c("study", "locus_id", "SNP"), all = FALSE)
  m_uv_kg_full    <- m_uv_kg_full[!is.na(PIP_ukbb) & !is.na(PIP_1kg)]
  m_uv_topld_full <- m_uv_topld_full[!is.na(PIP_ukbb) & !is.na(PIP_topld)]
  m_kg_topld_full <- m_kg_topld_full[!is.na(PIP_1kg)  & !is.na(PIP_topld)]
  cat(sprintf("  2-way paired variants — UKBB∩1KG=%d  UKBB∩TOPLD=%d  1KG∩TOPLD=%d\n",
              nrow(m_uv_kg_full), nrow(m_uv_topld_full), nrow(m_kg_topld_full)))

  # PIP > 0.1 filter: keep variants where PIP > 0.1 in either panel of the pair.
  m_uv_kg    <- m_uv_kg_full[PIP_ukbb > 0.1 | PIP_1kg   > 0.1]
  m_uv_topld <- m_uv_topld_full[PIP_ukbb > 0.1 | PIP_topld > 0.1]
  m_kg_topld <- m_kg_topld_full[PIP_1kg  > 0.1 | PIP_topld > 0.1]
  m_fm <- m_uv_kg  # for the if-exists guard below; not used downstream

  fm_long <- rbind(
    m_uv_kg[, .(comparison = sprintf("UKBB vs 1KG    (r=%.2f, n=%d)",
                                  cor(PIP_ukbb, PIP_1kg), .N),
             x = PIP_ukbb,  y = PIP_1kg)],
    m_uv_topld[, .(comparison = sprintf("UKBB vs TOP-LD (r=%.2f, n=%d)",
                                  cor(PIP_ukbb, PIP_topld), .N),
             x = PIP_ukbb,  y = PIP_topld)],
    m_kg_topld[, .(comparison = sprintf("1KG  vs TOP-LD (r=%.2f, n=%d)",
                                  cor(PIP_1kg, PIP_topld), .N),
             x = PIP_1kg,   y = PIP_topld)]
  )
  cat(sprintf("  PIP > 0.1 (either-panel) variant counts: UKBB-vs-1KG=%d  UKBB-vs-TOPLD=%d  1KG-vs-TOPLD=%d\n",
              nrow(m_uv_kg), nrow(m_uv_topld), nrow(m_kg_topld)))

  p_fm_scatter <- ggplot(fm_long, aes(x = x, y = y)) +
    geom_abline(slope = 1, intercept = 0,
                color = "#B0B0B0", linetype = "dashed", linewidth = 0.8) +
    rasterise(geom_point(alpha = 0.55, size = 1.4, color = "#0D47A1"), dpi = 220) +
    geom_hline(yintercept = 0.5, color = "grey60", linetype = "dotted", linewidth = 0.4) +
    geom_vline(xintercept = 0.5, color = "grey60", linetype = "dotted", linewidth = 0.4) +
    facet_wrap(~ comparison, ncol = 3) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "PIP", y = "PIP", title = "SuSiE PIP — 3 EUR LD panels (PIP > 0.1)") +
    theme_masld(base_size = 18) +
    theme(plot.title = element_text(size = 22, face = "plain"),
          strip.text = element_text(size = 16, face = "plain"),
          axis.title = element_text(size = 18),
          axis.text  = element_text(size = 14))
  ggsave(file.path(OUT_DIR, "fm/scatter.pdf"), p_fm_scatter,
         width = 20, height = 8, device = cairo_pdf)
  cat("  Wrote: fm_3way_scatter.pdf\n")

  # Per-panel RAW counts (no intersection — the honest number per panel).
  # The 3-way-intersected counts collapse because the strict variant×locus filter
  # drops strong-signal loci that fail to converge in 1+ panels.
  raw_v1 <- fm_v1[converged_flag == TRUE]
  raw_kg <- fm_kg[converged_flag == TRUE]
  raw_tl <- fm_topld[converged_flag == TRUE]
  hits_fm <- data.table(
    panel = factor(rep(c("UKBB v1", "1KG EUR", "TOP-LD EUR"), each = 2),
                   levels = c("UKBB v1", "1KG EUR", "TOP-LD EUR")),
    threshold = factor(rep(c("PIP > 0.5", "PIP > 0.9"), 3),
                       levels = c("PIP > 0.5", "PIP > 0.9")),
    n = c(
      sum(raw_v1$PIP > 0.5, na.rm = TRUE), sum(raw_v1$PIP > 0.9, na.rm = TRUE),
      sum(raw_kg$PIP > 0.5, na.rm = TRUE), sum(raw_kg$PIP > 0.9, na.rm = TRUE),
      sum(raw_tl$PIP > 0.5, na.rm = TRUE), sum(raw_tl$PIP > 0.9, na.rm = TRUE)
    )
  )
  p_fm_hits <- ggplot(hits_fm, aes(x = panel, y = n, fill = panel)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = n), vjust = -0.35, size = 7, fontface = "plain") +
    facet_wrap(~ threshold, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = c("UKBB v1" = "#0D47A1", "1KG EUR" = "#7B1FA2",
                                  "TOP-LD EUR" = "#C2185B"), guide = "none") +
    labs(x = NULL, y = "# variants",
         title = "High-PIP variant counts (per panel, raw)") +
    theme_masld(base_size = 18) +
    theme(plot.title = element_text(size = 26, face = "plain"),
          strip.text = element_text(size = 18, face = "plain"),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 18),
          panel.grid.major.x = element_blank()) +
    expand_limits(y = max(hits_fm$n) * 1.18)
  ggsave(file.path(OUT_DIR, "fm/pip_hits_bar.pdf"), p_fm_hits,
         width = 12, height = 7, device = cairo_pdf)
  cat("  Wrote: fm_3way_pip_hits_bar.pdf (per-panel raw)\n")

  # Also: count loci with at least one high-PIP variant per panel
  hi_loci <- data.table(
    panel = factor(c("UKBB v1", "1KG EUR", "TOP-LD EUR"),
                   levels = c("UKBB v1", "1KG EUR", "TOP-LD EUR")),
    n_loci = c(
      uniqueN(raw_v1[PIP > 0.5, .(study, locus_id)]),
      uniqueN(raw_kg[PIP > 0.5, .(study, locus_id)]),
      uniqueN(raw_tl[PIP > 0.5, .(study, locus_id)])
    )
  )
  p_fm_loci <- ggplot(hi_loci, aes(x = panel, y = n_loci, fill = panel)) +
    geom_col(width = 0.55) +
    geom_text(aes(label = n_loci), vjust = -0.3, size = 9, fontface = "plain") +
    scale_fill_manual(values = c("UKBB v1" = "#0D47A1", "1KG EUR" = "#7B1FA2",
                                  "TOP-LD EUR" = "#C2185B"), guide = "none") +
    labs(x = NULL, y = "# loci with ≥1 PIP>0.5 variant",
         title = "Loci with high-PIP variant (per panel)") +
    theme_masld(base_size = 18) +
    theme(plot.title = element_text(size = 26, face = "plain"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 20),
          panel.grid.major.x = element_blank()) +
    expand_limits(y = max(hi_loci$n_loci) * 1.18)
  ggsave(file.path(OUT_DIR, "fm/loci_with_high_pip.pdf"), p_fm_loci,
         width = 9, height = 7, device = cairo_pdf)
  cat("  Wrote: fm_3way_loci_with_high_pip.pdf\n")
} else {
  cat("\n  TOP-LD standalone FM not yet available — fm_3way figures skipped.\n")
  cat("  Re-run this script after fm_topld jobs complete.\n")
}

# ==========================================================================
# SuSiEX cross-ancestry joint fine-mapping — 3-way
#   - results/susiex/        (UKBB v1)
#   - results/susiex_1kg/    (1KG EUR)
#   - results/susiex_topld/  (TOP-LD EUR)
# Each panel's *.summary file gives one row per credible set with MAX_PIP_SNP,
# CS_LENGTH, MAX_PIP, etc. Aggregate to per-locus top PIP for cross-panel
# comparison.
# ==========================================================================
cat("\n=== Loading SuSiEX outputs ===\n")
load_susiex <- function(root) {
  if (!dir.exists(root)) return(NULL)
  files <- list.files(root, pattern = "\\.summary$", recursive = TRUE,
                      full.names = TRUE)
  if (length(files) == 0) return(NULL)
  rows <- lapply(files, function(f) {
    ln <- try(readLines(f), silent = TRUE)
    if (inherits(ln, "try-error") || length(ln) < 2) return(NULL)
    data_lines <- ln[!grepl("^#", ln) & nchar(ln) > 0]
    if (length(data_lines) < 2) return(NULL)
    dt <- try(fread(text = paste(data_lines, collapse = "\n")), silent = TRUE)
    if (inherits(dt, "try-error") || nrow(dt) == 0) return(NULL)
    dt[, locus_id := sub("\\.summary$", "", basename(f))]
    dt[, trait    := basename(dirname(f))]
    dt
  })
  rbindlist(rows, fill = TRUE)
}
sx_v1    <- load_susiex(file.path(FM_DIR, "results/susiex"))
sx_kg    <- load_susiex(file.path(FM_DIR, "results/susiex_1kg"))
sx_topld <- load_susiex(file.path(FM_DIR, "results/susiex_topld"))
cat(sprintf("  UKBB v1: %s\n", if (is.null(sx_v1))    "MISSING" else paste(nrow(sx_v1),    "CS rows across", uniqueN(sx_v1$locus_id), "loci")))
cat(sprintf("  1KG EUR: %s\n", if (is.null(sx_kg))    "MISSING" else paste(nrow(sx_kg),    "CS rows across", uniqueN(sx_kg$locus_id), "loci")))
cat(sprintf("  TOP-LD : %s\n", if (is.null(sx_topld)) "MISSING" else paste(nrow(sx_topld), "CS rows across", uniqueN(sx_topld$locus_id), "loci")))

dir.create(file.path(OUT_DIR, "susiex"), recursive = TRUE, showWarnings = FALSE)

if (!is.null(sx_topld) && nrow(sx_topld) > 0) {
  agg_sx <- function(dt) {
    dt[, .(top_pip = max(MAX_PIP, na.rm = TRUE),
           n_cs    = .N,
           top_snp = MAX_PIP_SNP[which.max(MAX_PIP)]),
       by = .(trait, locus_id)]
  }
  ag_v1    <- agg_sx(sx_v1)
  ag_kg    <- agg_sx(sx_kg)
  ag_topld <- agg_sx(sx_topld)

  m_sx <- Reduce(function(a, b) merge(a, b, by = c("trait", "locus_id"), all = FALSE),
    list(
      ag_v1[,    .(trait, locus_id, top_pip_ukbb  = top_pip, n_cs_ukbb  = n_cs, top_snp_ukbb  = top_snp)],
      ag_kg[,    .(trait, locus_id, top_pip_1kg   = top_pip, n_cs_1kg   = n_cs, top_snp_1kg   = top_snp)],
      ag_topld[, .(trait, locus_id, top_pip_topld = top_pip, n_cs_topld = n_cs, top_snp_topld = top_snp)]
    ))
  cat(sprintf("  3-way paired loci: %d\n", nrow(m_sx)))

  if (nrow(m_sx) > 0) {
    sx_long <- rbind(
      m_sx[, .(comparison = sprintf("UKBB vs 1KG    (r = %.3f)", cor(top_pip_ukbb, top_pip_1kg)),
               x = top_pip_ukbb,  y = top_pip_1kg,   trait)],
      m_sx[, .(comparison = sprintf("UKBB vs TOP-LD (r = %.3f)", cor(top_pip_ukbb, top_pip_topld)),
               x = top_pip_ukbb,  y = top_pip_topld, trait)],
      m_sx[, .(comparison = sprintf("1KG  vs TOP-LD (r = %.3f)", cor(top_pip_1kg, top_pip_topld)),
               x = top_pip_1kg,   y = top_pip_topld, trait)]
    )
    p_sx_scatter <- ggplot(sx_long, aes(x = x, y = y, color = trait)) +
      geom_abline(slope = 1, intercept = 0,
                  color = "#B0B0B0", linetype = "dashed", linewidth = 0.8) +
      geom_point(alpha = 0.7, size = 3) +
      scale_color_manual(values = c(ALT = "#0D47A1", AST = "#7B1FA2", GGT = "#C2185B")) +
      facet_wrap(~ comparison, ncol = 3) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(x = "Top PIP (panel on x)", y = "Top PIP (panel on y)",
           title = "SuSiEX per-locus top PIP across 3 LD panels") +
      theme_masld(base_size = 18) +
      theme(plot.title = element_text(size = 26, face = "plain"),
            strip.text = element_text(size = 18, face = "plain"),
            axis.title = element_text(size = 20),
            axis.text  = element_text(size = 16),
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 18))
    ggsave(file.path(OUT_DIR, "susiex/scatter.pdf"), p_sx_scatter,
           width = 21, height = 8, device = cairo_pdf)
    cat("  Wrote: susiex/scatter.pdf\n")

    # Per-panel CS-count + top-PIP>0.5 hit counts (raw, no intersection)
    hits_sx <- data.table(
      panel = factor(rep(c("UKBB v1", "1KG EUR", "TOP-LD EUR"), each = 2),
                     levels = c("UKBB v1", "1KG EUR", "TOP-LD EUR")),
      metric = factor(rep(c("Loci with top PIP > 0.5", "Total credible sets"), 3),
                      levels = c("Loci with top PIP > 0.5", "Total credible sets")),
      n = c(
        sum(ag_v1$top_pip    > 0.5, na.rm = TRUE), nrow(sx_v1),
        sum(ag_kg$top_pip    > 0.5, na.rm = TRUE), nrow(sx_kg),
        sum(ag_topld$top_pip > 0.5, na.rm = TRUE), nrow(sx_topld)
      )
    )
    p_sx_hits <- ggplot(hits_sx, aes(x = panel, y = n, fill = panel)) +
      geom_col(width = 0.7) +
      geom_text(aes(label = n), vjust = -0.35, size = 7, fontface = "plain") +
      facet_wrap(~ metric, ncol = 2, scales = "free_y") +
      scale_fill_manual(values = c("UKBB v1" = "#0D47A1", "1KG EUR" = "#7B1FA2",
                                    "TOP-LD EUR" = "#C2185B"), guide = "none") +
      labs(x = NULL, y = "count",
           title = "SuSiEX hit counts (per panel, raw)") +
      theme_masld(base_size = 18) +
      theme(plot.title = element_text(size = 26, face = "plain"),
            strip.text = element_text(size = 18, face = "plain"),
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 14),
            axis.title.y = element_text(size = 18),
            panel.grid.major.x = element_blank()) +
      expand_limits(y = max(hits_sx$n) * 1.18)
    ggsave(file.path(OUT_DIR, "susiex/hits_bar.pdf"), p_sx_hits,
           width = 14, height = 7, device = cairo_pdf)
    cat("  Wrote: susiex/hits_bar.pdf\n")

    # Top-PIP SNP agreement (3-way)
    same_all  <- sum(m_sx$top_snp_ukbb == m_sx$top_snp_1kg & m_sx$top_snp_1kg == m_sx$top_snp_topld, na.rm = TRUE)
    same_v1kg <- sum(m_sx$top_snp_ukbb == m_sx$top_snp_1kg, na.rm = TRUE)
    same_v1tl <- sum(m_sx$top_snp_ukbb == m_sx$top_snp_topld, na.rm = TRUE)
    same_kgtl <- sum(m_sx$top_snp_1kg  == m_sx$top_snp_topld, na.rm = TRUE)
    snp_agree <- data.table(
      pair = factor(c("All 3 agree", "UKBB ↔ 1KG", "UKBB ↔ TOP-LD", "1KG ↔ TOP-LD"),
                    levels = c("All 3 agree", "UKBB ↔ 1KG", "UKBB ↔ TOP-LD", "1KG ↔ TOP-LD")),
      n = c(same_all, same_v1kg, same_v1tl, same_kgtl),
      pct = c(same_all, same_v1kg, same_v1tl, same_kgtl) / nrow(m_sx) * 100
    )
    p_sx_snp <- ggplot(snp_agree, aes(x = pair, y = pct, fill = pair)) +
      geom_col(width = 0.55) +
      geom_text(aes(label = sprintf("%d (%.0f%%)", n, pct)),
                vjust = -0.3, size = 7, fontface = "plain") +
      scale_fill_manual(values = c("All 3 agree" = "#0D47A1", "UKBB ↔ 1KG" = "#1F77B4",
                                    "UKBB ↔ TOP-LD" = "#9467BD", "1KG ↔ TOP-LD" = "#FF7F0E"),
                        guide = "none") +
      labs(x = NULL, y = "% loci with same top-PIP SNP",
           title = sprintf("Top-PIP SNP agreement (n=%d loci)", nrow(m_sx))) +
      theme_masld(base_size = 18) +
      theme(plot.title = element_text(size = 24, face = "plain"),
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 14),
            axis.title.y = element_text(size = 18),
            panel.grid.major.x = element_blank()) +
      expand_limits(y = max(snp_agree$pct) * 1.18)
    ggsave(file.path(OUT_DIR, "susiex/top_snp_agreement.pdf"), p_sx_snp,
           width = 11, height = 7, device = cairo_pdf)
    cat("  Wrote: susiex/top_snp_agreement.pdf\n")
  }
} else {
  cat("\n  TOP-LD SuSiEX outputs not yet available — susiex 3-way figures skipped.\n")
  cat("  Re-run this script after susiex_topld jobs complete.\n")
}

# ==========================================================================
# Summary text
# ==========================================================================
sink(file.path(OUT_DIR, "summary.txt"))
cat("3-way LD-panel comparison: UKBB v1 (sghatan, N≈337K) vs 1KG EUR (N=379) vs TOP-LD EUR (N=13,160)\n")
cat("Generated: ", as.character(Sys.time()), "\n\n", sep = "")

cat("=== Gene-level COLOC ===\n")
cat(sprintf("Paired gene × GWAS rows: %d\n", nrow(m_coloc)))
for (p in panels) {
  r <- cor(m_coloc[[p[[2]]]], m_coloc[[p[[3]]]])
  rho <- cor(m_coloc[[p[[2]]]], m_coloc[[p[[3]]]], method = "spearman")
  cat(sprintf("  %-22s  r = %.4f   ρ = %.4f\n", p[[1]], r, rho))
}
cat("\nHit counts at PP.H4 thresholds (within 3-way paired set):\n")
print(hits_coloc, row.names = FALSE)
cat("\n3-way set agreement at PP.H4 > 0.5:\n")
print(sets, row.names = FALSE)
cat("\nHallmark gene head-to-head:\n")
print(hl_tbl, row.names = FALSE)

cat("\n=== Standalone fine-mapping (per-variant SuSiE PIP) ===\n")
if (exists("m_fm")) {
  cat(sprintf("Paired variants: %d\n", nrow(m_fm)))
  for (lab in list(c("UKBB vs 1KG","PIP_ukbb","PIP_1kg"),
                   c("UKBB vs TOP-LD","PIP_ukbb","PIP_topld"),
                   c("1KG  vs TOP-LD","PIP_1kg","PIP_topld"))) {
    r   <- cor(m_fm[[lab[2]]], m_fm[[lab[3]]])
    rho <- cor(m_fm[[lab[2]]], m_fm[[lab[3]]], method = "spearman")
    cat(sprintf("  %-22s  r = %.4f   ρ = %.4f\n", lab[1], r, rho))
  }
  if (exists("hits_fm")) {
    cat("\nHigh-PIP counts:\n")
    print(hits_fm, row.names = FALSE)
  }
} else {
  cat("TOP-LD standalone FM not yet available (m_fm has 0 paired rows or doesn't exist).\n")
}

cat("\n=== SuSiEX cross-ancestry (per-locus top PIP) ===\n")
if (exists("m_sx")) {
  cat(sprintf("Paired loci: %d\n", nrow(m_sx)))
  for (lab in list(c("UKBB vs 1KG","top_pip_ukbb","top_pip_1kg"),
                   c("UKBB vs TOP-LD","top_pip_ukbb","top_pip_topld"),
                   c("1KG  vs TOP-LD","top_pip_1kg","top_pip_topld"))) {
    r   <- cor(m_sx[[lab[2]]], m_sx[[lab[3]]])
    rho <- cor(m_sx[[lab[2]]], m_sx[[lab[3]]], method = "spearman")
    cat(sprintf("  %-22s  r = %.4f   ρ = %.4f\n", lab[1], r, rho))
  }
  if (exists("hits_sx")) {
    cat("\nHit counts (raw, per panel):\n")
    print(hits_sx, row.names = FALSE)
  }
  if (exists("snp_agree")) {
    cat("\nTop-PIP SNP agreement (3-way):\n")
    print(snp_agree, row.names = FALSE)
  }
} else {
  cat("TOP-LD SuSiEX outputs not yet available.\n")
}
sink()

cat("\nAll outputs in: ", OUT_DIR, "\n", sep = "")
