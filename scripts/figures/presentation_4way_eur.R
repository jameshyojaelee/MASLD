#!/usr/bin/env Rscript
# NOTE (Phase 10 swap, 2026-05-06): This script intentionally loads the FROZEN
# UKBB v1 sghatan snapshot as a sensitivity-comparison panel. v1 is no longer
# the production EUR LD reference — PolyFun replaced it (concordance r=0.986).
# v1 archive: archive/ld_panel_v1_sghatan_2026-05-06/
#
# presentation_4way_eur.R
# 4-way EUR LD-panel comparison for COLOC (gene-level) and standalone fine-mapping
# (per-variant PIP). Adds PolyFun UKBB EUR (~337K, Weissbrod 2020) as a 4th panel
# alongside UKBB v1 (sghatan, ~337K), 1KG EUR (n=379), TOP-LD EUR (n=13,160).
#
# Output: figures/supplementary/figS04_coloc/ld_panel_4way_eur/
#   coloc/scatter_pairwise.pdf      — 6 pairwise scatters (2×3 grid)
#   coloc/hits_bar.pdf              — Hit counts per panel per threshold (4 bars × 3 thresholds)
#   coloc/set_agreement_4way.pdf    — 4-way set agreement
#   coloc/hallmarks.pdf             — Hallmark genes 4-column table
#   fm/scatter_pairwise.pdf         — Per-variant SuSiE PIP, 6 facets
#   fm/pip_hits_bar.pdf             — PIP>0.5 hit counts per panel
#   summary.txt                     — text summary of all numbers

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

FM_DIR  <- file.path(PROJ, "GWAS/finemapping")
OUT_DIR <- file.path(PROJ, "figures/supplementary/figS04_coloc/ld_panel_4way_eur")
for (sub in c("coloc", "fm")) dir.create(file.path(OUT_DIR, sub),
                                          recursive = TRUE, showWarnings = FALSE)

PANEL_COLORS <- c("UKBB v1" = "#0D47A1", "1KG EUR" = "#7B1FA2",
                  "TOP-LD EUR" = "#C2185B", "PolyFun" = "#42A5F5")

# ==========================================================================
# COLOC (gene-level) 4-way
# ==========================================================================
load_coloc <- function(subdir) {
  root <- file.path(FM_DIR, "results", subdir)
  if (!dir.exists(root)) {
    cat("  WARNING: ", root, " missing — skipping\n")
    return(NULL)
  }
  files <- list.files(root, pattern = "susie_coloc_chr[0-9]+\\.csv$",
                      recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) return(NULL)
  rows <- lapply(files, function(f) {
    dt <- try(fread(f), silent = TRUE)
    if (inherits(dt, "try-error") || nrow(dt) == 0) return(NULL)
    dt
  })
  dt <- rbindlist(rows, fill = TRUE)
# PP.H4.susie is the LD-DEPENDENT comparable metric across panels.
  # PP.H4.abf is LD-independent (coloc.abf takes only marginal sumstats),
  # so cross-panel ABF comparisons are tautologically identical.
  dt[, max_pp4 := PP.H4.susie]
  dt[is.infinite(max_pp4), max_pp4 := NA]
  # Dedupe duplicate (gwas_name, ensembl) keys (~0.05% of rows from genes
  # appearing in multiple chr CSVs). Without this, the multi-way merge
  # downstream produces Cartesian-product mispairings.
  dt <- dt[order(-max_pp4, na.last = TRUE)]
  dt <- unique(dt, by = c("gwas_name", "ensembl"))
  dt[]
}

cat("\n=== Loading COLOC outputs (4 panels) ===\n")
v1      <- load_coloc("../../archive/ld_panel_v1_sghatan_2026-05-06/susie_coloc_v1_2026-04-21")
kg      <- load_coloc("susie_coloc_1kg")
topld   <- load_coloc("susie_coloc_topld")
polyfun <- load_coloc("susie_coloc_polyfun")
for (nm in c("v1","kg","topld","polyfun"))
  cat(sprintf("  %-8s : %s rows\n", nm, if (is.null(get(nm))) "MISSING" else format(nrow(get(nm)))))

if (is.null(v1) || is.null(kg) || is.null(topld) || is.null(polyfun)) {
  cat("\n[WARN] Missing one of v1/1kg/topld/polyfun COLOC outputs — skipping COLOC comparison.\n")
} else {
  # Restrict to EUR GWAS (1KG panel is EUR-only by construction)
  eur_gwas <- unique(kg$gwas_name)
  v1      <- v1[gwas_name %in% eur_gwas]
  topld   <- topld[gwas_name %in% eur_gwas]
  polyfun <- polyfun[gwas_name %in% eur_gwas]

  m_coloc <- Reduce(function(a,b) merge(a, b, by = c("gwas_name","ensembl"), all = FALSE),
    list(
      v1[,      .(gwas_name, ensembl, gene, pp4_v1      = max_pp4)],
      kg[,      .(gwas_name, ensembl,       pp4_1kg     = max_pp4)],
      topld[,   .(gwas_name, ensembl,       pp4_topld   = max_pp4)],
      polyfun[, .(gwas_name, ensembl,       pp4_polyfun = max_pp4)]
    ))
  m_coloc <- m_coloc[!is.na(pp4_v1) & !is.na(pp4_1kg) & !is.na(pp4_topld) & !is.na(pp4_polyfun)]
  cat(sprintf("  4-way paired rows (PP.H4.susie non-NA in all): %d\n", nrow(m_coloc)))

  pair_specs <- list(
    list("UKBB v1 vs 1KG",     "pp4_v1",    "pp4_1kg"),
    list("UKBB v1 vs TOP-LD",  "pp4_v1",    "pp4_topld"),
    list("UKBB v1 vs PolyFun", "pp4_v1",    "pp4_polyfun"),
    list("1KG vs TOP-LD",      "pp4_1kg",   "pp4_topld"),
    list("1KG vs PolyFun",     "pp4_1kg",   "pp4_polyfun"),
    list("TOP-LD vs PolyFun",  "pp4_topld", "pp4_polyfun")
  )
  for (p in pair_specs) {
    r  <- cor(m_coloc[[p[[2]]]], m_coloc[[p[[3]]]])
    rho <- cor(m_coloc[[p[[2]]]], m_coloc[[p[[3]]]], method = "spearman")
    cat(sprintf("  %-25s r = %.3f   rho = %.3f\n", p[[1]], r, rho))
  }

  # Scatter labels: three canonical MASLD loci (THRB, GCKR, ADH4) plus
  # disagreement-quadrant exemplars — NEAT1 / TCEA2 illustrate top-left
  # (panel-specific high signal), IRF5 / LEAP2 illustrate bottom-right
  # (claimed by one panel, lost in others).
  scatter_hallmark <- c("THRB", "GCKR", "ADH4", "RORA", "HKDC1",
                        "TCEA2", "IRF5", "LEAP2")
  hallmark <- c("RORA","THRB","HKDC1","MTTP","GCKR","CFLAR","EPHA2",
                "TNFSF10","SLC39A8","ADH4","DGAT2","LPL","PNPLA3","TM6SF2",
                "APOE","HSD17B13","MARC1")

  # ---- 6-panel scatter (2 rows × 3 cols)
  panel_long <- function(name, x_col, y_col) {
    d <- m_coloc[, .(gene, x = get(x_col), y = get(y_col))]
    r <- cor(d$x, d$y)
    d[, comparison := sprintf("%s   (r = %.3f)", name, r)]
    d
  }
  long_coloc <- rbindlist(lapply(pair_specs, function(p) panel_long(p[[1]], p[[2]], p[[3]])))
  long_coloc[, comparison := factor(comparison, levels = unique(comparison))]
  lbl_long <- long_coloc[gene %in% scatter_hallmark & pmax(x, y, na.rm = TRUE) > 0.5]
  lbl_long <- unique(lbl_long, by = c("gene", "comparison"))

  p_coloc_scatter <- ggplot(long_coloc, aes(x = x, y = y)) +
    geom_abline(slope = 1, intercept = 0,
                color = "#B0B0B0", linetype = "dashed", linewidth = 0.8) +
    rasterise(geom_point(alpha = 0.25, size = 0.7, color = "#1565C0"), dpi = 220) +
    geom_hline(yintercept = 0.5, color = "grey60", linetype = "dotted", linewidth = 0.4) +
    geom_vline(xintercept = 0.5, color = "grey60", linetype = "dotted", linewidth = 0.4) +
    geom_point(data = lbl_long, color = "black", size = 2.2) +
    geom_text_repel(data = lbl_long, aes(label = gene),
                    size = 4.2, color = "black", fontface = "bold",
                    box.padding = 0.4, max.overlaps = 25,
                    segment.color = "grey30", segment.size = 0.3) +
    facet_wrap(~ comparison, ncol = 3) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "PP.H4.susie (panel on x)", y = "PP.H4.susie (panel on y)",
         title = "Gene-level SuSiE-COLOC concordance — 4 EUR LD panels (6 pairwise)") +
    theme_masld(base_size = 16) +
    theme(plot.title = element_text(size = 22, face = "bold"),
          strip.text = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 18),
          axis.text  = element_text(size = 14))
  ggsave(file.path(OUT_DIR, "coloc/scatter_pairwise.pdf"), p_coloc_scatter,
         width = 18, height = 12, device = cairo_pdf)
  cat("  Wrote: coloc/scatter_pairwise.pdf\n")

  # ---- Hits bar (4 panels × 3 thresholds, SuSiE PP.H4)
  hits_coloc <- data.table(
    panel = factor(rep(c("UKBB v1","1KG EUR","TOP-LD EUR","PolyFun"), each = 3),
                   levels = c("UKBB v1","1KG EUR","TOP-LD EUR","PolyFun")),
    threshold = factor(rep(c("PP.H4.susie > 0.5","PP.H4.susie > 0.8","PP.H4.susie > 0.9"), 4),
                       levels = c("PP.H4.susie > 0.5","PP.H4.susie > 0.8","PP.H4.susie > 0.9")),
    n = c(
      sum(m_coloc$pp4_v1      > 0.5), sum(m_coloc$pp4_v1      > 0.8), sum(m_coloc$pp4_v1      > 0.9),
      sum(m_coloc$pp4_1kg     > 0.5), sum(m_coloc$pp4_1kg     > 0.8), sum(m_coloc$pp4_1kg     > 0.9),
      sum(m_coloc$pp4_topld   > 0.5), sum(m_coloc$pp4_topld   > 0.8), sum(m_coloc$pp4_topld   > 0.9),
      sum(m_coloc$pp4_polyfun > 0.5), sum(m_coloc$pp4_polyfun > 0.8), sum(m_coloc$pp4_polyfun > 0.9)
    )
  )
  p_hits <- ggplot(hits_coloc, aes(x = panel, y = n, fill = panel)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = n), vjust = -0.35, size = 5.5, fontface = "bold") +
    facet_wrap(~ threshold, ncol = 3, scales = "free_y") +
    scale_fill_manual(values = PANEL_COLORS, guide = "none") +
    labs(x = NULL, y = "# gene × GWAS",
         title = "SuSiE-COLOC hit counts at increasing PP.H4 thresholds") +
    theme_masld(base_size = 16) +
    theme(plot.title = element_text(size = 22, face = "bold"),
          strip.text = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size = 13, angle = 25, hjust = 1),
          axis.text.y = element_text(size = 13)) +
    expand_limits(y = max(hits_coloc$n) * 1.18)
  ggsave(file.path(OUT_DIR, "coloc/hits_bar.pdf"), p_hits,
         width = 16, height = 7, device = cairo_pdf)
  cat("  Wrote: coloc/hits_bar.pdf\n")

  # ---- 4-way set agreement at SuSiE PP.H4>0.5
  hi <- m_coloc[, .(v = pp4_v1 > 0.5, k = pp4_1kg > 0.5,
                    t = pp4_topld > 0.5, p = pp4_polyfun > 0.5)]
  set_label <- function(v,k,t,p) {
    s <- c(if (v) "V" else NULL, if (k) "K" else NULL,
           if (t) "T" else NULL, if (p) "P" else NULL)
    paste(s, collapse = "+")
  }
  hi[, lbl := mapply(set_label, v, k, t, p)]
  set_counts <- hi[lbl != "", .N, by = lbl][order(-N)]
  set_counts[, lbl := factor(lbl, levels = lbl)]
  p_set <- ggplot(set_counts, aes(x = lbl, y = N, fill = nchar(as.character(lbl)))) +
    geom_col(width = 0.65) +
    geom_text(aes(label = N), vjust = -0.3, size = 5, fontface = "bold") +
    scale_fill_gradient(low = "#90CAF9", high = "#0D47A1", guide = "none") +
    labs(x = "Subset (V=UKBB v1, K=1KG, T=TOP-LD, P=PolyFun) at SuSiE PP.H4>0.5",
         y = "# gene × GWAS",
         title = "4-way set agreement at SuSiE PP.H4 > 0.5") +
    theme_masld(base_size = 16) +
    theme(plot.title = element_text(size = 22, face = "bold"),
          axis.text.x = element_text(size = 11, angle = 30, hjust = 1)) +
    expand_limits(y = max(set_counts$N) * 1.18)
  ggsave(file.path(OUT_DIR, "coloc/set_agreement_4way.pdf"), p_set,
         width = 13, height = 7, device = cairo_pdf)
  cat("  Wrote: coloc/set_agreement_4way.pdf\n")

  # ---- Hallmark genes 4-column table
  hl_tbl <- m_coloc[gene %in% hallmark][order(-pmax(pp4_v1, pp4_1kg, pp4_topld, pp4_polyfun))]
  hl_tbl <- unique(hl_tbl, by = "gene")[1:min(15, .N)]
  hl_tbl <- hl_tbl[!is.na(gene),
                   .(gene,
                     UKBB_v1 = sprintf("%.3f", pp4_v1),
                     `1KG`   = sprintf("%.3f", pp4_1kg),
                     TOPLD   = sprintf("%.3f", pp4_topld),
                     PolyFun = sprintf("%.3f", pp4_polyfun))]
  N_HL <- nrow(hl_tbl)
  hl_tbl[, row_y := rev(seq_len(N_HL))]
  p_hl <- ggplot(hl_tbl) +
    geom_text(aes(x = 0,    y = row_y, label = gene),    hjust = 0, size = 6, fontface = "bold") +
    geom_text(aes(x = 1.4,  y = row_y, label = UKBB_v1), hjust = 1, size = 6) +
    geom_text(aes(x = 2.6,  y = row_y, label = `1KG`),   hjust = 1, size = 6) +
    geom_text(aes(x = 3.8,  y = row_y, label = TOPLD),   hjust = 1, size = 6) +
    geom_text(aes(x = 5.0,  y = row_y, label = PolyFun), hjust = 1, size = 6) +
    annotate("text", x = c(0, 1.4, 2.6, 3.8, 5.0), y = N_HL + 0.8,
             label = c("Gene","UKBB v1","1KG","TOPLD","PolyFun"),
             hjust = c(0,1,1,1,1), size = 7, fontface = "bold") +
    scale_x_continuous(limits = c(-0.1, 5.2)) +
    scale_y_continuous(limits = c(0.5, N_HL + 1.3)) +
    labs(title = "Hallmark genes — SuSiE PP.H4 head-to-head (4 panels)") +
    theme_void(base_size = 18) +
    theme(plot.title = element_text(size = 22, face = "bold"),
          plot.margin = margin(15,15,15,15))
  ggsave(file.path(OUT_DIR, "coloc/hallmarks.pdf"), p_hl,
         width = 12, height = 9, device = cairo_pdf)
  cat("  Wrote: coloc/hallmarks.pdf\n")
}

# ==========================================================================
# Standalone fine-mapping (per-variant PIP) 4-way
# ==========================================================================
cat("\n=== Loading standalone FM outputs (4 panels) ===\n")
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
    dt[, study := basename(dirname(dirname(dirname(f))))]
    dt[, locus_id := gsub("_cov0\\.95_.*", "", basename(f))]
    dt[, converged_flag := !grepl("_notconverged", basename(f))]
    # Allele orientation differs across panels (e.g. v1 emits "19:18207397:A:C"
    # while 1KG/TOP-LD emit "19:18207397:C:A" for the same variant). Build a
    # strand-agnostic key by sorting the two alleles. Without this, the 4-way
    # SNP merge drops virtually every high-PIP variant.
    parts <- tstrsplit(dt$SNP, ":", fixed = TRUE, type.convert = FALSE)
    a_lex <- pmin(parts[[3]], parts[[4]])
    b_lex <- pmax(parts[[3]], parts[[4]])
    dt[, snp_key := paste(parts[[1]], parts[[2]], a_lex, b_lex, sep = ":")]
    dt
  })
  rbindlist(rows, fill = TRUE)
}

fm_v1      <- load_fm("EUR_0.5Mb")
fm_kg      <- load_fm("1kg_eur_0.5Mb")
fm_topld   <- load_fm("topld_eur_0.5Mb")
fm_polyfun <- load_fm("polyfun_eur_0.5Mb")

for (nm in c("fm_v1","fm_kg","fm_topld","fm_polyfun"))
  cat(sprintf("  %-12s : %s rows\n", nm, if (is.null(get(nm))) "MISSING" else format(nrow(get(nm)))))

if (is.null(fm_polyfun)) {
  cat("\n[WARN] PolyFun FM results missing — skipping FM 4-way (run 8e first).\n")
} else if (!is.null(fm_v1) && !is.null(fm_kg) && !is.null(fm_topld)) {
  # Within each panel, dedupe to the max-PIP row per (study, locus_id, snp_key)
  # before merging — sorted-allele keys can collapse multi-allelic rows.
  dedupe_fm <- function(dt) {
    dt <- dt[converged_flag == TRUE, .(study, locus_id, snp_key, PIP)]
    dt <- dt[order(-PIP)]
    unique(dt, by = c("study","locus_id","snp_key"))
  }
  m_fm <- Reduce(function(a, b) merge(a, b, by = c("study","locus_id","snp_key"), all = FALSE),
    list(
      setnames(dedupe_fm(fm_v1),      "PIP", "PIP_v1"),
      setnames(dedupe_fm(fm_kg),      "PIP", "PIP_1kg"),
      setnames(dedupe_fm(fm_topld),   "PIP", "PIP_topld"),
      setnames(dedupe_fm(fm_polyfun), "PIP", "PIP_polyfun")
    ))
  m_fm <- m_fm[!is.na(PIP_v1) & !is.na(PIP_1kg) & !is.na(PIP_topld) & !is.na(PIP_polyfun)]
  cat(sprintf("  4-way paired variants: %d\n", nrow(m_fm)))

  fm_pairs <- list(
    list("UKBB v1 vs 1KG",     "PIP_v1",    "PIP_1kg"),
    list("UKBB v1 vs TOP-LD",  "PIP_v1",    "PIP_topld"),
    list("UKBB v1 vs PolyFun", "PIP_v1",    "PIP_polyfun"),
    list("1KG vs TOP-LD",      "PIP_1kg",   "PIP_topld"),
    list("1KG vs PolyFun",     "PIP_1kg",   "PIP_polyfun"),
    list("TOP-LD vs PolyFun",  "PIP_topld", "PIP_polyfun")
  )
  fm_long <- rbindlist(lapply(fm_pairs, function(p) {
    r <- cor(m_fm[[p[[2]]]], m_fm[[p[[3]]]])
    data.table(comparison = sprintf("%s   (r = %.3f)", p[[1]], r),
               x = m_fm[[p[[2]]]], y = m_fm[[p[[3]]]])
  }))
  fm_long[, comparison := factor(comparison, levels = unique(comparison))]

  p_fm <- ggplot(fm_long, aes(x = x, y = y)) +
    geom_abline(slope = 1, intercept = 0,
                color = "#B0B0B0", linetype = "dashed", linewidth = 0.8) +
    rasterise(geom_point(alpha = 0.25, size = 0.7, color = "#1565C0"), dpi = 220) +
    geom_hline(yintercept = 0.5, color = "grey60", linetype = "dotted") +
    geom_vline(xintercept = 0.5, color = "grey60", linetype = "dotted") +
    facet_wrap(~ comparison, ncol = 3) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "PIP (panel on x)", y = "PIP (panel on y)",
         title = "Per-variant SuSiE PIP concordance — 4 EUR LD panels (6 pairwise)") +
    theme_masld(base_size = 16) +
    theme(plot.title = element_text(size = 22, face = "bold"),
          strip.text = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 18))
  ggsave(file.path(OUT_DIR, "fm/scatter_pairwise.pdf"), p_fm,
         width = 18, height = 12, device = cairo_pdf)
  cat("  Wrote: fm/scatter_pairwise.pdf\n")

  # PIP > 0.5 hits bar
  hits_fm <- data.table(
    panel = factor(c("UKBB v1","1KG EUR","TOP-LD EUR","PolyFun"),
                   levels = c("UKBB v1","1KG EUR","TOP-LD EUR","PolyFun")),
    n = c(sum(m_fm$PIP_v1 > 0.5), sum(m_fm$PIP_1kg > 0.5),
          sum(m_fm$PIP_topld > 0.5), sum(m_fm$PIP_polyfun > 0.5))
  )
  p_fm_hits <- ggplot(hits_fm, aes(x = panel, y = n, fill = panel)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = n), vjust = -0.35, size = 6, fontface = "bold") +
    scale_fill_manual(values = PANEL_COLORS, guide = "none") +
    labs(x = NULL, y = "# variants with PIP > 0.5",
         title = "PIP > 0.5 variant counts per LD panel") +
    theme_masld(base_size = 18) +
    theme(plot.title = element_text(size = 22, face = "bold")) +
    expand_limits(y = max(hits_fm$n) * 1.18)
  ggsave(file.path(OUT_DIR, "fm/pip_hits_bar.pdf"), p_fm_hits,
         width = 10, height = 7, device = cairo_pdf)
  cat("  Wrote: fm/pip_hits_bar.pdf\n")
}

# ==========================================================================
# Summary text
# ==========================================================================
sink(file.path(OUT_DIR, "summary.txt"))
cat(sprintf("4-way EUR LD-panel comparison summary  (%s)\n\n", Sys.time()))
if (exists("m_coloc")) {
  cat("=== SuSiE-COLOC (gene-level, LD-DEPENDENT metric) ===\n")
  cat(sprintf("4-way paired rows (PP.H4.susie non-NA in v1, 1KG, TOP-LD, PolyFun): %d\n\n", nrow(m_coloc)))
  for (p in pair_specs) {
    r <- cor(m_coloc[[p[[2]]]], m_coloc[[p[[3]]]])
    rho <- cor(m_coloc[[p[[2]]]], m_coloc[[p[[3]]]], method = "spearman")
    cat(sprintf("  %-25s r = %.3f   rho = %.3f\n", p[[1]], r, rho))
  }
  cat("\nHit counts (SuSiE PP.H4 > 0.5 / >0.8 / >0.9):\n")
  cat(sprintf("  UKBB v1:  %5d / %5d / %5d\n", sum(m_coloc$pp4_v1 > 0.5), sum(m_coloc$pp4_v1 > 0.8), sum(m_coloc$pp4_v1 > 0.9)))
  cat(sprintf("  1KG EUR:  %5d / %5d / %5d\n", sum(m_coloc$pp4_1kg > 0.5), sum(m_coloc$pp4_1kg > 0.8), sum(m_coloc$pp4_1kg > 0.9)))
  cat(sprintf("  TOP-LD :  %5d / %5d / %5d\n", sum(m_coloc$pp4_topld > 0.5), sum(m_coloc$pp4_topld > 0.8), sum(m_coloc$pp4_topld > 0.9)))
  cat(sprintf("  PolyFun:  %5d / %5d / %5d\n", sum(m_coloc$pp4_polyfun > 0.5), sum(m_coloc$pp4_polyfun > 0.8), sum(m_coloc$pp4_polyfun > 0.9)))
}
if (exists("m_fm")) {
  cat("\n=== Per-variant FM ===\n")
  cat(sprintf("4-way paired variants: %d\n", nrow(m_fm)))
  for (p in fm_pairs) {
    r <- cor(m_fm[[p[[2]]]], m_fm[[p[[3]]]])
    cat(sprintf("  %-25s r = %.3f\n", p[[1]], r))
  }
}
sink()
cat("\nWrote: summary.txt\n")
cat("\nALL DONE. Outputs in:", OUT_DIR, "\n")
