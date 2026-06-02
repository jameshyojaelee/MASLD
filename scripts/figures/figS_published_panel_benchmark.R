#!/usr/bin/env Rscript
# =============================================================================
# Gap #9: Head-to-head benchmark of our atlas against published MASLD gene panels
#
# Compares recovery, direction concordance, atlas rank, and unique discoveries
# for Govaere 25, SteatoSITE 15, Piras 685, and Feng 2026 panels
#
# Output:
#   figures/supplementary/figS_sensitivity/figS_published_panel_benchmark.pdf
#   figures/supplementary/figS_sensitivity/published_panel_benchmark.csv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

cat("=== Gap #9: Published Panel Benchmark ===\n")

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
dream <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results_ashr.csv"))
atlas <- fread(file.path(BASE,
  "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"))

cat("Dream results:", nrow(dream), "genes\n")
cat("Atlas:", nrow(atlas), "genes x", ncol(atlas), "columns\n")

# Use symbol from dream; ensure unique
dream[, symbol := as.character(symbol)]
dream <- dream[!is.na(symbol) & symbol != ""]
# Keep one row per symbol (most significant)
dream <- dream[order(padj)][!duplicated(symbol)]

# Atlas uses human_symbol
atlas[, human_symbol := as.character(human_symbol)]
atlas <- atlas[!is.na(human_symbol) & human_symbol != ""]

# Merge dream into atlas-level table for comprehensive annotation
# Atlas already has dream columns, so just use atlas directly
# But dream has the authoritative DEG calls
setkey(dream, symbol)

# DEG definition: padj < 0.05, |logFC| > 0.5
dream[, is_deg := !is.na(padj) & padj < 0.05 & abs(logFC) > 0.5]
n_degs <- sum(dream$is_deg, na.rm = TRUE)
cat("DEGs (padj<0.05, |logFC|>0.5):", n_degs, "\n")

# Rank genes by absolute logFC (higher = better)
dream[, abs_logFC := abs(logFC)]
dream[, logFC_rank := frank(-abs_logFC, ties.method = "average")]
dream[, logFC_percentile := 100 * (1 - logFC_rank / .N)]

# ---------------------------------------------------------------------------
# 2. Define published panels
# ---------------------------------------------------------------------------

# Govaere 2020 — 25-gene fibrosis progression signature (from Script 73)
govaere_25 <- c("AKR1B10", "DUSP6", "GDF15", "THBS2", "A2M", "CDH2",
                "COL1A1", "COL3A1", "COL4A1", "COL4A2", "COL6A3", "DCN",
                "FBN1", "FSTL1", "IGFBP7", "LUM", "MFAP4", "MMP2",
                "POSTN", "SPARC", "SPP1", "TAGLN", "THY1", "TIMP1", "VCAN")

# SteatoSITE 15-gene endotype classifier (Teufel/Raverdy 2024)
steatosite_15 <- c("MT1F", "KCNH7", "COL25A1", "RASD2", "CTGF", "STC1",
                   "GDNF", "PRRX1", "FGF7", "LCNL1", "DPEP1", "CHRDL2",
                   "LHX6", "POU4F1", "CDH16")

# Piras & DiStefano 2024 — key MASLD genes (GWAS-validated + top DEGs)
# Full 685-DEG list not publicly available; use their highlighted genes + GWAS overlap
# We include their top reported genes from the paper
piras_key <- c("HSD17B13", "PNPLA3", "TM6SF2", "MARC1", "MBOAT7",
               "CIDEB", "GPAM", "APOH", "SERPINA6", "GHR",
               "CYP7A1", "ABCB4", "SLC22A1", "SLC27A5", "ACSM3",
               "ADH4", "ADH1B", "CYP2E1", "CYP4A11", "CES1",
               "AKR1D1", "SLC10A1", "HGF", "CRP", "LCAT",
               "CETP", "ALDOB", "PCK1", "GLS2", "HAO1")

# Feng et al. 2026 — top-ranked genes from 6-layer scoring + validated targets
# EFHD1 ranked #1, MLIP validated. Plus their 27-gene core progression signature
feng_top <- c("EFHD1", "MLIP", "TREM2", "SPP1", "GPNMB", "CCL2",
              "CCL20", "CXCL1", "CXCL6", "IL1B", "IL32",
              "COL1A1", "COL1A2", "COL3A1", "FN1", "LOX", "LOXL2",
              "ACTA2", "PDGFRB", "TGFB1", "SERPINE1", "MMP9",
              "AKR1B10", "GDF15", "THY1", "THBS2", "LUM")

# User-specified Govaere panel from prompt (fibrosis-focused, overlapping but different)
govaere_fibrosis <- c("ACTA2", "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2",
                      "CTGF", "CXCL8", "DCN", "FN1", "IL1B", "IL6", "LOX", "LOXL2",
                      "MMP2", "MMP9", "PDGFRB", "SERPINE1", "TGFB1", "TGFB2",
                      "TIMP1", "TIMP2", "TNF", "VIM", "WISP1")

panels <- list(
  "Govaere 25\n(Sci Transl Med 2020)"    = govaere_25,
  "Govaere Fibrosis 25\n(Cell Metab 2020)" = govaere_fibrosis,
  "SteatoSITE 15\n(Nat Med 2024)"         = steatosite_15,
  "Piras Key 30\n(Life Sci Alliance 2024)" = piras_key,
  "Feng Top 27\n(bioRxiv 2026)"           = feng_top
)

# Alias resolution: CTGF -> CCN2, WISP1 -> CCN4 (HGNC updates)
resolve_alias <- function(genes, available) {
  alias_map <- c("CTGF" = "CCN2", "WISP1" = "CCN4", "CXCL8" = "IL8")
  resolved <- ifelse(genes %in% available, genes,
                     ifelse(genes %in% names(alias_map) & alias_map[genes] %in% available,
                            alias_map[genes], genes))
  # Try reverse alias too
  rev_map <- setNames(names(alias_map), alias_map)
  resolved <- ifelse(resolved %in% available, resolved,
                     ifelse(resolved %in% names(rev_map) & rev_map[resolved] %in% available,
                            rev_map[resolved], resolved))
  resolved
}

available_genes <- dream$symbol
for (nm in names(panels)) {
  panels[[nm]] <- resolve_alias(panels[[nm]], available_genes)
}

# ---------------------------------------------------------------------------
# 3. Compute benchmark metrics per panel
# ---------------------------------------------------------------------------
benchmark_rows <- list()

for (panel_name in names(panels)) {
  genes <- panels[[panel_name]]
  n_total <- length(genes)

  # Which genes are in our dream results?
  in_dream <- genes[genes %in% dream$symbol]
  n_found <- length(in_dream)
  missing <- setdiff(genes, dream$symbol)

  if (n_found == 0) {
    benchmark_rows[[panel_name]] <- data.table(
      panel = panel_name, n_panel = n_total, n_found = 0, n_missing = length(missing),
      missing_genes = paste(missing, collapse = ";"),
      n_deg = 0, recovery_pct = 0, direction_concordance = NA_real_,
      median_percentile = NA_real_, mean_percentile = NA_real_,
      n_coloc = 0, n_conserved = 0, n_any_atlas_support = 0
    )
    next
  }

  sub <- dream[symbol %in% in_dream]

  # Recovery: % that are DEGs
  n_deg <- sum(sub$is_deg, na.rm = TRUE)
  recovery_pct <- 100 * n_deg / n_total


  # Direction concordance among those that are DEG:
  # For published MASLD panels, expected direction is UP (fibrosis/inflammation genes)
  # We check if logFC > 0 for DEGs
  deg_sub <- sub[is_deg == TRUE]
  if (nrow(deg_sub) > 0) {
    direction_concordance <- 100 * sum(deg_sub$logFC > 0) / nrow(deg_sub)
  } else {
    direction_concordance <- NA_real_
  }

  # Percentile rank
  median_pctl <- median(sub$logFC_percentile, na.rm = TRUE)
  mean_pctl <- mean(sub$logFC_percentile, na.rm = TRUE)

  # Atlas evidence layers
  atlas_sub <- atlas[human_symbol %in% in_dream]

  # COLOC support (any coloc_susie_best_pp4 > 0.5)
  n_coloc <- 0
  if ("coloc_susie_best_pp4" %in% names(atlas_sub) && nrow(atlas_sub) > 0) {
    n_coloc <- sum(!is.na(atlas_sub$coloc_susie_best_pp4) &
                     atlas_sub$coloc_susie_best_pp4 > 0.5, na.rm = TRUE)
  }

  # Conserved
  n_cc <- 0
  if ("is_conserved" %in% names(atlas_sub) && nrow(atlas_sub) > 0) {
    n_cc <- sum(atlas_sub$is_conserved == TRUE, na.rm = TRUE)
  }

  # Any atlas support (sources_active > 0)
  n_any <- 0
  if ("sources_active" %in% names(atlas_sub) && nrow(atlas_sub) > 0) {
    n_any <- sum(!is.na(atlas_sub$sources_active) &
                   atlas_sub$sources_active > 0, na.rm = TRUE)
  }

  benchmark_rows[[panel_name]] <- data.table(
    panel = panel_name, n_panel = n_total, n_found = n_found,
    n_missing = length(missing),
    missing_genes = paste(missing, collapse = ";"),
    n_deg = n_deg, recovery_pct = recovery_pct,
    direction_concordance = direction_concordance,
    median_percentile = median_pctl, mean_percentile = mean_pctl,
    n_coloc = n_coloc, n_conserved = n_cc,
    n_any_atlas_support = n_any
  )
}

benchmark <- rbindlist(benchmark_rows)
cat("\n=== Panel Benchmark Summary ===\n")
print(benchmark[, .(panel, n_panel, n_found, n_deg, recovery_pct,
                    direction_concordance, median_percentile,
                    n_coloc, n_conserved)])

# ---------------------------------------------------------------------------
# 4. Build per-gene detail table for figure data
# ---------------------------------------------------------------------------
all_panel_genes <- unique(unlist(panels))
gene_detail <- dream[symbol %in% all_panel_genes,
                     .(symbol, logFC, padj, is_deg, logFC_percentile)]

# Add atlas columns
atlas_slim <- atlas[human_symbol %in% all_panel_genes,
                    .(human_symbol, is_conserved, coloc_susie_best_pp4,
                      sources_active, n_coloc_sources)]
gene_detail <- merge(gene_detail, atlas_slim,
                     by.x = "symbol", by.y = "human_symbol", all.x = TRUE)

# Add panel membership columns
for (nm in names(panels)) {
  short_name <- gsub("\n.*", "", nm)
  gene_detail[, (short_name) := symbol %in% panels[[nm]]]
}

# Unique discoveries: our DEGs NOT in any published panel
our_degs <- dream[is_deg == TRUE]$symbol
all_published <- unique(unlist(panels))
unique_discoveries <- setdiff(our_degs, all_published)
cat("\nOur DEGs:", length(our_degs), "\n")
cat("Genes in any published panel:", length(all_published), "\n")
cat("Unique discoveries (our DEGs not in any panel):", length(unique_discoveries), "\n")

# ---------------------------------------------------------------------------
# 5. Panel (a): Recovery heatmap — panels x evidence layers
# ---------------------------------------------------------------------------

# Build matrix: rows = panels, cols = evidence layers
evidence_layers <- c("DEG\n(padj<0.05,|LFC|>0.5)", "COLOC\n(PP4>0.5)",
                     "Conserved\nCore", "Any Atlas\nSupport")

heatmap_data <- data.table()
for (panel_name in names(panels)) {
  genes <- panels[[panel_name]]
  in_dream_genes <- genes[genes %in% dream$symbol]
  n_total <- length(genes)

  # DEG recovery
  deg_n <- sum(dream[symbol %in% in_dream_genes]$is_deg, na.rm = TRUE)

  # COLOC
  a_sub <- atlas[human_symbol %in% in_dream_genes]
  coloc_n <- sum(!is.na(a_sub$coloc_susie_best_pp4) &
                   a_sub$coloc_susie_best_pp4 > 0.5, na.rm = TRUE)

  # CC
  cc_n <- sum(a_sub$is_conserved == TRUE, na.rm = TRUE)

  # Any support
  any_n <- sum(!is.na(a_sub$sources_active) & a_sub$sources_active > 0, na.rm = TRUE)

  heatmap_data <- rbind(heatmap_data, data.table(
    panel = panel_name,
    layer = evidence_layers,
    n = c(deg_n, coloc_n, cc_n, any_n),
    pct = 100 * c(deg_n, coloc_n, cc_n, any_n) / n_total
  ))
}

# Short panel labels for plotting
heatmap_data[, panel_short := gsub("\n.*", "", panel)]
heatmap_data[, layer := factor(layer, levels = evidence_layers)]
heatmap_data[, panel_short := factor(panel_short,
  levels = rev(c("Govaere 25", "Govaere Fibrosis 25", "SteatoSITE 15",
                 "Piras Key 30", "Feng Top 27")))]

p_a <- ggplot(heatmap_data, aes(x = layer, y = panel_short, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%d%%\n(%d/%d)",
                                round(pct), n,
                                as.integer(round(n * 100 / pct)))),
            size = 2, color = "black", lineheight = 0.85) +
  scale_fill_gradient(low = "white", high = masld_colors$up,
                      limits = c(0, 100), name = "Recovery (%)") +
  labs(x = NULL, y = NULL,
       title = "Recovery of published panels in our atlas") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 6),
        legend.position = "right")

# Fix the geom_text labels for 0% (avoid division by zero)
heatmap_data[, n_panel := as.integer(round(n * 100 / ifelse(pct == 0, 1, pct)))]

p_a <- ggplot(heatmap_data, aes(x = layer, y = panel_short, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%d%%\n(%d/%d)", round(pct), n, n_panel)),
            size = 2, color = ifelse(heatmap_data$pct > 60, "white", "black"),
            lineheight = 0.85) +
  scale_fill_gradient(low = "white", high = masld_colors$up,
                      limits = c(0, 100), name = "Recovery (%)") +
  labs(x = NULL, y = NULL,
       title = "Recovery of published panels in our atlas") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 6),
        legend.position = "right")

# ---------------------------------------------------------------------------
# 6. Panel (b): Rank distribution — violin/box of atlas percentiles
# ---------------------------------------------------------------------------

# Build data for each panel + random background
rank_data <- data.table()
for (panel_name in names(panels)) {
  genes <- panels[[panel_name]]
  in_dream_genes <- genes[genes %in% dream$symbol]
  sub <- dream[symbol %in% in_dream_genes]
  if (nrow(sub) > 0) {
    rank_data <- rbind(rank_data, data.table(
      panel = gsub("\n.*", "", panel_name),
      percentile = sub$logFC_percentile
    ))
  }
}

# Random background: 1000 random genes sampled 100 times
set.seed(42)
random_pctls <- numeric()
for (i in 1:100) {
  samp <- dream[sample(.N, 30)]  # match typical panel size
  random_pctls <- c(random_pctls, samp$logFC_percentile)
}
rank_data <- rbind(rank_data, data.table(
  panel = "Random\nBackground",
  percentile = random_pctls
))

rank_data[, panel := factor(panel,
  levels = c("Govaere 25", "Govaere Fibrosis 25", "SteatoSITE 15",
             "Piras Key 30", "Feng Top 27", "Random\nBackground"))]

# Compute medians for annotation
medians <- rank_data[, .(med = median(percentile)), by = panel]

# Wilcoxon test vs random
random_vals <- rank_data[panel == "Random\nBackground"]$percentile
p_vals <- rank_data[panel != "Random\nBackground",
                    .(p = wilcox.test(percentile, random_vals,
                                     alternative = "greater")$p.value),
                    by = panel]

p_b <- ggplot(rank_data, aes(x = panel, y = percentile, fill = panel)) +
  geom_violin(alpha = 0.6, scale = "width", linewidth = 0.3) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5, linewidth = 0.3) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_fill_manual(values = c(
    "Govaere 25" = masld_colors$up,
    "Govaere Fibrosis 25" = masld_colors$fibrosis,
    "SteatoSITE 15" = masld_colors$twas,
    "Piras Key 30" = masld_colors$down,
    "Feng Top 27" = "#7B1FA2",
    "Random\nBackground" = masld_colors$ns
  )) +
  # Add p-value annotations
  geom_text(data = p_vals, aes(x = panel, y = 102,
            label = ifelse(p < 0.001, "***",
                    ifelse(p < 0.01, "**",
                    ifelse(p < 0.05, "*", "ns")))),
            inherit.aes = FALSE, size = 2.5, vjust = 0) +
  labs(x = NULL, y = "|logFC| percentile rank",
       title = "Atlas rank distribution of published panel genes") +
  theme_masld() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1, size = 6))

# ---------------------------------------------------------------------------
# 7. Panel (c): UpSet-style overlap — atlas DEGs vs all panels
# ---------------------------------------------------------------------------

# For UpSet, use simple membership columns
upset_genes <- unique(c(our_degs, all_published))
upset_dt <- data.table(gene = upset_genes)
upset_dt[, Our_DEGs := gene %in% our_degs]
for (nm in names(panels)) {
  short <- gsub("\n.*", "", nm)
  upset_dt[, (short) := gene %in% panels[[nm]]]
}

# Since UpSet plots need UpSetR or ComplexHeatmap, and we want to keep
# dependencies minimal, use a grouped bar chart of overlap counts instead
# showing: total in panel, recovered as DEG, unique to us

overlap_summary <- data.table()
for (nm in names(panels)) {
  genes <- panels[[nm]]
  short <- gsub("\n.*", "", nm)
  in_both <- sum(genes %in% our_degs)
  in_panel_only <- sum(!genes %in% our_degs & genes %in% dream$symbol)
  not_found <- sum(!genes %in% dream$symbol)
  overlap_summary <- rbind(overlap_summary, data.table(
    panel = short,
    category = c("Recovered as DEG", "In atlas (not DEG)", "Not in atlas"),
    count = c(in_both, in_panel_only, not_found)
  ))
}

# Add our unique count
overlap_summary <- rbind(overlap_summary, data.table(
  panel = "Our Atlas",
  category = c("Unique DEGs", "Shared with panels", ""),
  count = c(length(unique_discoveries),
            length(intersect(our_degs, all_published)), 0)
))

overlap_summary[, panel := factor(panel,
  levels = c("Govaere 25", "Govaere Fibrosis 25", "SteatoSITE 15",
             "Piras Key 30", "Feng Top 27", "Our Atlas"))]
overlap_summary <- overlap_summary[count > 0]
overlap_summary[, category := factor(category,
  levels = c("Recovered as DEG", "In atlas (not DEG)", "Not in atlas",
             "Unique DEGs", "Shared with panels"))]

cat_colors <- c(
  "Recovered as DEG"   = masld_colors$up,
  "In atlas (not DEG)" = masld_colors$ns,
  "Not in atlas"       = "#E0E0E0",
  "Unique DEGs"        = masld_colors$down,
  "Shared with panels" = masld_colors$twas
)

p_c <- ggplot(overlap_summary, aes(x = panel, y = count, fill = category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_text(aes(label = ifelse(count > 0, count, "")),
            position = position_stack(vjust = 0.5), size = 2, color = "white") +
  scale_fill_manual(values = cat_colors, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Gene count",
       title = "Panel overlap with our atlas DEGs") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
        legend.position = "bottom",
        legend.text = element_text(size = 5))

# ---------------------------------------------------------------------------
# 8. Assemble and save
# ---------------------------------------------------------------------------
fig <- p_a / (p_b | p_c) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

out_dir <- FIGS_SENS_DIR
dir.create(file.path(out_dir, "panels"), showWarnings = FALSE, recursive = TRUE)

out_pdf <- file.path(out_dir, "figS_published_panel_benchmark.pdf")
save_fig_tall(fig, out_pdf, width = fig_full_width, height = 7)
cat("\nSaved figure:", out_pdf, "\n")

# Save companion CSV
out_csv <- file.path(out_dir, "published_panel_benchmark.csv")
# Clean panel names for CSV
benchmark[, panel := gsub("\n", " ", panel)]
fwrite(benchmark, out_csv)
cat("Saved CSV:", out_csv, "\n")

# Also save per-gene detail
out_gene <- file.path(out_dir, "published_panel_benchmark_genes.csv")
fwrite(gene_detail, out_gene)
cat("Saved per-gene detail:", out_gene, "\n")

# ---------------------------------------------------------------------------
# 9. Summary statistics for manuscript
# ---------------------------------------------------------------------------
cat("\n=== KEY NUMBERS FOR MANUSCRIPT ===\n")
for (i in seq_len(nrow(benchmark))) {
  row <- benchmark[i]
  cat(sprintf("%s: %d/%d recovered (%.0f%%), median pctl=%.1f, %d COLOC, %d CC\n",
              row$panel, row$n_deg, row$n_panel, row$recovery_pct,
              row$median_percentile, row$n_coloc, row$n_conserved))
}
cat(sprintf("\nUnique discoveries (our DEGs not in any panel): %d / %d (%.1f%%)\n",
            length(unique_discoveries), length(our_degs),
            100 * length(unique_discoveries) / length(our_degs)))

# Enrichment test: are published panel genes enriched among our DEGs?
all_in_dream <- nrow(dream)
n_published_in_dream <- sum(all_published %in% dream$symbol)
n_published_deg <- sum(all_published %in% our_degs)
expected <- n_published_in_dream * n_degs / all_in_dream
or <- (n_published_deg / (n_published_in_dream - n_published_deg)) /
      (n_degs / (all_in_dream - n_degs))
ft <- fisher.test(matrix(c(
  n_published_deg,
  n_published_in_dream - n_published_deg,
  n_degs - n_published_deg,
  all_in_dream - n_published_in_dream - n_degs + n_published_deg
), nrow = 2))
cat(sprintf("\nFisher enrichment of published panel genes among our DEGs:\n"))
cat(sprintf("  Observed: %d/%d (%.1f%%), Expected: %.0f, OR=%.2f, p=%.2e\n",
            n_published_deg, n_published_in_dream,
            100 * n_published_deg / n_published_in_dream,
            expected, ft$estimate, ft$p.value))

cat("\n=== DONE ===\n")
