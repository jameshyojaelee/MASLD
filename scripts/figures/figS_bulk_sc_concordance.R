#!/usr/bin/env Rscript
# figS_bulk_sc_concordance.R
#
# Four-panel figure demonstrating bulk RNA-seq vs scRNA-seq pseudobulk concordance:
#   A: logFC scatter, all shared genes (baseline ρ — the "low" number)
#   B: logFC scatter, hepatocyte-intrinsic DEGs only (deconvolution-attributed subset)
#   C: Direction concordance vs bulk DEG stringency threshold
#   D: fgsea — bulk up/down DEG gene sets enriched in sc hepatocyte t-stat ranking
#
# Output: figures/supplementary/figS03_deconvolution/figS_bulk_sc_concordance.pdf
# Environment: rnaseq

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(fgsea)
  library(scales)
})

source("scripts/figures/publication_theme.R")
source("scripts/figures/load_figure_data.R")

SC_DE_DIR   <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudobulk_de")
DECONV_PATH <- file.path(CAUSAL, "deconv_attribution_scores.csv")
DREAM_PATH  <- file.path(INT_RESULTS, "dream_results_ashr.csv")

# ---------------------------------------------------------------------------
# 1. Load & harmonise gene IDs
# ---------------------------------------------------------------------------
message("Loading bulk dream results...")
bulk <- fread(DREAM_PATH, select = c("gene", "logFC", "t", "padj", "symbol"))
bulk[, ensg_base := sub("\\.\\d+$", "", gene)]
setnames(bulk, c("logFC", "t", "padj"), c("bulk_lfc", "bulk_t", "bulk_padj"))

message("Loading sc hepatocyte pseudobulk DE...")
sc <- fread(file.path(SC_DE_DIR, "Hepatocytes_de.csv"),
            select = c("gene", "logFC", "t_stat", "padj"))
sc[, ensg_base := sub("\\.\\d+$", "", gene)]
setnames(sc, c("logFC", "t_stat", "padj"), c("sc_lfc", "sc_t", "sc_padj"))

merged <- merge(
  bulk[, .(ensg_base, symbol, bulk_lfc, bulk_t, bulk_padj)],
  sc[,   .(ensg_base, sc_lfc, sc_t, sc_padj)],
  by = "ensg_base"
)
merged <- merged[!is.na(bulk_lfc) & !is.na(sc_lfc)]
message(sprintf("  Shared genes: %d", nrow(merged)))

message("Loading deconvolution attribution...")
deconv <- fread(DECONV_PATH)
deconv[, ensg_base := sub("\\.\\d+$", "", gene)]
hep_intrinsic_ids <- deconv[category == "Hepatocyte_intrinsic", ensg_base]
message(sprintf("  Hepatocyte-intrinsic: %d", length(hep_intrinsic_ids)))

# ---------------------------------------------------------------------------
# 2. Panel A — all genes scatter
# ---------------------------------------------------------------------------
rho_all <- round(cor(merged$bulk_lfc, merged$sc_lfc,
                     method = "spearman", use = "complete.obs"), 3)
n_all   <- nrow(merged)

# Subsample for rendering speed; keep all for ρ computation
set.seed(42)
plt_all <- merged[sample(.N, min(.N, 8000))]

pA <- ggplot(plt_all, aes(bulk_lfc, sc_lfc)) +
  geom_point(alpha = 0.12, size = 0.35, color = "#4472C4") +
  geom_smooth(method = "lm", color = "#C0392B", linewidth = 0.7, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey50") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
           label = sprintf("rho = %.3f\nn = %s", rho_all, format(n_all, big.mark = ",")),
           size = 2.5, fontface = "italic") +
  labs(x = "Bulk RNA-seq logFC\n(MASLD vs Healthy)",
       y = "sc Hepatocyte pseudobulk logFC\n(MASLD vs Healthy)",
       title = "All shared genes") +
  theme_masld()

# ---------------------------------------------------------------------------
# 3. Panel B — hepatocyte-intrinsic subset
# ---------------------------------------------------------------------------
merged_hi <- merged[ensg_base %in% hep_intrinsic_ids]
rho_hi    <- round(cor(merged_hi$bulk_lfc, merged_hi$sc_lfc,
                       method = "spearman", use = "complete.obs"), 3)
n_hi      <- nrow(merged_hi)
message(sprintf("  Hepatocyte-intrinsic shared: %d  ρ=%.3f", n_hi, rho_hi))

pB <- ggplot(merged_hi, aes(bulk_lfc, sc_lfc)) +
  geom_point(alpha = 0.4, size = 0.8, color = "#27AE60") +
  geom_smooth(method = "lm", color = "#C0392B", linewidth = 0.7, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey50") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
           label = sprintf("rho = %.3f\nn = %s", rho_hi, format(n_hi, big.mark = ",")),
           size = 2.5, fontface = "italic") +
  labs(x = "Bulk RNA-seq logFC\n(MASLD vs Healthy)",
       y = "sc Hepatocyte pseudobulk logFC\n(MASLD vs Healthy)",
       title = "Hepatocyte-intrinsic genes\n(deconvolution-attributed)") +
  theme_masld()

# ---------------------------------------------------------------------------
# 4. Panel C — direction concordance vs stringency
# ---------------------------------------------------------------------------
thresholds <- list(
  list(label = "All tested\ngenes", data = merged),
  list(label = "Bulk\npadj<0.20", data = merged[bulk_padj < 0.20]),
  list(label = "Bulk\npadj<0.10", data = merged[bulk_padj < 0.10]),
  list(label = "Bulk DEGs\n(padj<0.05\n|LFC|>0.5)", data = merged[bulk_padj < 0.05 & abs(bulk_lfc) > 0.5])
)

conc_dt <- rbindlist(lapply(thresholds, function(x) {
  sub <- x$data
  conc <- mean(sign(sub$bulk_lfc) == sign(sub$sc_lfc), na.rm = TRUE)
  data.table(label = x$label, concordance = conc, n = nrow(sub))
}))
conc_dt[, label := factor(label, levels = label)]

pC <- ggplot(conc_dt, aes(label, concordance)) +
  geom_col(fill = "#4472C4", width = 0.55, color = "white") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%s)", concordance * 100,
                                format(n, big.mark = ","))),
            vjust = -0.2, size = 2.2) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.05), expand = c(0, 0)) +
  labs(x = NULL,
       y = "Direction concordance\n(bulk vs sc hepatocyte)",
       title = "Direction agreement\nby stringency") +
  theme_masld() +
  theme(axis.text.x = element_text(size = 6, lineheight = 1.1))

# ---------------------------------------------------------------------------
# 5. Panel D — fgsea enrichment of bulk DEG sets in sc ranking
# ---------------------------------------------------------------------------
# sc ranked by t-stat (most upregulated in sc → highest rank)
sc_for_rank <- merge(
  sc[!is.na(sc_t), .(ensg_base, sc_t)],
  bulk[, .(ensg_base, symbol)],
  by = "ensg_base"
)[!is.na(symbol) & symbol != ""]

# Deduplicate symbols (keep highest |t|)
sc_for_rank <- sc_for_rank[order(-abs(sc_t))][!duplicated(symbol)]
ranks <- setNames(sc_for_rank$sc_t, sc_for_rank$symbol)

# Build gene sets from bulk DEGs
bulk_up <- merged[bulk_padj < 0.05 & bulk_lfc > 0.5,  unique(symbol)]
bulk_dn <- merged[bulk_padj < 0.05 & bulk_lfc < -0.5, unique(symbol)]
gene_sets <- list(
  "Bulk up-DEGs"   = intersect(bulk_up, names(ranks)),
  "Bulk down-DEGs" = intersect(bulk_dn, names(ranks))
)
message(sprintf("  fgsea gene sets: up=%d, down=%d",
                length(gene_sets[["Bulk up-DEGs"]]),
                length(gene_sets[["Bulk down-DEGs"]])))

set.seed(42)
fres <- fgsea(pathways = gene_sets, stats = ranks, eps = 0, nPermSimple = 2000)
fres[, pval_label := ifelse(pval < 0.001,
                             sprintf("p<0.001"),
                             sprintf("p=%.3f", pval))]
fres[, nes_label := sprintf("NES=%.2f\n%s", NES, pval_label)]
fres[, pathway := factor(pathway, levels = rev(c("Bulk up-DEGs", "Bulk down-DEGs")))]

pD <- ggplot(fres, aes(NES, pathway, fill = NES > 0)) +
  geom_col(width = 0.45, show.legend = FALSE) +
  geom_text(aes(label = nes_label,
                x = ifelse(NES > 0, NES + 0.05, NES - 0.05),
                hjust = ifelse(NES > 0, 0, 1)),
            size = 2.4) +
  scale_fill_manual(values = c("TRUE" = "#C0392B", "FALSE" = "#2980B9")) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_x_continuous(expand = expansion(mult = 0.4)) +
  labs(x = "Normalized Enrichment Score",
       y = NULL,
       title = "Bulk DEG sets enriched\nin sc hepatocyte ranking") +
  theme_masld()

# ---------------------------------------------------------------------------
# 6. Panel E — ρ comparison: all-cell vs hepatocyte-only vs intrinsic-only
#    Uses the 3-way comparison table from 314_allcell_pseudobulk_de.R so the
#    gene set is identical across all three bars.
# ---------------------------------------------------------------------------
comp_path <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/pseudobulk_de/allcell_vs_bulk_comparison.csv")

pE <- NULL
if (file.exists(comp_path)) {
  comp <- fread(comp_path)

  deconv2 <- fread(DECONV_PATH)
  deconv2[, ensg_base := sub("\\.\\d+$", "", gene)]
  hi_ids2  <- deconv2[category == "Hepatocyte_intrinsic", ensg_base]

  rho_ac  <- cor(comp$bulk_lfc, comp$allcell_lfc, method = "spearman", use = "complete.obs")
  rho_hep <- cor(comp$bulk_lfc, comp$hep_lfc,     method = "spearman", use = "complete.obs")
  comp_hi <- comp[ensg_base %in% hi_ids2]
  rho_hi2 <- cor(comp_hi$bulk_lfc, comp_hi$hep_lfc, method = "spearman", use = "complete.obs")
  n_shared <- nrow(comp)

  rho_dt <- data.table(
    comparison = factor(
      c("All-cell\npseudopulk", "Hepatocyte\npseudopulk", "Hepatocyte-\nintrinsic only"),
      levels = c("All-cell\npseudopulk", "Hepatocyte\npseudopulk", "Hepatocyte-\nintrinsic only")
    ),
    rho = c(rho_ac, rho_hep, rho_hi2),
    fill_col = c("wrong", "right", "right")
  )
  message(sprintf("Panel E rho: all-cell=%.3f  hep=%.3f  hep-intrinsic=%.3f  (n_base=%d)",
                  rho_ac, rho_hep, rho_hi2, n_shared))

  pE <- ggplot(rho_dt, aes(comparison, rho, fill = fill_col)) +
    geom_col(width = 0.55, show.legend = FALSE) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    geom_text(aes(label = sprintf("rho=%.3f", rho),
                  y = ifelse(rho >= 0, rho + 0.015, rho - 0.015),
                  vjust = ifelse(rho >= 0, 0, 1)),
              size = 2.4) +
    scale_fill_manual(values = c("wrong" = "#E74C3C", "right" = "#2980B9")) +
    scale_y_continuous(limits = c(min(rho_dt$rho) - 0.08, max(rho_dt$rho) + 0.08)) +
    labs(x = NULL, y = "Spearman rho vs bulk",
         title = sprintf("Comparison level determines\nbulk-sc concordance (n=%d genes)", n_shared),
         subtitle = "Red = wrong comparison; Blue = cell-type matched") +
    theme_masld() +
    theme(axis.text.x = element_text(size = 6.5, lineheight = 1.1),
          plot.subtitle = element_text(size = 5.5, color = "grey40"))
} else {
  message("WARNING: allcell_vs_bulk_comparison.csv not found; skipping Panel E")
}

# ---------------------------------------------------------------------------
# 7. Assemble and save
# ---------------------------------------------------------------------------
if (!is.null(pE)) {
  fig <- (pA + pB + pE) / (pC + pD + plot_spacer()) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 8, face = "bold"))
  fig_width <- 13; fig_height <- 8
} else {
  fig <- (pA + pB) / (pC + pD) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 8, face = "bold"))
  fig_width <- 10; fig_height <- 8
}

out_path <- file.path(FIGS03_DIR, "figS_bulk_sc_concordance.pdf")
ggsave(out_path, fig, width = fig_width, height = fig_height)
message("Saved: ", out_path)

# ---------------------------------------------------------------------------
# 8. Write summary stats
# ---------------------------------------------------------------------------
summary_lines <- c(
  sprintf("Shared genes (all):              %d   rho=%.3f", n_all, rho_all),
  sprintf("Hepatocyte-intrinsic shared:     %d   rho=%.3f", n_hi, rho_hi),
  sprintf("Direction concordance — all:     %.1f%%", conc_dt[1, concordance] * 100),
  sprintf("Direction concordance — DEGs:    %.1f%%", conc_dt[4, concordance] * 100),
  sprintf("fgsea up NES=%.2f  p=%.3f", fres[pathway == "Bulk up-DEGs", NES],
          fres[pathway == "Bulk up-DEGs", pval]),
  sprintf("fgsea dn NES=%.2f  p=%.3f", fres[pathway == "Bulk down-DEGs", NES],
          fres[pathway == "Bulk down-DEGs", pval])
)
writeLines(summary_lines)
writeLines(summary_lines,
           file.path(FIGS03_DIR, "figS_bulk_sc_concordance_stats.txt"))
message("Done.")
