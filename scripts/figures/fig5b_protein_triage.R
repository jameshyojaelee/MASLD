#!/usr/bin/env Rscript
# ==============================================================================
# Protein-layer triage of the prioritized MASLD target set.
# Across the WHOLE prioritized universe measured in an independent liver proteomics
# cohort (PXD051911 liver DIA-MS, n=58), mRNA rank only weakly predicts protein rank.
# This panel shows the prioritized-gene mRNA-vs-protein log2FC distribution. The adjacent
# bar partitions only the protein-significant subset by direction.
#   scatter: x = mRNA log2FC (pooled bulk RNA-seq), y = protein log2FC (PXD051911)
#            colour = 3-way significance/direction class within prioritized genes
#   bar    : protein-significant prioritized genes partitioned by direction
# NOT a pass/fail filter: the majority "not yet protein-corroborated" class is buffering
# AND/OR below DIA-MS detection depth, NOT evidence the target is wrong.
# Output: figures/main/fig5_molecular_context/panels/fig5b_protein_triage.pdf
# Env: rnaseq
# ==============================================================================
# KEY MESSAGE: Adjusted liver DIA-MS resolves a small protein-confirmed core while most measured prioritized genes remain below the protein significance gate.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
grDevices::pdf.options(useDingbats = FALSE)
set.seed(42)
FAM  <- "Helvetica"
BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
CANDIDATE_ROOT <- Sys.getenv("FIGURE_CANDIDATE_ROOT", "")
OUTPUT_DIR <- if (nzchar(CANDIDATE_ROOT)) {
  file.path(CANDIDATE_ROOT, "figure5", "panels")
} else {
  file.path(FIG5_CONTEXT_DIR, "panels")
}

# ── data: current adjusted liver DIA-MS model + canonical bulk atlas ──────────
# Protein estimates use quantile-normalized PXD051911 liver DIA-MS with
# acquisition-batch, age, BMI, and sex adjustment. The atlas provides one
# canonical bulk effect per HGNC symbol, avoiding ambiguous transcript collapse.
prot <- fread(file.path(
  BASE, "Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv"
))[, .(gene, protein_logFC = logFC, protein_padj = padj)]
bulk <- fread(
  file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
  select = c("human_symbol", "bulk_logFC", "bulk_padj")
)[, .(gene = human_symbol, bulk_logFC, bulk_padj)]
stopifnot(!anyDuplicated(prot$gene), !anyDuplicated(bulk$gene))
con <- merge(prot, bulk, by = "gene", all = FALSE, sort = FALSE)
con <- con[is.finite(bulk_logFC) & is.finite(protein_logFC)]
uni <- trimws(readLines(file.path(BASE, "Analysis/Spatial/results/universe_validation/prioritized_universe_FINAL.txt")))
uni <- uni[uni != ""]
con[, in_uni := gene %in% uni]

# whole-proteome vs prioritized-universe correlation (Spearman) + direction concordance
sgn_conc <- function(a, b) mean(sign(a) == sign(b))
rho_full  <- cor(con$bulk_logFC, con$protein_logFC, method = "spearman")
conc_full <- sgn_conc(con$bulk_logFC, con$protein_logFC)
du <- con[in_uni == TRUE]
rho_uni   <- cor(du$bulk_logFC, du$protein_logFC, method = "spearman")
conc_uni  <- sgn_conc(du$bulk_logFC, du$protein_logFC)
fit_uni <- lm(protein_logFC ~ bulk_logFC, data = du)
fit_intercept <- unname(coef(fit_uni)[[1]])
fit_slope <- unname(coef(fit_uni)[[2]])
fit_r2 <- summary(fit_uni)$r.squared

# ── 3-way classification of the PRIORITIZED measured set ──────────────────────
con[, cls := fifelse(!in_uni, "rest",
              fifelse(sign(bulk_logFC) == sign(protein_logFC) & protein_padj < 0.05, "confirmed",
              fifelse(sign(bulk_logFC) != sign(protein_logFC) & protein_padj < 0.05, "discordant",
              "not_corr")))]
n_uni   <- nrow(du)
n_conf  <- sum(con$cls == "confirmed")
n_disc  <- sum(con$cls == "discordant")
n_notc  <- sum(con$cls == "not_corr")
n_rest  <- sum(con$cls == "rest")
cat(sprintf("[protein-triage] measured=%d (universe=%d, rest=%d) | full rho=%.3f conc=%.1f%% | uni rho=%.3f conc=%.1f%%\n",
            nrow(con), n_uni, n_rest, rho_full, 100*conc_full, rho_uni, 100*conc_uni))
cat(sprintf("[protein-triage] 3-way (universe): confirmed=%d (%.1f%%) not_corroborated=%d (%.1f%%) discordant=%d (%.1f%%)\n",
            n_conf, 100*n_conf/n_uni, n_notc, 100*n_notc/n_uni, n_disc, 100*n_disc/n_uni))
cat(sprintf("[protein-triage] prioritized-gene OLS: protein_logFC = %.4f + %.4f * bulk_logFC; R2=%.4f\n",
            fit_intercept, fit_slope, fit_r2))
summary_out <- data.table(
  protein_model = "PXD051911_quantile_batch_age_bmi_sex_adjusted",
  bulk_source = "multi_evidence_atlas_unique_hgnc_symbol",
  n_full_measured = nrow(con),
  n_prioritized_measured = n_uni,
  rho_full = rho_full,
  direction_full = conc_full,
  rho_prioritized = rho_uni,
  direction_prioritized = conc_uni,
  ols_intercept_prioritized = fit_intercept,
  ols_slope_prioritized = fit_slope,
  ols_r2_prioritized = fit_r2,
  n_significant_concordant = n_conf,
  n_significant_discordant = n_disc,
  n_not_significant = n_notc,
  n_protein_significant = n_conf + n_disc,
  pct_concordant_within_significant = 100 * n_conf / (n_conf + n_disc),
  pct_discordant_within_significant = 100 * n_disc / (n_conf + n_disc)
)
summary_path <- file.path(OUTPUT_DIR, "data", "fig5b_protein_triage_summary.tsv")
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
fwrite(summary_out, summary_path, sep = "\t", quote = FALSE)

# colours / draw order / point styling per class. The muted blue–rose pair is
# shared with Figure 5C's program palette; gray remains the unsupported/background
# state. This replaces the saturated teal–orange pair formerly used in 5B.
COL <- c(not_corr = "#9E9E9E", confirmed = "#3B6EA5", discordant = "#B05A5A")
SZ  <- c(not_corr = 0.7,       confirmed = 1.15,      discordant = 1.15)
AL  <- c(not_corr = 0.55,      confirmed = 0.95,      discordant = 0.95)
lab_map <- c(not_corr = "Not significant", confirmed = "Significant, concordant",
             discordant = "Discordant")
plot_dt <- con[in_uni == TRUE]
plot_dt[, cls := factor(cls, levels = c("not_corr", "confirmed", "discordant"))]
setorder(plot_dt, cls)   # draw not significant -> concordant -> discordant

# ── scatter ───────────────────────────────────────────────────────────────────
# single essential data annotation: the headline correlation for the prioritized set
# (full-proteome rho and cross-cohort/buffering narrative live in the message() CAPTION).
ann <- sprintf("ρ = %.2f · %.0f%% direction-concordant", rho_uni, 100*conc_uni)
xr <- range(plot_dt$bulk_logFC); yr <- range(plot_dt$protein_logFC)
sc <- ggplot(plot_dt, aes(bulk_logFC, protein_logFC, colour = cls, size = cls, alpha = cls)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.22, colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.22, colour = "grey70") +
  rasterize_layer(geom_point(stroke = 0), dpi = 400) +
  # Explicit OLS best-fit over the prioritized universe. The coefficients are
  # written to the source sidecar so the plotted line can be rederived exactly.
  geom_abline(intercept = fit_intercept, slope = fit_slope,
              linetype = "solid", linewidth = 0.4, colour = house_ink) +
  scale_colour_manual(values = COL, labels = lab_map, name = NULL,
                      guide = guide_legend(override.aes = list(size = 1.4, alpha = 1), ncol = 1)) +
  scale_size_manual(values = SZ, guide = "none") +
  scale_alpha_manual(values = AL, guide = "none") +
  annotate("text", x = xr[1], y = yr[2], hjust = 0, vjust = 1, size = 6/.pt, family = FAM,
           label = ann) +
  labs(x = expression("mRNA "*log[2]*"FC"),
       y = expression("protein "*log[2]*"FC")) +
  theme_masld(base_size = 6) +
  theme(legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
        legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
        legend.key.size = unit(0.28, "lines"), legend.text = element_text(size = 6),
        legend.margin = margin(1, 2, 1, 2), plot.margin = margin(2, 1, 1, 2))

# ── protein-significant subset bar ───────────────────────────────────────────
# The scatter retains every measured prioritized gene. The adjacent bar uses
# only protein-significant genes, with percentages defined within that subset.
n_sig <- n_conf + n_disc
bar <- data.table(
  cls = factor(c("confirmed", "discordant"),
               levels = c("confirmed", "discordant")),
  n   = c(n_conf, n_disc))
setorder(bar, cls)
bar[, ymax := cumsum(n)][, ymin := ymax - n][, mid := (ymin + ymax)/2]
bar[, pct := 100 * n / n_sig]
bar_lab <- c(
  confirmed = sprintf("Concordant\n%d (%.1f%%)", n_conf, 100 * n_conf / n_sig),
  discordant = sprintf("Discordant\n%d (%.1f%%)", n_disc, 100 * n_disc / n_sig)
)
bp <- ggplot(bar) +
  geom_rect(aes(xmin = 0, xmax = 0.22, ymin = ymin, ymax = ymax, fill = cls),
            colour = "white", linewidth = 0.3) +
  geom_text(aes(x = 0.29, y = mid, label = bar_lab[as.character(cls)]),
            hjust = 0, size = 6/.pt, family = FAM, colour = house_ink, lineheight = 0.9) +
  scale_fill_manual(values = COL[c("confirmed", "discordant")], guide = "none") +
  scale_x_continuous(limits = c(0, 1.0), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, n_sig), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(y = sprintf("protein-significant\n(n = %s)", format(n_sig, big.mark = ",")),
       x = NULL) +
  theme_void() +
  theme(text = element_text(family = FAM, face = "plain"),
        axis.title.y = element_text(size = 6, family = FAM, angle = 90, margin = margin(r = 2)),
        plot.margin = margin(t = 2, r = 1, b = 1, l = 1))

# ── compose (scatter left, vertical bar right) ─────────────────────────────────
p <- (sc | bp) + plot_layout(widths = c(4.25, 1.10))
out <- file.path(OUTPUT_DIR, "fig5b_protein_triage.pdf")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
ggsave(out, p, width = 3.20, height = 1.95, device = grDevices::cairo_pdf)
cat("[protein-triage] saved:", out, "\n")

caption_panel <- "Fig 5B"
target_panel <- "Fig 5"

message(sprintf(paste0(
  "CAPTION (", caption_panel, "): Protein-layer triage of the prioritized MASLD target set in an independent liver ",
  "proteomics cohort (PXD051911 liver DIA-MS, n=58 patients). Protein effects use quantile-normalized intensities ",
  "with acquisition-batch, age, BMI, and sex adjustment. This is a different cohort from the bulk RNA-seq, so the ",
  "mRNA-vs-protein log2FC correlation is attenuated by cross-cohort measurement noise ON TOP OF any real ",
  "post-transcriptional regulation). Each point is a prioritized gene measured at both layers; x = mRNA log2FC ",
  "(pooled bulk RNA-seq), y = protein log2FC (PXD051911). Among %s measured genes from the %s-gene shared ", target_panel, " ",
  "target set, Spearman rho=%.2f and %.0f%% agree in direction. mRNA rank therefore only partially predicts ",
  "protein rank across cohorts. ",
  "The adjacent bar is restricted to the %d protein-significant prioritized genes: ",
  "direction-concordant (n=%d, %.1f%% of significant) and discordant (n=%d, %.1f%% of significant). ",
  "The solid line is the ordinary least-squares fit across the prioritized genes only ",
  "(protein log2FC = %.3f + %.3f × mRNA log2FC; R2=%.3f); the displayed rho is Spearman rank correlation. ",
  "The remaining %s of %s measured prioritized genes were not protein-significant and remain visible as gray ",
  "scatter points; they may reflect buffering or DIA-MS sensitivity and are not treated as negative evidence. ",
  "Fig 5C provides a separate ",
  "selection-conditioned description of protein abundance, histology, and mRNA/protein effects in PXD051911."),
  format(n_uni, big.mark=","), format(length(uni), big.mark=","), rho_uni, 100*conc_uni,
  n_sig, n_conf, 100*n_conf/n_sig, n_disc, 100*n_disc/n_sig,
  fit_intercept, fit_slope, fit_r2,
  format(n_notc, big.mark=","), format(n_uni, big.mark=",")))
