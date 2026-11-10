#!/usr/bin/env Rscript
# KEY MESSAGE: Among 3,769 genes evaluated by multi-signal COLOC, the 426 liver-eQTL
# colocalized genes and the 220 disease vs control DEGs share 26 genes, inside the
# within-AveExpr-decile null (expected 21.2, central 95% 13 to 30; MH OR 1.29,
# 0.84 to 2.01; exact P 0.28).
#
# Copied from figures/candidates/review-r3-fixes-20261015T112148Z/F1/repair/fig1b_U_eval_repair.R
# (the producer of figures/main/fig1_atlas_overview/fig1b_complementary_maps.pdf, sha256 43f87b54...)
# and redrawn to the 1B redraw spec of Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/PANEL_SPECS.md.
# Changes: the two 100% stacked bars (93.9% / 88.2% non-overlap) are replaced by one count axis
# "Genes in both lists" with the exact within-decile null distribution, the 10-bin expectation,
# the unstratified expectation, the observed count and the matched control count; the MH OR and
# its 95% CI are drawn on a log axis; list sizes, the evaluated universe and the DEGs outside the
# comparison are printed. Black text only, 6 pt, no class colours.
#
# Values: every plotted value is recomputed here from the upstream bulk fit and COLOC r2 table
# (light steps) and asserted equal to the deep-review recomputation
# (Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig1/recompute_1b.json and
# recompute_1b_null_distribution_10bins.tsv). The control count 13 and its P are read from the
# prespecified H57 run (nc1_lists.tsv).
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales); library(jsonlite); library(patchwork)
})
grDevices::pdf.options(useDingbats = FALSE)
set.seed(20261017)  # no randomness used; set by project rule
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
setwd(root)
source("scripts/figures/publication_theme.R")
out <- file.path(root, "figures/candidates/panel-deep-fixes-20261017T035633Z/fig1")
pdf_out <- file.path(out, "fig1b_complementary_maps.pdf")
src_out <- file.path(out, "source_tables/fig1b.tsv")
deg_out <- file.path(out, "source_tables/fig1b_degs_outside_comparison.tsv")
stopifnot(!file.exists(pdf_out), !file.exists(src_out), !file.exists(deg_out))

deep <- file.path(root, "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig1")
f_deg <- "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
f_col <- "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
f_gtf <- "data/gencode_v49_gene_metadata.tsv.gz"
f_nc1 <- "Analysis/Hypotheses/H57-20261015T112421Z/run-20261015T113107Z/nc1_lists.tsv"
rj <- fromJSON(gsub("\\bNaN\\b", "null", paste(readLines(file.path(deep, "recompute_1b.json")), collapse = "\n")))  # Python json writes NaN
null_deep <- fread(file.path(deep, "recompute_1b_null_distribution_10bins.tsv"))

# ---- light recomputation from upstream ----
deg <- fread(f_deg)
deg[, ens := sub("\\.[0-9]+$", "", gene)]
stopifnot(!anyDuplicated(deg$ens))
deg[, is_deg := padj < 0.05 & abs(logFC) > 0.5]
col <- fread(f_col)
stopifnot(!anyDuplicated(col$ensembl))
full <- merge(deg, col, by.x = "ens", by.y = "ensembl")
u <- full[susie_state_t12 %in% c("positive", "tested_le_0.5")]
u[, gen := susie_state_t12 == "positive"]
N <- nrow(u); G <- sum(u$gen); D <- sum(u$is_deg); O <- sum(u$gen & u$is_deg)
stopifnot(N == rj$N, G == rj$G, D == rj$D, O == rj$O)
# pandas qcut on rank(method="first") after sorting by Ensembl ID (deep-review recompute_1b.py)
setorder(u, ens)
u[, r := rank(AveExpr, ties.method = "first")]
edges <- quantile(u$r, probs = seq(0, 1, 0.1), type = 7)
u[, bin := cut(r, edges, include.lowest = TRUE, right = TRUE, labels = FALSE)]
bt <- u[, .(n = .N, G = sum(gen), D = sum(is_deg), O = sum(gen & is_deg)), by = bin][order(bin)]
bins_deep <- fread(file.path(deep, "recompute_1b_bins10.tsv"))
stopifnot(all(bt$n == bins_deep$n), all(bt$G == bins_deep$G), all(bt$D == bins_deep$D),
          all(bt$O == bins_deep$O))
E10 <- bt[, sum(G * D / n)]
stopifnot(abs(E10 - rj$strat_10bins$E) < 1e-9)
dist <- 1
for (i in seq_len(nrow(bt))) {
  ks <- 0:min(bt$G[i], bt$D[i])
  dist <- convolve(dist, rev(dhyper(ks, bt$D[i], bt$n[i] - bt$D[i], bt$G[i])), type = "open")
}
xs <- seq_along(dist) - 1L
dist[dist < 0] <- 0
stopifnot(max(abs(dist[seq_len(nrow(null_deep))] - null_deep$prob)) < 1e-9)
p_exact <- sum(dist[abs(xs - E10) >= abs(O - E10) - 1e-9])
stopifnot(abs(p_exact - rj$strat_10bins$p_exact) < 1e-6)
cdf <- cumsum(dist)
q025 <- xs[which(cdf >= 0.025)[1]]; q975 <- xs[which(cdf >= 0.975)[1]]
med <- xs[which(cdf >= 0.5)[1]]
stopifnot(q025 == rj$strat_10bins$null_q025, q975 == rj$strat_10bins$null_q975, med == 21L)
p_ge_obs <- sum(dist[xs >= O])
tabs <- array(0, dim = c(2, 2, nrow(bt)))
for (i in seq_len(nrow(bt))) {
  tabs[, , i] <- matrix(c(bt$O[i], bt$D[i] - bt$O[i], bt$G[i] - bt$O[i],
                          bt$n[i] - bt$G[i] - bt$D[i] + bt$O[i]), 2)
}
mh <- mantelhaen.test(tabs, exact = FALSE)
mh_or <- unname(mh$estimate); mh_lo <- mh$conf.int[1]; mh_hi <- mh$conf.int[2]
stopifnot(abs(mh_or - rj$strat_10bins$MH_OR) < 1e-6, abs(mh_lo - rj$strat_10bins$lo) < 1e-4,
          abs(mh_hi - rj$strat_10bins$hi) < 1e-4)
E_un <- G * D / N
stopifnot(abs(E_un - rj$E_unstrat) < 1e-9)
nc1 <- fread(f_nc1)[universe == "U_eval"]
stopifnot(nrow(nc1) == 1000L, uniqueN(nc1$O) == 1L)
O_ctrl <- nc1$O[1]; p_ctrl <- nc1$p_exact[1]
stopifnot(O_ctrl == rj$nc1_O_min, O_ctrl == 13L)

# DEGs outside the comparison
n_deg_all <- sum(deg$is_deg)
n_deg_nores <- full[is_deg == TRUE & susie_state_t12 == "no_susie_result", .N]
deg_norec <- deg[is_deg == TRUE & !(ens %in% col$ensembl)]
stopifnot(n_deg_all == rj$n_bulk_deg_total, n_deg_nores == rj$deg_no_susie_result_in_U_full,
          nrow(deg_norec) == n_deg_all - rj$D_full)
gtf <- fread(f_gtf)[, .(ensembl_base, gene_biotype)]
stopifnot(!anyDuplicated(gtf$ensembl_base))
deg_bt <- merge(deg[is_deg == TRUE, .(ens, symbol)], gtf, by.x = "ens", by.y = "ensembl_base", all.x = TRUE)
deg_bt[is.na(gene_biotype), gene_biotype := "not_in_gencode_v49"]
deg_bt[, coloc_record := ens %in% col$ensembl]
deg_bt[, in_U_eval := ens %in% u$ens]
deg_bt[, state := fifelse(in_U_eval, "in_comparison",
                   fifelse(coloc_record, "coloc_record_no_multisignal_result_or_untestable", "no_coloc_record"))]
bt_tab <- deg_bt[, .N, by = .(state, gene_biotype)][order(state, -N)]
n_norec <- deg_bt[state == "no_coloc_record", .N]
n_norec_lnc <- deg_bt[state == "no_coloc_record" & gene_biotype == "lncRNA", .N]
n_lnc_deg <- deg_bt[gene_biotype == "lncRNA", .N]
n_lnc_rec <- deg_bt[gene_biotype == "lncRNA" & coloc_record == TRUE, .N]
fwrite(bt_tab, deg_out, sep = "\t")

# ---- source table: one row per mark ----
pmf <- data.table(part = "A_null_distribution", mark = "null_pmf_bar", x = xs, y = dist,
                  in_central_95 = xs >= q025 & xs <= q975)
src <- rbind(
  pmf[x <= 40, .(part, mark, x, y, x_lo = NA_real_, x_hi = NA_real_, in_central_95,
                 printed_label = NA_character_,
                 source = "exact convolution of 10 within-decile hypergeometrics (recomputed; equals recompute_1b_null_distribution_10bins.tsv)")],
  data.table(part = c("A_null_distribution", "B_counts", "B_counts", "B_counts", "B_counts",
                      "A_null_distribution", "C_odds_ratio", "C_odds_ratio"),
             mark = c("expected_10bin_dashed_line", "expected_10bin_dashed_line", "observed_point",
                      "unstratified_expected_open_circle", "matched_control_gray_diamond",
                      "central_95_label", "mh_or_point_with_95ci", "reference_line_or_1"),
             x = c(E10, E10, O, E_un, O_ctrl, NA, mh_or, 1),
             y = NA_real_, x_lo = c(NA, NA, NA, NA, NA, q025, mh_lo, NA),
             x_hi = c(NA, NA, NA, NA, NA, q975, mh_hi, NA), in_central_95 = NA,
             printed_label = c(sprintf("%.1f", E10), NA, sprintf("%d", O), sprintf("%.1f", E_un),
                               sprintf("%d", O_ctrl), sprintf("%d to %d", q025, q975),
                               sprintf("%.2f (%.2f to %.2f) / exact P = %.2f", mh_or, mh_lo, mh_hi, p_exact), NA),
             source = c("sum over 10 AveExpr bins of G_k*D_k/n_k", "same", "genes positive and DEG in U_eval",
                        "G*D/N", "H57 nc1_lists.tsv (U_eval, all 1,000 lists O = 13)",
                        "2.5% and 97.5% quantiles of the exact null", "mantelhaen.test, RBG 95% CI; exact two-sided P from the null",
                        "OR = 1"))
)
meta <- data.table(part = "header_counts", mark = c("N_evaluated", "G_colocalized", "D_degs", "deg_all",
  "deg_no_multisignal_result", "deg_no_coloc_record", "deg_no_coloc_record_lncRNA", "null_median",
  "P_null_ge_observed", "P_matched_control", "lncRNA_degs_total", "lncRNA_degs_with_coloc_record"),
  x = c(N, G, D, n_deg_all, n_deg_nores, n_norec, n_norec_lnc, med, p_ge_obs, p_ctrl, n_lnc_deg, n_lnc_rec),
  y = NA_real_, x_lo = NA_real_, x_hi = NA_real_, in_central_95 = NA, printed_label = NA_character_,
  source = c(rep("deg_results.csv x gene_level_coloc_tier12.csv", 3), "deg_results.csv",
             "deg_results.csv x gene_level_coloc_tier12.csv", "deg_results.csv not in gene_level_coloc_tier12.csv",
             "same, gene_biotype from data/gencode_v49_gene_metadata.tsv.gz", "exact null",
             "exact null", "H57 nc1_lists.tsv p_exact", "gencode v49", "gencode v49"))
src <- rbind(src, meta)
src[mark == "N_evaluated", printed_label := comma(N)]
src[mark == "G_colocalized", printed_label := comma(G)]
src[mark == "D_degs", printed_label := comma(D)]
src[mark == "deg_no_multisignal_result", printed_label := comma(n_deg_nores)]
src[mark == "deg_no_coloc_record", printed_label := comma(n_norec)]
src[mark == "deg_no_coloc_record_lncRNA", printed_label := comma(n_norec_lnc)]
fwrite(src, src_out, sep = "\t")

# ---- plot ----
xlim <- c(0, 40)
gA <- ggplot(pmf[x <= 40], aes(x, y)) +
  geom_col(aes(fill = in_central_95), width = 0.85, show.legend = FALSE) +
  scale_fill_manual(values = c(`TRUE` = "#BDBDBD", `FALSE` = "#E6E6E6")) +
  geom_vline(xintercept = E10, linetype = "22", linewidth = 0.4, colour = "black") +
  annotate("text", x = 0.5, y = 0.112, hjust = 0, vjust = 1, size = 6 / .pt, colour = "black",
           label = sprintf("dashed line: expected within\nAveExpr deciles, %.1f", E10), lineheight = 0.95) +
  annotate("text", x = q975 + 1.2, y = 0.050, hjust = 0, vjust = 0.5, size = 6 / .pt, colour = "black",
           label = sprintf("central 95%% of\nthe null: %d to %d", q025, q975), lineheight = 0.95) +
  scale_x_continuous(limits = xlim + c(-0.5, 0.5), breaks = seq(0, 40, 10), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.115), breaks = c(0, 0.05, 0.10), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Null probability") +
  theme_masld(base_size = 6) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())

rows <- data.table(
  row = factor(c("Observed", "Expected without stratification",
                 "Expression-matched non-colocalized\ngenes (control)"),
               levels = rev(c("Observed", "Expected without stratification",
                              "Expression-matched non-colocalized\ngenes (control)"))),
  x = c(O, E_un, O_ctrl), kind = c("obs", "unstrat", "ctrl"),
  lab = c(sprintf("%d", O), sprintf("%.1f", E_un), sprintf("%d", O_ctrl)))
gB <- ggplot(rows, aes(x, row)) +
  geom_vline(xintercept = E10, linetype = "22", linewidth = 0.4, colour = "black") +
  geom_point(aes(shape = kind, fill = kind), size = 1.8, stroke = 0.4, colour = "black",
             show.legend = FALSE) +
  scale_shape_manual(values = c(obs = 21, unstrat = 21, ctrl = 23)) +
  scale_fill_manual(values = c(obs = "black", unstrat = "white", ctrl = "#9E9E9E")) +
  geom_text(aes(label = lab), nudge_x = 0.9, hjust = 0, size = 6 / .pt, colour = "black") +
  scale_x_continuous(limits = xlim + c(-0.5, 0.5), breaks = seq(0, 40, 10), expand = c(0, 0)) +
  labs(x = "Genes in both lists", y = NULL) +
  theme_masld(base_size = 6) + theme(axis.ticks.y = element_blank())

gC <- ggplot(data.table(row = "Mantel-Haenszel OR\n(95% CI)", x = mh_or, lo = mh_lo, hi = mh_hi),
             aes(x, row)) +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = "black") +
  geom_segment(aes(x = lo, xend = hi, yend = row), linewidth = 0.4, colour = "black") +
  geom_point(shape = 21, fill = "black", colour = "black", size = 1.8, stroke = 0.4) +
  annotate("text", x = mh_hi * 1.08, y = 1, hjust = 0, size = 6 / .pt, colour = "black", lineheight = 0.95,
           label = sprintf("%.2f (%.2f to %.2f)\nexact P = %.2f", mh_or, mh_lo, mh_hi, p_exact)) +
  scale_x_log10(limits = c(0.4, 8), breaks = c(0.5, 1, 2, 4), labels = c("0.5", "1", "2", "4"),
                expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "Odds ratio within AveExpr deciles (log scale)", y = NULL) +
  theme_masld(base_size = 6) + theme(axis.ticks.y = element_blank())

hdr <- sprintf(paste0("%s genes evaluated by multi-signal COLOC (unit: gene)\n",
                      "%d colocalized with liver eQTL (SuSiE PP.H4 > 0.5)\n",
                      "%d disease vs control DEGs (BH q < 0.05, abs log2FC > 0.5; 774 genetic-identity clusters)"),
               comma(N), G, D)
ftr <- sprintf("Not in this comparison: %d DEGs without a multi-signal COLOC result;\n%d DEGs without a COLOC record (%d lncRNA)",
               n_deg_nores, n_norec, n_norec_lnc)
plot <- (gA / gB / gC) + plot_layout(heights = c(1.15, 0.85, 0.32)) +
  plot_annotation(title = hdr, caption = ftr,
                  theme = theme(plot.title = element_text(size = 6, colour = "black", family = "Helvetica",
                                                          face = "plain", hjust = 0, lineheight = 1.05),
                                plot.caption = element_text(size = 6, colour = "black", family = "Helvetica",
                                                            face = "plain", hjust = 0, lineheight = 1.05),
                                plot.title.position = "plot", plot.caption.position = "plot"))
ggsave(pdf_out, plot, width = 3.9, height = 3.0, device = cairo_pdf, bg = "white")
writeLines(capture.output(sessionInfo()), file.path(out, "logs/sessionInfo_fig1b.txt"))
message(sprintf("Caption: Among %s genes evaluated by multi-signal COLOC, %d colocalized with liver eQTL and %d disease vs control DEGs share %d; %.2f expected within 10 AveExpr bins (central 95%% of the exact null %d to %d; median %d); unstratified %.2f; MH OR %.3f (%.3f to %.3f); exact two-sided P %.4f; P(X>=%d) %.3f; matched control list %d DEGs (P %.3f). Outside: %d DEGs no multi-signal result; %d DEGs no COLOC record (%d lncRNA); lncRNA DEGs %d, %d with a record.",
  comma(N), G, D, O, E10, q025, q975, med, E_un, mh_or, mh_lo, mh_hi, p_exact, O, p_ge_obs, O_ctrl, p_ctrl,
  n_deg_nores, n_norec, n_norec_lnc, n_lnc_deg, n_lnc_rec))
cat("FIG1B_CANDIDATE\tPASS\n")
