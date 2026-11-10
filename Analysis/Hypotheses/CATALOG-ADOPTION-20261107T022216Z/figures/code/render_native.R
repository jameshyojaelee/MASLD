#!/usr/bin/env Rscript
# KEY MESSAGE: Genetic and disease-state overlap remains compatible with its
# expression-stratified null after excluding untestable SuSiE genes.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(scales); library(patchwork)})
grDevices::pdf.options(useDingbats = FALSE)
setDTthreads(1); set.seed(20261106)
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
out <- Sys.getenv("CATALOG_FIGURE_OUT")
stopifnot(nzchar(out))
source(file.path(root, "scripts/figures/publication_theme.R"))
run <- file.path(root, "Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN1/gene_level_20261021T165252Z")
native_path <- file.path(run, "h57_universe_rule_a_sens.tsv")
bulk_path <- file.path(root, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
u <- fread(native_path)
u[, `:=`(gen = as.character(gen) %in% c("True", "TRUE", "1"),
          is_deg = as.character(deg) %in% c("True", "TRUE", "1"))]
setorder(u, g)
bulk <- fread(bulk_path)
bulk[, g := sub("\\.[0-9]+$", "", gene)]
stopifnot(!anyDuplicated(u$g), !anyDuplicated(bulk$g))
stopifnot(nrow(u) == 3589L, sum(u$gen) == 426L, sum(u$is_deg) == 209L,
          sum(u$gen & u$is_deg) == 26L)
u <- merge(u, bulk[, .(g, symbol, logFC, padj)], by = "g", sort = TRUE)
stopifnot(nrow(u) == 3589L, all(u$is_deg == (u$padj < .05 & abs(u$logFC) > .5)))
u[, r := rank(AveExpr, ties.method = "first")]
edges <- quantile(u$r, probs = seq(0, 1, .1), type = 7)
u[, bin := cut(r, edges, include.lowest = TRUE, right = TRUE, labels = FALSE)]
bt <- u[, .(n = .N, G = sum(gen), D = sum(is_deg), O = sum(gen & is_deg)), by = bin][order(bin)]
N <- nrow(u); G <- sum(u$gen); D <- sum(u$is_deg); O <- sum(u$gen & u$is_deg)
states <- fread(file.path(run, "gene_states.tsv"), select = c("ensembl", "state_a_sens"))
bulk[, is_deg := padj < .05 & abs(logFC) > .5]
outside <- merge(bulk[is_deg == TRUE, .(g, symbol)], states,
                 by.x = "g", by.y = "ensembl", all.x = TRUE)
gtf <- fread(file.path(root, "data/gencode_v49_gene_metadata.tsv.gz"), select = c("ensembl_base", "gene_biotype"))
outside <- merge(outside, gtf, by.x = "g", by.y = "ensembl_base", all.x = TRUE)
outside[, in_comparison := g %in% u$g]
n_no_result <- outside[state_a_sens == "no_susie_result", .N]
n_untestable <- outside[state_a_sens == "untestable", .N]
n_no_record <- outside[is.na(state_a_sens), .N]
n_no_record_lnc <- outside[is.na(state_a_sens) & gene_biotype == "lncRNA", .N]
stopifnot(sum(bulk$is_deg) == 1854L, n_no_result == 935L, n_untestable == 11L,
          n_no_record == 699L, n_no_record_lnc == 503L)
fwrite(outside, file.path(out, "source_tables/fig1b_degs_outside_comparison.tsv"), sep = "\t")
E10 <- bt[, sum(G * D / n)]
dist <- 1
for (i in seq_len(nrow(bt))) {
  ks <- 0:min(bt$G[i], bt$D[i])
  dist <- convolve(dist, rev(dhyper(ks, bt$D[i], bt$n[i] - bt$D[i], bt$G[i])), type = "open")
}
xs <- seq_along(dist) - 1L; dist[dist < 0] <- 0
p_exact <- sum(dist[abs(xs - E10) >= abs(O - E10) - 1e-9])
cdf <- cumsum(dist)
q025 <- xs[which(cdf >= .025)[1]]; q975 <- xs[which(cdf >= .975)[1]]
tabs <- array(0, dim = c(2, 2, nrow(bt)))
for (i in seq_len(nrow(bt))) {
  tabs[, , i] <- matrix(c(bt$O[i], bt$D[i] - bt$O[i], bt$G[i] - bt$O[i],
    bt$n[i] - bt$G[i] - bt$D[i] + bt$O[i]), 2)
}
mh <- mantelhaen.test(tabs, exact = FALSE)
mh_or <- unname(mh$estimate); mh_lo <- mh$conf.int[1]; mh_hi <- mh$conf.int[2]
rho <- cor.test(u$logFC, u$susie, method = "spearman", exact = FALSE)
E_un <- G * D / N
stopifnot(abs(E10 - 21.396150075473457) < 1e-10,
          abs(p_exact - .27698706604034284) < 1e-10,
          abs(mh_or - 1.2815058404412019) < 1e-9,
          abs(mh_lo - .8261843360996866) < 1e-8,
          abs(mh_hi - 1.9877612626230643) < 1e-8)
catalog_matrix <- fread(file.path(root, "RNA-seq/results/evidence_passports/releases/gene-catalog-eligibility-20261106/call_matrix.tsv"))
catalog_jointly_eligible <- catalog_matrix[genetics_call != "untestable", sum(indeterminate + supported)]
native_interface <- length(intersect(bulk$g, states$ensembl))
stopifnot(catalog_jointly_eligible == 14502L, native_interface == 14704L)
summary <- data.table(N, G, D, O, E10, E_un, mh_or, mh_lo, mh_hi, p_exact,
  null_q025 = q025, null_q975 = q975, rho = unname(rho$estimate), rho_p = rho$p.value,
  native_interface, catalog_jointly_eligible,
  deg_no_susie_result = n_no_result, deg_untestable_susie = n_untestable,
  deg_no_coloc_record = n_no_record, deg_no_coloc_record_lncRNA = n_no_record_lnc,
  source = native_path, eligibility_rule = "positive > untestable > tested below threshold > no SuSiE result")
fwrite(summary, file.path(out, "source_tables/native_overlap_summary.tsv"), sep = "\t")
fwrite(bt, file.path(out, "source_tables/native_overlap_bins.tsv"), sep = "\t")
pmf <- data.table(x = xs, y = dist, in_central_95 = xs >= q025 & xs <= q975)
fwrite(pmf, file.path(out, "source_tables/native_overlap_null.tsv"), sep = "\t")
u[, display_class := fifelse(gen & is_deg, "Both", fifelse(gen, "Genetics", fifelse(is_deg, "Transcriptomics", "Neither")))]
fwrite(u, file.path(out, "source_tables/fig3g.tsv"), sep = "\t")
fwrite(u[, .N, by = display_class], file.path(out, "source_tables/fig3g_class_counts.tsv"), sep = "\t")

# Preserve the reviewed count/null/odds-ratio arrangement; the old matched-control
# list was defined in the former universe and is not carried into this figure.
xlim <- c(0, 40)
gA <- ggplot(pmf[x <= 40], aes(x, y)) +
  geom_col(aes(fill = in_central_95), width = .85, show.legend = FALSE) +
  scale_fill_manual(values = c(`TRUE` = "#9E9E9E", `FALSE` = "#E6E6E6")) +
  geom_vline(xintercept = E10, linetype = "22", linewidth = .4, colour = "black") +
  annotate("text", x = .5, y = .112, hjust = 0, vjust = 1, size = 6 / .pt,
    label = sprintf("dashed line: expected within\nAveExpr deciles, %.1f", E10), lineheight = .95) +
  annotate("text", x = q975 + 1.2, y = .050, hjust = 0, size = 6 / .pt,
    label = sprintf("central 95%% of\nthe null: %d to %d", q025, q975), lineheight = .95) +
  scale_x_continuous(limits = xlim + c(-.5, .5), breaks = seq(0, 40, 10), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, .115), breaks = c(0, .05, .10), expand = c(0, 0)) +
  coord_cartesian(clip = "off") + labs(x = NULL, y = "Null probability") + theme_masld(base_size = 6) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
rows <- data.table(row = factor(c("Observed", "Expected without stratification"),
  levels = c("Expected without stratification", "Observed")), x = c(O, E_un),
  kind = c("obs", "unstrat"), lab = c(sprintf("%d", O), sprintf("%.1f", E_un)))
fwrite(rows, file.path(out, "source_tables/fig1b_count_marks.tsv"), sep = "\t")
gB <- ggplot(rows, aes(x, row)) +
  geom_vline(xintercept = E10, linetype = "22", linewidth = .4, colour = "black") +
  geom_point(aes(fill = kind), shape = 21, size = 1.8, stroke = .4, colour = "black", show.legend = FALSE) +
  scale_fill_manual(values = c(obs = "black", unstrat = "white")) +
  geom_text(aes(label = lab), nudge_x = .9, hjust = 0, size = 6 / .pt, colour = "black") +
  scale_x_continuous(limits = xlim + c(-.5, .5), breaks = seq(0, 40, 10), expand = c(0, 0)) +
  labs(x = "Genes in both lists", y = NULL) + theme_masld(base_size = 6) + theme(axis.ticks.y = element_blank())
gC <- ggplot(data.table(row = "Mantel-Haenszel OR\n(95% CI)", x = mh_or, lo = mh_lo, hi = mh_hi), aes(x, row)) +
  geom_vline(xintercept = 1, linewidth = .3, colour = "black") +
  geom_segment(aes(x = lo, xend = hi, yend = row), linewidth = .4, colour = "black") +
  geom_point(shape = 21, fill = "black", size = 1.8, stroke = .4) +
  annotate("text", x = mh_hi * 1.08, y = 1, hjust = 0, size = 6 / .pt, lineheight = .95,
    label = sprintf("%.2f (%.2f to %.2f)\nexact P = %.2f", mh_or, mh_lo, mh_hi, p_exact)) +
  scale_x_log10(limits = c(.4, 8), breaks = c(.5, 1, 2, 4), expand = c(0, 0)) +
  coord_cartesian(clip = "off") + labs(x = "Odds ratio within AveExpr deciles (log scale)", y = NULL) +
  theme_masld(base_size = 6) + theme(axis.ticks.y = element_blank())
hdr <- sprintf("%s genes evaluated by multi-signal COLOC (unit: gene)\n%d colocalized with liver eQTL (SuSiE PP.H4 > 0.5)\n%d disease vs control DEGs (BH q < 0.05, abs log2FC > 0.5; 774 genetic-identity clusters)", comma(N), G, D)
ftr <- sprintf("Not compared: %d DEGs without a multi-signal result; %d with untestable SuSiE\n%d DEGs without a COLOC record (%d lncRNA)",
               n_no_result, n_untestable, n_no_record, n_no_record_lnc)
plot <- (gA / gB / gC) + plot_layout(heights = c(1.15, .85, .32)) +
  plot_annotation(title = hdr, caption = ftr, theme = theme(
    plot.title = element_text(size = 6, family = "Helvetica", face = "plain", hjust = 0, lineheight = 1.05),
    plot.caption = element_text(size = 6, family = "Helvetica", face = "plain", hjust = 0, lineheight = 1.05),
    plot.title.position = "plot", plot.caption.position = "plot"))
ggsave(file.path(out, "panels/fig1b_complementary_maps.pdf"), plot, width = 3.9, height = 3, device = cairo_pdf, bg = "white")

# Figure 3G preserves its native 170 x 210 pt layout and four-class legend.
txt6 <- element_text(size = 6, face = "plain", family = "Helvetica", color = "black")
lab_of <- c(Neither = "Neither", Transcriptomics = "Transcriptomics sig.", Genetics = "Genetics sig.", Both = "Both")
u[, state := factor(lab_of[display_class], levels = lab_of)]
n_by <- u[, .N, by = state]
leg <- setNames(sprintf("%s\n(%s)", n_by$state, prettyNum(n_by$N, big.mark = ",")), n_by$state)
cols <- c(Neither = "#9E9E9E", `Transcriptomics sig.` = "#C9265E", `Genetics sig.` = "#1565C0", Both = "#6A3D9A")
p_g <- ggplot(u[order(state)], aes(logFC, susie, color = state)) +
  geom_point(size = .45, alpha = .7, stroke = 0) + scale_color_manual(values = cols, labels = leg[names(cols)], name = NULL) +
  scale_y_continuous(breaks = seq(0, 1, .2), limits = c(0, 1), expand = expansion(mult = .02)) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 1.2, alpha = 1))) +
  labs(x = "Disease-state log2FC", y = "SuSiE PP.H4") + theme_masld(base_size = 6) +
  theme(legend.position = "bottom", axis.text = txt6, axis.title = txt6, legend.text = txt6,
    legend.key.size = unit(.25, "cm"), legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"), panel.grid = element_blank())
ggsave(file.path(out, "panels/fig3g_genetic_state_complementarity.pdf"), p_g,
  width = 170 / 72, height = 210 / 72, device = cairo_pdf, family = "Helvetica")
writeLines(capture.output(sessionInfo()), file.path(out, "logs/sessionInfo_native.txt"))
print(summary)
