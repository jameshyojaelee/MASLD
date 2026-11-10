# KEY MESSAGE: in 493 SteatoSITE donors, discovery stage effects recur with the same sign (90-96% of genes) at every
# contrast; only the F3F4 vs F0F1 slope is powered by the prespecified rule, and it replicates (b 1.46, Holm p 0.002).
# Reads finished S3a outputs only. Unit: donor. Model ~ batch + inferred_sex + stage group (voomWithQualityWeights).
suppressPackageStartupMessages({library(ggplot2); library(data.table)})
set.seed(20261028)
SS <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z"
A <- file.path(SS, "out_s3a-20261026T234211Z"); OUT <- commandArgs(trailingOnly = TRUE)[1]
FS <- 6; GT <- FS / .pt; UPC <- "#C9265E"; NSC <- "#9E9E9E"
theme_ss <- function() theme_classic(base_size = FS, base_family = "Helvetica") %+replace% theme(
  text = element_text(size = FS, face = "plain", colour = "black"),
  plot.title = element_text(size = FS, face = "plain", hjust = 0, margin = margin(b = 3)),
  axis.text = element_text(size = FS, colour = "black"), axis.title = element_text(size = FS),
  legend.text = element_text(size = FS), legend.title = element_text(size = FS),
  axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
  legend.key.size = unit(0.25, "cm"), plot.margin = margin(4, 6, 4, 4))
save_pdf <- function(p, name, w, h) { cairo_pdf(file.path(OUT, paste0(name, ".pdf")), width = w, height = h, family = "Helvetica"); print(p); invisible(dev.off()) }

# Panel 1: replication slope per contrast
p1 <- fread(file.path(A, "F_SS1_primary.tsv")); p2 <- fread(file.path(A, "F_SS2_steps.tsv"))
S <- rbind(p1[, .(contrast = c(F3F4 = "F3–F4 vs F0–F1", F2 = "F2 vs F0–F1")[contrast], group = "Stage vs F0–F1",
                  n_genes, slope, lo = ci_norm_lo, hi = ci_norm_hi, conc_obs, conc_null_mean, MDE_slope, testable, p_holm)],
           p2[, .(contrast = c(S2mS1 = "F0–F1 → F2", S3mS2 = "F2 → F3", S4mS3 = "F3 → F4")[step], group = "Adjacent step",
                  n_genes, slope, lo = ci_norm_lo, hi = ci_norm_hi, conc_obs, conc_null_mean, MDE_slope, testable, p_holm)])
S[, contrast := factor(contrast, levels = rev(contrast))]
S[, status := ifelse(testable, "Powered (MDE ≤ 0.30)", "Not powered (descriptive)")]
S[, lab := sprintf("%s genes, %.1f%% same sign", format(n_genes, big.mark = ","), 100 * conc_obs)]
fwrite(S, file.path(OUT, "p5_replication_slopes_source.tsv"), sep = "\t")
g1 <- ggplot(S, aes(slope, contrast)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = NSC) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = "dashed", colour = NSC) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.4) +
  geom_point(aes(shape = status), size = 1.4, colour = UPC, fill = UPC, stroke = 0.5) +
  geom_text(aes(x = 1.7, label = lab), hjust = 0, size = GT) +
  scale_shape_manual(values = c("Powered (MDE ≤ 0.30)" = 21, "Not powered (descriptive)" = 1), name = NULL) +
  scale_x_continuous(limits = c(-0.1, 3.7), breaks = c(0, 0.5, 1, 1.5)) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Replication slope, SteatoSITE logFC on discovery logFC (95% CI)", y = NULL,
       title = "SteatoSITE, 493 donors; discovery-called genes per contrast") +
  theme_ss() + theme(legend.position = "bottom", strip.background = element_blank(),
                     strip.text.y = element_text(size = FS, angle = 0, hjust = 0), plot.margin = margin(4, 10, 4, 4))
save_pdf(g1, "p5_replication_slopes", 4.0, 2.0)

# Panel 2: F3F4 gene-level scatter (supplement)
G <- fread(file.path(A, "per_gene_states.tsv.gz"))[contrast == "F3F4" & !is.na(G_log2FC) & !is.na(D0_log2FC)]
b <- coef(lm(G_log2FC ~ D0_log2FC, G))
fwrite(G[, .(gene, gene_name, D0_log2FC, G_log2FC, state)], file.path(OUT, "p6_F3F4_scatter_source.tsv"), sep = "\t")
writeLines(sprintf("n %d; OLS slope %.4f; intercept %.4f; same sign %.4f", nrow(G), b[2], b[1],
                   mean(sign(G$G_log2FC) == sign(G$D0_log2FC))), file.path(OUT, "p6_check.txt"))
lim <- quantile(abs(c(G$D0_log2FC, G$G_log2FC)), 0.999)
g2 <- ggplot(G, aes(D0_log2FC, G_log2FC)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = NSC) + geom_vline(xintercept = 0, linewidth = 0.3, colour = NSC) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = "dashed", colour = NSC) +
  geom_point(size = 0.25, alpha = 0.25, stroke = 0, colour = "black") +
  geom_abline(slope = b[2], intercept = b[1], linewidth = 0.4, colour = UPC) +
  annotate("text", x = -lim, y = lim, hjust = 0, vjust = 1, size = GT,
           label = sprintf("n %s genes\nslope %.2f\n%.1f%% same sign", format(nrow(G), big.mark = ","), b[2],
                           100 * mean(sign(G$G_log2FC) == sign(G$D0_log2FC)))) +
  coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(x = "Discovery log2FC, F3–F4 vs F0–F1", y = "SteatoSITE log2FC, F3–F4 vs F0–F1",
       title = "Discovery F3–F4 genes in SteatoSITE") + theme_ss()
save_pdf(g2, "p6_F3F4_scatter", 2.3, 2.3)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
