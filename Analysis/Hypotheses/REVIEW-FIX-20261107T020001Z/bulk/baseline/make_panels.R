# Candidate SSX panels (not adopted). Reads finished SSX outputs only.
# A: GTEx liver eQTL direction vs SteatoSITE allelic shift (50 genes).
# B: intronic vs exonic F3-F4 change (12,010 genes): stage change is mirrored in pre-mRNA.
# C: Fig 3E NMF axis stage profiles, discovery vs SteatoSITE.
suppressPackageStartupMessages({library(ggplot2); library(data.table)})
X <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z"
OUT <- commandArgs(trailingOnly = TRUE)[1]
FS <- 6; GT <- FS / .pt; MAG <- "#C9265E"; BLU <- "#1565C0"; GRY <- "#9E9E9E"
theme_ss <- theme_classic(base_size = FS, base_family = "Helvetica") %+replace% theme(
  text = element_text(size = FS, colour = "black"), axis.text = element_text(size = FS, colour = "black"),
  axis.title = element_text(size = FS), plot.title = element_text(size = FS, hjust = 0, margin = margin(b = 3)),
  strip.text = element_text(size = FS), strip.background = element_blank(), legend.text = element_text(size = FS),
  legend.title = element_blank(), legend.key.size = unit(6, "pt"),
  axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3), plot.margin = margin(4, 6, 4, 4))
save <- function(g, f, w, h) { cairo_pdf(file.path(OUT, f), width = w, height = h, family = "Helvetica"); print(g); invisible(dev.off()) }

# A
A <- fread(file.path(X, "out_a-20261102T133643Z/primary_genes.tsv"))
A[, dir := factor(fifelse(agree, "Same direction as GTEx", "Opposite direction"), levels = c("Same direction as GTEx", "Opposite direction"))]
fwrite(A[, .(gene_id = id, variant_id, gtex_liver_slope = slope, n_het_donors = n_het, pooled_alt_fraction = pooled_f,
             null_median_logit = null_med_logit, allelic_shift_logit = shift, same_direction = agree)], file.path(OUT, "ssx_a_source.tsv"), sep = "\t")
gA <- ggplot(A, aes(slope, shift, colour = dir)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey60") + geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey60") +
  geom_point(size = 1, stroke = 0) + scale_colour_manual(values = c(MAG, GRY)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2, size = GT,
           label = sprintf("%d of %d genes (%.0f%%)\nrandom-sign null 50%%, p = 0.008", sum(A$agree), nrow(A), 100 * mean(A$agree))) +
  labs(x = "GTEx liver eQTL slope (ALT allele)", y = "ALT allele read share in SteatoSITE\nheterozygotes (logit, null-centred)",
       title = "Liver eQTL direction in diseased FFPE liver (515 donors)") + theme_ss + theme(legend.position = "bottom")
save(gA, "ssx_a_eqtl_direction.pdf", 2.6, 2.6)

# B
B <- fread(file.path(X, "out_b-20261102T144813Z/eisa_all_filtered_genes.tsv.gz"))
named <- c("GNMT", "MAT1A", "COL1A1", "LUM", "IGFBP7", "TNFSF10")
B[, grp := fifelse(eisa_FDR < 0.05 & abs(dEmI) >= 0.5, "Exon/intron change >= 1.41-fold (q < 0.05)", "Other genes")]
fwrite(B[, .(gene, symbol, dE_log2 = dE, dI_log2 = dI, dE_minus_dI = dEmI, eisa_BH_q = eisa_FDR, class)], file.path(OUT, "ssx_b_source.tsv"), sep = "\t")
rho <- cor(B$dE, B$dI, method = "spearman")
gB <- ggplot(B, aes(dE, dI)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.25, colour = "grey60", linetype = 2) +
  geom_point(data = B[grp == "Other genes"], colour = GRY, size = 0.3, stroke = 0, alpha = 0.5) +
  geom_point(data = B[grp != "Other genes"], colour = MAG, size = 0.6, stroke = 0) +
  geom_point(data = B[symbol %in% named], colour = BLU, size = 1, stroke = 0) +
  geom_text(data = B[symbol %in% named][, `:=`(hj = fifelse(symbol == "COL1A1", 1.15, -0.15), vj = fifelse(symbol == "COL1A1", -0.3, 1.1))],
            aes(label = symbol, hjust = hj, vjust = vj), colour = BLU, size = GT) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2, size = GT,
           label = sprintf("%s genes, Spearman %.2f\nmagenta: %d genes with exon/intron\nchange >= 1.41-fold, q < 0.05",
                           format(nrow(B), big.mark = ","), rho, sum(B$grp != "Other genes"))) +
  coord_equal() +
  labs(x = "Exonic change, F3-F4 vs F0-F1 (log2)", y = "Intronic change, F3-F4 vs F0-F1 (log2)",
       title = "Stage change is mirrored in pre-mRNA (493 donors)") + theme_ss
save(gB, "ssx_b_exon_vs_intron.pdf", 2.6, 2.7)

# C
D <- fread(file.path(X, "out_c-20261102T134043Z/discovery_vs_steatosite.tsv"))
C <- rbind(D[, .(program, label, cohort = "Discovery (3 cohorts, 387 units)", F0F1 = 0, F2 = F2_disc, F3 = F3_disc, F4 = F4_disc)],
           D[, .(program, label, cohort = "SteatoSITE (493 donors)", F0F1 = 0, F2 = F2_ss, F3 = F3_ss, F4 = F4_ss)])
C <- melt(C, id.vars = c("program", "label", "cohort"), variable.name = "stage", value.name = "activity")
C[, stage := factor(sub("F0F1", "F0-F1", stage), levels = c("F0-F1", "F2", "F3", "F4"))]
C[, panel := factor(paste0(program, " ", label), levels = unique(paste0(D$program, " ", D$label)))]
fwrite(C, file.path(OUT, "ssx_c_source.tsv"), sep = "\t")
gC <- ggplot(C, aes(stage, activity, colour = cohort, group = cohort)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey60") +
  geom_line(linewidth = 0.4) + geom_point(size = 0.9, stroke = 0) +
  scale_colour_manual(values = c(BLU, MAG)) + facet_wrap(~ panel, nrow = 1, labeller = label_wrap_gen(16)) +
  labs(x = "Fibrosis stage", y = "Axis activity vs F0-F1 (SD)", title = "Figure 3E NMF axes by stage: discovery and SteatoSITE") +
  theme_ss + theme(legend.position = "bottom", axis.text.x = element_text(size = FS, angle = 45, hjust = 1))
save(gC, "ssx_c_nmf_stage_profiles.pdf", 4.8, 2.1)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
