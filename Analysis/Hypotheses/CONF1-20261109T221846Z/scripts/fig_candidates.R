# CONF1 candidate panels (not adopted). Individual PDFs, 6 pt, no titles.
# Fig 2 candidate A  KEY MESSAGE: many colocalized liver-enzyme loci raise ALT without raising MASLD risk,
#                    unlike PNPLA3 and TM6SF2.
# Fig 2 candidate B  KEY MESSAGE: enzyme-only loci also carry little FinnGen MASLD or liver-fat signal.
# Fig 3 candidate A  KEY MESSAGE: of the replicated stage genes, 1,155 change inside hepatocytes; most of the rest
#                    cannot be separated from cell-mix change in bulk.
# Fig 3 candidate B  KEY MESSAGE: within-hepatocyte change follows the bulk direction for these genes.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
args <- commandArgs(trailingOnly = TRUE); OUT <- args[1]
if (is.na(OUT) || dir.exists(OUT)) stop("output directory missing or already exists")
dir.create(OUT)
pdf.options(useDingbats = FALSE)
HYP <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses"
FS <- 6; PT <- FS / .pt
up <- "#C9265E"; down <- "#1565C0"; ns <- "#9E9E9E"
th <- theme_classic(base_size = FS, base_family = "Helvetica") +
  theme(axis.text = element_text(size = FS, colour = "black"), axis.title = element_text(size = FS),
        legend.text = element_text(size = FS), legend.title = element_blank(), strip.text = element_text(size = FS),
        strip.background = element_blank(), axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25),
        legend.key.size = unit(6, "pt"))
save <- function(p, f, w, h) ggsave(file.path(OUT, f), p, width = w, height = h, units = "in", device = cairo_pdf)

# ---- Fig 2 candidates (DIG1 job 73 loci; MVP EUR z oriented to ALT-raising allele)
L <- fread(file.path(HYP, "DIG1-20261024T113951Z/out_release_vs_disease_loci/loci.tsv"))
# benchmarks from job 73 log 22402952 (not in loci.tsv)
Bm <- data.table(genes = c("PNPLA3", "TM6SF2"), z_MVP_ALT_EUR = c(49.5, 25.3), z_MVP_NAFLD_EUR = c(42.7, 23.3),
                 z_FinnGen_NAFLD = c(18.2, 9.5), z_Pazoki_PDFF = c(27.1, 22.3))
K <- 0.891
CLS <- c(release = "Enzyme only", intermediate = "Intermediate", disease = "Disease liability", opposite = "Opposite")
L[, class := factor(CLS[cls], levels = CLS)]
L[, lab := fifelse(genes %in% c("TNFSF10", "PANX1"), genes, NA_character_)]
colc <- c("Enzyme only" = down, "Intermediate" = ns, "Disease liability" = up, "Opposite" = "#F4A674")
pA <- ggplot(L, aes(z_MVP_ALT_EUR, z_MVP_NAFLD_EUR)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey70") +
  geom_abline(slope = K, intercept = 0, linewidth = 0.25, linetype = "dashed", colour = "black") +
  geom_point(aes(colour = class), size = 0.9, alpha = 0.85, stroke = 0) +
  geom_point(data = Bm, colour = "black", size = 1.1, shape = 17) +
  geom_text(data = Bm, aes(label = genes), size = PT, fontface = "italic", hjust = 1.15, vjust = -0.3) +
  geom_text(aes(label = lab), size = PT, fontface = "italic", hjust = -0.15, na.rm = TRUE) +
  scale_colour_manual(values = colc, drop = TRUE) +
  labs(x = "ALT z (MVP)", y = "MASLD z (MVP)") + th + theme(legend.position = c(0.25, 0.85))
save(pA, "fig2cand_A_alt_vs_masld_z.pdf", 2.6, 2.4)

# out of sample: observed z against the z expected if the locus behaved like PNPLA3/TM6SF2 (k_trait x MVP ALT z)
kF <- mean(Bm$z_FinnGen_NAFLD / Bm$z_MVP_ALT_EUR); kP <- mean(Bm$z_Pazoki_PDFF / Bm$z_MVP_ALT_EUR)
long <- rbind(L[, .(locus, class, trait = "MASLD (FinnGen)", pred = kF * z_MVP_ALT_EUR, z = z_FinnGen_NAFLD)],
              L[, .(locus, class, trait = "Liver fat (PDFF)", pred = kP * z_MVP_ALT_EUR, z = z_Pazoki_PDFF)])
long <- long[is.finite(z) & class != "Opposite"]
long[, trait := factor(trait, levels = c("MASLD (FinnGen)", "Liver fat (PDFF)"))]
pB <- ggplot(long, aes(pred, z)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.25, linetype = "dashed") +
  geom_point(data = long[class == "Intermediate"], colour = ns, size = 0.6, alpha = 0.6, stroke = 0) +
  geom_point(data = long[class != "Intermediate"], aes(colour = class), size = 0.8, alpha = 0.9, stroke = 0) +
  facet_wrap(~trait, nrow = 1, scales = "free") + scale_colour_manual(values = colc, guide = "none") +
  labs(x = "Expected z from ALT effect", y = "Observed z") + th
save(pB, "fig2cand_B_outofsample_expected_vs_observed.pdf", 2.8, 1.5)
fwrite(rbind(long[, .(n = .N, sum_obs_over_sum_expected = sum(z) / sum(pred), median_z = median(z)), by = .(trait, class)],
             data.table(trait = c("k_FinnGen", "k_PDFF"), class = "benchmark", n = 2, sum_obs_over_sum_expected = c(kF, kP), median_z = NA),
             fill = TRUE), file.path(OUT, "fig2cand_B_values.tsv"), sep = "\t")
fwrite(L[, .N, by = class], file.path(OUT, "fig2cand_A_values.tsv"), sep = "\t")

# ---- Fig 3 candidates (C3/C3b post hoc partition)
P <- fread(file.path(HYP, "CONF1-20261109T221846Z/out_c3b_ambient/final_partition_posthoc.tsv.gz"))
FL <- c(within_hepatocyte = "Inside hepatocytes", composition_associated = "Not separable from cell mix",
        bulk_composition_independent = "Beyond cell mix (bulk)", unresolved = "Unresolved",
        within_other_lineage_uncorrected = NA)
P <- P[final != "within_other_lineage_uncorrected"]
P[, cls := factor(FL[final], levels = rev(FL[!is.na(FL)]))]
cnt <- P[, .N, by = .(cls, dir)][, n_signed := fifelse(dir == "up", N, -N)]
fwrite(cnt, file.path(OUT, "fig3cand_A_values.tsv"), sep = "\t")
pC <- ggplot(cnt, aes(n_signed, cls, fill = dir)) + geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  scale_fill_manual(values = c(up = up, down = down), labels = c(up = "Higher in F3-F4", down = "Lower in F3-F4")) +
  scale_x_continuous(labels = abs) + labs(x = "Genes", y = NULL) + th + theme(legend.position = "bottom")
save(pC, "fig3cand_A_partition_counts.pdf", 2.8, 1.6)

H <- P[is.finite(corr_logFC)]
H[, hl := final == "within_hepatocyte"]
NAMED <- c("CCL20", "AKR1B10", "SOX9", "TRIM31", "CXCL10", "CDKN1A", "HKDC1", "CYP2C19", "GNMT", "MAT1A")
pD <- ggplot(H, aes(bulk_log2fc, corr_logFC)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey70") + geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey70") +
  geom_point(data = H[hl == FALSE], colour = ns, size = 0.3, alpha = 0.3, stroke = 0) +
  geom_point(data = H[hl == TRUE], aes(colour = dir), size = 0.45, alpha = 0.7, stroke = 0) +
  geom_point(data = H[symbol %in% NAMED], colour = "black", size = 0.7, stroke = 0) +
  ggrepel::geom_text_repel(data = H[symbol %in% NAMED], aes(label = symbol), size = PT, fontface = "italic", segment.size = 0.15,
                           min.segment.length = 0, box.padding = 0.2, seed = 20261109) +
  scale_colour_manual(values = c(up = up, down = down), guide = "none") +
  labs(x = "Bulk log2FC, F3-F4 vs F0-F1", y = "Hepatocyte log2FC, corrected") + th
save(pD, "fig3cand_B_bulk_vs_hepatocyte.pdf", 2.6, 2.4)
fwrite(H[symbol %in% NAMED, .(symbol, final, bulk_log2fc, corr_logFC, corr_q, ambient_share)], file.path(OUT, "fig3cand_B_named.tsv"), sep = "\t")
cat("spearman within_hepatocyte", cor(H[hl == TRUE, bulk_log2fc], H[hl == TRUE, corr_logFC], method = "spearman"),
    "sign agreement", mean(sign(H[hl == TRUE, bulk_log2fc]) == sign(H[hl == TRUE, corr_logFC])), "\n")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
