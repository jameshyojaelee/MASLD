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

L <- fread(file.path(HYP, "EQTL-DEDUP-PARITY-20261110T041858Z/fig2fg/out_loci_23160547/loci.tsv")); stopifnot(nrow(L) == 167L)
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
save(pA, "Fig2F_enzyme_vs_disease_liability.pdf", 2.6, 2.4)

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
save(pB, "Fig2G_enzyme_loci_out_of_sample.pdf", 2.8, 1.5)
fwrite(rbind(long[, .(n = .N, sum_obs_over_sum_expected = sum(z) / sum(pred), median_z = median(z)), by = .(trait, class)],
             data.table(trait = c("k_FinnGen", "k_PDFF"), class = "benchmark", n = 2, sum_obs_over_sum_expected = c(kF, kP), median_z = NA),
             fill = TRUE), file.path(OUT, "Fig2G_enzyme_loci_out_of_sample_values.tsv"), sep = "\t")
fwrite(L[, .N, by = class], file.path(OUT, "Fig2F_enzyme_vs_disease_liability_values.tsv"), sep = "\t")
