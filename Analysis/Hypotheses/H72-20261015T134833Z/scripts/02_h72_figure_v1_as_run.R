# H72 supplementary figure: female against male per-allele effect at the 974 clumped leads and the
# 281 gene-arm SNPs, raw (U/L) and irnt (sex-specific SD) scales, one PDF per panel.
# KEY MESSAGE: how far female per-allele effects at liver-enzyme loci sit from the male effects,
# on each scale, relative to equality and (raw scale) to the female/male trait SD ratio.
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
OUT <- args[1]
FIG <- file.path(OUT, "figures")
dir.create(FIG, showWarnings = FALSE)
grDevices::pdf.options(useDingbats = FALSE)
FS <- 6
GRAY <- "#9E9E9E"; MAGENTA <- "#C9265E"
th <- theme_classic(base_size = FS) +
  theme(text = element_text(size = FS, colour = "black", face = "plain"),
        axis.text = element_text(size = FS, colour = "black"),
        axis.title = element_text(size = FS, colour = "black"),
        plot.title = element_text(size = FS, colour = "black", face = "plain", hjust = 0),
        legend.text = element_text(size = FS, colour = "black"),
        legend.title = element_blank(), legend.position = "bottom",
        legend.key.size = unit(6, "pt"), line = element_line(linewidth = 0.25))

loc <- fread(file.path(OUT, "locus_sex_effect.tsv"))
gen <- fread(file.path(OUT, "gene_arm_tests.tsv"))
f4 <- fread(file.path(OUT, "set_level_ratio.tsv"))
loc[, set := "leads"]; gen[, set := "gene"]
units <- c(irnt = "sex-specific SD per allele", raw = "U/L per allele")
for (tr in c("ALT", "AST", "GGT")) {
  sdr <- unique(f4[trait == tr & scale == "raw", sd_ratio_raw_PC3])
  for (st in c("leads", "gene")) {
    d0 <- if (st == "leads") loc[trait == tr] else gen[trait == tr]
    for (sc in c("irnt", "raw")) {
      bF <- d0[[paste0("beta_F_", sc)]]; bM <- d0[[paste0("beta_M_", sc)]]; o <- d0[[paste0("orient_", sc)]]
      d <- data.table(x = o * bM, y = o * bF,
                      grp = ifelse(d0$state == "heterogeneous", "different on both scales (BH q < 0.05)", "other states"))
      d <- d[is.finite(x) & is.finite(y)]
      if (nrow(d) == 0) next
      d[, grp := factor(grp, levels = c("other states", "different on both scales (BH q < 0.05)"))]
      setorder(d, grp)
      lim <- range(c(0, d$x, d$y))
      p <- ggplot(d, aes(x, y)) +
        geom_abline(slope = 1, intercept = 0, linewidth = 0.25, colour = "black", linetype = "dashed")
      if (sc == "raw" && length(sdr) == 1 && is.finite(sdr))
        p <- p + geom_abline(slope = sdr, intercept = 0, linewidth = 0.25, colour = "black", linetype = "dotted")
      p <- p + geom_point(aes(colour = grp, shape = grp), size = 0.5, stroke = 0) +
        scale_colour_manual(values = c("other states" = GRAY, "different on both scales (BH q < 0.05)" = MAGENTA), drop = FALSE) +
        scale_shape_manual(values = c("other states" = 16, "different on both scales (BH q < 0.05)" = 17), drop = FALSE) +
        coord_equal(xlim = lim, ylim = lim) +
        labs(x = paste0("Male effect (", units[[sc]], ")"), y = paste0("Female effect (", units[[sc]], ")"),
             title = paste0(tr, ", ", ifelse(st == "leads", "clumped leads", "gene-arm SNPs"), ", ", sc, " (n = ", nrow(d), ")")) + th
      fn <- file.path(FIG, sprintf("H72_%s_%s_%s.pdf", tr, st, sc))
      ggsave(fn, p, width = 2.4, height = 2.7, units = "in", useDingbats = FALSE)
      message(sprintf("%s: %s %s %s; n = %d; heterogeneous = %d; dashed line female = male%s. Alleles oriented so the inverse-variance combined effect is positive.",
                      basename(fn), tr, st, sc, nrow(d), sum(d$grp != "other states"),
                      ifelse(sc == "raw", sprintf("; dotted line slope SD_F/SD_M = %.3f", sdr), "")))
    }
  }
}
writeLines(capture.output(sessionInfo()), file.path(FIG, "sessionInfo_R.txt"))
