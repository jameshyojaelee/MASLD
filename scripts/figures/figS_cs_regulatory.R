#!/usr/bin/env Rscript
# figS_cs_regulatory.R  (2026-06-17)  — Fig 2 supplement (mechanism for para 7)
# Colocalizing credible-set LEAD variants are overwhelmingly non-coding
# regulatory (intron / intergenic / promoter), consistent with effector loci
# acting through cis-regulation rather than protein change. VEP consequence
# classification (not TSS-distance binning). All numbers from disk.
#
# Out: figures/main/fig2_genetics/panels/FigS2K_cs_regulatory_composition.pdf (+ source CSV)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- file.path(FIG3_DIR, "panels")

la <- fread(file.path(BASE, "RNA-seq/results/coloc_variant_classes/lead_causal_annotation.csv"))
la <- la[!is.na(fine_class) & fine_class != ""]
N <- nrow(la)
lev <- c("coding", "fiveUTR", "threeUTR", "promoter", "intron", "intergenic")
labmap <- c(coding = "Coding", fiveUTR = "5' UTR", threeUTR = "3' UTR",
            promoter = "Promoter", intron = "Intron", intergenic = "Intergenic")
cls_cols <- c("Coding" = "#C9265E", "5' UTR" = "#9575CD", "3' UTR" = "#7E57C2",
              "Promoter" = "#26A69A", "Intron" = "#1565C0", "Intergenic" = "#90A4AE")

comp <- la[, .(n = .N), by = fine_class]
comp[, cls := factor(labmap[fine_class], levels = rev(labmap[lev]))]
comp[, pct := 100 * n / N]
pct_noncoding <- 100 * sum(comp[fine_class != "coding", n]) / N

p <- ggplot(comp, aes(pct, cls, fill = cls)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = sprintf("%.0f%%  (%d)", pct, n)), hjust = -0.1, size = GEOM_TEXT_6PT, color = "black") +
  scale_fill_manual(values = cls_cols, guide = "none") +
  scale_x_continuous(limits = c(0, 52), expand = expansion(mult = c(0, 0.04))) +
  labs(x = "% of colocalizing credible-set lead variants", y = NULL) +
  theme_masld(base_size = 6) +
  theme(axis.text.y = element_text(size = 6))

message(sprintf("[caption] Colocalizing variants are non-coding regulatory: %d credible-set leads (VEP consequence); %.1f%% non-coding", N, pct_noncoding))
save_fig(p, file.path(PANEL_DIR, "FigS2K_cs_regulatory_composition.pdf"),
         width = fig_col_width * 1.05, height = 2.3)

fwrite(comp[order(-n), .(fine_class, n, pct = round(pct, 1))],
       file.path(PANEL_DIR, "FigS2K_cs_regulatory_composition_source.csv"))
cat(sprintf("[figS regulatory] N=%d  %.1f%% non-coding regulatory\n", N, pct_noncoding))
