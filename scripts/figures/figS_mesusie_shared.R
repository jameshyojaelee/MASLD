#!/usr/bin/env Rscript
# figS_mesusie_shared.R  (2026-06-17)  — Fig 2 supplement (defends para 4)
# Multi-ancestry joint fine-mapping (meSuSiE) distinguishes causal credible sets
# SHARED across EUR-EAS from ancestry-SPECIFIC ones. Most loci carry a shared
# signal, but a substantial minority are ancestry-restricted — discoverable only
# with multi-ancestry inference. All numbers from disk.
#
# Out: figures/main/fig2_genetics/panels/mesusie_shared_specific.pdf (+ source CSV)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- file.path(FIG3_DIR, "panels")

ml <- fread(file.path(BASE, "GWAS/finemapping/results/mesusie/mesusie_locus_summary.csv"))
sh <- sum(ml$n_cs_shared, na.rm = TRUE); eu <- sum(ml$n_cs_eur, na.rm = TRUE); ea <- sum(ml$n_cs_eas, na.rm = TRUE)
tot <- sh + eu + ea
nloci <- nrow(ml); nshared_loci <- sum(ml$n_cs_shared > 0, na.rm = TRUE)

dec <- data.table(type = c("Shared EUR-EAS", "EUR-only", "EAS-only"), n = c(sh, eu, ea))
dec[, type := factor(type, levels = c("EAS-only", "EUR-only", "Shared EUR-EAS"))]
dec[, pct := 100 * n / tot]
dec_cols <- c("Shared EUR-EAS" = "#00695C", "EUR-only" = "#1565C0", "EAS-only" = "#E68A2E")

p <- ggplot(dec, aes(n, type, fill = type)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = sprintf("%d  (%.0f%%)", n, pct)), hjust = -0.12, size = 2.8, color = "grey15") +
  scale_fill_manual(values = dec_cols, guide = "none") +
  scale_x_continuous(limits = c(0, 1.18 * sh), expand = expansion(mult = c(0, 0))) +
  labs(x = "meSuSiE credible sets", y = NULL,
       title = "Shared vs ancestry-specific causal signals",
       subtitle = sprintf("%d EUR-EAS loci jointly fine-mapped; %d (%.0f%%) carry >=1 shared credible set",
                          nloci, nshared_loci, 100 * nshared_loci / nloci)) +
  theme_masld(base_size = 9) +
  theme(plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 5.8, color = "grey35"),
        axis.text.y = element_text(size = 8, face = "bold"))

save_fig(p, file.path(PANEL_DIR, "mesusie_shared_specific.pdf"),
         width = fig_col_width * 1.05, height = 2.2)

fwrite(dec[order(-n), .(type, n, pct = round(pct, 1))],
       file.path(PANEL_DIR, "mesusie_shared_specific_source.csv"))
cat(sprintf("[figS meSuSiE] shared=%d eur=%d eas=%d (%.0f%% shared); %d/%d loci shared\n",
            sh, eu, ea, 100 * sh / tot, nshared_loci, nloci))
