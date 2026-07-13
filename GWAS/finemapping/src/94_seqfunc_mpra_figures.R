#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
IN <- file.path(ROOT, "GWAS/finemapping/results/seqfunc/mpra_benchmark")
OUT <- file.path(ROOT, "figures/supplementary/seqfunc")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
d <- fread(file.path(IN, "mpra_model_benchmark.tsv"))
d <- d[feature == "cbp_abs_logfc" & !is.na(auroc)]
d[, context := factor(context, levels = c("HepG2_ctrl", "HepG2_PAOA", "LX2_ctrl", "LX2_TGFb"))]
p <- ggplot(d, aes(context, auroc, color = cell_model)) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "grey55") +
  geom_errorbar(aes(ymin = auroc_lo, ymax = auroc_hi), width = 0.12) +
  geom_point(size = 2.4) +
  scale_color_manual(values = c(HepG2 = "#2166AC", LX2 = "#B2182B")) +
  coord_cartesian(ylim = c(0.25, 1)) +
  labs(x = NULL, y = "ChromBPNet AUROC\n(1-Mb block bootstrap)", color = NULL) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "top")
ggsave(file.path(OUT, "mpra_chrombpnet_context_transfer.pdf"), p,
       width = 3.4, height = 2.8, device = cairo_pdf)

o <- fread(file.path(IN, "mpra_overlap_summary.tsv"))
keep <- c("published_DAV_positions_overlapping_full_seqfunc_substrate",
          "published_DAV_overlaps_maxPIP_ge_0.1",
          "published_DAV_overlaps_maxPIP_ge_0.5",
          "published_DAV_overlaps_maxPIP_ge_0.9",
          "published_DAV_overlaps_retained_as_final_nomination_lead")
o <- o[metric %in% keep]
o[, label := factor(c("Any DAV", "PIP >= 0.1", "PIP >= 0.5", "PIP >= 0.9",
                      "Final lead")[match(metric, keep)], levels = rev(c(
                        "Any DAV", "PIP >= 0.1", "PIP >= 0.5", "PIP >= 0.9", "Final lead")))]
p2 <- ggplot(o, aes(value, label)) + geom_col(fill = "#4D4D4D", width = .68) +
  geom_text(aes(label = value), hjust = -0.2, size = 2.6) +
  scale_x_continuous(expand = expansion(mult = c(0, .12))) +
  labs(x = "Official Hu DAV positions", y = NULL) + theme_classic(base_size = 8)
ggsave(file.path(OUT, "mpra_nomination_overlap_audit.pdf"), p2,
       width = 3.2, height = 2.4, device = cairo_pdf)
message("[94] wrote MPRA supplementary PDFs to ", OUT)
