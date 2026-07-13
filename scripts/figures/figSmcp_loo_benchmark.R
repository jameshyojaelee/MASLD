#!/usr/bin/env Rscript
# figSmcp_loo_benchmark.R
# Leave-one-dataset-out (LOO) program stability and plain-NMF benchmark panels.

suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT <- FIGS_MCP_DIR
RD <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/mcp/reviewer_defense")

# Panel: LOO stability
loo <- fread(file.path(RD, "loo_stability_GSE244832_k16.tsv"))
loo[, full_program := as.character(full_program)]
loo[, pass05 := jaccard_top100 >= 0.5]
p_loo <- ggplot(loo[order(-jaccard_top100)],
                aes(reorder(full_program, jaccard_top100), jaccard_top100)) +
  geom_col(aes(fill = pass05)) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "gray40") +
  geom_hline(yintercept = 0.3, linetype = 3, color = "gray60") +
  scale_fill_manual(values = c(`TRUE` = masld_colors$up, `FALSE` = "gray70")) +
  labs(x = "cNMF program (full-atlas)", y = "Top-100 Jaccard vs LOO fit",
       fill = "Jaccard ≥ 0.5") +
  coord_flip() +
  theme_minimal(base_size = 6)
ggsave(file.path(OUT, "figSmcp_f_loo_stability.pdf"), p_loo, width = 6, height = 4, device = cairo_pdf)
message("[caption] LOO stability: GSE244832 excluded (62.5% of cells)")
cat("[figS] panel F (LOO) written\n")

# Panel: Plain NMF vs cNMF
pn <- fread(file.path(BASE, "Analysis/SingleCell/results_gpu_v2/mcp/benchmarks/plain_nmf_vs_cnmf.tsv"))
cophenetics <- data.table(
  method = c("cNMF (consensus)", "Plain NMF"),
  cophenetic = c(0.951, pn$cophenetic[1])
)
p_pn <- ggplot(cophenetics, aes(method, cophenetic)) +
  geom_col(fill = masld_colors$up, alpha = 0.7) +
  geom_hline(yintercept = 0.95, linetype = 2, color = "gray40") +
  coord_cartesian(ylim = c(0.9, 1.0)) +
  labs(x = NULL, y = "Cophenetic correlation (k=16)") +
  theme_minimal(base_size = 6)
ggsave(file.path(OUT, "figSmcp_g_cnmf_vs_plain.pdf"), p_pn, width = 4, height = 3, device = cairo_pdf)
message("[caption] cNMF vs plain NMF at 20 replicates")
cat("[figS] panel G (cNMF vs plain NMF) written\n")
cat("[figSmcp_loo_benchmark] DONE\n")
