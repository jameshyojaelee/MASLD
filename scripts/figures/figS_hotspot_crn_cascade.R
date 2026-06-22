#!/usr/bin/env Rscript
# ============================================================================
# figS_hotspot_crn_cascade.R
#
# KEY MESSAGE: Top hotspot genes for the 5 canonical hepatocyte progression
# modules shown in Fig 2 (Hep-19, 20, 24, 26, 27).
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- FIGS_HOTSPOT_PANELS_DIR
DATA_DIR  <- FIGS_HOTSPOT_DATA_DIR
PREFIX    <- "crn_"

HS_RES <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")

read_topgenes <- function(ct, mod, n = 10) {
  f <- file.path(HS_RES, ct, "module_genes.tsv")
  if (!file.exists(f)) { message("WARN: ", f, " not found"); return(NULL) }
  d <- fread(f)
  d[module == mod][order(-weight)][1:min(n, .N)]
}

# The 5 canonical hepatocyte programs from Fig 2, in display order
FIVE_MODS <- list(
  list(ct = "hepatocytes", mod = 19,
       label = "Hep-19\nHNF4A-AS1 · EHHADH · BDH1\n(HNF4A identity, loss)"),
  list(ct = "hepatocytes", mod = 27,
       label = "Hep-27\nCFH · NR1H4 · C8A\n(complement/FXR, loss)"),
  list(ct = "hepatocytes", mod = 20,
       label = "Hep-20\nBICC1 · GLS · ITGAV\n(glutamine/TGFβ, gain)"),
  list(ct = "hepatocytes", mod = 24,
       label = "Hep-24\nSOD2 · HKDC1 · SQSTM1\n(NRF2 antioxidant, gain)"),
  list(ct = "hepatocytes", mod = 26,
       label = "Hep-26\nJUN · ATF3 · SERPINE1\n(AP-1 injury, gain)")
)

gene_q <- rbindlist(lapply(FIVE_MODS, function(x) {
  d <- read_topgenes(x$ct, x$mod, n = 10)
  if (is.null(d)) return(NULL)
  d[, mod_label := x$label]
  d[, gene_label := gene]
  d
}))

gene_q[, mod_label := factor(mod_label, levels = sapply(FIVE_MODS, `[[`, "label"))]
gene_q[, gene_label := factor(gene_label,
                              levels = unique(gene_q[order(mod_label, weight), gene_label]))]

HEP_COLOR <- ct_palette[["Hepatocytes"]]

pC <- ggplot(gene_q, aes(weight, gene_label)) +
  geom_col(width = 0.7, fill = HEP_COLOR, color = "black", linewidth = 0.15) +
  facet_wrap(~ mod_label, nrow = 2, scales = "free_y") +
  labs(x = "Hotspot module weight", y = NULL) +
  theme_masld() + theme_pub() +
  theme(axis.text.y = element_text(size = 5),
        strip.text = element_text(size = PUB_AXIS_TITLE, face = "bold"))

# ----------------------------------------------------------------------------
# Save
# ----------------------------------------------------------------------------
save_fig(pC, file.path(PANEL_DIR, paste0(PREFIX, "C_five_hep_programs.pdf")),
         width = fig_full_width, height = 5.0)
save_fig(pC,
         file.path(FIGS_HOTSPOT_DIR, "figS_hotspot_crn_cascade.pdf"),
         width = fig_full_width, height = 5.0)

fwrite(gene_q, file.path(DATA_DIR, paste0(PREFIX, "C_five_hep_programs.csv")))

cat("Wrote figS_hotspot_crn_cascade: 5-program hepatocyte top-gene panel\n")
