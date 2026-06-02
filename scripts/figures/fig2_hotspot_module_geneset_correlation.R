#!/usr/bin/env Rscript
# fig2_hotspot_module_geneset_correlation.R
#
# fgsea NES heatmap: 12 passing hepatocyte Hotspot modules vs Hallmark gene sets.
# Module gene weights used as ranking statistic.
#
# Output: figures/main/fig2_progression_sex/panels/fig2_hotspot_geneset_correlation.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(msigdbr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

HS_RES    <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")
PANEL_DIR <- file.path(FIG2_DIR, "panels")
OUT_PDF   <- file.path(PANEL_DIR, "fig2_hotspot_geneset_correlation.pdf")

# ── passing modules ───────────────────────────────────────────────────────────
mods <- fread(file.path(HS_RES, "all_modules.tsv"))
passing <- mods[cell_type == "hepatocytes" &
                progression_module == TRUE &
                stability_fail == FALSE,
                .(module, disease_stage_beta, disease_stage_q)]
setorder(passing, disease_stage_q)
passing[, module_label := paste0("Hep-", module,
                                  ifelse(disease_stage_beta > 0, " (+)", " (−)"))]
cat(sprintf("[modules] %d passing\n", nrow(passing)))

# ── gene universe + module membership ────────────────────────────────────────
autocorr <- fread(file.path(HS_RES, "hepatocytes/autocorr.tsv"))
autocorr  <- autocorr[!grepl("^ENSG", gene)]   # keep symbol-named genes only
universe  <- autocorr$gene
cat(sprintf("[universe] %d genes\n", length(universe)))

mod_genes <- fread(file.path(HS_RES, "hepatocytes/module_genes.tsv"))
mod_genes <- mod_genes[!grepl("^ENSG", gene)]

# ── Hallmark gene sets (intersected with universe) ────────────────────────────
cat("[msigdbr] loading Hallmark...\n")
h         <- msigdbr(species = "Homo sapiens", collection = "H")
pathways  <- split(h$gene_symbol, h$gs_name)
pathways  <- lapply(pathways, function(g) intersect(g, universe))
pathways  <- pathways[lengths(pathways) >= 10]
cat(sprintf("[msigdbr] %d Hallmark gene sets after universe filter\n", length(pathways)))

# ── overrepresentation analysis (Fisher's exact) per module × gene set ────────
N <- length(universe)

results <- rbindlist(lapply(passing$module, function(mod) {
  mod_list <- mod_genes[module == mod, unique(gene)]
  k <- length(mod_list)

  rbindlist(lapply(names(pathways), function(gs) {
    gs_genes <- pathways[[gs]]
    M  <- length(gs_genes)
    x  <- length(intersect(mod_list, gs_genes))
    # Fisher's exact: 2x2 table
    mat <- matrix(c(x, k - x, M - x, N - k - M + x), nrow = 2)
    p   <- fisher.test(mat, alternative = "greater")$p.value
    or  <- (x / (k - x)) / ((M - x) / (N - k - M + x))
    data.table(module = mod, pathway = gs, n_overlap = x,
               n_module = k, n_gs = M, odds_ratio = or, p = p)
  }))
}))

results[, padj := p.adjust(p, method = "BH"), by = module]
results[, log2or := log2(pmax(odds_ratio, 1e-3))]
cat(sprintf("[ORA] done — %d rows\n", nrow(results)))

# ── filter: padj < 0.05 in at least 1 module, min 2 overlapping genes ────────
sig_paths <- results[padj < 0.05 & n_overlap >= 2, unique(pathway)]
cat(sprintf("[filter] %d significant pathways\n", length(sig_paths)))

if (length(sig_paths) == 0) {
  sig_paths <- results[p < 0.05 & n_overlap >= 2, unique(pathway)]
  cat(sprintf("[filter] relaxed to nominal p<0.05: %d pathways\n", length(sig_paths)))
}

plot_dt <- results[pathway %in% sig_paths]
plot_dt <- merge(plot_dt, passing[, .(module, module_label)], by = "module")

plot_dt[, pathway_clean := gsub("HALLMARK_", "", pathway)]
plot_dt[, pathway_clean := gsub("_", " ", pathway_clean)]
plot_dt[, pathway_clean := tools::toTitleCase(tolower(pathway_clean))]

path_order <- plot_dt[, .(mean_lor = mean(log2or, na.rm = TRUE)), by = pathway_clean
                      ][order(mean_lor), pathway_clean]
plot_dt[, pathway_clean := factor(pathway_clean, levels = path_order)]
plot_dt[, module_label  := factor(module_label, levels = passing$module_label)]

plot_dt[, sig_mark := fcase(padj < 0.001, "***",
                             padj < 0.01,  "**",
                             padj < 0.05,  "*",
                             default = "")]

# ── plot ──────────────────────────────────────────────────────────────────────
lim <- ceiling(max(abs(plot_dt$log2or), na.rm = TRUE) * 10) / 10

p <- ggplot(plot_dt, aes(x = module_label, y = pathway_clean, fill = log2or)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sig_mark), size = 2.5, color = "black", vjust = 0.75) +
  scale_fill_gradient2(low      = masld_colors$down,
                       mid      = "white",
                       high     = masld_colors$up,
                       midpoint = 0,
                       limits   = c(0, lim),
                       name     = "log₂ OR") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title    = "Hepatocyte Hotspot modules — Hallmark enrichment",
       subtitle = "Fisher's exact log₂OR  |  * padj<0.05  ** padj<0.01  *** padj<0.001") +
  theme_masld(base_size = 10) +
  theme(axis.text.x      = element_text(size = 9, angle = 40, hjust = 1),
        axis.text.y      = element_text(size = 8.5),
        legend.key.height = unit(1.2, "cm"),
        legend.key.width  = unit(0.35, "cm"),
        panel.grid        = element_blank(),
        axis.line         = element_blank(),
        axis.ticks        = element_blank(),
        plot.title        = element_text(size = 11, face = "bold"),
        plot.subtitle     = element_text(size = 8, color = "gray40"))

n_gs  <- length(unique(plot_dt$pathway_clean))
n_mod <- nrow(passing)
ggsave(OUT_PDF, p,
       width  = 3 + n_mod * 0.6,
       height = 2 + n_gs  * 0.3,
       device = cairo_pdf)
cat(sprintf("Saved: %s  (%d modules × %d pathways)\n", OUT_PDF, n_mod, n_gs))
