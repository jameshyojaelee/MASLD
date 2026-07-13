#!/usr/bin/env Rscript
# hotspot_module_concordance_plane.R  — CANDIDATE panel (evidence-plane grammar; N=117 = sparse cloud)
# All 117 Hotspot co-expression modules (5 cell types) on a single-cell vs bulk
# concordance plane: x = single-cell disease-stage slope, y = bulk mean log2FC of the
# module's genes. Extends the hepatocyte-only N=30 scatter to every cell type.
# Original hep_module_sc_vs_bulk_concordance_scatter.R UNTOUCHED.
#
# Output: FIG2_DIR/panels/hotspot_module_concordance.pdf (exploratory; no Fig-3 letter)

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
set.seed(42)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- file.path(FIG2_DIR, "panels"); DATA_DIR <- file.path(PANEL_DIR, "data")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
lab_size <- 6 / ggplot2::.pt

d <- fread(file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules/all_modules.tsv"),
  select = c("cell_type","module","disease_stage_beta","disease_stage_q","mean_bulk_logFC",
             "bulk_replicated","module_name"))
d <- d[!is.na(disease_stage_beta) & !is.na(mean_bulk_logFC)]
disp <- c(hepatocytes="Hepatocytes", macrophages="Macrophages", fibroblasts="Fibroblasts",
          cholangiocytes="Cholangiocytes", tcells="T cells")
d[, ct := disp[cell_type]]
d[, progressing := !is.na(disease_stage_q) & disease_stage_q < 0.05]
# italic gene-symbol label = the marker gene(s) in the module_name parentheses
d[, gene_lab := ifelse(grepl("\\(", module_name), sub(".*\\(([^)]*)\\).*", "\\1", module_name), module_name)]
rho <- cor(d$disease_stage_beta, d$mean_bulk_logFC, method = "spearman")
cat(sprintf("[hotspot concordance] %d modules, Spearman rho(sc-beta, bulk-logFC)=%.3f\n", nrow(d), rho))

# anchors = top |beta| modules (the disease-associated heroes)
lab <- d[order(-abs(disease_stage_beta))][1:11]

CT_COL <- c(Hepatocytes="#0D47A1", Macrophages="#C2185B", Fibroblasts="#F57F17",
            Cholangiocytes="#1565C0", `T cells`="#7B1FA2")
p <- ggplot(d, aes(disease_stage_beta, mean_bulk_logFC)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "gray85") +
  geom_vline(xintercept = 0, linewidth = 0.2, colour = "gray85") +
  geom_point(aes(colour = ct, alpha = progressing, size = progressing), shape = 16) +
  geom_text_repel(data = lab, aes(label = gene_lab, colour = ct),
    colour = "black", size = lab_size, fontface = "italic",
    box.padding = 0.45, point.padding = 0.3, segment.size = 0.2, segment.color = "gray60",
    min.segment.length = 0, max.overlaps = Inf, seed = 42, force = 9,
    bg.color = "white", bg.r = 0.12) +
  scale_colour_manual(values = CT_COL, name = NULL) +
  scale_alpha_manual(values = c(`TRUE` = 0.95, `FALSE` = 0.35), guide = "none") +
  scale_size_manual(values = c(`TRUE` = 1.7, `FALSE` = 1.0), guide = "none") +
  guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1))) +
  labs(x = "Single-cell disease-stage slope (Hotspot module)",
       y = expression("Bulk mean " * log[2] * "FC of module genes")) +
  theme_masld_compact() + theme(legend.position = "right", legend.key.size = unit(0.28, "cm"))

message(sprintf(paste0("CAPTION (Hotspot module concordance): all %d Hotspot modules across 5 cell ",
  "types; x = single-cell disease-stage slope, y = bulk mean log2FC of the module's genes. Filled/",
  "larger = disease-progressing (q<0.05). Concordant modules (rise/fall together) sit in the upper-",
  "right / lower-left quadrants (Spearman rho=%.2f). Extends the hepatocyte-only N=30 view to every ",
  "cell type; N=117 reads as a sparse cloud."), nrow(d), rho))

out <- file.path(PANEL_DIR, "hotspot_module_concordance.pdf")
save_fig(p, out, width = fig_half_width + 0.6, height = 3.0)
cat("[hotspot concordance] Saved:", out, "\n")
fwrite(d[order(-abs(disease_stage_beta))][, .(cell_type, module, module_name,
        disease_stage_beta = round(disease_stage_beta,3), mean_bulk_logFC = round(mean_bulk_logFC,3),
        progressing, bulk_replicated)], file.path(DATA_DIR, "hotspot_module_concordance.csv"))
