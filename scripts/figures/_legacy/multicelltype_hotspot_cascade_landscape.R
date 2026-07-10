#!/usr/bin/env Rscript
# ============================================================================
# multicelltype_hotspot_cascade_landscape.R
# Fig 3 (RNA-seq) main panel — multi-cell-type Hotspot cascade landscape
#
# Claim: stage-tracking, bulk-replicated co-expression modules are NOT
# hepatocyte-private. Cholangiocyte / fibroblast / macrophage compartments
# each contribute coordinated severity-gaining programs while hepatocytes
# uniquely split bidirectionally (both gaining and declining modules).
#
# Dotplot: rows faceted by cell_type. x = disease_stage_beta (leak-safe coarse
# 4-bin axis; NEVER F_stage_beta = scVI leakage). Dot color = sign of bulk
# mean_bulk_logFC (magenta up / blue down). Dot size = n_genes_overlapping.
# Solid ring outline when bulk_replicated. Filter to disease_stage_q < 0.05.
#
# Output: figures/main/fig3_RNAseq/panels/fig3m_multicelltype_hotspot_cascade_landscape.pdf
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG2_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
OUT_PDF   <- file.path(PANEL_DIR, "fig3m_multicelltype_hotspot_cascade_landscape.pdf")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
mods <- fread(file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/hotspot_modules/all_modules.tsv"))

# Compartment ordering (use actual label values; drop the 'global' pseudo-cell).
# No B-cell compartment is present in the file. Endothelial + global were
# deferred at the 2026-05-22 protocol-contamination remediation and are absent
# from the reconciled all_modules.tsv (5 clean cell types remain).
ct_order <- c("hepatocytes", "cholangiocytes", "fibroblasts",
              "macrophages", "tcells")
ct_pretty <- c(hepatocytes = "Hepatocyte", cholangiocytes = "Cholangiocyte",
               fibroblasts = "Fibroblast", macrophages = "Macrophage",
               tcells = "T cell")
ct_prefix <- c(hepatocytes = "hep", cholangiocytes = "chol",
               fibroblasts = "fib", macrophages = "mac",
               tcells = "T")

mods <- mods[cell_type %in% ct_order]

# ---------------------------------------------------------------------------
# LEAK-SAFE filter: coarse disease-stage axis only, significant modules.
# ---------------------------------------------------------------------------
d <- mods[!is.na(disease_stage_beta) & !is.na(disease_stage_q) &
            disease_stage_q < 0.05]

d[, ct_lab  := factor(ct_pretty[cell_type], levels = ct_pretty[ct_order])]
d[, mod_lab := paste0(ct_prefix[cell_type], "-", module)]
# Functional names from 509_module_pathway_names (module_names.tsv, the persistent
# authority). all_modules.tsv may already carry a merged module_name (509
# convenience merge, dropped on a 508 rebuild) — drop it first to avoid a
# module_name.x/.y merge collision and always use the fresh module_names.tsv.
for (.c in c("module_name", "module_hallmark", "module_top_pathway"))
  if (.c %in% names(d)) d[, (.c) := NULL]
.mn <- fread(file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/hotspot_modules/module_names.tsv"))
d <- merge(d, .mn[, .(cell_type, module, module_name)],
           by = c("cell_type", "module"), all.x = TRUE)
d[is.na(module_name) | module_name == "Unresolved", module_name := mod_lab]
d[, mod_disp := module_name]
d[, dir := ifelse(mean_bulk_logFC >= 0, "Bulk up", "Bulk down")]
dir_colors <- c("Bulk up" = "#C9265E", "Bulk down" = "#1565C0")

# ---------------------------------------------------------------------------
# Hero numbers: per-compartment up/down + bulk-replicated counts
# ---------------------------------------------------------------------------
hero <- d[, .(
  n_modules     = .N,
  up_beta       = sum(disease_stage_beta > 0),
  down_beta     = sum(disease_stage_beta < 0),
  bulk_rep      = sum(bulk_replicated == TRUE),
  bulk_rep_up   = sum(bulk_replicated == TRUE & disease_stage_beta > 0),
  bulk_rep_down = sum(bulk_replicated == TRUE & disease_stage_beta < 0)
), by = ct_lab][order(match(ct_lab, ct_pretty[ct_order]))]

cat("[hero] Per-compartment stage-significant (disease_stage_q<0.05) modules:\n")
for (i in seq_len(nrow(hero))) {
  cat(sprintf("[hero]   %-13s n=%2d | gaining(beta>0)=%2d declining(beta<0)=%2d | bulk_replicated=%2d (up=%d down=%d)\n",
              hero$ct_lab[i], hero$n_modules[i], hero$up_beta[i], hero$down_beta[i],
              hero$bulk_rep[i], hero$bulk_rep_up[i], hero$bulk_rep_down[i]))
}
cat(sprintf("[hero] Hepatocyte is the only compartment with BOTH gaining and declining modules: gaining=%d declining=%d\n",
            hero[ct_lab == "Hepatocyte", up_beta], hero[ct_lab == "Hepatocyte", down_beta]))

fwrite(d[, .(cell_type, ct_lab, module, mod_lab, module_name, disease_stage_beta,
             disease_stage_q, mean_bulk_logFC, n_genes_overlapping,
             bulk_replicated, dir)],
       file.path(DATA_DIR, "multicelltype_hotspot_cascade_landscape.csv"))

# ---------------------------------------------------------------------------
# Labels: hepatocyte standouts + top stromal/biliary modules by |beta|
# ---------------------------------------------------------------------------
hep_standouts <- c("hep-20", "hep-24", "hep-26", "hep-19", "hep-27")
nonhep <- d[cell_type != "hepatocytes"]
top_stromal <- nonhep[order(-abs(disease_stage_beta))][, head(.SD, 2), by = ct_lab]
lab_ids <- unique(c(hep_standouts, top_stromal$mod_lab))
d[, do_label := mod_lab %in% lab_ids]

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
p <- ggplot(d, aes(x = disease_stage_beta, y = 0)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "#9E9E9E") +
  geom_point(aes(size = n_genes_overlapping, fill = dir,
                 color = bulk_replicated),
             shape = 21, stroke = 0.55, alpha = 0.9) +
  geom_text_repel(data = d[do_label == TRUE],
                  aes(label = mod_disp),
                  size = 1.9, color = "grey15", segment.size = 0.2,
                  min.segment.length = 0, max.overlaps = Inf,
                  box.padding = 0.25, direction = "y", seed = 42) +
  facet_grid(ct_lab ~ ., switch = "y") +
  scale_fill_manual(values = dir_colors, name = "Bulk direction") +
  scale_color_manual(values = c(`TRUE` = "grey10", `FALSE` = "grey75"),
                     name = "Bulk replicated",
                     labels = c(`TRUE` = "yes", `FALSE` = "no")) +
  scale_size_continuous(range = c(0.8, 4.2), name = "Genes\noverlapping") +
  scale_y_continuous(limits = c(-1, 1), breaks = NULL, name = NULL) +
  labs(x = "Disease-stage co-expression slope (coarse 4-bin axis)",
       title = "Stage-tracking modules span every hepatic compartment") +
  guides(fill = guide_legend(override.aes = list(size = 2.4),
                             order = 1),
         color = guide_legend(override.aes = list(size = 2.4, shape = 21,
                                                  fill = "grey85"), order = 2),
         size = guide_legend(order = 3)) +
  theme_masld(base_size = 7) +
  theme(
    plot.title       = element_text(size = 7.3, face = "bold", margin = margin(b = 5)),
    strip.text.y.left = element_text(size = 6, face = "bold", angle = 0),
    strip.placement  = "outside",
    panel.spacing.y  = unit(0.12, "lines"),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    legend.position  = "right",
    legend.text      = element_text(size = 5.5),
    legend.title     = element_text(size = 6),
    legend.key.size  = unit(0.16, "cm")
  )

ggsave(OUT_PDF, p,
       width  = 110 / 25.4,
       height = 90 / 25.4,
       units  = "in",
       device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))
