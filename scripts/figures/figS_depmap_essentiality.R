#!/usr/bin/env Rscript
# Supplementary Figure: DepMap Cas9 Essentiality (Liver)
# 2 panels: (a) Chronos score distribution, (b) Essentiality vs DEG status

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT <- file.path(FIGS08_DIR, "figS_depmap_essentiality.pdf")
dir.create(file.path(FIGS08_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)

# --- Load data ----------------------------------------------------------------
cat("Loading atlas...\n")
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"))

# If atlas hasn't been rebuilt yet, compute essentiality inline
if (!"essentiality_chronos" %in% names(atlas) ||
    sum(!is.na(atlas$essentiality_chronos)) < 1000) {
  cat("Atlas essentiality column missing or sparse — computing from raw DepMap...\n")
  models <- fread(file.path(BASE, "Analysis/downstream_analysis/essentiality/Model.csv"))
  liver_ids <- models[OncotreeLineage == "Liver", ModelID]

  depmap <- fread(file.path(BASE, "Analysis/downstream_analysis/essentiality/CRISPRGeneEffect.csv"))
  setnames(depmap, names(depmap)[1], "ModelID")
  liver_ids <- intersect(liver_ids, depmap$ModelID)

  depmap_liver <- depmap[ModelID %in% liver_ids]
  gene_cols <- setdiff(names(depmap_liver), "ModelID")
  gene_symbols <- sub(" \\(.*", "", gene_cols)

  scores <- depmap_liver[, lapply(.SD, function(x) mean(x, na.rm = TRUE)),
                         .SDcols = gene_cols]

  ess_dt <- data.table(
    human_symbol         = gene_symbols,
    essentiality_chronos = as.numeric(scores[1, ])
  )
  ess_dt[, is_essential := essentiality_chronos < -0.5]
  ess_dt <- ess_dt[!duplicated(human_symbol)]

  atlas <- merge(atlas, ess_dt, by = "human_symbol", all.x = TRUE, suffixes = c(".old", ""))
  # Drop old columns if they exist
  old_cols <- grep("\\.old$", names(atlas), value = TRUE)
  if (length(old_cols)) atlas[, (old_cols) := NULL]

  rm(depmap, depmap_liver, models, scores); gc(verbose = FALSE)
}

# Classify DEG status using padj < 0.05 + |logFC| > 0.5
atlas[, is_deg_flag := is_dream_deg(atlas)]
atlas[, deg_status := fifelse(
  is_deg_flag & !is.na(bulk_logFC) & bulk_logFC > 0, "Upregulated",
  fifelse(
    is_deg_flag & !is.na(bulk_logFC) & bulk_logFC < 0, "Downregulated",
    "Not significant"
  )
)]
atlas[, is_deg_flag := NULL]

# Subset to genes with essentiality data
ess_genes <- atlas[!is.na(essentiality_chronos)]
cat("Genes with essentiality data:", nrow(ess_genes), "\n")
cat("Essential (Chronos < -0.5):", sum(ess_genes$is_essential, na.rm = TRUE), "\n")

# --- Panel a: Chronos score distribution --------------------------------------
pa <- ggplot(ess_genes, aes(x = essentiality_chronos)) +
  geom_histogram(aes(y = after_stat(density)), bins = 80,
                 fill = masld_colors$ns, color = NA, alpha = 0.8) +
  geom_density(linewidth = 0.4, color = "black") +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = masld_colors$up,
             linewidth = 0.4) +
  annotate("text", x = -0.55, y = Inf, label = "Essential\nthreshold",
           hjust = 1, vjust = 1.2, size = 2, color = masld_colors$up) +
  annotate("text", x = -4.3, y = Inf,
           label = sprintf("n = %s genes\n%s essential (%.1f%%)",
                           format(nrow(ess_genes), big.mark = ","),
                           format(sum(ess_genes$is_essential, na.rm = TRUE), big.mark = ","),
                           100 * mean(ess_genes$is_essential, na.rm = TRUE)),
           hjust = 0, vjust = 1.2, size = 2, color = "gray30") +
  scale_x_continuous(limits = c(-4.5, 1), breaks = seq(-4, 1, 1)) +
  labs(x = "DepMap Chronos score (liver cell lines)",
       y = "Density",
       tag = "a") +
  theme_masld()

# --- Panel b: Essentiality by DEG status --------------------------------------
# Compute counts for a stacked bar
deg_ess_counts <- ess_genes[deg_status != "Not significant",
                            .(n_essential = sum(is_essential, na.rm = TRUE),
                              n_total = .N),
                            by = deg_status]
deg_ess_counts[, n_non_essential := n_total - n_essential]
deg_ess_counts[, pct_essential := 100 * n_essential / n_total]

# Background rate
bg_rate <- 100 * mean(ess_genes$is_essential, na.rm = TRUE)

# Long format for plotting
plot_dt <- melt(deg_ess_counts, id.vars = "deg_status",
                measure.vars = c("n_essential", "n_non_essential"),
                variable.name = "ess_class", value.name = "count")
plot_dt[, ess_class := fifelse(ess_class == "n_essential", "Essential", "Non-essential")]
plot_dt[, ess_class := factor(ess_class, levels = c("Non-essential", "Essential"))]
plot_dt[, deg_status := factor(deg_status, levels = c("Upregulated", "Downregulated"))]

pb <- ggplot(plot_dt, aes(x = deg_status, y = count, fill = ess_class)) +
  geom_col(width = 0.6) +
  annotate("text",
           x = c(1, 2),
           y = deg_ess_counts[match(c("Upregulated", "Downregulated"), deg_status)]$n_total,
           label = sprintf("%.1f%%", deg_ess_counts[match(c("Upregulated", "Downregulated"), deg_status)]$pct_essential),
           vjust = -0.5, hjust = 0.5, size = 2, color = masld_colors$up) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c("Non-essential" = masld_colors$ns,
                               "Essential" = masld_colors$up),
                    name = "Essentiality") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "MASLD DEG status (padj < 0.05, |LFC| > 0.5)",
       y = "Number of genes",
       tag = "b") +
  theme_masld() +
  theme(legend.position = c(0.8, 0.85),
        legend.background = element_blank())

# --- Combine and save ---------------------------------------------------------
fig <- pa | pb
fig <- fig + plot_layout(widths = c(1.2, 1))

dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
ggsave(OUT, fig, width = fig_full_width, height = 3, device = pdf_device)
cat("Saved:", OUT, "\n")
