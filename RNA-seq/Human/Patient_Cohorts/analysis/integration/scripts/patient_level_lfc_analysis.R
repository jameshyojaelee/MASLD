#!/usr/bin/env Rscript
# patient_level_lfc_analysis.R
# ---------------------------------------------------------------------------
# Per-patient LFC analysis: For each gene, compute log2FC of each MASLD
# patient relative to within-dataset healthy control mean. Then count how
# many genes are consistently up/downregulated across X% of patients at
# various LFC cutoffs.
#
# Approach:
#   1. Load merged_dge.rds (RLE-normalized DGEList, 1,444 QC-passing samples)
#   2. Compute logCPM per sample (using edgeR::cpm with library sizes + norm factors)
#   3. For each dataset, compute mean logCPM of Control samples
#   4. For each Disease sample, compute per-gene LFC = logCPM_patient - mean_logCPM_controls
#   5. Exclude GSE167523 (same as dream)
#   6. Sweep LFC cutoffs × patient-percentage thresholds
#
# Output:
#   - results/integration/patient_lfc_sweep.csv (summary table)
#   - results/integration/patient_lfc_sweep.pdf (visualization)
#   - results/integration/patient_lfc_matrix.csv.gz (full per-patient LFC matrix)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(edgeR)
  library(ggplot2)
  library(patchwork)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT  <- file.path(BASE, "analysis/integration")
RDIR <- file.path(INT, "results/integration")

# --- Load merged DGE (same object used by dream) ---
cat("Loading merged_dge.rds...\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
cat("  Samples:", ncol(dge), "  Genes:", nrow(dge), "\n")

# --- Mega-analysis cohort selection: enforce config/human_datasets.yaml include_in_mega ---
ycfg <- yaml::read_yaml(file.path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design", "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
keep <- dge$samples$dataset %in% mega_cohorts
dge <- dge[, keep]
cat("Mega cohorts (yaml):", paste(mega_cohorts, collapse = ", "), ":", ncol(dge), "samples\n")

# --- Compute logCPM (uses library size + RLE norm factors) ---
cat("Computing logCPM...\n")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)  # log2(CPM + prior)

# --- Build sample info ---
info <- data.table(
  sample_id = colnames(dge),
  dataset   = dge$samples$dataset,
  group     = as.character(dge$samples$group_binary)
)

cat("\nDataset × Group distribution:\n")
print(dcast(info, dataset ~ group, fun.aggregate = length))

datasets_with_controls <- info[, .(n_ctrl = sum(group == "Control"),
                                    n_dis  = sum(group == "Disease")), by = dataset]
cat("\nDatasets with controls:\n")
print(datasets_with_controls)

# Only use datasets that have both Control and Disease samples
valid_datasets <- datasets_with_controls[n_ctrl > 0 & n_dis > 0, dataset]
cat("\nValid datasets (have both groups):", paste(valid_datasets, collapse = ", "), "\n")

# --- Compute per-patient LFC relative to within-dataset control mean ---
cat("\nComputing per-patient LFC (within-dataset normalization)...\n")

disease_samples <- info[group == "Disease" & dataset %in% valid_datasets, sample_id]
cat("  Disease samples:", length(disease_samples), "\n")

# Pre-compute per-dataset control means
ctrl_means <- list()
for (ds in valid_datasets) {
  ctrl_ids <- info[dataset == ds & group == "Control", sample_id]
  if (length(ctrl_ids) == 1) {
    ctrl_means[[ds]] <- logcpm[, ctrl_ids]
  } else {
    ctrl_means[[ds]] <- rowMeans(logcpm[, ctrl_ids, drop = FALSE])
  }
}

# Compute LFC matrix: genes × disease samples
patient_lfc <- matrix(NA_real_, nrow = nrow(logcpm), ncol = length(disease_samples),
                      dimnames = list(rownames(logcpm), disease_samples))

for (ds in valid_datasets) {
  ds_disease <- info[dataset == ds & group == "Disease", sample_id]
  for (sid in ds_disease) {
    patient_lfc[, sid] <- logcpm[, sid] - ctrl_means[[ds]]
  }
}

cat("  Patient LFC matrix:", nrow(patient_lfc), "genes ×", ncol(patient_lfc), "patients\n")

# --- Also load dream results for significance filter ---
dream <- fread(file.path(RDIR, "dream_results.csv"))
setnames(dream, "adj.P.Val", "padj", skip_absent = TRUE)

# --- Sweep: LFC cutoffs × patient percentage thresholds ---
lfc_cutoffs <- c(0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0)
pct_thresholds <- c(50, 60, 70, 80, 90)

n_patients <- ncol(patient_lfc)
cat("\nTotal disease patients:", n_patients, "\n")

results <- list()
for (lfc_cut in lfc_cutoffs) {
  # For each gene, count patients with |LFC| > cutoff and in correct direction
  n_up   <- rowSums(patient_lfc > lfc_cut)
  n_down <- rowSums(patient_lfc < -lfc_cut)

  for (pct in pct_thresholds) {
    min_patients <- ceiling(n_patients * pct / 100)

    # Genes upregulated in >= pct% of patients
    genes_up <- sum(n_up >= min_patients)
    # Genes downregulated in >= pct% of patients
    genes_down <- sum(n_down >= min_patients)

    # Also count among dream-significant genes only (padj < 0.1)
    sig_genes <- dream[padj < 0.1, gene]
    sig_idx <- rownames(patient_lfc) %in% sig_genes
    genes_up_sig <- sum(n_up[sig_idx] >= min_patients)
    genes_down_sig <- sum(n_down[sig_idx] >= min_patients)

    results[[length(results) + 1]] <- data.table(
      lfc_cutoff = lfc_cut,
      pct_threshold = pct,
      min_patients = min_patients,
      n_up_all = genes_up,
      n_down_all = genes_down,
      n_total_all = genes_up + genes_down,
      n_up_sig = genes_up_sig,
      n_down_sig = genes_down_sig,
      n_total_sig = genes_up_sig + genes_down_sig
    )
  }
}
sweep <- rbindlist(results)

cat("\n===== SWEEP RESULTS =====\n")
print(sweep)

# --- Save sweep table ---
fwrite(sweep, file.path(RDIR, "patient_lfc_sweep.csv"))
cat("\nSaved:", file.path(RDIR, "patient_lfc_sweep.csv"), "\n")

# --- Save compressed patient LFC matrix ---
cat("Saving patient LFC matrix (compressed)...\n")
lfc_dt <- as.data.table(patient_lfc, keep.rownames = "gene")
fwrite(lfc_dt, file.path(RDIR, "patient_lfc_matrix.csv.gz"))
cat("Saved:", file.path(RDIR, "patient_lfc_matrix.csv.gz"), "\n")

# --- Per-gene summary stats ---
gene_summary <- data.table(
  gene = rownames(patient_lfc),
  mean_lfc = rowMeans(patient_lfc),
  median_lfc = apply(patient_lfc, 1, median),
  sd_lfc = apply(patient_lfc, 1, sd),
  pct_up_any = rowSums(patient_lfc > 0) / n_patients * 100,
  pct_up_0.5 = rowSums(patient_lfc > 0.5) / n_patients * 100,
  pct_up_1.0 = rowSums(patient_lfc > 1.0) / n_patients * 100,
  pct_down_any = rowSums(patient_lfc < 0) / n_patients * 100,
  pct_down_0.5 = rowSums(patient_lfc < -0.5) / n_patients * 100,
  pct_down_1.0 = rowSums(patient_lfc < -1.0) / n_patients * 100
)
# Merge dream stats
gene_summary <- merge(gene_summary, dream[, .(gene, dream_logFC = logFC, dream_padj = padj)],
                      by = "gene", all.x = TRUE)
fwrite(gene_summary, file.path(RDIR, "patient_lfc_gene_summary.csv"))
cat("Saved:", file.path(RDIR, "patient_lfc_gene_summary.csv"), "\n")

# --- Visualization ---
cat("\nGenerating figures...\n")

# Theme — compact canvas with larger relative text
theme_pub <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.position = "right",
    plot.margin = margin(4, 4, 4, 4)
  )

# ---- Panel A: Sweep heatmap (all genes) ----
sweep_long <- melt(sweep[, .(lfc_cutoff, pct_threshold, n_up_all, n_down_all)],
                   id.vars = c("lfc_cutoff", "pct_threshold"),
                   variable.name = "direction", value.name = "n_genes")
sweep_long[, direction := fifelse(direction == "n_up_all", "Upregulated", "Downregulated")]
sweep_long[, pct_label := paste0(pct_threshold, "%")]
sweep_long[, lfc_label := as.character(lfc_cutoff)]

p_heatmap_all <- ggplot(sweep_long, aes(x = factor(lfc_cutoff), y = factor(pct_threshold),
                                         fill = n_genes)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = format(n_genes, big.mark = ",")), size = 3.5, color = "black") +
  facet_wrap(~ direction) +
  scale_fill_viridis_c(option = "plasma", trans = "log1p", name = "Genes",
                        labels = scales::comma) +
  labs(x = expression("Log"[2]*"FC cutoff"), y = "% of patients",
       title = "All genes: consistently dysregulated across patients") +
  theme_pub

# ---- Panel B: Sweep heatmap (dream-significant genes only) ----
sweep_long_sig <- melt(sweep[, .(lfc_cutoff, pct_threshold, n_up_sig, n_down_sig)],
                       id.vars = c("lfc_cutoff", "pct_threshold"),
                       variable.name = "direction", value.name = "n_genes")
sweep_long_sig[, direction := fifelse(direction == "n_up_sig", "Upregulated", "Downregulated")]

p_heatmap_sig <- ggplot(sweep_long_sig, aes(x = factor(lfc_cutoff), y = factor(pct_threshold),
                                              fill = n_genes)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = format(n_genes, big.mark = ",")), size = 3.5, color = "black") +
  facet_wrap(~ direction) +
  scale_fill_viridis_c(option = "plasma", trans = "log1p", name = "Genes",
                        labels = scales::comma) +
  labs(x = expression("Log"[2]*"FC cutoff"), y = "% of patients",
       title = "Integrated-significant genes (padj < 0.1): consistently dysregulated") +
  theme_pub

# ---- Panel C: Line plot — gene count vs LFC cutoff, colored by patient % ----
sweep_line <- melt(sweep[, .(lfc_cutoff, pct_threshold, n_up_all, n_down_all)],
                   id.vars = c("lfc_cutoff", "pct_threshold"),
                   variable.name = "direction", value.name = "n_genes")
sweep_line[, direction := fifelse(direction == "n_up_all", "Upregulated", "Downregulated")]
sweep_line[, pct_label := factor(paste0(pct_threshold, "%"),
                                  levels = paste0(sort(pct_thresholds), "%"))]

p_line <- ggplot(sweep_line, aes(x = lfc_cutoff, y = n_genes, color = pct_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ direction, scales = "free_y") +
  scale_color_brewer(palette = "Set1", name = "Patient %") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = expression("|Log"[2]*"FC| cutoff"),
       y = "Number of genes",
       title = "Gene count vs fold-change cutoff by patient consistency") +
  theme_pub

# ---- Panel D: Distribution of patient-level LFC for top DEGs ----
# Pick top 20 dream DEGs by absolute logFC
top_genes <- dream[padj < 0.1][order(-abs(logFC))][1:20, gene]
top_lfc <- patient_lfc[top_genes, , drop = FALSE]
top_lfc_long <- melt(as.data.table(top_lfc, keep.rownames = "gene"),
                     id.vars = "gene", variable.name = "sample", value.name = "lfc")
# Add dream direction
top_lfc_long <- merge(top_lfc_long, dream[, .(gene, dream_logFC = logFC)], by = "gene")
top_lfc_long[, gene_label := factor(gene, levels = top_genes)]

p_violin <- ggplot(top_lfc_long, aes(x = reorder(gene, -abs(dream_logFC)), y = lfc,
                                      fill = ifelse(dream_logFC > 0, "Up", "Down"))) +
  geom_violin(scale = "width", alpha = 0.7, linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8"), name = "Integrated direction") +
  labs(x = NULL, y = expression("Patient-level log"[2]*"FC"),
       title = "Per-patient fold-change distribution (top 20 integrated DEGs by |LFC|)") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# ---- Panel E: Scatter — dream LFC vs median patient LFC ----
gene_summary_sig <- gene_summary[gene %in% dream[padj < 0.1, gene]]

p_scatter <- ggplot(gene_summary_sig, aes(x = dream_logFC, y = median_lfc)) +
  geom_hex(bins = 80) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_fill_viridis_c(option = "inferno", trans = "log1p", name = "Genes") +
  labs(x = expression("Integrated mega-analysis log"[2]*"FC"),
       y = expression("Median patient-level log"[2]*"FC"),
       title = "Population vs individual effect sizes (integrated-significant genes)") +
  coord_equal() +
  theme_pub

# ---- Panel F: Histogram of patient consistency (pct_up for upregulated DEGs) ----
up_degs <- gene_summary[dream_padj < 0.1 & dream_logFC > 0]
down_degs <- gene_summary[dream_padj < 0.1 & dream_logFC < 0]

hist_data <- rbind(
  up_degs[, .(gene, pct = pct_up_any, direction = "Upregulated DEGs\n(% patients with LFC > 0)")],
  down_degs[, .(gene, pct = pct_down_any, direction = "Downregulated DEGs\n(% patients with LFC < 0)")]
)

p_hist <- ggplot(hist_data, aes(x = pct, fill = direction)) +
  geom_histogram(binwidth = 2, alpha = 0.8, color = "white", linewidth = 0.2) +
  facet_wrap(~ direction, scales = "free_y") +
  scale_fill_manual(values = c("Upregulated DEGs\n(% patients with LFC > 0)" = "#E41A1C",
                                "Downregulated DEGs\n(% patients with LFC < 0)" = "#377EB8")) +
  labs(x = "% of patients with concordant direction",
       y = "Number of genes",
       title = "How consistently are DEGs dysregulated across patients?") +
  theme_pub +
  theme(legend.position = "none")

# ---- Save individual panels ----
FIGDIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/supplementary/figS02_progression"
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
PANELS_DIR <- file.path(FIGDIR, "panels")
dir.create(PANELS_DIR, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(PANELS_DIR, "heatmap_all.pdf"), p_heatmap_all,
       width = 8, height = 4, device = cairo_pdf)
ggsave(file.path(PANELS_DIR, "heatmap_sig.pdf"), p_heatmap_sig,
       width = 8, height = 4, device = cairo_pdf)
ggsave(file.path(PANELS_DIR, "line_plot.pdf"), p_line,
       width = 8, height = 4, device = cairo_pdf)
ggsave(file.path(PANELS_DIR, "violin_top_degs.pdf"), p_violin,
       width = 8, height = 4.5, device = cairo_pdf)
ggsave(file.path(PANELS_DIR, "scatter_dream_vs_patient.pdf"), p_scatter,
       width = 5.5, height = 5, device = cairo_pdf)
ggsave(file.path(PANELS_DIR, "histogram_consistency.pdf"), p_hist,
       width = 8, height = 4, device = cairo_pdf)

# ---- Combined multi-panel (compact) ----
combined <- (p_heatmap_all | p_line) / (p_heatmap_sig | p_hist) / (p_scatter | p_violin) +
  plot_annotation(
    title = "Patient-level fold-change analysis: MASLD vs healthy controls",
    subtitle = sprintf("%s disease patients across %d datasets | %s genes tested",
                       format(n_patients, big.mark = ","),
                       length(valid_datasets),
                       format(nrow(patient_lfc), big.mark = ",")),
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 12, color = "grey30"))
  )

ggsave(file.path(FIGDIR, "patient_lfc_sweep.pdf"), combined,
       width = 14, height = 14, device = cairo_pdf)
cat("Saved:", file.path(FIGDIR, "patient_lfc_sweep.pdf"), "\n")

# --- Print key highlights ---
cat("\n===== KEY HIGHLIGHTS =====\n")
cat(sprintf("Disease patients: %d across %d datasets\n", n_patients, length(valid_datasets)))
cat(sprintf("Genes tested: %s\n", format(nrow(patient_lfc), big.mark = ",")))

for (pct in c(50, 70, 90)) {
  row_lfc0 <- sweep[pct_threshold == pct & lfc_cutoff == 0]
  row_lfc05 <- sweep[pct_threshold == pct & lfc_cutoff == 0.5]
  row_lfc1 <- sweep[pct_threshold == pct & lfc_cutoff == 1.0]
  cat(sprintf("\n>=%d%% of patients:\n", pct))
  cat(sprintf("  |LFC| > 0:   %s up, %s down\n", format(row_lfc0$n_up_all, big.mark = ","), format(row_lfc0$n_down_all, big.mark = ",")))
  cat(sprintf("  |LFC| > 0.5: %s up, %s down\n", format(row_lfc05$n_up_all, big.mark = ","), format(row_lfc05$n_down_all, big.mark = ",")))
  cat(sprintf("  |LFC| > 1.0: %s up, %s down\n", format(row_lfc1$n_up_all, big.mark = ","), format(row_lfc1$n_down_all, big.mark = ",")))
}

# Correlation between dream LFC and median patient LFC
r <- cor(gene_summary$dream_logFC, gene_summary$median_lfc, use = "complete.obs", method = "spearman")
cat(sprintf("\nSpearman rho (dream LFC vs median patient LFC): %.3f\n", r))

cat("\nDone.\n")
