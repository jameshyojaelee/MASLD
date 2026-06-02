#!/usr/bin/env Rscript
# figS_control_variability.R
# Diagnostic panels: within-control expression variability for canonical DEGs.
# Motivates the question of whether an expression filter reduces control noise.
#
# Panels:
#   V1  Within-control SD distribution per dataset (violin)
#   V2  |LFC| distribution: control vs disease (density overlay)
#   V3  Dream |logFC| vs within-control SD per gene (scatter)
#   V4  Mean expression (logCPM) vs within-control SD per gene (scatter)
#
# Output: figures/supplementary/figS_lfc_sensitivity/

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(scales)
  library(edgeR)
  library(yaml)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

OUTDIR  <- file.path(BASE, "figures/supplementary/figS_lfc_sensitivity")
INT     <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")

theme_pub <- theme_minimal(base_size = 11) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.3),
        axis.ticks = element_line(colour = "black", linewidth = 0.3),
        legend.background = element_blank(), legend.key = element_blank(),
        strip.background = element_blank(),
        strip.text  = element_text(face = "bold", size = 10),
        plot.title  = element_text(face = "bold", size = 12),
        axis.title  = element_text(size = 10),
        axis.text   = element_text(size = 9),
        legend.text = element_text(size = 9),
        plot.margin = margin(8, 10, 8, 8))
theme_set(theme_pub)

col_dis  <- "#00695C"   # teal — disease
col_ctrl <- "#757575"   # gray — control

# ── Load data ─────────────────────────────────────────────────────────────────
message("Loading dream STAR results...")
dream <- fread(file.path(INT, "dream_results_ashr.csv"))
can_up <- dream[padj < 0.05 & logFC >  0.5, gene]
can_dn <- dream[padj < 0.05 & logFC < -0.5, gene]
can_degs <- c(can_up, can_dn)
dream_fc <- setNames(dream$logFC, dream$gene)
message(sprintf("  Canonical DEGs: %d up + %d down = %d total",
                length(can_up), length(can_dn), length(can_degs)))

message("Loading STAR DGE...")
dge <- readRDS(file.path(INT, "merged_dge.rds"))
ycfg <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
mega_c <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
dge  <- dge[, dge$samples$dataset %in% mega_c]
lcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

# ── TPM: counts / length_kb, then per-million normalise ──────────────────────
message("Computing TPM...")
gene_len <- fread(
  file.path(BASE, "RNA-seq/Human/Patient_Cohorts/results/GSE126848/counts/featurecounts/gene_counts.txt"),
  skip = 1, select = c("Geneid", "Length"), col.names = c("gene", "length_bp"))
gl <- gene_len$length_bp[match(rownames(dge$counts), gene_len$gene)]
gl[is.na(gl)] <- median(gl, na.rm = TRUE)
rpk     <- sweep(dge$counts, 1, gl / 1e3, FUN = "/")
tpm_mat <- sweep(rpk, 2, colSums(rpk) / 1e6, FUN = "/")

sinfo <- data.table(sample_id = colnames(dge),
                    dataset   = dge$samples$dataset,
                    group     = as.character(dge$samples$group_binary))
valid_ds <- sinfo[, .(nc = sum(group == "Control"), nd = sum(group == "Disease")),
                  by = dataset][nc > 0 & nd > 0, dataset]

# Within-dataset control means
cmeans <- lapply(setNames(valid_ds, valid_ds), function(ds) {
  ids <- sinfo[dataset == ds & group == "Control", sample_id]
  if (length(ids) == 1) lcpm[, ids] else rowMeans(lcpm[, ids, drop = FALSE])
})

# Control LFC matrix (genes × ctrl samples)
ctrl_sids <- sinfo[group == "Control" & dataset %in% valid_ds, sample_id]
ctrl_lfc  <- matrix(NA_real_, nrow = nrow(lcpm), ncol = length(ctrl_sids),
                    dimnames = list(rownames(lcpm), ctrl_sids))
for (ds in valid_ds) {
  ids <- sinfo[dataset == ds & group == "Control", sample_id]
  for (sid in ids) ctrl_lfc[, sid] <- lcpm[, sid] - cmeans[[ds]]
}

# Disease LFC matrix
message("Loading patient LFC matrix...")
lfc_dis_dt <- fread(file.path(INT, "patient_lfc_matrix.csv.gz"))
setkey(lfc_dis_dt, gene)
dis_sids <- setdiff(colnames(lfc_dis_dt), "gene")
dis_lfc  <- as.matrix(lfc_dis_dt[gene %in% can_degs, dis_sids, with = FALSE])
rownames(dis_lfc) <- lfc_dis_dt[gene %in% can_degs, gene]

# Subset both to canonical DEGs present in ctrl_lfc
can_common <- intersect(can_degs, rownames(ctrl_lfc))
ctrl_lfc_can <- ctrl_lfc[can_common, , drop = FALSE]  # DEGs × ctrl

# ── Per-gene within-control SD ────────────────────────────────────────────────
ctrl_sd   <- apply(ctrl_lfc_can, 1, sd)
mean_expr <- rowMeans(lcpm[can_common, , drop = FALSE])
mean_tpm  <- rowMeans(tpm_mat[can_common, , drop = FALSE])

gene_stats <- data.table(
  gene          = can_common,
  ctrl_sd       = ctrl_sd,
  dream_abs_lfc = abs(dream_fc[can_common]),
  mean_logcpm   = mean_expr,
  mean_log2tpm1 = log2(mean_tpm + 1),
  direction = ifelse(can_common %in% can_up, "Up", "Down")
)

# ── Panel V1: Within-control SD per dataset (violin) ─────────────────────────
sd_per_ds <- rbindlist(lapply(valid_ds, function(ds) {
  ids  <- sinfo[dataset == ds & group == "Control", sample_id]
  mat  <- ctrl_lfc_can[, ids[ids %in% colnames(ctrl_lfc_can)], drop = FALSE]
  if (ncol(mat) < 2) return(NULL)
  data.table(dataset = ds, ctrl_sd = apply(mat, 1, sd),
             n_ctrl = ncol(mat))
}))
# Add an "All" aggregate
sd_all <- data.table(dataset = "All cohorts", ctrl_sd = ctrl_sd, n_ctrl = length(ctrl_sids))
sd_plot <- rbind(sd_per_ds, sd_all)
sd_plot[, ds_label := sprintf("%s\n(n=%d)", dataset, n_ctrl)]

# Compute per-group median for ordering
med_order <- sd_plot[, .(med = median(ctrl_sd)), by = ds_label][order(med)]
sd_plot[, ds_label := factor(ds_label, levels = med_order$ds_label)]

pV1 <- ggplot(sd_plot, aes(x = ds_label, y = ctrl_sd)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey60", linewidth = 0.35) +
  geom_violin(fill = col_ctrl, color = col_ctrl, alpha = 0.5, linewidth = 0.3,
              scale = "width") +
  geom_boxplot(width = 0.12, outlier.size = 0.4, outlier.alpha = 0.3,
               fill = "white", color = col_ctrl, linewidth = 0.4) +
  scale_y_continuous(breaks = seq(0, 4, 0.5)) +
  labs(x = NULL,
       y = "Within-control SD (log\u2082CPM)",
       title = "V1  Within-control variability of canonical DEGs") +
  theme(axis.text.x = element_text(size = 8))

ggsave(file.path(OUTDIR, "ctrl_variability_per_dataset.pdf"),
       pV1, width = 6.5, height = 4, device = cairo_pdf)
message("Saved: V1_ctrl_sd_per_dataset.pdf")

# ── Panel V2: |LFC| distribution — control vs disease ────────────────────────
# Sample for density (avoid overplotting with millions of points)
set.seed(42)
v2_ctrl <- data.table(
  abs_lfc = as.vector(abs(ctrl_lfc_can)),
  group   = "Control"
)
dis_can <- dis_lfc[rownames(dis_lfc) %in% can_common, , drop = FALSE]
v2_dis  <- data.table(
  abs_lfc = as.vector(abs(dis_can)),
  group   = "Disease"
)
v2 <- rbind(v2_ctrl[abs_lfc <= 4], v2_dis[abs_lfc <= 4])
v2[, group := factor(group, levels = c("Disease", "Control"))]

group_colors <- c("Disease" = col_dis, "Control" = col_ctrl)

pV2 <- ggplot(v2, aes(x = abs_lfc, color = group, fill = group)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey55", linewidth = 0.35) +
  geom_density(alpha = 0.25, linewidth = 0.55, adjust = 0.7) +
  scale_color_manual(values = group_colors, name = NULL) +
  scale_fill_manual(values  = group_colors, name = NULL) +
  scale_x_continuous(breaks = seq(0, 4, 0.5),
                     labels = function(x) sprintf("%.1f", x)) +
  annotate("text", x = 0.52, y = Inf, label = "|LFC| = 0.5",
           hjust = 0, vjust = 1.5, size = 2.8, color = "grey40") +
  labs(x = "|Patient LFC| (log\u2082CPM, relative to control mean)",
       y = "Density",
       title = "V2  |LFC| distribution: disease patients vs controls") +
  theme(legend.position = c(0.82, 0.82))

ggsave(file.path(OUTDIR, "abs_lfc_disease_vs_control.pdf"),
       pV2, width = 5.5, height = 4, device = cairo_pdf)
message("Saved: V2_abs_lfc_distribution.pdf")

# ── Panel V3: Dream |logFC| vs within-control SD (per gene) ──────────────────
pV3 <- ggplot(gene_stats, aes(x = dream_abs_lfc, y = ctrl_sd, color = direction)) +
  geom_point(size = 0.6, alpha = 0.4) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.7,
              color = "black", fill = "grey80") +
  scale_color_manual(values = c("Up" = "#C9265E", "Down" = "#1565C0"), name = NULL) +
  scale_x_continuous(breaks = seq(0, 5, 0.5)) +
  scale_y_continuous(breaks = seq(0, 4, 0.5)) +
  labs(x = "Dream |log\u2082FC| (population effect size)",
       y = "Within-control SD (log\u2082CPM)",
       title = "V3  Larger population effect = more variable in controls") +
  theme(legend.position = c(0.88, 0.15))

ggsave(file.path(OUTDIR, "dream_lfc_vs_ctrl_variability.pdf"),
       pV3, width = 5.5, height = 4, device = cairo_pdf)
message("Saved: V3_dream_lfc_vs_ctrl_sd.pdf")

# ── Panel V4: Mean TPM vs within-control SD ──────────────────────────────────
tpm1_vline <- log2(1 + 1)   # log2(TPM=1 + 1)

pV4 <- ggplot(gene_stats, aes(x = mean_log2tpm1, y = ctrl_sd, color = direction)) +
  geom_point(size = 0.6, alpha = 0.4) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.7,
              color = "black", fill = "grey80") +
  geom_vline(xintercept = tpm1_vline, linetype = "dashed",
             color = "grey55", linewidth = 0.35) +
  scale_color_manual(values = c("Up" = "#C9265E", "Down" = "#1565C0"), name = NULL) +
  scale_x_continuous(breaks = seq(0, 14, 2)) +
  scale_y_continuous(breaks = seq(0, 4, 0.5)) +
  annotate("text", x = tpm1_vline + 0.1, y = Inf, label = "TPM = 1",
           hjust = 0, vjust = 1.5, size = 2.8, color = "grey40") +
  labs(x = "Mean log\u2082(TPM + 1) across all samples",
       y = "Within-control SD (log\u2082CPM)",
       title = "Expression level vs within-control variability") +
  theme(legend.position = c(0.88, 0.85))

ggsave(file.path(OUTDIR, "expression_vs_ctrl_variability.pdf"),
       pV4, width = 5.5, height = 4, device = cairo_pdf)
message("Saved: V4_expr_vs_ctrl_sd.pdf")

# ── Quick summary numbers ─────────────────────────────────────────────────────
message("\n===== CONTROL VARIABILITY SUMMARY =====")
message(sprintf("Canonical DEGs:  %d", length(can_common)))
message(sprintf("Within-ctrl SD:  median=%.3f  mean=%.3f  pct>0.5: %.1f%%  pct>1.0: %.1f%%",
                median(ctrl_sd), mean(ctrl_sd),
                mean(ctrl_sd > 0.5) * 100, mean(ctrl_sd > 1.0) * 100))
message(sprintf("Corr(dream|LFC|, ctrl_SD):  r=%.3f",
                cor(gene_stats$dream_abs_lfc, gene_stats$ctrl_sd, method = "pearson")))
message(sprintf("Corr(mean_logCPM, ctrl_SD): r=%.3f",
                cor(gene_stats$mean_logcpm,   gene_stats$ctrl_sd, method = "pearson")))
