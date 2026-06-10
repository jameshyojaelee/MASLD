#!/usr/bin/env Rscript
# figS_deg_expression.R
# Expression distribution of canonical MASLD DEGs in disease patients.
# Shows both logCPM (library-size normalised) and TPM (length + library-size
# normalised) so the TPM ≥ 1 filter threshold is visible on the right scale.
#
# Outputs (figures/supplementary/figS_methods_validation/lfc_sensitivity/):
#   deg_expression_distribution.pdf  — logCPM | log2(TPM+1) side-by-side
#   tpm_filtered_degs.csv            — canonical DEGs passing median TPM ≥ 1
#                                      across disease patients (used by
#                                      figS_patient_concordance.R)

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(scales)
  library(edgeR)
  library(yaml)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

# 2026-05-28: repointed to canonical STAR -s 2 (was figures_STAR_s0 + _star inputs)
OUTDIR <- file.path(BASE, "figures/supplementary/figS_methods_validation/lfc_sensitivity")
INT    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")

theme_pub <- theme_minimal(base_size = 11) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line  = element_line(colour = "black", linewidth = 0.3),
        axis.ticks = element_line(colour = "black", linewidth = 0.3),
        legend.background = element_blank(), legend.key = element_blank(),
        strip.background  = element_blank(),
        strip.text  = element_text(face = "bold", size = 10),
        plot.title  = element_text(face = "bold", size = 12),
        axis.title  = element_text(size = 10),
        axis.text   = element_text(size = 9),
        legend.text = element_text(size = 9),
        plot.margin = margin(8, 10, 8, 8))
theme_set(theme_pub)

col_up   <- "#C9265E"
col_down <- "#1565C0"
dir_colors <- c("Up" = col_up, "Down" = col_down)

# ── Canonical DEGs ────────────────────────────────────────────────────────────
message("Loading dream STAR results...")
dream <- fread(file.path(INT, "dream_results_ashr.csv"))
can_up <- dream[padj < 0.05 & logFC >  0.5, gene]
can_dn <- dream[padj < 0.05 & logFC < -0.5, gene]
can_degs <- c(can_up, can_dn)
message(sprintf("  Canonical DEGs: %d up + %d down = %d", length(can_up), length(can_dn), length(can_degs)))

# ── Gene lengths from featureCounts ──────────────────────────────────────────
message("Reading gene lengths from featureCounts output...")
fc_file <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/results/GSE126848/counts/featurecounts/gene_counts.txt")
fc_hdr <- fread(fc_file, skip = 1, nrows = 0)   # get column names
gene_len <- fread(fc_file, skip = 1,
                  select = c("Geneid", "Length"),
                  col.names = c("gene", "length_bp"))
message(sprintf("  Gene lengths: %d genes", nrow(gene_len)))

# ── Raw counts + mega-cohort sample selection ─────────────────────────────────
message("Loading STAR DGE...")
dge <- readRDS(file.path(INT, "merged_dge.rds"))
ycfg   <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
mega_c <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
dge    <- dge[, dge$samples$dataset %in% mega_c]

sinfo <- data.table(sample_id = colnames(dge),
                    dataset   = dge$samples$dataset,
                    group     = as.character(dge$samples$group_binary))
dis_sids <- sinfo[group == "Disease", sample_id]

# ── Compute logCPM ────────────────────────────────────────────────────────────
message("Computing logCPM...")
lcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

# ── Compute TPM ───────────────────────────────────────────────────────────────
message("Computing TPM (counts / length → per-million normalisation)...")
raw_counts <- dge$counts   # genes × samples (integer counts)

# Align gene lengths to DGE row order
genes_in_dge <- rownames(raw_counts)
gl <- gene_len[match(genes_in_dge, gene), length_bp]
# Any missing lengths → impute with median (shouldn't happen with GENCODE v49)
gl[is.na(gl)] <- median(gl, na.rm = TRUE)

# TPM = (count / length_kb) / sum(count / length_kb) × 1e6
rpk     <- sweep(raw_counts, 1, gl / 1e3, FUN = "/")   # reads per kilobase
tpm_mat <- sweep(rpk, 2, colSums(rpk) / 1e6, FUN = "/")
message(sprintf("  TPM matrix: %d genes × %d samples", nrow(tpm_mat), ncol(tpm_mat)))

# ── Disease patients — canonical DEGs ─────────────────────────────────────────
can_in_dge <- intersect(can_degs, rownames(lcpm))

# logCPM for disease patients, canonical DEGs
lcpm_dis <- lcpm[can_in_dge, dis_sids, drop = FALSE]

# TPM for disease patients, canonical DEGs
tpm_dis  <- tpm_mat[can_in_dge, dis_sids, drop = FALSE]

# ── TPM ≥ 1 filter: median TPM across disease patients ────────────────────────
median_tpm <- apply(tpm_dis, 1, median)
names(median_tpm) <- can_in_dge

tpm_pass <- names(median_tpm)[median_tpm >= 1]
message(sprintf("  DEGs with median TPM ≥ 1: %d / %d (%.1f%%)",
                length(tpm_pass), length(can_in_dge),
                length(tpm_pass) / length(can_in_dge) * 100))

# Save filtered gene list for use in concordance analysis
tpm_filter_tbl <- data.table(
  gene      = can_in_dge,
  direction = ifelse(can_in_dge %in% can_up, "Up", "Down"),
  dream_logFC = dream$logFC[match(can_in_dge, dream$gene)],
  median_tpm  = round(median_tpm, 3),
  tpm_pass    = can_in_dge %in% tpm_pass
)
fwrite(tpm_filter_tbl, file.path(OUTDIR, "tpm_filtered_degs.csv"))
message("Saved: tpm_filtered_degs.csv")

# ── Build long tables for plotting ────────────────────────────────────────────
make_long <- function(mat, value_col) {
  dt <- as.data.table(mat, keep.rownames = "gene")
  dt_long <- melt(dt, id.vars = "gene", variable.name = "sample_id",
                  value.name = value_col)
  dt_long[, direction := ifelse(gene %in% can_up, "Up", "Down")]
  dt_long
}

lcpm_long <- make_long(lcpm_dis, "logcpm")
tpm_long  <- make_long(log2(tpm_dis + 1), "log2tpm1")   # log2(TPM+1)

# ── Panel: logCPM distribution ────────────────────────────────────────────────
p_lcpm <- ggplot(lcpm_long, aes(x = logcpm, fill = direction, color = direction)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60,
                 alpha = 0.55, position = "identity", linewidth = 0.2) +
  scale_fill_manual(values  = dir_colors, name = NULL) +
  scale_color_manual(values = dir_colors, name = NULL) +
  scale_x_continuous(breaks = seq(-4, 14, 2)) +
  labs(x = "log\u2082CPM (per disease patient, per gene)",
       y = "Density",
       title = "log\u2082CPM") +
  theme(legend.position = "top",
        legend.key.size = unit(0.3, "cm"))

# ── Panel: log2(TPM+1) distribution with TPM=1 vline ─────────────────────────
p_tpm <- ggplot(tpm_long, aes(x = log2tpm1, fill = direction, color = direction)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60,
                 alpha = 0.55, position = "identity", linewidth = 0.2) +
  geom_vline(xintercept = log2(1 + 1), linetype = "dashed",
             color = "grey40", linewidth = 0.45) +
  annotate("text", x = log2(2) + 0.1, y = Inf,
           label = "TPM = 1", hjust = 0, vjust = 1.5,
           size = 2.8, color = "grey40") +
  scale_fill_manual(values  = dir_colors, name = NULL) +
  scale_color_manual(values = dir_colors, name = NULL) +
  scale_x_continuous(breaks = seq(0, 16, 2)) +
  labs(x = "log\u2082(TPM + 1) (per disease patient, per gene)",
       y = "Density",
       title = "log\u2082(TPM + 1)") +
  theme(legend.position = "none")

p_combined <- p_lcpm + p_tpm +
  plot_annotation(
    title    = "Expression distribution of canonical MASLD DEGs — disease patients",
    subtitle = sprintf("%d DEGs (padj < 0.05, |log\u2082FC| > 0.5)  \u00B7  %d disease patients",
                       length(can_in_dge), length(dis_sids)),
    theme = theme(plot.title    = element_text(face = "bold", size = 11),
                  plot.subtitle = element_text(size = 9, color = "grey40"))
  )

ggsave(file.path(OUTDIR, "deg_expression_distribution.pdf"),
       p_combined, width = 9, height = 4, device = cairo_pdf)
message("Saved: deg_expression_distribution.pdf")

message(sprintf("\nTPM filter summary: %d / %d DEGs pass median TPM >= 1",
                length(tpm_pass), length(can_in_dge)))
message(sprintf("  Up passing:   %d / %d (%.1f%%)",
                sum(tpm_pass %in% can_up), length(can_up),
                sum(tpm_pass %in% can_up) / length(can_up) * 100))
message(sprintf("  Down passing: %d / %d (%.1f%%)",
                sum(tpm_pass %in% can_dn), length(can_dn),
                sum(tpm_pass %in% can_dn) / length(can_dn) * 100))
