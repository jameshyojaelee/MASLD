#!/usr/bin/env Rscript
# figS_multimethod_batch_pca.R
# Panel J (batch-model PCA) for the multi-method comparison: PCA of the pooled
# 5-cohort mega samples under each method's ACTUAL batch model, so the reader sees
# how each DE method treats the dataset (cohort) axis. (panelA-I now belong to the
# degx battery; this bespoke batch-handling diagnostic is panelJ.)
#   panelJ_batch_model_pca.pdf
#   Row 1 (coloured by COHORT)  : Raw | Fixed-effect (DESeq2 ~dataset) | Random-effect (dream (1|dataset))
#   Row 2 (coloured by DISEASE) : same three corrections -> biology preserved
#   Row 3 (metafor)             : 5 per-cohort PCAs -> metafor never pools (no joint matrix)
# Fixed = limma::removeBatchEffect (OLS per-dataset offsets). Random = per-gene
# lme4 fit y ~ group + (1|dataset), subtract the SHRUNKEN dataset BLUPs (dream's
# partial pooling). metafor analyses each cohort separately -> shown per-cohort.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(edgeR); library(limma); library(matrixStats); library(lme4)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
OUT <- file.path(BASE, "figures/supplementary/figS_methods_validation/multimethod_validation/panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CTRL <- "#9E9E9E"

MEGA <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
cohort_short <- c(GSE126848 = "Suppli", GSE130970 = "Hoang", GSE135251 = "Govaere",
                  GSE162694 = "Bril", GSE213621 = "Chen")
cohort_pal <- c(Suppli = "#1F77B4", Hoang = "#FF7F0E", Govaere = "#2CA02C",
                Bril = "#D62728", Chen = "#9467BD")

# --- pooled counts: 5 mega cohorts only -------------------------------------
dge <- load_merged_dge(); stopifnot(!is.null(dge))
samp <- as.data.table(dge$samples, keep.rownames = "sample_id")
keep_s <- samp$dataset %in% MEGA
dge <- dge[, keep_s]; samp <- samp[keep_s]
cat(sprintf("Pooled mega samples: %d  (cohorts: %s)\n",
            ncol(dge), paste(sort(unique(samp$dataset)), collapse = ", ")))

logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
rv <- matrixStats::rowVars(logcpm)
top <- order(rv, decreasing = TRUE)[seq_len(min(2000, length(rv)))]
X   <- logcpm[top, ]                                  # 2000 x n, shared across corrections

grp <- factor(samp$group_binary, levels = c("Control", "Disease"))
ds  <- factor(samp$dataset)

# --- (1) Fixed-effect correction = DESeq2's ~dataset (OLS offsets) -----------
X_fixed <- limma::removeBatchEffect(X, batch = ds, design = model.matrix(~ grp))

# --- (2) Random-effect correction = dream's (1|dataset) (shrunken BLUPs) -----
# Per gene: y ~ group + (1|dataset); subtract the (shrunk) per-dataset random
# intercept from each sample. Falls back to the fixed offset if lmer fails.
cat("Fitting per-gene mixed models for the random-effect (dream) correction...\n")
ctrl_lmer <- lme4::lmerControl(optimizer = "nloptwrap",
                               check.conv.singular = "ignore")
ds_chr <- as.character(ds)
X_rand <- X
for (i in seq_len(nrow(X))) {
  y <- X[i, ]
  b <- tryCatch({
    fit <- suppressWarnings(suppressMessages(
      lme4::lmer(y ~ grp + (1 | ds), control = ctrl_lmer)))
    re <- lme4::ranef(fit)$ds[, 1]; names(re) <- rownames(lme4::ranef(fit)$ds)
    re[ds_chr]
  }, error = function(e) {                       # fallback: OLS dataset means (centered)
    m <- tapply(y, ds_chr, mean); (m - mean(m))[ds_chr]
  })
  X_rand[i, ] <- y - as.numeric(b)
  if (i %% 400 == 0) cat("  ", i, "/", nrow(X), "\n")
}

# --- PCA helper -------------------------------------------------------------
do_pca <- function(mat, label) {
  pr <- prcomp(t(mat), center = TRUE, scale. = FALSE)
  pve <- 100 * pr$sdev^2 / sum(pr$sdev^2)
  data.table(sample_id = colnames(mat), PC1 = pr$x[, 1], PC2 = pr$x[, 2],
             model = label, pc1 = pve[1], pc2 = pve[2])
}
MODELS <- c("Raw (no correction)",
            "Fixed-effect  ~dataset\n(DESeq2 batch model)",
            "Random-effect  (1|dataset)\n(dream batch model)")
pca <- rbindlist(list(do_pca(X, MODELS[1]),
                      do_pca(X_fixed, MODELS[2]),
                      do_pca(X_rand, MODELS[3])))
pca <- merge(pca, samp[, .(sample_id, dataset, group_binary)], by = "sample_id")
pca[, model   := factor(model, levels = MODELS)]
pca[, cohort  := factor(cohort_short[dataset], levels = unname(cohort_short))]
pca[, disease := factor(group_binary, levels = c("Control", "Disease"))]
ann <- pca[, .(lab = sprintf("PC1=%.0f%%  PC2=%.0f%%", unique(pc1), unique(pc2))), by = model]

base_pca <- function() theme_masld(base_size = 7) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_text(size = 6), panel.grid = element_blank(),
        strip.text = element_text(size = 6.6, face = "bold", lineheight = 0.9),
        strip.background = element_blank(),
        legend.position = "right", legend.title = element_text(size = 6.5, face = "bold"),
        legend.text = element_text(size = 6), legend.key.size = unit(0.25, "cm"),
        plot.title = element_text(size = 8, face = "bold"))

row_cohort <- ggplot(pca, aes(PC1, PC2, colour = cohort)) +
  geom_point(size = 0.4, alpha = 0.7) +
  facet_wrap(~ model, nrow = 1, scales = "free") +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
            hjust = -0.06, vjust = 1.4, size = 2.1, colour = "grey30") +
  scale_colour_manual(values = cohort_pal, name = "Cohort") +
  guides(colour = guide_legend(override.aes = list(size = 1.6, alpha = 1))) +
  labs(x = "PC1", y = "PC2",
       title = "Dataset (batch) axis: raw -> fixed-effect -> random-effect correction") +
  base_pca()

row_dis <- ggplot(pca, aes(PC1, PC2, colour = disease)) +
  geom_point(size = 0.4, alpha = 0.7) +
  facet_wrap(~ model, nrow = 1, scales = "free") +
  scale_colour_manual(values = c(Control = CTRL, Disease = masld_colors$nash), name = "Disease") +
  guides(colour = guide_legend(override.aes = list(size = 1.6, alpha = 1))) +
  labs(x = "PC1", y = "PC2",
       title = "Same projections coloured by disease (a minor axis ~0.9% of variance) - not expected to separate in PCA") +
  base_pca()

# --- (3) metafor: per-cohort PCAs (no pooled matrix) ------------------------
percoh <- rbindlist(lapply(MEGA, function(d) {
  sub <- which(samp$dataset == d)
  lc  <- edgeR::cpm(dge[, sub], log = TRUE, prior.count = 1)
  rvc <- matrixStats::rowVars(lc)
  tc  <- order(rvc, decreasing = TRUE)[seq_len(min(2000, length(rvc)))]
  pr  <- prcomp(t(lc[tc, ]), center = TRUE, scale. = FALSE)
  pve <- 100 * pr$sdev^2 / sum(pr$sdev^2)
  data.table(PC1 = pr$x[, 1], PC2 = pr$x[, 2],
             cohort = cohort_short[d],
             disease = factor(samp$group_binary[sub], levels = c("Control", "Disease")),
             lab = sprintf("PC1=%.0f%%", pve[1]))
}))
percoh[, cohort := factor(cohort, levels = unname(cohort_short))]
pc_ann <- percoh[, .(lab = unique(lab)), by = cohort]
row_meta <- ggplot(percoh, aes(PC1, PC2, colour = disease)) +
  geom_point(size = 0.4, alpha = 0.7) +
  facet_wrap(~ cohort, nrow = 1, scales = "free") +
  geom_text(data = pc_ann, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
            hjust = -0.08, vjust = 1.4, size = 2.0, colour = "grey30") +
  scale_colour_manual(values = c(Control = CTRL, Disease = masld_colors$nash), name = "Disease") +
  guides(colour = guide_legend(override.aes = list(size = 1.6, alpha = 1))) +
  labs(x = "PC1", y = "PC2",
       title = "metafor: RAW within-cohort PCA (no cross-cohort batch to remove) - it never pools the samples, so there is no joint corrected matrix") +
  base_pca()

fig <- (row_cohort / row_dis / row_meta) +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(
    title = "Batch handling of the pooled mega samples, per DE method",
    subtitle = sprintf("Top-2,000 most-variable log2-CPM; %d samples, 5 control-bearing mega cohorts. Fixed-effect = limma::removeBatchEffect (DESeq2 ~dataset); random-effect = per-gene lme4 (1|dataset) BLUP removal (dream). Disease is ~0.9%% of expression variance (vs dataset ~23%%) - never a dominant PC, so it does not separate in PCA by ANY method (PCA = batch diagnostic, not a disease classifier; for supervised disease separation see panelK).",
                       ncol(dge)),
    theme = theme(plot.title = element_text(size = 9.5, face = "bold"),
                  plot.subtitle = element_text(size = 6.2, colour = "grey35")))

ggsave(file.path(OUT, "panelJ_batch_model_pca.pdf"), fig,
       width = 9.0, height = 8.2, device = cairo_pdf)
fwrite(pca[, .(sample_id, dataset, cohort, disease, model, PC1, PC2)],
       file.path(OUT, "panelJ_batch_model_pca_data.csv"))
cat("Wrote panelJ_batch_model_pca.pdf\n")
print(ann)
