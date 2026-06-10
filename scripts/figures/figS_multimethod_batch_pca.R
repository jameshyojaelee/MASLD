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
  library(msigdbr)
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
# HVG selection on within-cohort centered matrix so cohort mean-shifts don't
# inflate variance and bias gene selection toward cohort-discriminating genes.
logcpm_wc <- logcpm
for (d in unique(samp$dataset)) {
  idx <- which(samp$dataset == d)
  logcpm_wc[, idx] <- logcpm[, idx] - rowMeans(logcpm[, idx, drop = FALSE])
}
rv <- matrixStats::rowVars(logcpm_wc)
top <- order(rv, decreasing = TRUE)[seq_len(min(2000, length(rv)))]
X   <- logcpm[top, ]                                  # 2000 x n, shared across corrections
cat(sprintf("HVG selection: top 2000 by within-cohort variance (was raw variance)\n"))

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

# =============================================================================
# Loadings analysis: what genes drive PC1/PC2 in the corrected space?
# Use fixed-effect corrected matrix (rand is r=0.9999 identical, no need to redo)
# =============================================================================
pr_load <- prcomp(t(X_fixed), center = TRUE, scale. = FALSE)
pve_load <- 100 * pr_load$sdev^2 / sum(pr_load$sdev^2)
loads <- as.data.table(pr_load$rotation[, 1:4], keep.rownames = "gene")
cat(sprintf("PC1=%.1f%%  PC2=%.1f%%  PC3=%.1f%%  PC4=%.1f%%\n",
            pve_load[1], pve_load[2], pve_load[3], pve_load[4]))

# map versioned Ensembl IDs -> gene symbols for enrichment + barplot labels
meta <- fread(file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz"),
              select = c("gene_id", "gene_name"))
meta[, ensembl_base := sub("\\..*", "", gene_id)]
loads[, ensembl_base := sub("\\..*", "", gene)]
loads <- merge(loads, meta[, .(ensembl_base, gene_name)], by = "ensembl_base", all.x = TRUE)
loads[is.na(gene_name), gene_name := gene]   # fallback to Ensembl if no match
cat(sprintf("Symbol mapping: %d / %d genes matched\n",
            sum(!is.na(loads$gene_name) & loads$gene_name != loads$gene), nrow(loads)))

# top 30 genes by |loading| for PC1 and PC2
top_pc1 <- loads[order(-abs(PC1))][1:30, .(gene = gene_name, loading = PC1, pc = "PC1 (batch-corrected)")]
top_pc2 <- loads[order(-abs(PC2))][1:30, .(gene = gene_name, loading = PC2, pc = "PC2 (batch-corrected)")]
top_loads <- rbind(top_pc1, top_pc2)
top_loads[, gene_fac := factor(gene, levels = rev(unique(gene[order(pc, loading)])))]

pLoad <- ggplot(top_loads, aes(loading, gene_fac, fill = loading > 0)) +
  geom_col(width = 0.72, linewidth = 0) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  facet_wrap(~ pc, scales = "free") +
  scale_fill_manual(values = c("TRUE" = masld_colors$nash, "FALSE" = CTRL), guide = "none") +
  labs(x = "PC loading", y = NULL,
       title = sprintf("Top 30 genes driving PC1 (%.1f%% var) and PC2 (%.1f%% var) after batch correction",
                       pve_load[1], pve_load[2]),
       subtitle = "Positive loading = pushes samples rightward/upward on that PC; negative = opposite") +
  theme_masld(base_size = 7) +
  theme(axis.text.y = element_text(size = 5.5), strip.text = element_text(size = 7, face = "bold"))

# --- pathway enrichment: hypergeometric against MSigDB Hallmark --------------
hall <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))[, .(gs_name, gene_symbol)]
universe_sym <- loads$gene_name
n_univ       <- length(universe_sym)

enrich_res <- rbindlist(lapply(c("PC1", "PC2"), function(pc_col) {
  lo <- loads[[pc_col]]; names(lo) <- loads$gene_name
  dirs <- list(positive = names(sort(lo, decreasing = TRUE))[1:200],
               negative = names(sort(lo, decreasing = FALSE))[1:200])
  rbindlist(lapply(names(dirs), function(dir_name) {
    qgenes <- dirs[[dir_name]]
    rbindlist(lapply(split(hall, hall$gs_name), function(gs) {
      q <- sum(qgenes %in% gs$gene_symbol)
      if (q == 0) return(NULL)
      K <- sum(universe_sym %in% gs$gene_symbol)
      data.table(pc = pc_col, direction = dir_name,
                 pathway = sub("HALLMARK_", "", unique(gs$gs_name)),
                 n_overlap = q, n_pathway = K,
                 p_hyper = phyper(q - 1, K, n_univ - K, length(qgenes), lower.tail = FALSE))
    }))
  }))
}))
if (nrow(enrich_res) == 0) {
  cat("No Hallmark overlaps found — check gene symbol mapping\n")
  sig_enrich <- data.table()
} else {
  enrich_res[, padj := p.adjust(p_hyper, method = "BH"), by = .(pc, direction)]
  sig_enrich <- enrich_res[padj < 0.25][order(pc, direction, padj)]
}
cat("\nSignificant Hallmark pathways (FDR<0.25) driving top PC loadings:\n")
print(sig_enrich[, .(pc, direction, pathway, n_overlap, n_pathway, padj)])

# dotplot of top enrichments (up to 8 per PC × direction)
plot_enrich <- enrich_res[padj < 0.25][order(padj)][
  , head(.SD, 8), by = .(pc, direction)]
if (nrow(plot_enrich) > 0) {
  plot_enrich[, label := sprintf("%s (%s)", pathway, direction)]
  plot_enrich[, neg_log10p := -log10(padj + 1e-10)]
  pEnrich <- ggplot(plot_enrich, aes(neg_log10p, reorder(label, neg_log10p),
                                     colour = direction, size = n_overlap)) +
    geom_point() +
    facet_wrap(~ pc, scales = "free_y", ncol = 1) +
    scale_colour_manual(values = c(positive = masld_colors$nash, negative = CTRL)) +
    scale_size_continuous(range = c(1.5, 4), name = "n genes") +
    labs(x = "-log10(FDR)", y = NULL, colour = "Loading direction",
         title = "Hallmark pathways enriched in top-200 PC loadings (FDR < 0.25)") +
    theme_masld(base_size = 7) +
    theme(axis.text.y = element_text(size = 5.5), strip.text = element_text(size = 7, face = "bold"))
} else {
  pEnrich <- ggplot() + annotate("text", 0, 0, label = "No Hallmark pathways at FDR<0.25") +
    theme_void()
  cat("No significant Hallmark enrichments — top PC loadings are not pathway-structured\n")
}

fig_loads <- (pLoad / pEnrich) + plot_layout(heights = c(1.2, 1)) +
  plot_annotation(title = "PC loadings: what genes drive the top PCs after batch correction",
    theme = theme(plot.title = element_text(size = 9, face = "bold")))
ggsave(file.path(OUT, "panelJ_loadings.pdf"), fig_loads,
       width = 9.0, height = 10.0, device = cairo_pdf)
fwrite(loads[, .(gene, gene_name, PC1, PC2, PC3, PC4)],
       file.path(OUT, "panelJ_loadings_data.csv"))
fwrite(sig_enrich, file.path(OUT, "panelJ_loadings_enrichment.csv"))
cat("Wrote panelJ_loadings.pdf\n")

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
