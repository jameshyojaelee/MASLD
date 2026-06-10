#!/usr/bin/env Rscript
# figS_covariate_pc_association.R
# Covariate-PC association heatmap (Govaere 2020 style).
# Shows all 4 correction levels: Raw | Batch | Batch+Sex | Batch+Sex+SVA
# Categorical -> Kruskal-Wallis; Continuous/ordinal -> Spearman rho
# Output: panelJ_covariate_pc_association.pdf
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(edgeR); library(limma); library(matrixStats); library(sva)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
OUT <- file.path(BASE,
  "figures/supplementary/figS_methods_validation/multimethod_validation/panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

MEGA  <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
N_PCS <- 7L; N_HVG <- 2000L

# ---------------------------------------------------------------------------
# 1. Load data + infer sex (XIST/DDX3Y k-means for missing annotation)
# ---------------------------------------------------------------------------
dge  <- load_merged_dge(); stopifnot(!is.null(dge))
samp <- as.data.table(dge$samples, keep.rownames = "sample_id")
keep <- samp$dataset %in% MEGA
dge  <- dge[, keep]; samp <- samp[keep]
cat(sprintf("Samples: %d\n", ncol(dge)))

grp <- factor(samp$group_binary, levels = c("Control", "Disease"))
ds  <- factor(samp$dataset)
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

# sex inference (fills NA for Govaere + Chen)
xist <- rownames(dge)[grep("^ENSG00000229807", rownames(dge))]
ddy  <- rownames(dge)[grep("^ENSG00000067048", rownames(dge))]
sex_inferred <- samp$sex
if (length(xist) > 0 && length(ddy) > 0) {
  expr_sex <- t(logcpm[c(xist[1], ddy[1]), , drop = FALSE])
  km <- kmeans(expr_sex, centers = 2, nstart = 20, iter.max = 50)
  female_cluster <- as.integer(names(which.max(
    tapply(expr_sex[, 1], km$cluster, mean))))
  inferred <- ifelse(km$cluster == female_cluster, "F", "M")
  missing <- is.na(samp$sex) | samp$sex == ""
  sex_inferred[missing] <- inferred[missing]
  cat(sprintf("Sex: %d annotated + %d inferred. Final F=%d M=%d\n",
              sum(!missing), sum(missing),
              sum(sex_inferred == "F", na.rm = TRUE),
              sum(sex_inferred == "M", na.rm = TRUE)))
}
sex_inferred[is.na(sex_inferred) | sex_inferred == ""] <- "Unknown"

# ---------------------------------------------------------------------------
# 2. HVG selection (within-cohort centred) + 4 correction levels
# ---------------------------------------------------------------------------
logcpm_wc <- logcpm
for (d in unique(samp$dataset)) {
  idx <- which(samp$dataset == d)
  logcpm_wc[, idx] <- logcpm[, idx] - rowMeans(logcpm[, idx, drop = FALSE])
}
rv  <- matrixStats::rowVars(logcpm_wc)
top <- order(rv, decreasing = TRUE)[seq_len(N_HVG)]
X   <- logcpm[top, ]

des  <- model.matrix(~ grp)
sex_cov <- as.numeric(sex_inferred == "F")
X_batch     <- limma::removeBatchEffect(X, batch = ds, design = des)
X_batch_sex <- limma::removeBatchEffect(X, batch = ds,
                                         covariates = sex_cov, design = des)
cat("Running SVA...\n")
mod0    <- model.matrix(~ ds + sex_cov)
mod1    <- model.matrix(~ grp + ds + sex_cov)
n_sv      <- sva::num.sv(X, mod1, method = "leek")
sva_label <- sprintf("Batch + sex + SVA (%d SV)", n_sv)
sva_fit   <- sva::sva(X, mod1, mod0, n.sv = n_sv)
X_sva   <- limma::removeBatchEffect(X, batch = ds,
                                     covariates = cbind(sex_cov, sva_fit$sv),
                                     design = des)
cat(sprintf("SVA: %d surrogate variables\n", n_sv))

# ---------------------------------------------------------------------------
# 3. PCA (7 PCs each)
# ---------------------------------------------------------------------------
run_pca <- function(mat) {
  pr  <- prcomp(t(mat), center = TRUE, scale. = FALSE)
  pve <- round(100 * pr$sdev^2 / sum(pr$sdev^2), 1)
  list(scores = pr$x[, seq_len(N_PCS)], pve = pve[seq_len(N_PCS)])
}
pca_raw <- run_pca(X);           cat("Raw          PVE:", pca_raw$pve, "\n")
pca_bat <- run_pca(X_batch);     cat("Batch        PVE:", pca_bat$pve, "\n")
pca_sex <- run_pca(X_batch_sex); cat("Batch+sex    PVE:", pca_sex$pve, "\n")
pca_sva <- run_pca(X_sva);       cat("Batch+sex+SVA PVE:", pca_sva$pve,"\n")

# merge clinical metadata AFTER matrix ops (preserves column order)
meta <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"),
  select = c("sample_id", "fibrosis_stage", "nas_score", "diagnosis_harmonized"))
setkey(samp, sample_id); setkey(meta, sample_id)
samp <- meta[samp]   # right join preserving samp rows
samp <- samp[colnames(logcpm)]  # restore original column order

# ---------------------------------------------------------------------------
# 4. Covariate definitions
# ---------------------------------------------------------------------------
covars <- list(
  "Disease Status" = list(vals = samp$group_binary,         type = "cat"),
  "Cohort (Batch)" = list(vals = samp$dataset,              type = "cat"),
  "Sex"            = list(vals = sex_inferred,              type = "cat"),
  "Diagnosis"      = list(vals = samp$diagnosis_harmonized, type = "cat"),
  "Fibrosis Stage" = list(vals = samp$fibrosis_stage,       type = "cont"),
  "NAS Score"      = list(vals = samp$nas_score,            type = "cont"),
  "Library Size"   = list(vals = log10(samp$lib.size + 1),  type = "cont")
)

# ---------------------------------------------------------------------------
# 5. Association tests
# ---------------------------------------------------------------------------
assoc_pval <- function(pc_scores, covars) {
  pc_names  <- colnames(pc_scores)
  cov_names <- names(covars)
  res <- data.table(
    PC        = rep(pc_names, each = length(cov_names)),
    covariate = rep(cov_names, times = length(pc_names)),
    neglog10p = NA_real_)
  for (pc in pc_names) {
    y <- pc_scores[, pc]
    for (cv in cov_names) {
      x <- covars[[cv]]$vals; type <- covars[[cv]]$type
      ok <- !is.na(x) & !is.na(y) & (type == "cont" | x != "")
      if (sum(ok) < 10) next
      p <- tryCatch(
        if (type == "cat") kruskal.test(y[ok] ~ factor(x[ok]))$p.value
        else cor.test(y[ok], as.numeric(x[ok]),
                      method = "spearman", exact = FALSE)$p.value,
        error = function(e) NA_real_)
      res[PC == pc & covariate == cv, neglog10p := -log10(p + 1e-300)]
    }
  }
  res
}

res_raw <- assoc_pval(pca_raw$scores, covars); res_raw[, label := "Raw"]
res_bat <- assoc_pval(pca_bat$scores, covars); res_bat[, label := "Batch"]
res_sex <- assoc_pval(pca_sex$scores, covars); res_sex[, label := "Batch + sex"]
res_sva <- assoc_pval(pca_sva$scores, covars)
res_sva[, label := sva_label]
res_all <- rbind(res_raw, res_bat, res_sex, res_sva)

lbl_order <- c("Raw", "Batch", "Batch + sex", sva_label)
pc_order  <- paste0("PC", seq_len(N_PCS))
cov_order <- names(covars)
res_all[, PC        := factor(PC, levels = pc_order)]
res_all[, covariate := factor(covariate, levels = cov_order)]
res_all[, label     := factor(label, levels = lbl_order)]

# PVE label per correction level (for y-axis)
pve_labs <- function(pve)
  setNames(sprintf("PC%d\n(%.1f%%)", seq_len(N_PCS), pve), pc_order)
pve_list  <- setNames(
  list(pve_labs(pca_raw$pve), pve_labs(pca_bat$pve),
       pve_labs(pca_sex$pve), pve_labs(pca_sva$pve)),
  c("Raw", "Batch", "Batch + sex", sva_label))

# ---------------------------------------------------------------------------
# 6. Plot: 2x2 grid of heatmaps
# ---------------------------------------------------------------------------
make_heatmap <- function(dat, lbl) {
  pve_labs_l <- pve_list[[lbl]]
  dat2 <- copy(dat)
  dat2[, PC_lab := factor(pve_labs_l[as.character(PC)],
                           levels = pve_labs_l[pc_order])]
  sig <- -log10(0.05)
  ggplot(dat2, aes(covariate, PC_lab, fill = neglog10p)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(aes(label = ifelse(!is.na(neglog10p) & neglog10p > sig,
                                 sprintf("%.0f", neglog10p), "")),
              size = 1.9, colour = "grey20") +
    scale_fill_gradient(low = "white", high = "#D62728",
                        name = expression(-log[10](p)),
                        na.value = "grey92", limits = c(0, 35),
                        oob = scales::squish) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL, title = lbl) +
    theme_masld(base_size = 6.5) +
    theme(axis.text.x  = element_text(angle = 40, hjust = 0, size = 6),
          axis.text.y  = element_text(size = 6),
          axis.ticks   = element_blank(),
          panel.grid   = element_blank(),
          legend.position = "none",
          plot.title   = element_text(size = 7.5, face = "bold"),
          plot.margin  = margin(2, 6, 2, 6))
}

panels <- lapply(lbl_order, function(lbl)
  make_heatmap(res_all[label == lbl], lbl))

# shared legend from first panel
leg_p <- ggplot(data.table(x=1,y=1,z=17.5), aes(x,y,fill=z)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="#D62728",
                      name=expression(-log[10](p)), limits=c(0,35)) +
  theme_void() +
  theme(legend.position="right",
        legend.title=element_text(size=6.5), legend.text=element_text(size=6),
        legend.key.height=unit(0.5,"cm"), legend.key.width=unit(0.25,"cm"))
leg <- cowplot::get_legend(leg_p)

library(cowplot)
grid <- plot_grid(plotlist = panels, nrow = 2, ncol = 2, align = "hv")
fig  <- plot_grid(grid, leg, nrow = 1, rel_widths = c(1, 0.07))

ggsave(file.path(OUT, "panelJ_covariate_pc_association.pdf"), fig,
       width = 10, height = 7.5, device = cairo_pdf)
fwrite(res_all, file.path(OUT, "panelJ_covariate_pc_association_data.csv"))
cat("Wrote panelJ_covariate_pc_association.pdf\n")

cat("\nTop associations per correction level:\n")
print(res_all[, .SD[order(-neglog10p)][1:3],
              by=label][, .(label, PC, covariate, neglog10p=round(neglog10p,1))])
