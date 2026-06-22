#!/usr/bin/env Rscript
# figS_supervised_disease_pca_top100DEG.R
# SUPERVISED disease-separation PCA: same 846 mega samples, but genes selected by
# canonical DEG |t| (top 200) instead of by variance. Batch+sex corrected logCPM
# (same recipe as the unsupervised panels) so cohort/sex axes are removed.
#
# CAVEAT (printed on the figure): DEGs were derived from THESE SAME samples, so the
# control/disease separation is partly circular ("double dipping") — descriptive only,
# NOT an out-of-sample classifier. For honest accuracy use the LOOCV/held-out panels.
#
# Outputs -> panels/pca_top100/
#   supervised_disease_top100DEG.pdf        (2D: cohort | disease | sex)
#   supervised_disease_top100DEG_3d.html    (3D interactive: cohort | disease | sex)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(edgeR); library(limma); library(plotly); library(htmlwidgets)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
OUT <- file.path(BASE,
  "figures/supplementary/figS_methods_validation/multimethod_validation/panels/pca_top100")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CTRL <- "#9E9E9E"
NDEG <- 100L

MEGA <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
cohort_short <- c(GSE126848 = "GSE126848", GSE130970 = "GSE130970", GSE135251 = "GSE135251",
                  GSE162694 = "GSE162694", GSE213621 = "GSE213621")
cohort_pal  <- c(GSE126848 = "#1F77B4", GSE130970 = "#FF7F0E", GSE135251 = "#2CA02C",
                 GSE162694 = "#D62728", GSE213621 = "#9467BD")
disease_pal <- c(Control = CTRL, Disease = "#C0392B")
sex_pal     <- c(F = "#C0392B", M = "#1F77B4", Unknown = "grey80")

# --- data: 5 mega cohorts ---------------------------------------------------
dge  <- load_merged_dge(); stopifnot(!is.null(dge))
samp <- as.data.table(dge$samples, keep.rownames = "sample_id")
keep <- samp$dataset %in% MEGA
dge  <- dge[, keep]; samp <- samp[keep]
cat(sprintf("n = %d samples, %d cohorts\n", ncol(dge), length(unique(samp$dataset))))
grp <- factor(samp$group_binary, levels = c("Control", "Disease"))
ds  <- factor(samp$dataset)

# --- sex inference (XIST/DDX3Y k-means, same as unsupervised panels) ---------
logcpm_all <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
xist <- rownames(dge)[grep("^ENSG00000229807", rownames(dge))]
ddy  <- rownames(dge)[grep("^ENSG00000067048", rownames(dge))]
sex_inferred <- samp$sex
if (length(xist) > 0 && length(ddy) > 0) {
  es <- t(logcpm_all[c(xist[1], ddy[1]), , drop = FALSE])
  km <- kmeans(es, centers = 2, nstart = 20, iter.max = 50)
  fem <- as.integer(names(which.max(tapply(es[, 1], km$cluster, mean))))
  inf <- ifelse(km$cluster == fem, "F", "M")
  miss <- is.na(samp$sex) | samp$sex == ""
  sex_inferred[miss] <- inf[miss]
}
sex_inferred[is.na(sex_inferred) | sex_inferred == ""] <- "Unknown"

# --- SUPERVISED gene selection: top-200 by canonical DEG |t| -----------------
deg <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"))
deg[, base := sub("\\..*", "", gene)]
rn_base <- sub("\\..*", "", rownames(dge))
deg <- deg[base %in% rn_base & is.finite(t)]
deg[, abst := abs(t)]
setorder(deg, -abst)                                    # rank by |t|
sel_base <- head(deg$base, NDEG)
sel_idx  <- match(sel_base, rn_base)
sel_idx  <- sel_idx[!is.na(sel_idx)]
cat(sprintf("Selected %d DEGs by |t| (top %d requested)\n", length(sel_idx), NDEG))

# --- batch+sex corrected logCPM on the selected DEGs ------------------------
X       <- logcpm_all[sel_idx, ]
des_grp <- model.matrix(~ grp)
sex_cov <- as.numeric(sex_inferred == "F")
X_corr  <- limma::removeBatchEffect(X, batch = ds, covariates = sex_cov, design = des_grp)

pr  <- prcomp(t(X_corr), center = TRUE, scale. = FALSE)
pve <- round(100 * pr$sdev^2 / sum(pr$sdev^2), 1)
cat(sprintf("Supervised PCA PVE: PC1=%.1f%% PC2=%.1f%% PC3=%.1f%%\n", pve[1], pve[2], pve[3]))

dt <- data.table(sample_id = samp$sample_id,
                 cohort  = factor(cohort_short[samp$dataset], levels = unname(cohort_short)),
                 disease = factor(samp$group_binary, levels = c("Control", "Disease")),
                 sex     = sex_inferred,
                 PC1 = pr$x[, 1], PC2 = pr$x[, 2], PC3 = pr$x[, 3])

# orient PC1 so Disease is on the positive side (cosmetic, for readability)
if (mean(dt$PC1[dt$disease == "Disease"]) < mean(dt$PC1[dt$disease == "Control"])) {
  dt[, PC1 := -PC1]; pr$x[, 1] <- -pr$x[, 1]
}

# honest separation metric: AUC of PC1 separating Disease vs Control (Mann-Whitney)
auc_pc1 <- {
  r <- rank(dt$PC1); nD <- sum(dt$disease == "Disease"); nC <- sum(dt$disease == "Control")
  (sum(r[dt$disease == "Disease"]) - nD * (nD + 1) / 2) / (nD * nC)
}
cat(sprintf("PC1 Disease-vs-Control AUC (in-sample, double-dipped) = %.3f\n", auc_pc1))
fwrite(dt, file.path(OUT, "supervised_disease_top100DEG_data.csv"))

# --- 2D figure: cohort | disease | sex --------------------------------------
base_theme <- function() theme_masld(base_size = 7) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
        legend.position = "right", legend.title = element_text(size = 6.5, face = "bold"),
        legend.text = element_text(size = 6), legend.key.size = unit(0.28, "cm"),
        plot.title = element_text(size = 7.5, face = "bold"))
sc <- function(cby, pal, nm, ti) ggplot(dt, aes(PC1, PC2, colour = .data[[cby]])) +
  geom_point(size = 0.6, alpha = 0.7) + scale_colour_manual(values = pal, name = nm) +
  guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1))) +
  labs(x = sprintf("PC1 (%.1f%%)", pve[1]), y = sprintf("PC2 (%.1f%%)", pve[2]), title = ti) +
  base_theme()
fig2d <- (sc("cohort", cohort_pal, "Cohort", "by Cohort") |
          sc("disease", disease_pal, "Disease", "by Disease") |
          sc("sex", sex_pal, "Sex", "by Sex")) +
  plot_annotation(
    title = sprintf("Supervised PCA: top-%d DEGs by |t|, batch+sex corrected (n=%d)", length(sel_idx), ncol(dge)),
    subtitle = sprintf("PC1 separates disease (AUC=%.2f) BY CONSTRUCTION — DEGs were selected on these same samples (double-dipping). Descriptive only; for honest accuracy see the LOOCV/held-out panels.", auc_pc1),
    theme = theme(plot.title = element_text(size = 9, face = "bold"),
                  plot.subtitle = element_text(size = 6.2, colour = "grey35")))
ggsave(file.path(OUT, "supervised_disease_top100DEG.pdf"), fig2d, width = 10, height = 3.8, device = cairo_pdf)
cat("Wrote supervised_disease_top100DEG.pdf\n")

# --- 3D interactive ---------------------------------------------------------
mk <- function(cby, pal) plot_ly(dt, x = ~PC1, y = ~PC2, z = ~PC3, color = dt[[cby]], colors = pal,
  type = "scatter3d", mode = "markers",
  marker = list(size = 2.5, opacity = 0.72, line = list(width = 0))) |>
  layout(scene = list(
    xaxis = list(title = sprintf("PC1 (%.1f%%)", pve[1]), titlefont = list(size = 9), tickfont = list(size = 7)),
    yaxis = list(title = sprintf("PC2 (%.1f%%)", pve[2]), titlefont = list(size = 9), tickfont = list(size = 7)),
    zaxis = list(title = sprintf("PC3 (%.1f%%)", pve[3]), titlefont = list(size = 9), tickfont = list(size = 7)),
    camera = list(eye = list(x = 1.5, y = 1.5, z = 0.8))))
fig3d <- subplot(mk("cohort", cohort_pal), mk("disease", disease_pal), mk("sex", sex_pal),
                 nrows = 1, shareX = FALSE, shareY = FALSE, titleX = TRUE, titleY = TRUE) |>
  layout(title = list(text = sprintf("Supervised 3D PCA: top-%d DEGs by |t| (batch+sex corrected) — Cohort | Disease | Sex<br><sub>in-sample PC1 AUC=%.2f (double-dipped, descriptive)</sub>", length(sel_idx), auc_pc1),
                      font = list(size = 12)),
         margin = list(l = 0, r = 0, t = 60, b = 0))
htmlwidgets::saveWidget(fig3d, file.path(OUT, "supervised_disease_top100DEG_3d.html"),
                        selfcontained = TRUE, title = "supervised_disease_top100DEG_3d")
cat("Wrote supervised_disease_top100DEG_3d.html\nSUPERVISED_DONE\n")
