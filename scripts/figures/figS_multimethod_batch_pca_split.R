#!/usr/bin/env Rscript
# figS_multimethod_batch_pca_split.R
# PCA panels showing progressive batch correction.
# One combined PDF per correction level (cohort | disease | sex side by side).
# Sex inferred from XIST/DDX3Y k-means for samples with missing annotation.
#
# Outputs (all in panels/):
#   panelJ_pca_raw.pdf
#   panelJ_pca_batch.pdf
#   panelJ_pca_batch_sex.pdf
#   panelJ_pca_batch_sex_sva.pdf
#   panelJ_3d_raw.pdf
#   panelJ_3d_batch.pdf
#   panelJ_3d_batch_sex.pdf
#   panelJ_3d_batch_sex_sva.pdf
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(edgeR); library(limma); library(matrixStats); library(sva)
  library(reticulate); library(plotly)
})
reticulate::use_python(system("which python", intern = TRUE), required = TRUE)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
OUT <- file.path(BASE,
  "figures/supplementary/figS_methods_validation/multimethod_validation/panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CTRL <- "#9E9E9E"

MEGA <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
cohort_short <- c(GSE126848 = "Suppli", GSE130970 = "Hoang",
                  GSE135251 = "Govaere", GSE162694 = "Bril", GSE213621 = "Chen")
cohort_pal   <- c(Suppli = "#1F77B4", Hoang = "#FF7F0E", Govaere = "#2CA02C",
                  Bril = "#D62728", Chen = "#9467BD")
disease_pal  <- c(Control = CTRL, Disease = masld_colors$nash)
sex_pal      <- c(F = masld_colors$nash, M = "#1F77B4", Unknown = "grey80")

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
dge  <- load_merged_dge(); stopifnot(!is.null(dge))
samp <- as.data.table(dge$samples, keep.rownames = "sample_id")
keep <- samp$dataset %in% MEGA
dge  <- dge[, keep]; samp <- samp[keep]
cat(sprintf("n = %d samples\n", ncol(dge)))

grp <- factor(samp$group_binary, levels = c("Control", "Disease"))
ds  <- factor(samp$dataset)

# ---------------------------------------------------------------------------
# 2. Infer sex from XIST / DDX3Y for samples missing annotation
# ---------------------------------------------------------------------------
logcpm_all <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

xist  <- rownames(dge)[grep("^ENSG00000229807", rownames(dge))]  # XIST
ddy   <- rownames(dge)[grep("^ENSG00000067048", rownames(dge))]  # DDX3Y
cat(sprintf("Sex marker genes found: XIST=%d  DDX3Y=%d\n",
            length(xist), length(ddy)))

sex_inferred <- samp$sex
if (length(xist) > 0 && length(ddy) > 0) {
  expr_sex <- t(logcpm_all[c(xist[1], ddy[1]), , drop = FALSE])
  km <- kmeans(expr_sex, centers = 2, nstart = 20, iter.max = 50)
  # cluster with higher XIST = Female
  xist_mean <- tapply(expr_sex[, 1], km$cluster, mean)
  female_cluster <- as.integer(names(which.max(xist_mean)))
  inferred <- ifelse(km$cluster == female_cluster, "F", "M")
  # fill in only where annotation is missing
  missing <- is.na(samp$sex) | samp$sex == ""
  sex_inferred[missing] <- inferred[missing]
  cat(sprintf("Sex inference: %d annotated + %d inferred from XIST/DDX3Y\n",
              sum(!missing), sum(missing)))
}
cat(sprintf("Final sex: F=%d  M=%d  Unknown=%d\n",
            sum(sex_inferred == "F", na.rm = TRUE),
            sum(sex_inferred == "M", na.rm = TRUE),
            sum(is.na(sex_inferred) | sex_inferred == "")))
sex_inferred[is.na(sex_inferred) | sex_inferred == ""] <- "Unknown"

# ---------------------------------------------------------------------------
# 3. HVG selection (within-cohort centred) and corrections
# ---------------------------------------------------------------------------
logcpm <- logcpm_all
logcpm_wc <- logcpm
for (d in unique(samp$dataset)) {
  idx <- which(samp$dataset == d)
  logcpm_wc[, idx] <- logcpm[, idx] - rowMeans(logcpm[, idx, drop = FALSE])
}
rv  <- matrixStats::rowVars(logcpm_wc)
top <- order(rv, decreasing = TRUE)[seq_len(2000L)]
X   <- logcpm[top, ]

des_grp <- model.matrix(~ grp)
sex_cov <- as.numeric(sex_inferred == "F")

# batch only
X_batch <- limma::removeBatchEffect(X, batch = ds, design = des_grp)

# batch + sex
X_batch_sex <- limma::removeBatchEffect(X, batch = ds,
                                         covariates = sex_cov,
                                         design = des_grp)

# batch + sex + SVA
cat("Running SVA...\n")
mod0    <- model.matrix(~ ds + sex_cov)
mod1    <- model.matrix(~ grp + ds + sex_cov)
n_sv    <- sva::num.sv(X, mod1, method = "leek")
cat(sprintf("SVA: %d surrogate variables\n", n_sv))
sva_fit <- sva::sva(X, mod1, mod0, n.sv = n_sv)
X_sva   <- limma::removeBatchEffect(X, batch = ds,
                                     covariates = cbind(sex_cov, sva_fit$sv),
                                     design = des_grp)

# ---------------------------------------------------------------------------
# 4. PCA (3 components each)
# ---------------------------------------------------------------------------
do_pca3 <- function(mat) {
  pr  <- prcomp(t(mat), center = TRUE, scale. = FALSE)
  pve <- round(100 * pr$sdev^2 / sum(pr$sdev^2), 1)
  list(scores = pr$x[, 1:3], pve = pve[1:3])
}
pca_raw <- do_pca3(X);           cat("Raw          PVE:", pca_raw$pve,    "\n")
pca_bat <- do_pca3(X_batch);     cat("Batch        PVE:", pca_bat$pve,    "\n")
pca_sex <- do_pca3(X_batch_sex); cat("Batch+sex    PVE:", pca_sex$pve,    "\n")
pca_sva <- do_pca3(X_sva);       cat("Batch+sex+SVA PVE:", pca_sva$pve,  "\n")

make_dt <- function(scores, pve) {
  data.table(samp[, .(sample_id, dataset, group_binary)],
             sex     = sex_inferred,
             cohort  = factor(cohort_short[samp$dataset], levels = unname(cohort_short)),
             disease = factor(samp$group_binary, levels = c("Control", "Disease")),
             PC1 = scores[, 1], PC2 = scores[, 2], PC3 = scores[, 3],
             pve1 = pve[1], pve2 = pve[2], pve3 = pve[3])
}
dt_raw <- make_dt(pca_raw$scores, pca_raw$pve)
dt_bat <- make_dt(pca_bat$scores, pca_bat$pve)
dt_sex <- make_dt(pca_sex$scores, pca_sex$pve)
dt_sva <- make_dt(pca_sva$scores, pca_sva$pve)

# ---------------------------------------------------------------------------
# 5. Plot helpers
# ---------------------------------------------------------------------------
base_theme <- function() {
  theme_masld(base_size = 7) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = "right",
          legend.title = element_text(size = 6.5, face = "bold"),
          legend.text  = element_text(size = 6),
          legend.key.size = unit(0.28, "cm"),
          plot.title  = element_text(size = 7.5, face = "bold"))
}

scatter2d <- function(dt, colour_by, pal, leg_name, title) {
  ggplot(dt, aes(PC1, PC2, colour = .data[[colour_by]])) +
    geom_point(size = 0.5, alpha = 0.65) +
    scale_colour_manual(values = pal, name = leg_name) +
    guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1))) +
    labs(x = sprintf("PC1 (%.1f%%)", dt$pve1[1]),
         y = sprintf("PC2 (%.1f%%)", dt$pve2[1]),
         title = title) +
    base_theme()
}

# combined figure: cohort | disease | sex  (3 panels per correction level)
make_combined <- function(dt, subtitle) {
  pC <- scatter2d(dt, "cohort",  cohort_pal,  "Cohort",  "by Cohort")
  pD <- scatter2d(dt, "disease", disease_pal, "Disease", "by Disease")
  pS <- scatter2d(dt, "sex",     sex_pal,     "Sex",     "by Sex")
  (pC | pD | pS) +
    plot_annotation(subtitle = subtitle,
      theme = theme(plot.subtitle = element_text(size = 6.5, colour = "grey35")))
}

# ---------------------------------------------------------------------------
# 6. Save 2D combined figures
# ---------------------------------------------------------------------------
save_combined <- function(fig, filename, title, width = 10, height = 3.8) {
  fig2 <- fig + plot_annotation(
    title = title,
    theme = theme(plot.title = element_text(size = 9, face = "bold")))
  ggsave(file.path(OUT, filename), fig2, width = width, height = height,
         device = cairo_pdf)
  cat("Wrote", filename, "\n")
}

save_combined(
  make_combined(dt_raw, sprintf("Top-2000 HVGs, no correction. PC1=%.1f%% PC2=%.1f%%",
                                 dt_raw$pve1[1], dt_raw$pve2[1])),
  "panelJ_pca_raw.pdf",
  "Raw logCPM — cohort dominates, disease invisible, sex axis visible")

save_combined(
  make_combined(dt_bat, sprintf("Batch (cohort) corrected. PC1=%.1f%% PC2=%.1f%%",
                                 dt_bat$pve1[1], dt_bat$pve2[1])),
  "panelJ_pca_batch.pdf",
  "Batch-corrected — sex becomes the dominant binary structure within cohorts")

save_combined(
  make_combined(dt_sex, sprintf("Batch + sex corrected. PC1=%.1f%% PC2=%.1f%%",
                                 dt_sex$pve1[1], dt_sex$pve2[1])),
  "panelJ_pca_batch_sex.pdf",
  "Batch + sex corrected — diagonal residual from unmeasured confounders (BMI, age, RIN)")

save_combined(
  make_combined(dt_sva, sprintf("Batch + sex + SVA (%d SVs). PC1=%.1f%% PC2=%.1f%%",
                                 n_sv, dt_sva$pve1[1], dt_sva$pve2[1])),
  "panelJ_pca_batch_sex_sva.pdf",
  sprintf("Batch + sex + SVA (%d SVs) — unmeasured confounders removed; disease still not dominant in PCA", n_sv))

# ---------------------------------------------------------------------------
# 7. Save 3D figures (plotly + kaleido)
# ---------------------------------------------------------------------------
save_3d <- function(dt, subtitle_str, filename) {
  mk <- function(dt2, colour_by, pal, title_str) {
    plot_ly(dt2, x = ~PC1, y = ~PC2, z = ~PC3,
            color = dt2[[colour_by]], colors = pal,
            type = "scatter3d", mode = "markers",
            marker = list(size = 2, opacity = 0.7,
                          line = list(width = 0))) |>
      layout(
        title = list(text = title_str, font = list(size = 9)),
        scene = list(
          xaxis = list(title = sprintf("PC1 (%.1f%%)", dt2$pve1[1]),
                       titlefont = list(size = 8), tickfont = list(size = 6)),
          yaxis = list(title = sprintf("PC2 (%.1f%%)", dt2$pve2[1]),
                       titlefont = list(size = 8), tickfont = list(size = 6)),
          zaxis = list(title = sprintf("PC3 (%.1f%%)", dt2$pve3[1]),
                       titlefont = list(size = 8), tickfont = list(size = 6)),
          camera = list(eye = list(x = 1.5, y = 1.5, z = 0.8))
        ),
        legend = list(font = list(size = 7)),
        margin = list(l = 0, r = 0, t = 35, b = 0))
  }
  fig <- subplot(
    mk(dt, "cohort",  cohort_pal,  "by Cohort"),
    mk(dt, "disease", disease_pal, "by Disease"),
    mk(dt, "sex",     sex_pal,     "by Sex"),
    nrows = 1, shareX = FALSE, shareY = FALSE, titleX = TRUE, titleY = TRUE
  ) |> layout(title = list(text = subtitle_str, font = list(size = 10)))

  out_path <- file.path(OUT, filename)
  tryCatch(
    plotly::save_image(fig, out_path, width = 1400, height = 520, scale = 2),
    error = function(e) {
      htmlwidgets::saveWidget(fig, sub("\\.pdf$", ".html", out_path),
                              selfcontained = TRUE)
      cat("  kaleido fallback -> HTML\n")
    }
  )
  cat("Wrote", filename, "\n")
}

save_3d(dt_raw, "Raw logCPM",         "panelJ_3d_raw.pdf")
save_3d(dt_bat, "Batch-corrected",    "panelJ_3d_batch.pdf")
save_3d(dt_sex, "Batch+sex corrected","panelJ_3d_batch_sex.pdf")
save_3d(dt_sva, sprintf("Batch+sex+SVA (%d SVs)", n_sv), "panelJ_3d_batch_sex_sva.pdf")

# ---------------------------------------------------------------------------
# 8. Raw vs. SVA comparison: 2-row x 3-col 3D figure (6 scenes)
# Uses plotly multi-scene layout so both rows share the same figure.
# Row 1 = Raw, Row 2 = Batch+sex+SVA
# Cols = Cohort | Disease | Sex
# ---------------------------------------------------------------------------
scene_spec <- function(dt, col_x, row_y, w, h) {
  list(
    domain = list(x = c(col_x, col_x + w), y = c(row_y, row_y + h)),
    xaxis  = list(title = sprintf("PC1 %.1f%%", dt$pve1[1]),
                  titlefont = list(size = 7), tickfont = list(size = 5.5)),
    yaxis  = list(title = sprintf("PC2 %.1f%%", dt$pve2[1]),
                  titlefont = list(size = 7), tickfont = list(size = 5.5)),
    zaxis  = list(title = sprintf("PC3 %.1f%%", dt$pve3[1]),
                  titlefont = list(size = 7), tickfont = list(size = 5.5)),
    camera = list(eye = list(x = 1.5, y = 1.5, z = 0.8))
  )
}

# add all traces for one coloring into the figure, assigning a specific scene
add_group_traces <- function(fig, dt, colour_by, pal, scene_id, show_legend) {
  for (grp in unique(dt[[colour_by]])) {
    sub <- dt[dt[[colour_by]] == grp, ]
    fig <- fig |> add_trace(
      x = sub$PC1, y = sub$PC2, z = sub$PC3,
      type = "scatter3d", mode = "markers", name = as.character(grp),
      marker = list(color = pal[as.character(grp)], size = 2.2,
                    opacity = 0.72, line = list(width = 0)),
      scene = scene_id, showlegend = show_legend,
      legendgroup = as.character(grp))
  }
  fig
}

gap <- 0.03; w <- (1 - 2*gap)/3; h <- 0.46
y_top <- 0.52; y_bot <- 0.01

fig_cmp <- plot_ly()
fig_cmp <- add_group_traces(fig_cmp, dt_raw, "cohort",  cohort_pal,  "scene",  TRUE)
fig_cmp <- add_group_traces(fig_cmp, dt_raw, "disease", disease_pal, "scene2", TRUE)
fig_cmp <- add_group_traces(fig_cmp, dt_raw, "sex",     sex_pal,     "scene3", TRUE)
fig_cmp <- add_group_traces(fig_cmp, dt_sva, "cohort",  cohort_pal,  "scene4", FALSE)
fig_cmp <- add_group_traces(fig_cmp, dt_sva, "disease", disease_pal, "scene5", FALSE)
fig_cmp <- add_group_traces(fig_cmp, dt_sva, "sex",     sex_pal,     "scene6", FALSE)
fig_cmp <- fig_cmp |> layout(
    scene  = scene_spec(dt_raw, 0,           y_top, w, h),
    scene2 = scene_spec(dt_raw, w + gap,     y_top, w, h),
    scene3 = scene_spec(dt_raw, 2*(w+gap),   y_top, w, h),
    scene4 = scene_spec(dt_sva, 0,           y_bot, w, h),
    scene5 = scene_spec(dt_sva, w + gap,     y_bot, w, h),
    scene6 = scene_spec(dt_sva, 2*(w+gap),   y_bot, w, h),
    annotations = list(
      list(text = "RAW", x = 0.01, y = 0.99,
           xref="paper", yref="paper", showarrow=FALSE,
           font=list(size=12, color="grey25"), xanchor="left"),
      list(text = sprintf("BATCH + SEX + SVA (%d SV)", n_sv),
           x = 0.01, y = 0.48,
           xref="paper", yref="paper", showarrow=FALSE,
           font=list(size=12, color="grey25"), xanchor="left"),
      list(text="Cohort",  x=w/2,           y=1.02, xref="paper",
           yref="paper", showarrow=FALSE, font=list(size=9), xanchor="center"),
      list(text="Disease", x=w+gap+w/2,     y=1.02, xref="paper",
           yref="paper", showarrow=FALSE, font=list(size=9), xanchor="center"),
      list(text="Sex",     x=2*(w+gap)+w/2, y=1.02, xref="paper",
           yref="paper", showarrow=FALSE, font=list(size=9), xanchor="center")
    ),
    legend = list(font=list(size=7)),
    margin = list(l=0, r=0, t=45, b=0)
  )

out_cmp <- file.path(OUT, "panelJ_3d_comparison_raw_vs_sva.pdf")
tryCatch(
  plotly::save_image(fig_cmp, out_cmp, width = 1400, height = 900, scale = 2),
  error = function(e) {
    htmlwidgets::saveWidget(fig_cmp,
      sub("\\.pdf$", ".html", out_cmp), selfcontained = TRUE)
    cat("  kaleido fallback -> HTML\n")
  }
)
cat("Wrote panelJ_3d_comparison_raw_vs_sva.pdf\n")

cat("\nDone. All panels in", OUT, "\n")
