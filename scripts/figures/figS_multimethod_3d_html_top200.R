#!/usr/bin/env Rscript
# figS_multimethod_3d_html.R
# Generates interactive HTML versions of the 3D PCA panels showing progressive
# batch correction. Mirrors figS_multimethod_batch_pca_split.R but:
#   - skips the per-gene lme4 loop (only needed for the 2D batch-model-PCA PDFs)
#   - saves HTML via htmlwidgets::saveWidget (selfcontained = TRUE)
#
# Outputs (all in panels/):
#   3d_raw.html
#   3d_batch.html
#   3d_batch_sex.html
#   3d_batch_sex_sva.html
#   3d_comparison_raw_vs_sva.html
suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma)
  library(matrixStats); library(sva); library(plotly); library(htmlwidgets)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
OUT <- file.path(BASE,
  "figures/supplementary/figS_methods_validation/multimethod_validation/panels/pca_top200")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

MEGA        <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
cohort_short <- c(GSE126848 = "GSE126848", GSE130970 = "GSE130970",
                  GSE135251 = "GSE135251", GSE162694 = "GSE162694", GSE213621 = "GSE213621")
cohort_pal   <- c(GSE126848 = "#1F77B4", GSE130970 = "#FF7F0E", GSE135251 = "#2CA02C",
                  GSE162694 = "#D62728", GSE213621 = "#9467BD")
CTRL         <- "#9E9E9E"
disease_pal  <- c(Control = CTRL, Disease = "#C0392B")
sex_pal      <- c(F = "#C0392B", M = "#1F77B4", Unknown = "grey80")

# ---------------------------------------------------------------------------
# 1. Load and filter to mega cohorts
# ---------------------------------------------------------------------------
dge  <- load_merged_dge(); stopifnot(!is.null(dge))
samp <- as.data.table(dge$samples, keep.rownames = "sample_id")
keep <- samp$dataset %in% MEGA
dge  <- dge[, keep]; samp <- samp[keep]
cat(sprintf("n = %d samples across %d cohorts\n", ncol(dge),
            length(unique(samp$dataset))))

grp <- factor(samp$group_binary, levels = c("Control", "Disease"))
ds  <- factor(samp$dataset)

# ---------------------------------------------------------------------------
# 2. Sex inference from XIST / DDX3Y
# ---------------------------------------------------------------------------
logcpm_all <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
xist <- rownames(dge)[grep("^ENSG00000229807", rownames(dge))]
ddy  <- rownames(dge)[grep("^ENSG00000067048", rownames(dge))]
cat(sprintf("Sex markers: XIST=%d  DDX3Y=%d\n", length(xist), length(ddy)))

sex_inferred <- samp$sex
if (length(xist) > 0 && length(ddy) > 0) {
  expr_sex <- t(logcpm_all[c(xist[1], ddy[1]), , drop = FALSE])
  km <- kmeans(expr_sex, centers = 2, nstart = 20, iter.max = 50)
  xist_mean     <- tapply(expr_sex[, 1], km$cluster, mean)
  female_cluster <- as.integer(names(which.max(xist_mean)))
  inferred  <- ifelse(km$cluster == female_cluster, "F", "M")
  missing   <- is.na(samp$sex) | samp$sex == ""
  sex_inferred[missing] <- inferred[missing]
  cat(sprintf("Sex: %d annotated + %d inferred\n", sum(!missing), sum(missing)))
}
sex_inferred[is.na(sex_inferred) | sex_inferred == ""] <- "Unknown"

# ---------------------------------------------------------------------------
# 3. HVG selection (within-cohort centred) and batch corrections
# ---------------------------------------------------------------------------
logcpm <- logcpm_all
logcpm_wc <- logcpm
for (d in unique(samp$dataset)) {
  idx <- which(samp$dataset == d)
  logcpm_wc[, idx] <- logcpm[, idx] - rowMeans(logcpm[, idx, drop = FALSE])
}
rv  <- matrixStats::rowVars(logcpm_wc)
top <- order(rv, decreasing = TRUE)[seq_len(200L)]
X   <- logcpm[top, ]

des_grp <- model.matrix(~ grp)
sex_cov <- as.numeric(sex_inferred == "F")

X_batch     <- limma::removeBatchEffect(X, batch = ds, design = des_grp)
X_batch_sex <- limma::removeBatchEffect(X, batch = ds,
                                         covariates = sex_cov, design = des_grp)

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
do_pca3 <- function(mat, label) {
  pr  <- prcomp(t(mat), center = TRUE, scale. = FALSE)
  pve <- round(100 * pr$sdev^2 / sum(pr$sdev^2), 1)
  cat(sprintf("%-22s PVE: PC1=%.1f%% PC2=%.1f%% PC3=%.1f%%\n",
              label, pve[1], pve[2], pve[3]))
  list(scores = pr$x[, 1:3], pve = pve[1:3])
}

pca_raw <- do_pca3(X,           "Raw")
pca_bat <- do_pca3(X_batch,     "Batch")
pca_sex <- do_pca3(X_batch_sex, "Batch+sex")
pca_sva <- do_pca3(X_sva,       "Batch+sex+SVA")

make_dt <- function(scores, pve) {
  data.table(
    sample_id = samp$sample_id,
    cohort    = factor(cohort_short[samp$dataset], levels = unname(cohort_short)),
    disease   = factor(samp$group_binary, levels = c("Control", "Disease")),
    sex       = sex_inferred,
    PC1 = scores[, 1], PC2 = scores[, 2], PC3 = scores[, 3],
    pve1 = pve[1], pve2 = pve[2], pve3 = pve[3])
}
dt_raw <- make_dt(pca_raw$scores, pca_raw$pve)
dt_bat <- make_dt(pca_bat$scores, pca_bat$pve)
dt_sex <- make_dt(pca_sex$scores, pca_sex$pve)
dt_sva <- make_dt(pca_sva$scores, pca_sva$pve)

# save PCA coordinates for future lightweight re-use
fwrite(dt_raw, file.path(OUT, "3d_pca_raw.csv"))
fwrite(dt_bat, file.path(OUT, "3d_pca_batch.csv"))
fwrite(dt_sex, file.path(OUT, "3d_pca_batch_sex.csv"))
fwrite(dt_sva, file.path(OUT, "3d_pca_batch_sex_sva.csv"))
cat("PCA data CSVs written\n")

# ---------------------------------------------------------------------------
# 5. Plotly helpers
# ---------------------------------------------------------------------------
mk_trace <- function(dt, colour_by, pal, axis_labels) {
  plot_ly(dt,
    x = ~PC1, y = ~PC2, z = ~PC3,
    color  = dt[[colour_by]], colors = pal,
    type   = "scatter3d", mode = "markers",
    marker = list(size = 2.5, opacity = 0.72, line = list(width = 0))) |>
    layout(scene = list(
      xaxis = list(title = axis_labels[1], titlefont = list(size = 9),
                   tickfont = list(size = 7)),
      yaxis = list(title = axis_labels[2], titlefont = list(size = 9),
                   tickfont = list(size = 7)),
      zaxis = list(title = axis_labels[3], titlefont = list(size = 9),
                   tickfont = list(size = 7)),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 0.8))
    ),
    legend = list(font = list(size = 9)),
    margin = list(l = 0, r = 0, t = 45, b = 0))
}

ax_labels <- function(dt)
  sprintf(c("PC1 (%.1f%%)", "PC2 (%.1f%%)", "PC3 (%.1f%%)"),
          dt$pve1[1], dt$pve2[1], dt$pve3[1])

make_3panel <- function(dt, title_str) {
  subplot(
    mk_trace(dt, "cohort",  cohort_pal,  ax_labels(dt)),
    mk_trace(dt, "disease", disease_pal, ax_labels(dt)),
    mk_trace(dt, "sex",     sex_pal,     ax_labels(dt)),
    nrows = 1, shareX = FALSE, shareY = FALSE,
    titleX = TRUE, titleY = TRUE
  ) |>
    layout(
      title  = list(text = title_str, font = list(size = 13, color = "#333")),
      legend = list(font = list(size = 9)),
      margin = list(l = 10, r = 10, t = 60, b = 10),
      annotations = list(
        list(text = "Cohort",  x = 0.13, y = 1.06, xref = "paper",
             yref = "paper", showarrow = FALSE, font = list(size = 11)),
        list(text = "Disease", x = 0.50, y = 1.06, xref = "paper",
             yref = "paper", showarrow = FALSE, font = list(size = 11)),
        list(text = "Sex",     x = 0.87, y = 1.06, xref = "paper",
             yref = "paper", showarrow = FALSE, font = list(size = 11))
      )
    )
}

save_html <- function(fig, filename) {
  path <- file.path(OUT, filename)
  htmlwidgets::saveWidget(fig, path, selfcontained = TRUE,
                          title = sub("\\.html$", "", filename))
  cat("Wrote", filename, "\n")
}

# ---------------------------------------------------------------------------
# 6. Save the four single-correction HTML files
# ---------------------------------------------------------------------------
save_html(
  make_3panel(dt_raw, "Raw logCPM — 3D PCA (no correction)<br><sub>Cohort | Disease | Sex</sub>"),
  "3d_raw.html")

save_html(
  make_3panel(dt_bat, "Batch-corrected — 3D PCA<br><sub>Cohort | Disease | Sex</sub>"),
  "3d_batch.html")

save_html(
  make_3panel(dt_sex, "Batch + Sex corrected — 3D PCA<br><sub>Cohort | Disease | Sex</sub>"),
  "3d_batch_sex.html")

save_html(
  make_3panel(dt_sva,
    sprintf("Batch + Sex + SVA (%d SVs) — 3D PCA<br><sub>Cohort | Disease | Sex</sub>", n_sv)),
  "3d_batch_sex_sva.html")

# ---------------------------------------------------------------------------
# 7. Raw vs. SVA comparison: 2-row × 3-col (6 scenes)
# ---------------------------------------------------------------------------
gap <- 0.03; w <- (1 - 2 * gap) / 3; h <- 0.44
y_top <- 0.53; y_bot <- 0.02

scene_spec <- function(dt, col_x, row_y) {
  list(
    domain = list(x = c(col_x, col_x + w), y = c(row_y, row_y + h)),
    xaxis  = list(title = sprintf("PC1 %.1f%%", dt$pve1[1]),
                  titlefont = list(size = 8), tickfont = list(size = 6)),
    yaxis  = list(title = sprintf("PC2 %.1f%%", dt$pve2[1]),
                  titlefont = list(size = 8), tickfont = list(size = 6)),
    zaxis  = list(title = sprintf("PC3 %.1f%%", dt$pve3[1]),
                  titlefont = list(size = 8), tickfont = list(size = 6)),
    camera = list(eye = list(x = 1.5, y = 1.5, z = 0.8))
  )
}

add_traces <- function(fig, dt, colour_by, pal, scene_id, show_legend) {
  for (g in unique(dt[[colour_by]])) {
    sub <- dt[dt[[colour_by]] == g, ]
    fig <- fig |> add_trace(
      x = sub$PC1, y = sub$PC2, z = sub$PC3,
      type = "scatter3d", mode = "markers",
      name = as.character(g),
      marker = list(color = pal[as.character(g)], size = 2.5,
                    opacity = 0.72, line = list(width = 0)),
      scene = scene_id, showlegend = show_legend,
      legendgroup = as.character(g))
  }
  fig
}

fig_cmp <- plot_ly()
fig_cmp <- add_traces(fig_cmp, dt_raw, "cohort",  cohort_pal,  "scene",  TRUE)
fig_cmp <- add_traces(fig_cmp, dt_raw, "disease", disease_pal, "scene2", TRUE)
fig_cmp <- add_traces(fig_cmp, dt_raw, "sex",     sex_pal,     "scene3", TRUE)
fig_cmp <- add_traces(fig_cmp, dt_sva, "cohort",  cohort_pal,  "scene4", FALSE)
fig_cmp <- add_traces(fig_cmp, dt_sva, "disease", disease_pal, "scene5", FALSE)
fig_cmp <- add_traces(fig_cmp, dt_sva, "sex",     sex_pal,     "scene6", FALSE)
fig_cmp <- fig_cmp |> layout(
  scene  = scene_spec(dt_raw, 0,           y_top),
  scene2 = scene_spec(dt_raw, w + gap,     y_top),
  scene3 = scene_spec(dt_raw, 2*(w+gap),   y_top),
  scene4 = scene_spec(dt_sva, 0,           y_bot),
  scene5 = scene_spec(dt_sva, w + gap,     y_bot),
  scene6 = scene_spec(dt_sva, 2*(w+gap),   y_bot),
  annotations = list(
    list(text = "RAW", x = 0.01, y = 0.99,
         xref = "paper", yref = "paper", showarrow = FALSE,
         font = list(size = 13, color = "grey25"), xanchor = "left"),
    list(text = sprintf("BATCH + SEX + SVA (%d SV)", n_sv),
         x = 0.01, y = 0.47,
         xref = "paper", yref = "paper", showarrow = FALSE,
         font = list(size = 13, color = "grey25"), xanchor = "left"),
    list(text = "Cohort",  x = w/2,           y = 1.02, xref = "paper",
         yref = "paper", showarrow = FALSE, font = list(size = 10), xanchor = "center"),
    list(text = "Disease", x = w+gap+w/2,     y = 1.02, xref = "paper",
         yref = "paper", showarrow = FALSE, font = list(size = 10), xanchor = "center"),
    list(text = "Sex",     x = 2*(w+gap)+w/2, y = 1.02, xref = "paper",
         yref = "paper", showarrow = FALSE, font = list(size = 10), xanchor = "center")
  ),
  legend = list(font = list(size = 8)),
  margin = list(l = 0, r = 0, t = 55, b = 0)
)

save_html(fig_cmp, "3d_comparison_raw_vs_sva.html")

cat("\nDone. All 5 HTML panels in", OUT, "\n")
