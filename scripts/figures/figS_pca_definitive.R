#!/usr/bin/env Rscript
# figS_pca_definitive.R
# Definitive disease-vs-control PCA, with intermediates projected onto the
# definitive-derived axes.
#
# Motivation: "Control" is not one phenotype across the 5 mega cohorts (normal
# histology / fibrosis-F0 / clinical label), and disease spans NAFL -> NASH ->
# fibrosis. The standard PCA overlays everything, blurring the disease axis.
# Here we keep only the unambiguous extremes per cohort, fit the PCA on them,
# and PROJECT the eliminated (intermediate) samples to see where they fall.
#
# Definitive disease = NAS>=5 OR fibrosis>=3 (biopsy) | NASH label (Suppli).
# Definitive control  = NAS==0 & F0 (biopsy) | Control label (Suppli/Chen), unless
#                       already disease.  Everything else = intermediate.
#
# Outputs (main Fig 3 — figures/main/fig3_RNAseq/panels/, Fig S3U):
#   figs3u_pca_definitive_control_vs_disease.pdf  (2x3: raw | corrected x Disease/NAS/Fibrosis)
#   figs3u_pca_definitive_raw_nas_fib.pdf         (compact 1x3: raw Disease | corrected NAS | corrected Fibrosis)
#   data/pca_definitive_counts.csv                (per-cohort class counts)
# (The intermediate-projection figure figs3u_pca_definitive_projection.pdf was
#  retired 2026-06-23 at user request.)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(edgeR); library(limma); library(matrixStats)
})
grDevices::pdf.options(useDingbats = FALSE)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
# Output to main Fig 3 (RNA-seq) panels dir (moved here 2026-06-18; was the
# multimethod_validation/panels/pca QC dir). FIG2_DIR resolves to fig3_RNAseq.
OUT      <- file.path(FIG2_DIR, "panels")
DATA_DIR <- file.path(OUT, "data")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

MEGA <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
cohort_short <- c(GSE126848 = "GSE126848", GSE130970 = "GSE130970", GSE135251 = "GSE135251",
                  GSE162694 = "GSE162694", GSE213621 = "GSE213621")
CTRL <- "#9E9E9E"   # Liang control gray (invariant)
DIS  <- "#C9265E"   # Liang disease magenta

# ---------------------------------------------------------------------------
# 1. Load + subset to mega cohorts; merge clinical fields
# ---------------------------------------------------------------------------
dge  <- load_merged_dge(); stopifnot(!is.null(dge))
samp <- as.data.table(dge$samples, keep.rownames = "sample_id")
keep <- samp$dataset %in% MEGA; dge <- dge[, keep]; samp <- samp[keep]

meta <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"),
  select = c("sample_id", "condition", "diagnosis_harmonized", "fibrosis_stage", "nas_score"))
samp[meta, on = "sample_id",
     `:=`(condition = i.condition, diagnosis_harmonized = i.diagnosis_harmonized,
          fibrosis_stage = i.fibrosis_stage, nas_score = i.nas_score)]

# ---------------------------------------------------------------------------
# 2. Sex inference (XIST / DDX3Y k-means), numeric covariate (F vs not-F)
# ---------------------------------------------------------------------------
logcpm_all <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
xist <- rownames(dge)[grep("^ENSG00000229807", rownames(dge))]
ddy  <- rownames(dge)[grep("^ENSG00000067048", rownames(dge))]
sx <- samp$sex
if (length(xist) > 0 && length(ddy) > 0) {
  es <- t(logcpm_all[c(xist[1], ddy[1]), , drop = FALSE])
  km <- kmeans(es, centers = 2, nstart = 20, iter.max = 50)
  fem <- as.integer(names(which.max(tapply(es[, 1], km$cluster, mean))))
  inf <- ifelse(km$cluster == fem, "F", "M")
  mi <- is.na(samp$sex) | samp$sex == ""; sx[mi] <- inf[mi]
}
sex_cov <- as.numeric(sx == "F")

# ---------------------------------------------------------------------------
# 3. Definitive classifier  (disease takes precedence on conflict)
# ---------------------------------------------------------------------------
nas <- as.numeric(samp$nas_score); fib <- as.numeric(samp$fibrosis_stage)
cond <- samp$condition; dsid <- samp$dataset
def_disease <- (!is.na(nas) & nas >= 5) | (!is.na(fib) & fib >= 3) |
               (dsid == "GSE126848" & cond == "NASH")
def_control <- !def_disease &
               ((!is.na(nas) & nas == 0 & (is.na(fib) | fib == 0)) |   # NAS0 AND no fibrosis
                (dsid == "GSE126848" & cond == "Control") |            # Suppli: no biopsy
                (dsid == "GSE213621" & cond == "Control"))             # Chen control label (= F0)
samp[, def_class := fifelse(def_disease, "Disease",
                     fifelse(def_control, "Control", "Intermediate"))]
samp[, def_class := factor(def_class, levels = c("Control", "Disease", "Intermediate"))]

# severity proxy (0-1) for colouring the projected intermediates
samp[, sev := fcase(
  !is.na(nas_score),        pmin(as.numeric(nas_score), 8) / 8,
  !is.na(fibrosis_stage),   pmin(as.numeric(fibrosis_stage), 4) / 4,
  condition == "NASH",          0.85,
  condition == "NAFL",          0.30,
  condition == "Control_Obese", 0.08,
  condition == "Control",       0.0,
  default = NA_real_)]

cnt <- dcast(samp[, .N, by = .(cohort = factor(cohort_short[dsid], levels = unname(cohort_short)),
                               def_class)],
             cohort ~ def_class, value.var = "N", fill = 0)
fwrite(cnt, file.path(DATA_DIR, "pca_definitive_counts.csv"))
cat("Per-cohort class counts:\n"); print(cnt)
cat(sprintf("\nTotals: Control=%d  Disease=%d  Intermediate=%d  (n=%d)\n",
            sum(samp$def_class == "Control"), sum(samp$def_class == "Disease"),
            sum(samp$def_class == "Intermediate"), nrow(samp)))
stopifnot(sum(samp$def_class == "Control") > 10, sum(samp$def_class == "Disease") > 10)

def <- which(samp$def_class %in% c("Control", "Disease"))

# ---------------------------------------------------------------------------
# 4. HVG (within-cohort variance) on DEFINITIVE samples only
# ---------------------------------------------------------------------------
wc <- logcpm_all[, def]; sd_ds <- samp$dataset[def]
for (d in unique(sd_ds)) {
  i <- which(sd_ds == d); wc[, i] <- wc[, i] - rowMeans(wc[, i, drop = FALSE])
}
G <- order(matrixStats::rowVars(wc), decreasing = TRUE)[seq_len(2000L)]

# ---------------------------------------------------------------------------
# 5. Unsupervised batch+sex correction: fit on definitive, apply to all
#    (group NOT protected -> honest separation test)
# ---------------------------------------------------------------------------
ds_fac <- factor(samp$dataset, levels = MEGA)
mm_all <- model.matrix(~ ds_fac + sex_cov)        # intercept + dataset + sex
mm_d   <- mm_all[def, , drop = FALSE]
Xall   <- logcpm_all[G, ]; Xd <- Xall[, def]
beta   <- limma::lmFit(Xd, mm_d)$coefficients     # genes x params
rmc    <- 2:ncol(mm_all)                          # dataset + sex cols (keep intercept)
Xcorr  <- Xall - t(mm_all[, rmc, drop = FALSE] %*% t(beta[, rmc, drop = FALSE]))

# ---------------------------------------------------------------------------
# 6. PCA fit on definitive (corrected); project all samples
# ---------------------------------------------------------------------------
pca  <- prcomp(t(Xcorr[, def]), center = TRUE, scale. = FALSE)
pve  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
R    <- pca$rotation[, 1:3]; ctr <- pca$center
sc   <- scale(t(Xcorr), center = ctr, scale = FALSE) %*% R     # all samples x 3
dmask <- samp$def_class == "Disease"; cmask <- samp$def_class == "Control"
if (mean(sc[dmask, 1]) < mean(sc[cmask, 1])) { sc[, 1] <- -sc[, 1]; R[, 1] <- -R[, 1] }

# raw (uncorrected) definitive PCA, same gene set
pr_raw <- prcomp(t(Xall[, def]), center = TRUE, scale. = FALSE)
pve_raw <- round(100 * pr_raw$sdev^2 / sum(pr_raw$sdev^2), 1)
xr <- pr_raw$x[, 1:2]
cls_def <- samp$def_class[def]
if (mean(xr[cls_def == "Disease", 1]) < mean(xr[cls_def == "Control", 1])) xr[, 1] <- -xr[, 1]

# group-PROTECTED definitive PCA: remove dataset + sex while CONDITIONING on
# Control/Disease -> identical covariate structure to the canonical DEG model
# (~ dataset + sex + group_binary). DEG-faithful view / upper bound. Definitive
# samples only (no projection needed for this figure).
Xd_prot <- limma::removeBatchEffect(Xall[, def], batch = factor(samp$dataset[def]),
                                    covariates = sex_cov[def],
                                    design = model.matrix(~ cls_def))
pr_prot <- prcomp(t(Xd_prot), center = TRUE, scale. = FALSE)
pve_prot <- round(100 * pr_prot$sdev^2 / sum(pr_prot$sdev^2), 1)
xp <- pr_prot$x[, 1:2]
if (mean(xp[cls_def == "Disease", 1]) < mean(xp[cls_def == "Control", 1])) xp[, 1] <- -xp[, 1]

# ---------------------------------------------------------------------------
# 7. AUC (Mann-Whitney) of Disease vs Control on a PC score
# ---------------------------------------------------------------------------
auc1 <- function(score, pos) {
  r <- rank(score); nP <- sum(pos); nN <- sum(!pos)
  (sum(r[pos]) - nP * (nP + 1) / 2) / (nP * nN)
}

nas_def  <- as.numeric(samp$nas_score)[def]
fib_def  <- as.numeric(samp$fibrosis_stage)[def]
shared_keep <- cls_def == "Control" | (cls_def == "Disease" & !is.na(nas_def) & nas_def > 0 & !is.na(fib_def))

pos_def <- cls_def[shared_keep] == "Disease"
auc_raw  <- auc1(xr[shared_keep, 1], pos_def)
auc_corr <- auc1(sc[def, 1][shared_keep], pos_def)
auc_best <- max(sapply(1:3, function(k) { a <- auc1(sc[def, k][shared_keep], pos_def); max(a, 1 - a) }))
auc_prot <- auc1(xp[shared_keep, 1], pos_def)
cat(sprintf("\nDefinitive PC1 AUC (metadata-complete subset)  raw=%.3f  group-blind(unsupervised)=%.3f  group-protected(DEG covariates)=%.3f  (group-blind best-of-PC1-3=%.3f)\n",
            auc_raw, auc_corr, auc_prot, auc_best))

# ---------------------------------------------------------------------------
# 8. Figure 1 — definitive only: raw | batch+sex, coloured Control/Disease
# ---------------------------------------------------------------------------
bt <- function() theme_masld(base_size = 6) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_text(size = 6, face = "plain"), panel.grid = element_blank(),
        legend.position = "right", legend.title = element_text(size = 6, face = "plain"),
        legend.text = element_text(size = 6, face = "plain"), legend.key.size = unit(0.22, "cm"),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
        plot.margin = margin(1, 1, 1, 1), plot.title = element_text(size = 6, face = "plain"))

# Fibrosis F0–F4 uses the canonical fibrosis_stage_colors gradient (blue ramp);
# NAS uses pink_gradient (named gradient, publication_color_themes.R) — a single-hue
# light->dark magenta ramp so HIGH NAS = DARK, matching the fibrosis high=dark
# convention and distinct from the fibrosis blue. Controls drawn separately as grey
# open circles (never on these gradients).
fibf_def <- factor(ifelse(is.na(fib_def), "Unknown", paste0("F", fib_def)),
                   levels = c("F0", "F1", "F2", "F3", "F4", "Unknown"))
dt_raw  <- data.table(PC1 = xr[, 1],    PC2 = xr[, 2],    def_class = cls_def, nas = nas_def, fibrosis = fibf_def)[shared_keep]
dt_corr <- data.table(PC1 = sc[def, 1], PC2 = sc[def, 2], def_class = cls_def, nas = nas_def, fibrosis = fibf_def)[shared_keep]
dt_prot <- data.table(PC1 = xp[, 1],    PC2 = xp[, 2],    def_class = cls_def, nas = nas_def, fibrosis = fibf_def)[shared_keep]

pc12_axis_spearman_label <- function(dt, score_col) {
  d <- copy(dt)
  d[, def_class := droplevels(def_class)]
  if (!all(c("Control", "Disease") %in% as.character(d$def_class))) return("PC1/2 rho = NA")

  c_ctrl <- colMeans(d[def_class == "Control", .(PC1, PC2)])
  c_dis  <- colMeans(d[def_class == "Disease", .(PC1, PC2)])
  axis_vec <- c_dis - c_ctrl
  axis_norm <- sqrt(sum(axis_vec^2))
  if (!is.finite(axis_norm) || axis_norm <= 0) return("PC1/2 rho = NA")
  axis_vec <- axis_vec / axis_norm

  x <- as.matrix(d[, .(PC1, PC2)])
  x <- sweep(x, 2, c_ctrl, "-")
  axis_score <- as.vector(x %*% axis_vec)
  y <- d[[score_col]]
  ok <- is.finite(axis_score) & is.finite(y)
  n_ok <- sum(ok)
  if (n_ok < 3L || length(unique(y[ok])) < 2L) return("PC1/2 rho = NA")
  ct <- suppressWarnings(cor.test(axis_score[ok], y[ok], method = "spearman", exact = FALSE))
  p <- ct$p.value
  pstr <- if (is.na(p)) {
    "p = NA"
  } else if (p < 0.001) {
    "p < 0.001"
  } else {
    sprintf("p = %.3f", p)
  }
  sprintf("PC1/2 rho = %.2f (n=%d)\n%s", unname(ct$estimate), n_ok, pstr)
}

# by Disease (Control/Disease colour)
pcvd <- function(dt, p1, p2, ttl) {
  dt_local <- copy(dt)
  dt_local[, def_class := droplevels(def_class)]
  dt_local[, disease_status := as.numeric(def_class == "Disease")]
  label_str <- pc12_axis_spearman_label(dt_local, "disease_status")
  x_pos <- min(dt_local$PC1) + (max(dt_local$PC1) - min(dt_local$PC1)) * 0.05
  y_pos <- max(dt_local$PC2) - (max(dt_local$PC2) - min(dt_local$PC2)) * 0.05

  ggplot(dt_local, aes(PC1, PC2, colour = def_class)) +
    geom_point(size = 0.85, alpha = 0.5) +
    scale_colour_manual(values = c(Control = CTRL, Disease = DIS), name = NULL) +
    annotate("text", x = x_pos, y = y_pos, label = label_str,
             size = 6 / ggplot2::.pt, colour = "black", lineheight = 0.85,
             hjust = 0, vjust = 1) +
    guides(colour = guide_legend(override.aes = list(size = 1.8))) +
    labs(x = sprintf("PC1 (%.1f%%)", p1), y = sprintf("PC2 (%.1f%%)", p2)) + bt()
}

# by NAS — controls = grey "Control" (open); disease coloured by NAS (>0; no Unknown)
pnas <- function(dt, p1, p2, ttl) {
  dt_local <- copy(dt)
  dt_local[, def_class := droplevels(def_class)]
  dt_local[, nas_for_rho := nas]
  dt_local[def_class == "Control", nas_for_rho := 0]
  x_pos <- min(dt_local$PC1) + (max(dt_local$PC1) - min(dt_local$PC1)) * 0.05
  y_pos <- max(dt_local$PC2) - (max(dt_local$PC2) - min(dt_local$PC2)) * 0.05
  
  ctl <- dt_local[def_class == "Control"]; dis <- dt_local[def_class == "Disease"]
  plot_dt <- rbindlist(list(ctl, dis), use.names = TRUE, fill = TRUE)
  label_str <- pc12_axis_spearman_label(plot_dt, "nas_for_rho")
  ggplot() +
    geom_point(data = ctl, aes(PC1, PC2, shape = "Control"),
               colour = "#9E9E9E", size = 0.85, stroke = 0.3, alpha = 0.5) +
    geom_point(data = dis, aes(PC1, PC2, colour = nas),
               shape = 16, size = 0.85, alpha = 0.5) +
    annotate("text", x = x_pos, y = y_pos, label = label_str,
             size = 6 / ggplot2::.pt, colour = "black", lineheight = 0.85,
             hjust = 0, vjust = 1) +
    scale_colour_gradientn(name = "NAS", colours = pink_gradient, limits = c(1, 8)) +
    scale_shape_manual(name = NULL, values = c(Control = 1)) +
    guides(colour = guide_colourbar(order = 1, barwidth = 0.4, barheight = 2.4),
           shape  = guide_legend(order = 2, override.aes = list(size = 1.8, colour = "grey30"))) +
    labs(x = sprintf("PC1 (%.1f%%)", p1), y = sprintf("PC2 (%.1f%%)", p2)) + bt()
}
# by Fibrosis — controls = grey OPEN circles; disease coloured by stage on a
# continuous blue gradient BAR (matches the NAS colourbar convention); no Unknown
pfib <- function(dt, p1, p2, ttl) {
  dt_local <- copy(dt)
  dt_local[, def_class := droplevels(def_class)]
  dt_local[, fib_num := NA_integer_]
  dt_local[grepl("^F[0-4]$", as.character(fibrosis)),
           fib_num := as.integer(sub("F", "", as.character(fibrosis)))]
  dt_local[, fib_for_rho := fib_num]
  dt_local[def_class == "Control", fib_for_rho := 0L]
  x_pos <- min(dt_local$PC1) + (max(dt_local$PC1) - min(dt_local$PC1)) * 0.05
  y_pos <- max(dt_local$PC2) - (max(dt_local$PC2) - min(dt_local$PC2)) * 0.05
  
  ctl <- dt_local[def_class == "Control"]
  dis <- copy(dt_local[def_class == "Disease"])
  plot_dt <- rbindlist(list(ctl, dis), use.names = TRUE, fill = TRUE)
  label_str <- pc12_axis_spearman_label(plot_dt, "fib_for_rho")
  ggplot() +
    geom_point(data = ctl, aes(PC1, PC2, shape = "Control"),
               colour = "#9E9E9E", size = 0.85, stroke = 0.3, alpha = 0.5) +
    geom_point(data = dis, aes(PC1, PC2, colour = fib_num),
               shape = 16, size = 0.85, alpha = 0.5) +
    annotate("text", x = x_pos, y = y_pos, label = label_str,
             size = 6 / ggplot2::.pt, colour = "black", lineheight = 0.85,
             hjust = 0, vjust = 1) +
    scale_colour_gradientn(name = "Fibrosis", colours = blue_gradient, limits = c(0, 4)) +
    scale_shape_manual(name = NULL, values = c(Control = 1)) +
    guides(colour = guide_colourbar(order = 1, barwidth = 0.4, barheight = 2.4),
           shape  = guide_legend(order = 2, override.aes = list(size = 1.8, colour = "grey30"))) +
    labs(x = sprintf("PC1 (%.1f%%)", p1), y = sprintf("PC2 (%.1f%%)", p2)) + bt()
}

# 2 columns: Raw | corrected (~dataset+sex+group). NAS/fibrosis rows keep only
# controls + score-bearing disease dots (no-score / zero-score disease dropped);
# the former "scored only" 3rd column is now redundant and was removed.
cT2 <- "Corrected (~dataset+sex+group)"
row1 <- (pcvd(dt_raw,            pve_raw[1],  pve_raw[2],  sprintf("Raw — Disease (AUC=%.2f)", auc_raw)) |
         pcvd(dt_prot,           pve_prot[1], pve_prot[2], sprintf("%s — Disease (AUC=%.2f)", cT2, auc_prot))) +
        plot_layout(guides = "collect")
row2 <- (pnas(dt_raw,            pve_raw[1],  pve_raw[2],  "Raw — NAS") |
         pnas(dt_prot,           pve_prot[1], pve_prot[2], paste(cT2, "— NAS"))) +
        plot_layout(guides = "collect")
row3 <- (pfib(dt_raw,            pve_raw[1],  pve_raw[2],  "Raw — Fibrosis") |
         pfib(dt_prot,           pve_prot[1], pve_prot[2], paste(cT2, "— Fibrosis"))) +
        plot_layout(guides = "collect")
fig1 <- (row1 / row2 / row3)
# Page-width, proportional to the other supplementary figures (was 11.5 x 7.2 in).
ggsave(file.path(OUT, "figs3u_pca_definitive_control_vs_disease.pdf"), fig1,
       width = fig_full_width, height = fig_full_width * 0.9, device = cairo_pdf)
message(sprintf(
  "[figs3u 6-panel] Definitive control (NAS0/F0) vs definitive disease (NASH/F3-F4), n=%d (%d ctrl / %d dis). Rows: Disease/NAS/Fibrosis; columns: raw (AUC=%.2f) | corrected ~dataset+sex+group (AUC=%.2f).",
  nrow(dt_raw), sum(dt_raw$def_class == "Control"), sum(dt_raw$def_class == "Disease"), auc_raw, auc_prot))
cat("Wrote figs3u_pca_definitive_control_vs_disease.pdf\n")

# ---------------------------------------------------------------------------
# 8b. MAIN Fig 3 panel 4 — distilled fibrosis-gradient scatter on the definitive
#     PCA (control NAS0/F0 vs disease NAS>=5|F>=3), covariate-adjusted
#     (~dataset+sex+group; group-PROTECTED `dt_prot`). Disease points are
#     coloured by fibrosis stage (F0->F4 blue ramp) so severity resolves along
#     the disease axis, foreshadowing the stage cascade in panel 5; controls are
#     grey open circles. NOTE: this PCA is a SEPARATION VISUAL on the clean
#     extremes and does NOT define the DEG set -- the canonical bulk DEG set
#     remains the POOLED all-sample analysis (canonical_deg_results.csv). Title
#     is labelled honestly as covariate-adjusted (PC1 is boosted by protecting
#     the binary Control/Disease label; the extremes-vs-pooled concordance is the
#     05i sensitivity arm). Individual panel, house style.
# ---------------------------------------------------------------------------
# Top PC1 marginal-density strip = the BINARY axis (Control grey vs Disease
# magenta, the Liang disease colour) on the SAME group-protected PC1, so it reads
# as a distinct layer from the blue fibrosis-severity gradient in the scatter
# below (binary shift on top; within-disease severity structure beneath).
pc1_rng  <- range(dt_prot$PC1, na.rm = TRUE)
x_expand <- ggplot2::expansion(mult = 0.03)
dprot_d  <- copy(dt_prot); dprot_d[, def_class := droplevels(def_class)]
top_dens <- ggplot(dprot_d, aes(PC1, colour = def_class, fill = def_class)) +
  geom_density(alpha = 0.25, linewidth = 0.4) +
  scale_colour_manual(values = c(Control = CTRL, Disease = DIS), guide = "none") +
  scale_fill_manual(values = c(Control = CTRL, Disease = DIS), guide = "none") +
  scale_x_continuous(limits = pc1_rng, expand = x_expand) +
  labs(x = NULL, y = NULL) +
  theme_masld(base_size = 6) +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), panel.grid = element_blank(),
        plot.margin = margin(1, 1, 0, 1))
fib_scatter <- pfib(dt_prot, pve_prot[1], pve_prot[2], "") +   # title lifted to composite top
  scale_x_continuous(limits = pc1_rng, expand = x_expand)
fig_fibgrad <- ((top_dens / fib_scatter) +
  plot_layout(heights = c(1, 4), guides = "collect")) &
  theme(legend.position = "right")
ggsave(file.path(OUT, "fig3d_pca_fibrosis_gradient.pdf"), fig_fibgrad,
       width = 2.42, height = 2.23, device = cairo_pdf)
message(sprintf(
  "[fig3 panel4 PCA] Definitive PCA, covariate-adjusted (~dataset+sex+group, group-PROTECTED). n=%d (%d ctrl / %d dis); PC1 %.1f%%, PC2 %.1f%%; Disease-vs-Control PC1 AUC=%.2f. Disease coloured by fibrosis stage; controls grey. DEG set = pooled canonical (unchanged); extremes-vs-pooled concordance lives in the 05i sensitivity arm.",
  nrow(dt_prot), sum(dt_prot$def_class == "Control"), sum(dt_prot$def_class == "Disease"),
  pve_prot[1], pve_prot[2], auc_prot))
cat("Wrote fig3d_pca_fibrosis_gradient.pdf\n")

# ---------------------------------------------------------------------------
# 9. Figure 2 — compact 3-panel variant, side-by-side:
#      raw (Disease/Control) | corrected NAS | corrected fibrosis
#    Reads left->right: raw definitive PCA is dominated by per-cohort batch
#    blobs; after ~dataset+sex+group correction the NAS and fibrosis severity
#    gradients align along PC1. (Replaces the retired intermediate-projection
#    figure, removed 2026-06-23 at user request.)
# ---------------------------------------------------------------------------
fig3 <- (pcvd(dt_raw,            pve_raw[1],  pve_raw[2],  sprintf("Raw — Disease (AUC=%.2f)", auc_raw)) |
         pnas(dt_prot,           pve_prot[1], pve_prot[2], "Corrected — NAS") |
         pfib(dt_prot,           pve_prot[1], pve_prot[2], "Corrected — Fibrosis"))
ggsave(file.path(OUT, "figs3u_pca_definitive_raw_nas_fib.pdf"), fig3,
       width = fig_full_width, height = 2.7, device = cairo_pdf)
# Descriptive caption -> stdout (house style: keep long descriptions off the panel).
message(sprintf(
  "[figs3u 3-panel] Definitive control (NAS0/F0) vs disease (NASH/F3-F4), n=%d (%d ctrl / %d dis). Raw definitive PCA reflects per-cohort batch structure (PC1 Disease AUC=%.2f); after ~dataset+sex+group correction the NAS and fibrosis severity gradients resolve along PC1 (AUC=%.2f).",
  nrow(dt_raw), sum(dt_raw$def_class == "Control"), sum(dt_raw$def_class == "Disease"), auc_raw, auc_prot))
cat("Wrote figs3u_pca_definitive_raw_nas_fib.pdf\nDEFINITIVE_DONE\n")
