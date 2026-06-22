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
#   figs3u_pca_definitive_control_vs_disease.pdf  (raw | batch+sex; definitive only)
#   figs3u_pca_definitive_projection.pdf          (definitive axes + intermediates projected)
#   data/pca_definitive_counts.csv                (per-cohort class counts)
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
pos_def <- cls_def == "Disease"
auc_raw  <- auc1(xr[, 1], pos_def)
auc_corr <- auc1(sc[def, 1], pos_def)
auc_best <- max(sapply(1:3, function(k) { a <- auc1(sc[def, k], pos_def); max(a, 1 - a) }))
auc_prot <- auc1(xp[, 1], pos_def)
cat(sprintf("\nDefinitive PC1 AUC  raw=%.3f  group-blind(unsupervised)=%.3f  group-protected(DEG covariates)=%.3f  (group-blind best-of-PC1-3=%.3f)\n",
            auc_raw, auc_corr, auc_prot, auc_best))

# ---------------------------------------------------------------------------
# 8. Figure 1 — definitive only: raw | batch+sex, coloured Control/Disease
# ---------------------------------------------------------------------------
bt <- function() theme_masld(base_size = 7) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_text(size = 6), panel.grid = element_blank(),
        legend.position = "right", legend.title = element_text(size = 6, face = "bold"),
        legend.text = element_text(size = 5.5), legend.key.size = unit(0.22, "cm"),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
        plot.margin = margin(1, 1, 1, 1), plot.title = element_text(size = 6.5, face = "bold"))

# Fibrosis F0–F4 uses the canonical fibrosis_stage_colors gradient (blue ramp);
# NAS uses pink_gradient (named gradient, publication_color_themes.R) — a single-hue
# light->dark magenta ramp so HIGH NAS = DARK, matching the fibrosis high=dark
# convention and distinct from the fibrosis blue. Controls drawn separately as grey
# open circles (never on these gradients).
nas_def  <- as.numeric(samp$nas_score)[def]
fib_def  <- as.numeric(samp$fibrosis_stage)[def]
fibf_def <- factor(ifelse(is.na(fib_def), "Unknown", paste0("F", fib_def)),
                   levels = c("F0", "F1", "F2", "F3", "F4", "Unknown"))
dt_raw  <- data.table(PC1 = xr[, 1],    PC2 = xr[, 2],    def_class = cls_def, nas = nas_def, fibrosis = fibf_def)
dt_corr <- data.table(PC1 = sc[def, 1], PC2 = sc[def, 2], def_class = cls_def, nas = nas_def, fibrosis = fibf_def)
dt_prot <- data.table(PC1 = xp[, 1],    PC2 = xp[, 2],    def_class = cls_def, nas = nas_def, fibrosis = fibf_def)
# NAS / fibrosis panels carry NO "Unknown" dots (request 2026-06-18). Controls
# are shown as their own "Control" category (recoloured in pnas/pfib, never
# Unknown); DISEASE dots are kept only if they actually carry the score (NAS > 0 /
# a fibrosis stage present) and are coloured by it.
nas_keep <- cls_def == "Control" | (cls_def == "Disease" & !is.na(nas_def) & nas_def > 0)
fib_keep <- cls_def == "Control" | (cls_def == "Disease" & !is.na(fib_def))

# by Disease (Control/Disease colour)
pcvd <- function(dt, p1, p2, ttl) {
  ggplot(dt, aes(PC1, PC2, colour = def_class)) +
    geom_point(size = 0.85, alpha = 0.8) +
    scale_colour_manual(values = c(Control = CTRL, Disease = DIS), name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 1.8))) +
    labs(x = sprintf("PC1 (%.1f%%)", p1), y = sprintf("PC2 (%.1f%%)", p2), title = ttl) + bt()
}
# by NAS — controls = grey "Control" (open); disease coloured by NAS (>0; no Unknown)
pnas <- function(dt, p1, p2, ttl) {
  ctl <- dt[def_class == "Control"]; dis <- dt[def_class == "Disease"]
  ggplot() +
    geom_point(data = ctl, aes(PC1, PC2, shape = "Control"),
               colour = "#9E9E9E", size = 0.85, stroke = 0.3, alpha = 0.85) +
    geom_point(data = dis, aes(PC1, PC2, colour = nas),
               shape = 16, size = 0.85, alpha = 0.85) +
    scale_colour_gradientn(name = "NAS", colours = pink_gradient, limits = c(1, 8)) +
    scale_shape_manual(name = NULL, values = c(Control = 1)) +
    guides(colour = guide_colourbar(order = 1, barwidth = 0.4, barheight = 2.4),
           shape  = guide_legend(order = 2, override.aes = list(size = 1.8, colour = "grey30"))) +
    labs(x = sprintf("PC1 (%.1f%%)", p1), y = sprintf("PC2 (%.1f%%)", p2), title = ttl) + bt()
}
# by Fibrosis — controls = grey OPEN circles; disease filled, coloured by stage; no Unknown
pfib <- function(dt, p1, p2, ttl) {
  ctl <- dt[def_class == "Control"]
  dis <- copy(dt[def_class == "Disease"])
  dis[, fib_disp := factor(as.character(fibrosis), levels = c("F0", "F1", "F2", "F3", "F4"))]
  ggplot() +
    geom_point(data = ctl, aes(PC1, PC2, shape = "Control"),
               colour = "#9E9E9E", size = 0.85, stroke = 0.3, alpha = 0.85) +
    geom_point(data = dis, aes(PC1, PC2, colour = fib_disp),
               shape = 16, size = 0.85, alpha = 0.85) +
    scale_colour_manual(name = "Fibrosis", values = fibrosis_stage_colors,
                        limits = c("F0", "F1", "F2", "F3", "F4"), drop = FALSE) +
    scale_shape_manual(name = NULL, values = c(Control = 1)) +
    guides(colour = guide_legend(order = 1, override.aes = list(size = 1.8)),
           shape  = guide_legend(order = 2, override.aes = list(size = 1.8, colour = "grey30"))) +
    labs(x = sprintf("PC1 (%.1f%%)", p1), y = sprintf("PC2 (%.1f%%)", p2), title = ttl) + bt()
}

# 2 columns: Raw | corrected (~dataset+sex+group). NAS/fibrosis rows keep only
# controls + score-bearing disease dots (no-score / zero-score disease dropped);
# the former "scored only" 3rd column is now redundant and was removed.
cT2 <- "Corrected (~dataset+sex+group)"
row1 <- (pcvd(dt_raw,            pve_raw[1],  pve_raw[2],  sprintf("Raw — Disease (AUC=%.2f)", auc_raw)) |
         pcvd(dt_prot,           pve_prot[1], pve_prot[2], sprintf("%s — Disease (AUC=%.2f)", cT2, auc_prot))) +
        plot_layout(guides = "collect")
row2 <- (pnas(dt_raw[nas_keep],  pve_raw[1],  pve_raw[2],  "Raw — NAS") |
         pnas(dt_prot[nas_keep], pve_prot[1], pve_prot[2], paste(cT2, "— NAS"))) +
        plot_layout(guides = "collect")
row3 <- (pfib(dt_raw[fib_keep],  pve_raw[1],  pve_raw[2],  "Raw — Fibrosis") |
         pfib(dt_prot[fib_keep], pve_prot[1], pve_prot[2], paste(cT2, "— Fibrosis"))) +
        plot_layout(guides = "collect")
fig1 <- (row1 / row2 / row3) +
  plot_annotation(
    title = sprintf("Definitive control (NAS0/F0) vs definitive disease (NASH/F3-F4), n=%d (%d ctrl / %d dis)",
                    length(def), sum(cls_def == "Control"), sum(cls_def == "Disease")),
    theme = theme(plot.title = element_text(size = 8.5, face = "bold")))
# Page-width, proportional to the other supplementary figures (was 11.5 x 7.2 in).
ggsave(file.path(OUT, "figs3u_pca_definitive_control_vs_disease.pdf"), fig1,
       width = fig_full_width, height = fig_full_width * 0.9, device = cairo_pdf)
cat("Wrote figs3u_pca_definitive_control_vs_disease.pdf\n")

# ---------------------------------------------------------------------------
# 9. Figure 2 — projection: definitive axes + intermediates projected on
# ---------------------------------------------------------------------------
proj <- data.table(PC1 = sc[, 1], PC2 = sc[, 2], def_class = samp$def_class, sev = samp$sev)
anc <- proj[def_class != "Intermediate"]
imd <- proj[def_class == "Intermediate"]
fig2 <- ggplot() +
  geom_point(data = imd, aes(PC1, PC2, fill = sev),
             shape = 21, colour = "grey55", stroke = 0.1, size = 0.9, alpha = 0.55) +
  geom_point(data = anc, aes(PC1, PC2, colour = def_class), size = 1.0, alpha = 0.9) +
  scale_fill_gradientn(name = "Intermediate\nseverity\n(NAS/fib)",
                       colours = c("#FEE08B", "#FDAE61", "#D73027"),
                       limits = c(0, 1), na.value = "grey80") +
  scale_colour_manual(name = "Definitive",
                      values = c(Control = CTRL, Disease = DIS)) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 2)),
         fill = guide_colourbar(order = 2, barwidth = 0.4, barheight = 3)) +
  labs(x = sprintf("PC1 (%.1f%%)", pve[1]), y = sprintf("PC2 (%.1f%%)", pve[2]),
       title = "Intermediate samples projected onto the definitive control–disease axis") +
  theme_masld(base_size = 7) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
        legend.position = "right", plot.title = element_text(size = 8.5, face = "bold"))
ggsave(file.path(OUT, "figs3u_pca_definitive_projection.pdf"), fig2,
       width = 5.6, height = 4, device = cairo_pdf)
cat("Wrote figs3u_pca_definitive_projection.pdf\nDEFINITIVE_DONE\n")
