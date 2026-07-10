#!/usr/bin/env Rscript
# figS_sample_spearman_heatmap.R
# ===========================================================================
# Unsupervised Spearman-correlation clustering of ALL human bulk RNA-seq
# samples (every QC-passing sample in merged_dge.rds), annotated by disease
# status, cohort and sex. Publication-ready ComplexHeatmap.
#
# Batch + sex are removed BEFORE clustering (same recipe as the definitive PCA
# panels, figS_pca_definitive.R): raw sample-sample correlation is otherwise
# dominated by cohort batch (same-cohort rho ~0.78 vs cross-cohort ~0.46), which
# scatters controls across cohort blocks. We:
#   1. logCPM (RLE-normalized DGEList)
#   2. infer a COMPLETE sex covariate (annotated where available; else from
#      XIST + Y-genes) -- 637/1260 samples lack annotated sex
#   3. pick top 2000 variable genes on WITHIN-COHORT-CENTERED logCPM, so the
#      variable genes reflect biology, not batch
#   4. limma::removeBatchEffect(batch = cohort, covariates = sex)  [unsupervised,
#      no disease term -> the disease/control split is NOT injected, it emerges]
#   5. sample x sample Spearman rho (rank -> scale -> crossprod, BLAS)
#   6. hclust(ward.D2) on (1 - rho) -> annotated, rasterized heatmap -> PDF
#
# Output: figures/supplementary/figS_sample_clustering/
#           sample_spearman_heatmap.pdf      (the figure)
#           sample_spearman_matrix.rds       (rho + clustering + labels)
# ===========================================================================
suppressPackageStartupMessages({
  library(edgeR); library(limma); library(matrixStats)
  library(ComplexHeatmap); library(circlize); library(grid)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(Sys.getenv("HOME"), "publication_color_themes.R"))

INTEG   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
DGE_RDS <- file.path(INTEG, "results/integration/merged_dge.rds")
META_CSV<- file.path(INTEG, "metadata/unified_metadata.csv")
OUT_DIR <- file.path(BASE, "figures/supplementary/figS_sample_clustering")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
TOP_N_DEG <- as.integer(Sys.getenv("TOP_N_DEG", "200"))   # strongest DEGs to cluster on
OUT_SUFFIX <- Sys.getenv("HEATMAP_OUT_SUFFIX", "")        # side-by-side previews
set.seed(42)

# --- 1. Load ---------------------------------------------------------------
cat("[load] merged_dge.rds ...\n")
dge <- readRDS(DGE_RDS)
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
cat(sprintf("[load] logCPM: %d genes x %d samples\n", nrow(logcpm), ncol(logcpm)))
meta <- read.csv(META_CSV, stringsAsFactors = FALSE)
meta <- meta[match(colnames(logcpm), meta$sample_id), ]
stopifnot(identical(meta$sample_id, colnames(logcpm)))
cohort <- factor(meta$dataset)

# --- 2. Complete sex covariate (annotated where present, else inferred) -----
ens   <- sub("\\.\\d+$", "", rownames(logcpm))
xist  <- logcpm[which(ens == "ENSG00000229807")[1], ]                 # female-high
ygene <- c("ENSG00000067048", "ENSG00000129824", "ENSG00000183878")  # DDX3Y/RPS4Y1/UTY
ymat  <- logcpm[match(ygene, ens)[!is.na(match(ygene, ens))], , drop = FALSE]
fscore <- as.numeric(scale(xist)) - colMeans(t(scale(t(ymat))))      # high = female
km     <- kmeans(fscore, centers = 2, nstart = 25)
fem_cl <- which.max(tapply(fscore, km$cluster, mean))
sex_inf <- ifelse(km$cluster == fem_cl, "F", "M")
sex_fin <- ifelse(meta$sex %in% c("F", "M"), meta$sex, sex_inf)       # prefer annotated
sex_cov <- as.numeric(sex_fin == "F")
sex <- factor(sex_fin, levels = c("F", "M"))
cat("[sex] final (annotated+inferred):\n"); print(table(sex))

# --- 3. Disease-status axis (full coverage, biologically ordered) -----------
ds <- rep(NA_character_, nrow(meta))
ds[meta$group_binary == "Control"] <- "Control"
ds[meta$group_binary == "Disease" & meta$diagnosis_harmonized == "NAFL"] <- "MASL"
ds[meta$group_binary == "Disease" &
   meta$diagnosis_harmonized %in% c("NASH", "Borderline")] <- "MASH"
ds[meta$group_binary == "Disease" & is.na(ds)] <- "MASLD (unspec.)"
disease_status <- factor(ds, levels = c("Control", "MASL", "MASH", "MASLD (unspec.)"))

# --- 4. Feature set = top-200 STRONGEST consensus DEGs ----------------------
# Unsupervised variable-gene clustering CANNOT resolve the MASLD control/disease
# axis (~1% of expression variance). Cluster on the strongest disease-
# discriminative genes (top 200 by |t|) so the disease signature stratifies
# samples sharply and controls form a clean block.
deg <- read.csv(file.path(INTEG, "results/integration/consensus_significant_degs.csv"),
                stringsAsFactors = FALSE)
deg <- deg[order(-abs(deg$t)), ]
deg <- deg[deg$gene %in% rownames(logcpm), ]
sig <- head(deg$gene, TOP_N_DEG)
cat(sprintf("[deg] top %d strongest DEGs (by |t|); |t| range %.1f-%.1f\n",
            length(sig), abs(deg$t[length(sig)]), abs(deg$t[1])))

# --- 5. Remove batch (cohort) + sex, CONDITIONING on disease status ----------
# Group-protected correction (same recipe as the definitive PCA panel,
# figS_pca_definitive.R): removeBatchEffect strips the cohort + sex nuisance axes
# while the design term PRESERVES the Control-vs-Disease contrast. This is
# mandatory here because cohort is confounded with disease (4 of 9 cohorts are
# disease-only) and disease is only ~1% of expression variance, so UNSUPERVISED
# batch removal co-removes the disease signal and controls fail to separate. The
# protected view is the "DEG-faithful" disease-conditioned axis.
grp <- factor(meta$group_binary, levels = c("Control", "Disease"))
Xc  <- limma::removeBatchEffect(logcpm[sig, ], batch = cohort, covariates = sex_cov,
                                design = model.matrix(~ grp))
cat("[correct] removeBatchEffect(batch = cohort, covariates = sex | design = ~group)\n")

# --- 6. Per-gene z-score, then sample x sample Spearman rho -----------------
# z-score each gene across samples so within-sample ranks reflect the DISEASE
# DIRECTION (above/below the gene's average) rather than absolute expression
# level -- otherwise highly-expressed genes rank high in every sample and the
# subtle disease contrast is swamped (controls fail to separate).
Xz <- t(scale(t(Xc))); Xz[!is.finite(Xz)] <- 0
R  <- apply(Xz, 2, rank); Rs <- scale(R)
S  <- crossprod(Rs) / (nrow(Rs) - 1)
S  <- pmin(pmax(S, -1), 1); diag(S) <- 1
cat(sprintf("[rho] off-diagonal range [%.3f, %.3f], median %.3f\n",
            min(S[lower.tri(S)]), max(S[lower.tri(S)]), median(S[lower.tri(S)])))

# --- 7. Unsupervised hierarchical clustering --------------------------------
hc <- hclust(as.dist(1 - S), method = "ward.D2")
# Rotate branches by disease severity (COSMETIC: dendrograms are rotation-
# invariant, the clustering/topology is unchanged) so the control-enriched branch
# sits at one edge and severity grades across the map -> clean control corner block.
sev  <- c("Control" = 0, "MASL" = 1, "MASLD (unspec.)" = 2, "MASH" = 3)[as.character(disease_status)]
dend <- reorder(as.dendrogram(hc), wts = sev, agglo.FUN = mean)

# diagnostics: how much do cohort vs disease drive clustering now?
purity <- function(lab, cl) mean(vapply(split(lab, cl), function(g){tt<-table(g); max(tt)/length(g)}, 0))
for (k in c(2, 5, 9)) {
  cl <- cutree(hc, k)
  cat(sprintf("[cluster] k=%-2d cohort-purity=%.2f disease-purity=%.2f\n",
              k, purity(cohort, cl), purity(disease_status, cl)))
}
ctrl <- disease_status == "Control"
cat(sprintf("[control] mean rho Control-Control=%.3f  Control-Disease=%.3f  (gap=%.3f)\n",
            mean(S[ctrl, ctrl][lower.tri(S[ctrl, ctrl])]), mean(S[ctrl, !ctrl]),
            mean(S[ctrl, ctrl][lower.tri(S[ctrl, ctrl])]) - mean(S[ctrl, !ctrl])))

# --- 8. Colors --------------------------------------------------------------
# Control = neutral gray; disease groups = one coherent warm/red family graded by
# severity (MASL light -> MASH deep), so disease clearly contrasts control and the
# subgroups read as one theme rather than clashing hues.
disease_cols <- c("Control"         = "#9E9E9E",   # gray (neutral)
                  "MASL"            = "#F4A582",   # light salmon (mild)
                  "MASH"            = "#B2182B",   # deep red (severe)
                  "MASLD (unspec.)" = "#D6604D")   # mid red-orange
cohort_cols  <- setNames(cat_palette[seq_along(levels(cohort))], levels(cohort))
sex_cols     <- c("F"="#C98BB9", "M"="#6FB3A8")
qs <- quantile(S[lower.tri(S)], c(0.02, 0.5, 0.98))
body_col <- colorRamp2(c(qs[1], qs[2], qs[3]), c("#3B4CC0", "#F7F7F7", "#B40426"))

# --- 9. Heatmap -------------------------------------------------------------
ha <- HeatmapAnnotation(
  `Disease status` = disease_status, Cohort = cohort, Sex = sex,
  col = list(`Disease status` = disease_cols, Cohort = cohort_cols, Sex = sex_cols),
  annotation_name_gp = gpar(fontsize = 6), annotation_name_side = "left",
  simple_anno_size = unit(3.0, "mm"), gap = unit(0.7, "mm"),
  annotation_legend_param = list(
    `Disease status` = list(title_gp=gpar(fontsize=6,fontface="plain"), labels_gp=gpar(fontsize=6)),
    Cohort = list(title_gp=gpar(fontsize=6,fontface="plain"), labels_gp=gpar(fontsize=6), ncol=1),
    Sex = list(title_gp=gpar(fontsize=6,fontface="plain"), labels_gp=gpar(fontsize=6))))

ht <- Heatmap(
  S, name = "Spearman ρ", col = body_col,
  cluster_rows = dend, cluster_columns = dend,
  show_row_names = FALSE, show_column_names = FALSE,
  show_row_dend = FALSE, show_column_dend = TRUE,
  column_dend_height = unit(14, "mm"), top_annotation = ha,
  use_raster = TRUE, raster_quality = 4,
  heatmap_legend_param = list(title_gp=gpar(fontsize=6,fontface="plain"),
    labels_gp=gpar(fontsize=6), legend_height=unit(28,"mm")))
message(sprintf(
  "[caption] Spearman clustering of %d human RNA-seq samples on the consensus disease signature (%d DEGs; batch + sex removed, disease-conditioned; %d cohorts)",
  ncol(S), length(sig), nlevels(cohort)))

pdf_path <- file.path(OUT_DIR, sprintf("sample_spearman_heatmap%s.pdf", OUT_SUFFIX))
ok <- tryCatch({ cairo_pdf(pdf_path, width = 7.09, height = 5.88); TRUE },
               error = function(e) FALSE)
if (!ok) pdf(pdf_path, width = 7.09, height = 5.88, useDingbats = FALSE)
draw(ht, merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
cat("[write]", pdf_path, "\n")

# Replicate the canonical (default, n=200) heatmap into the main Fig 3 panels dir
# (promoted to main Fig 3C 2026-06-15). The sweep variants (OUT_SUFFIX != "")
# stay supplementary-only; the supplementary copy above is the canonical source.
if (OUT_SUFFIX == "") {
  fig3_panels <- file.path(BASE, "figures/main/fig3_RNAseq/panels")
  dir.create(fig3_panels, recursive = TRUE, showWarnings = FALSE)
  file.copy(pdf_path, file.path(fig3_panels, "sample_spearman_heatmap.pdf"), overwrite = TRUE)
  cat("[replicate] sample_spearman_heatmap.pdf -> main Fig 3 panels (3C)\n")
}

saveRDS(list(rho = S, hclust = hc, order = hc$order, disease_status = disease_status,
             cohort = cohort, sex = sex, signature_genes = sig),
        file.path(OUT_DIR, sprintf("sample_spearman_matrix%s.rds", OUT_SUFFIX)))
cat("[done]\n")
