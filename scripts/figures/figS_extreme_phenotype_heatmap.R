#!/usr/bin/env Rscript
# figS_extreme_phenotype_heatmap.R
# ===========================================================================
# PAIRED heatmap to figS_sample_spearman_heatmap.R, restricted to the two
# histologic EXTREMES (definite-disease vs strict-control; 05i). Clustering the
# extreme-phenotype samples on the strongest extreme-contrast DEGs yields a
# near-BIMODAL separation -- the visual answer to "are the controls truly clean,
# or is the consensus signal diluted by ambiguous mild cases?"
#
# Same recipe as the all-sample heatmap:
#   1. logCPM (RLE DGEList), subset to the 419 extreme-phenotype samples (05i meta)
#   2. feature set = top-N strongest extreme DEGs (by |t|, 05i table)
#   3. removeBatchEffect(batch = cohort, covariates = sex, design = ~phenotype)
#      -- strips cohort + sex; the binary phenotype mean offset is protected
#   4. per-gene z-score -> sample x sample Spearman rho (rank -> scale -> crossprod)
#   5. hclust(ward.D2) on (1 - rho), branches rotated so StrictControl sits at one edge
#
# Output: figures/supplementary/figS_sample_clustering/
#           extreme_phenotype_heatmap[.SUFFIX].pdf
#           extreme_phenotype_matrix[.SUFFIX].rds
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
SDIR    <- file.path(INTEG, "results/integration/sensitivity")
DE_CSV  <- file.path(SDIR, "extreme_phenotype_definite_vs_strict.csv")
META_CSV<- file.path(SDIR, "extreme_phenotype_definite_vs_strict_meta.csv")
OUT_DIR <- file.path(BASE, "figures/supplementary/figS_sample_clustering")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
TOP_N_DEG  <- as.integer(Sys.getenv("TOP_N_DEG", "200"))
OUT_SUFFIX <- Sys.getenv("HEATMAP_OUT_SUFFIX", "")
set.seed(42)

# --- 1. Load + subset to the extreme-phenotype samples ----------------------
dge    <- readRDS(DGE_RDS)
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
emeta  <- read.csv(META_CSV, stringsAsFactors = FALSE)
keep   <- colnames(logcpm) %in% emeta$sample_id
logcpm <- logcpm[, keep]
emeta  <- emeta[match(colnames(logcpm), emeta$sample_id), ]
stopifnot(identical(emeta$sample_id, colnames(logcpm)))
cat(sprintf("[load] %d genes x %d extreme-phenotype samples\n", nrow(logcpm), ncol(logcpm)))

cohort    <- factor(emeta$dataset)
sex       <- factor(emeta$inferred_sex, levels = c("F", "M"))
sex_cov   <- as.numeric(sex == "F")
phenotype <- factor(emeta$phenotype, levels = c("StrictControl", "Definite"),
                    labels = c("Strict control", "Advanced disease"))
nas       <- suppressWarnings(as.numeric(emeta$nas_score))
fib       <- suppressWarnings(as.numeric(emeta$fibrosis_stage))
cat("[arms]\n"); print(table(phenotype, cohort))

# --- 2. Feature set = top-N strongest extreme DEGs (by |t|) ------------------
deg <- read.csv(DE_CSV, stringsAsFactors = FALSE)
deg <- deg[order(-abs(deg$t)), ]
deg <- deg[deg$gene %in% rownames(logcpm), ]
sig <- head(deg$gene, TOP_N_DEG)
cat(sprintf("[deg] top %d strongest extreme DEGs (by |t|); |t| range %.1f-%.1f\n",
            length(sig), abs(deg$t[length(sig)]), abs(deg$t[1])))

# --- 3. Remove cohort + sex, conditioning on phenotype ----------------------
Xc <- limma::removeBatchEffect(logcpm[sig, ], batch = cohort, covariates = sex_cov,
                               design = model.matrix(~ phenotype))
cat("[correct] removeBatchEffect(batch = cohort, covariates = sex | design = ~phenotype)\n")

# --- 4. Per-gene z-score -> sample x sample Spearman rho --------------------
Xz <- t(scale(t(Xc))); Xz[!is.finite(Xz)] <- 0
R  <- apply(Xz, 2, rank); Rs <- scale(R)
S  <- crossprod(Rs) / (nrow(Rs) - 1)
S  <- pmin(pmax(S, -1), 1); diag(S) <- 1

# --- 5. Cluster; rotate so StrictControl sits at one edge -------------------
hc   <- hclust(as.dist(1 - S), method = "ward.D2")
sev  <- c("Strict control" = 0, "Advanced disease" = 1)[as.character(phenotype)]
dend <- reorder(as.dendrogram(hc), wts = sev, agglo.FUN = mean)

ctrl <- phenotype == "Strict control"
cat(sprintf("[control] mean rho Strict-Strict=%.3f  Strict-Advanced=%.3f  (gap=%.3f)\n",
            mean(S[ctrl, ctrl][lower.tri(S[ctrl, ctrl])]), mean(S[ctrl, !ctrl]),
            mean(S[ctrl, ctrl][lower.tri(S[ctrl, ctrl])]) - mean(S[ctrl, !ctrl])))
for (k in c(2, 4)) {
  cl <- cutree(hc, k)
  pur <- function(lab) mean(vapply(split(lab, cl), function(g){tt<-table(g); max(tt)/length(g)}, 0))
  cat(sprintf("[cluster] k=%d cohort-purity=%.2f phenotype-purity=%.2f\n",
              k, pur(cohort), pur(phenotype)))
}

# --- 6. Colors --------------------------------------------------------------
pheno_cols  <- c("Strict control" = "#9E9E9E", "Advanced disease" = "#B2182B")
cohort_cols <- setNames(cat_palette[seq_along(levels(cohort))], levels(cohort))
sex_cols    <- c("F" = "#C98BB9", "M" = "#6FB3A8")
nas_col     <- colorRamp2(c(0, 4, 8), c("#FFF7EC", "#FC8D59", "#7F0000"))
fib_col     <- colorRamp2(c(0, 2, 4), c("#F7FCFD", "#8C96C6", "#4D004B"))
qs       <- quantile(S[lower.tri(S)], c(0.02, 0.5, 0.98))
body_col <- colorRamp2(c(qs[1], qs[2], qs[3]), c("#3B4CC0", "#F7F7F7", "#B40426"))

# --- 7. Heatmap -------------------------------------------------------------
ha <- HeatmapAnnotation(
  `Phenotype` = phenotype, NAS = nas, Fibrosis = fib, Cohort = cohort, Sex = sex,
  col = list(`Phenotype` = pheno_cols, NAS = nas_col, Fibrosis = fib_col,
             Cohort = cohort_cols, Sex = sex_cols),
  na_col = "grey92",
  annotation_name_gp = gpar(fontsize = 8), annotation_name_side = "left",
  simple_anno_size = unit(3.0, "mm"), gap = unit(0.7, "mm"),
  annotation_legend_param = list(
    `Phenotype` = list(title_gp = gpar(fontsize=8, fontface="bold"), labels_gp = gpar(fontsize=7)),
    Cohort = list(title_gp = gpar(fontsize=8, fontface="bold"), labels_gp = gpar(fontsize=7), ncol=1),
    Sex = list(title_gp = gpar(fontsize=8, fontface="bold"), labels_gp = gpar(fontsize=7))))

ht <- Heatmap(
  S, name = "Spearman ρ", col = body_col,
  cluster_rows = dend, cluster_columns = dend,
  show_row_names = FALSE, show_column_names = FALSE,
  show_row_dend = FALSE, show_column_dend = TRUE,
  column_dend_height = unit(14, "mm"), top_annotation = ha,
  use_raster = TRUE, raster_quality = 4,
  column_title = sprintf(
    "Sample clustering: advanced disease (NAS≥5 | F≥3) vs strict controls\n(%d samples; top %d DEGs; batch + sex removed; %d cohorts)",
    ncol(S), length(sig), nlevels(cohort)),
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  heatmap_legend_param = list(title_gp = gpar(fontsize=8, fontface="bold"),
    labels_gp = gpar(fontsize=7), legend_height = unit(28, "mm")))

pdf_path <- file.path(OUT_DIR, sprintf("extreme_phenotype_heatmap%s.pdf", OUT_SUFFIX))
ok <- tryCatch({ cairo_pdf(pdf_path, width = 9.4, height = 7.8); TRUE },
               error = function(e) FALSE)
if (!ok) pdf(pdf_path, width = 9.4, height = 7.8, useDingbats = FALSE)
draw(ht, merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
cat("[write]", pdf_path, "\n")

saveRDS(list(rho = S, hclust = hc, order = hc$order, phenotype = phenotype,
             cohort = cohort, sex = sex, nas = nas, fib = fib, signature_genes = sig),
        file.path(OUT_DIR, sprintf("extreme_phenotype_matrix%s.rds", OUT_SUFFIX)))
cat("[done]\n")
