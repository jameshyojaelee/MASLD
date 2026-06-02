#!/usr/bin/env Rscript
#SBATCH --partition=bigmem
#SBATCH --mem=240G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --job-name=net_281d_f2_deconv
#SBATCH --output=logs/net_281d_f2_deconv_%j.out
#SBATCH --error=logs/net_281d_f2_deconv_%j.err
# ===========================================================================
# Script 281d: Deconvolution-adjusted F2 delta-r for composition-driven flag.
# ---------------------------------------------------------------------------
# For every existing D-F2 edge (edges_d_f2_loco.csv), recompute per-stage
# Pearson correlation on residualized log-CPM (each gene residualized
# against the 5 major liver cell-type proportions). Flag an edge as
# `composition_driven` if the residual |delta_max| is less than 50% of
# the raw |delta_max| -- i.e., more than half the apparent stage-dependent
# correlation shift disappears once cell-composition is regressed out.
#
# Inputs:
#   RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_dge.rds
#   RNA-seq/results/celltype_attribution/persample_celltype_proportions.csv
#   RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv
#   RNA-seq/results/network/edges_d_f2_loco.csv  (source pairs to re-score)
#   RNA-seq/results/network/network_nodes.csv    (ENSG <-> symbol map)
#
# Output:
#   RNA-seq/results/network/edges_d_f2_deconv_flag.csv
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

t0 <- Sys.time()

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

NETDIR    <- file.path(PROJ, "RNA-seq/results/network")
DGE_PATH  <- file.path(PROJ,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_dge.rds")
PROP_PATH <- file.path(PROJ,
  "RNA-seq/results/celltype_attribution/persample_celltype_proportions.csv")
META_PATH <- file.path(PROJ,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
D_F2_PATH <- file.path(NETDIR, "edges_d_f2_loco.csv")
NODE_PATH <- file.path(NETDIR, "network_nodes.csv")
OUT       <- file.path(NETDIR, "edges_d_f2_deconv_flag.csv")

COMPOSITION_RATIO_THRESHOLD <- 0.5  # resid/raw < 0.5 -> flagged

# The 5 cell-type proportions we residualize against. Using the major
# liver cell types (parenchymal + non-parenchymal) that are reliably
# deconvolved by MuSiC.
COVARIATES <- c("Hepatocytes", "Macrophages", "Endothelial cells",
                "Fibroblasts", "Cholangiocytes")

message("[281d] Deconvolution-adjusted D-F2 starting at ", format(t0))

stopifnot(file.exists(DGE_PATH), file.exists(PROP_PATH),
          file.exists(META_PATH), file.exists(D_F2_PATH),
          file.exists(NODE_PATH))

nodes <- fread(NODE_PATH)
strip_version <- function(x) sub("\\..*$", "", x)
sym_map <- setNames(nodes$human_symbol, strip_version(nodes$ensembl_id))

message("[", Sys.time(), "] Loading DGE + metadata")
dge <- readRDS(DGE_PATH)
unified_meta <- fread(META_PATH)
prop <- fread(PROP_PATH)

# Align samples -> DGE column order
meta <- data.table(sample_id = colnames(dge))
meta <- merge(meta, unified_meta[, .(sample_id, dataset, fibrosis_stage, sex)],
              by = "sample_id", all.x = TRUE, sort = FALSE)
setnames(meta, c("dataset", "sex"), c("cohort", "inferred_sex"))
meta <- merge(meta,
              prop[, c("sample_id", COVARIATES), with = FALSE],
              by = "sample_id", all.x = TRUE, sort = FALSE)
stopifnot(all(meta$sample_id == colnames(dge)))

# Stage normalization (same as 281)
norm_stage <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("0", "F0", "1", "F1")] <- "F01"
  x[x %in% c("2", "F2")]            <- "F2"
  x[x %in% c("3", "F3", "4", "F4")] <- "F34"
  x[!(x %in% c("F01", "F2", "F34"))] <- NA_character_
  x
}
meta[, stage := norm_stage(fibrosis_stage)]

# Need: stage AND all covariates present
covariate_ok <- complete.cases(meta[, ..COVARIATES])
keep_s <- !is.na(meta$stage) & covariate_ok
message("  Samples with stage + all covariates: ", sum(keep_s),
        " / ", nrow(meta))
meta <- meta[keep_s]
dge <- dge[, keep_s, keep.lib.sizes = FALSE]

# log-CPM on filtered DGE
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
rn_sym <- sym_map[strip_version(rownames(dge))]
rownames(logcpm) <- rn_sym
logcpm <- logcpm[!is.na(rownames(logcpm)) & !duplicated(rownames(logcpm)), ,
                 drop = FALSE]

# Residualize each row of logcpm against covariates (pooled; same covariate
# model for all stages).
message("[", Sys.time(), "] Residualizing log-CPM against cell-type covariates")
X <- as.matrix(meta[, ..COVARIATES])
X <- cbind(1, X)   # intercept
# Y = logcpm^T (samples x genes)
Y <- t(logcpm)
qrX <- qr(X)
# residuals matrix (samples x genes)
resid_mat <- qr.resid(qrX, Y)
resid_logcpm <- t(resid_mat)
rownames(resid_logcpm) <- rownames(logcpm)
rm(Y, resid_mat); invisible(gc(verbose = FALSE))

# Edge set to re-evaluate
message("[", Sys.time(), "] Loading D-F2 edges to re-evaluate")
d_f2 <- fread(D_F2_PATH)
message("  D-F2 edges loaded: ", nrow(d_f2))

d_f2 <- d_f2[gene_a %in% rownames(resid_logcpm) &
             gene_b %in% rownames(resid_logcpm)]
message("  After restricting to genes in residualized matrix: ", nrow(d_f2))

# Per-stage sample indices
i01 <- which(meta$stage == "F01")
i2  <- which(meta$stage == "F2")
i34 <- which(meta$stage == "F34")
message("  Stage n: F01=", length(i01), " F2=", length(i2), " F34=", length(i34))

# Vectorized row-wise Pearson for a given column subset
row_cor <- function(A, B, idx) {
  a <- A[, idx, drop = FALSE]
  b <- B[, idx, drop = FALSE]
  am <- rowMeans(a); bm <- rowMeans(b)
  a <- a - am; b <- b - bm
  num <- rowSums(a * b)
  den <- sqrt(rowSums(a * a) * rowSums(b * b))
  ifelse(den > 0, num / den, NA_real_)
}

message("[", Sys.time(), "] Recomputing per-stage Pearson on residualized matrix")
A <- resid_logcpm[d_f2$gene_a, , drop = FALSE]
B <- resid_logcpm[d_f2$gene_b, , drop = FALSE]

r_F01_resid <- row_cor(A, B, i01)
r_F2_resid  <- row_cor(A, B, i2)
r_F34_resid <- row_cor(A, B, i34)
rm(A, B); invisible(gc(verbose = FALSE))

d_f2[, r_F01_resid := r_F01_resid]
d_f2[, r_F2_resid  := r_F2_resid]
d_f2[, r_F34_resid := r_F34_resid]
d_f2[, delta_F2_vs_F01_resid  := r_F2_resid  - r_F01_resid]
d_f2[, delta_F34_vs_F2_resid  := r_F34_resid - r_F2_resid]
d_f2[, delta_F34_vs_F01_resid := r_F34_resid - r_F01_resid]
d_f2[, delta_max_resid := pmax(abs(delta_F2_vs_F01_resid),
                               abs(delta_F34_vs_F01_resid),
                               abs(delta_F34_vs_F2_resid),
                               na.rm = TRUE)]

# Concordance ratio: resid vs raw delta_max. Small ratio => composition-driven.
d_f2[, deconv_concordance := ifelse(delta_max > 0,
                                    delta_max_resid / delta_max,
                                    NA_real_)]
d_f2[, composition_driven := !is.na(deconv_concordance) &
                             deconv_concordance < COMPOSITION_RATIO_THRESHOLD]

out <- d_f2[, .(gene_a, gene_b,
                delta_max,
                delta_max_resid,
                deconv_concordance,
                composition_driven)]

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
fwrite(out, OUT)
message("[", Sys.time(), "] Wrote ", nrow(out),
        " deconv-flag rows -> ", OUT)

# Diagnostics
message("  composition_driven flagged: ",
        sum(out$composition_driven, na.rm = TRUE),
        " / ", nrow(out),
        sprintf("  (%.1f%%)",
                100 * mean(out$composition_driven, na.rm = TRUE)))
message("  deconv_concordance quantiles:")
print(round(quantile(out$deconv_concordance,
                     c(0.1, 0.25, 0.5, 0.75, 0.9),
                     na.rm = TRUE), 3))

message("[281d] Elapsed: ",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
        " min")
