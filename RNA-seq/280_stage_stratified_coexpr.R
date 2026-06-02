#!/usr/bin/env Rscript
#SBATCH --partition=cpu
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --job-name=net_280_stage_coexpr
#SBATCH --output=logs/net_280_stage_coexpr_%j.out
#SBATCH --error=logs/net_280_stage_coexpr_%j.err
# ===========================================================================
# Script 280: Stage-stratified co-expression (F01 / F2 / F34)
# ---------------------------------------------------------------------------
# Computes per-fibrosis-stage Pearson correlation matrices on the MASLD
# atlas logCPM for downstream D-F2 edge construction (Architecture C, Phase 2).
# Housekeeping genes are excluded prior to correlation.
# Outputs upper-triangle edge tables per stage.
# Environment: micromamba activate rnaseq
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

t0 <- Sys.time()

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

NODE_PATH <- file.path(BASE, "RNA-seq/results/network/network_nodes.csv")
DGE_PATH  <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_dge.rds")

OUTDIR <- file.path(BASE, "RNA-seq/results/network/stage")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

MIN_STAGE_N <- 50

HOUSEKEEPING_PATTERNS <- c(
  "^MT-", "^RPL", "^RPS", "^RPSA", "^EEF", "^EIF", "^HSP",
  "^TUB", "^YWHA", "^SNRP", "^HNRNP", "^COX", "^NDUF",
  "^UQCR", "^ATP5", "^S100A", "^PSM[A-F]", "^MRPL", "^MRPS",
  "^UBE"
)
HOUSEKEEPING_EXACT <- c(
  "ACTB", "ACTG1", "GAPDH", "B2M", "UBA52", "UBB", "UBC", "FAU"
)

is_housekeeping <- function(sym) {
  hit <- rep(FALSE, length(sym))
  for (p in HOUSEKEEPING_PATTERNS) hit <- hit | grepl(p, sym)
  hit | sym %in% HOUSEKEEPING_EXACT
}

message("[", Sys.time(), "] Loading node set + merged DGEList")
nodes <- fread(NODE_PATH)
dge   <- readRDS(DGE_PATH)

strip_version <- function(x) sub("\\..*$", "", x)

sym_map <- setNames(nodes$human_symbol, strip_version(nodes$ensembl_id))

rn_stripped <- strip_version(rownames(dge))
keep_row <- rn_stripped %in% names(sym_map)
dge <- dge[keep_row, , keep.lib.sizes = FALSE]
rn_stripped <- strip_version(rownames(dge))
row_sym <- sym_map[rn_stripped]

message("[", Sys.time(), "] Computing logCPM for ",
        nrow(dge), " genes x ", ncol(dge), " samples")
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
rownames(logcpm) <- row_sym

dup <- duplicated(rownames(logcpm))
if (any(dup)) {
  message("Collapsing ", sum(dup), " duplicate symbols (keeping highest mean)")
  ord <- order(-rowMeans(logcpm))
  logcpm <- logcpm[ord, , drop = FALSE]
  logcpm <- logcpm[!duplicated(rownames(logcpm)), , drop = FALSE]
}

hk <- is_housekeeping(rownames(logcpm))
message("Excluding ", sum(hk), " housekeeping genes; retaining ", sum(!hk))
logcpm <- logcpm[!hk, , drop = FALSE]

UNIFIED_META <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
unified_meta <- fread(UNIFIED_META)
meta <- data.table(sample_id = colnames(dge))
meta <- merge(meta, unified_meta[, .(sample_id, dataset, fibrosis_stage, sex)], by = "sample_id", all.x = TRUE, sort = FALSE)
setnames(meta, "dataset", "cohort")
setnames(meta, "sex", "inferred_sex")
stopifnot(all(meta$sample_id == colnames(dge)))
stopifnot(all(colnames(logcpm) == colnames(dge)))

fs_raw <- as.character(meta$fibrosis_stage)
norm_stage <- function(x) {
  x <- toupper(trimws(x))
  x[x %in% c("0", "F0")] <- "F01"
  x[x %in% c("1", "F1")] <- "F01"
  x[x %in% c("2", "F2")] <- "F2"
  x[x %in% c("3", "F3")] <- "F34"
  x[x %in% c("4", "F4")] <- "F34"
  x[!(x %in% c("F01", "F2", "F34"))] <- NA_character_
  x
}
stage <- norm_stage(fs_raw)
message("Stage distribution: ",
        paste(paste0(names(table(stage, useNA = "ifany")), "=",
                     table(stage, useNA = "ifany")), collapse = ", "))

manifest <- meta[, .(sample_id, cohort, inferred_sex)][
  , stage := stage]
manifest <- manifest[!is.na(stage)]
manifest_summary <- manifest[, .(
  n_samples = .N,
  n_female  = sum(inferred_sex %in% c("F", "female", "Female"), na.rm = TRUE),
  n_male    = sum(inferred_sex %in% c("M", "male", "Male"), na.rm = TRUE)
), by = .(fibrosis_stage = stage, cohort)]
fwrite(manifest_summary,
       file.path(OUTDIR, "stage_sample_manifest.csv"))

compute_stage_edges <- function(stage_label) {
  idx <- which(stage == stage_label)
  if (length(idx) < MIN_STAGE_N) {
    message("Stage ", stage_label, " has only ", length(idx),
            " samples (< ", MIN_STAGE_N, "); skipping.")
    return(invisible(NULL))
  }
  t_stage <- Sys.time()
  message("[", t_stage, "] Stage ", stage_label, ": ",
          length(idx), " samples, ", nrow(logcpm), " genes")

  m <- logcpm[, idx, drop = FALSE]
  keep_var <- apply(m, 1, function(v) sd(v) > 0)
  m <- m[keep_var, , drop = FALSE]
  message("  Genes with non-zero variance: ", nrow(m))

  cmat <- cor(t(m), method = "pearson")
  message("  cor() done in ",
          round(as.numeric(difftime(Sys.time(), t_stage, units = "mins")), 2),
          " min")

  g <- rownames(cmat)
  ng <- length(g)
  ut <- which(upper.tri(cmat), arr.ind = TRUE)
  ord <- g[ut[, 1]] < g[ut[, 2]]
  ga <- ifelse(ord, g[ut[, 1]], g[ut[, 2]])
  gb <- ifelse(ord, g[ut[, 2]], g[ut[, 1]])

  edges <- data.table(
    gene_a    = ga,
    gene_b    = gb,
    r         = cmat[ut],
    n_samples = length(idx)
  )
  rm(cmat); gc()

  out_path <- file.path(OUTDIR,
                        paste0("edges_coexpr_", stage_label, ".csv"))
  fwrite(edges, out_path)
  message("  Wrote ", nrow(edges), " edges -> ", out_path,
          " (total ",
          round(as.numeric(difftime(Sys.time(), t_stage, units = "mins")), 2),
          " min)")
  invisible(NULL)
}

for (s in c("F01", "F2", "F34")) compute_stage_edges(s)

message("[", Sys.time(), "] Script 280 done in ",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
        " min")
