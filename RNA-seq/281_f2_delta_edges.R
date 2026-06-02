#!/usr/bin/env Rscript
#SBATCH --partition=cpu
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --job-name=net_281_f2_delta
#SBATCH --output=logs/net_281_f2_delta_%j.out
#SBATCH --error=logs/net_281_f2_delta_%j.err
# ===========================================================================
# Script 281: F2 Delta-r edges (D-F2 edge type)
# ---------------------------------------------------------------------------
# Merges per-stage correlation tables from Script 280, computes stage-to-stage
# delta_r, classifies emergence_stage, and assesses significance via a
# permutation null that shuffles fibrosis_stage labels within cohort x sex
# strata. Applies BH FDR on the permutation p-values and filters on
# |delta_max| >= 0.2. Output is edges_d_f2.csv.
# Environment: micromamba activate rnaseq
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(parallel)
})

t0 <- Sys.time()

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

NETDIR <- file.path(BASE, "RNA-seq/results/network")
STGDIR <- file.path(NETDIR, "stage")
DGE_PATH <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_dge.rds")
NODE_PATH <- file.path(BASE, "RNA-seq/results/network/network_nodes.csv")

N_PERM      <- 100
N_PAIR_NULL <- 10000L
FDR_THRESH  <- 0.05
DELTA_FLOOR <- 0.2
N_CORES     <- min(8L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",
                                             unset = "8")))

message("[", Sys.time(), "] Loading stage edge tables")
e01 <- fread(file.path(STGDIR, "edges_coexpr_F01.csv"))
e2  <- fread(file.path(STGDIR, "edges_coexpr_F2.csv"))
e34 <- fread(file.path(STGDIR, "edges_coexpr_F34.csv"))

setnames(e01, "r", "r_F01"); e01[, n_samples := NULL]
setnames(e2,  "r", "r_F2");  n_F2 <- e2$n_samples[1]; e2[, n_samples := NULL]
setnames(e34, "r", "r_F34"); e34[, n_samples := NULL]

message("  F01 edges=", nrow(e01), " F2 edges=", nrow(e2),
        " F34 edges=", nrow(e34))

merged <- merge(e01, e2,  by = c("gene_a", "gene_b"), all = TRUE)
merged <- merge(merged, e34, by = c("gene_a", "gene_b"), all = TRUE)
message("  Merged pairs: ", nrow(merged))

merged[, delta_F2_vs_F01  := r_F2  - r_F01]
merged[, delta_F34_vs_F2  := r_F34 - r_F2]
merged[, delta_F34_vs_F01 := r_F34 - r_F01]
merged[, delta_max := pmax(abs(delta_F2_vs_F01),
                           abs(delta_F34_vs_F01), na.rm = TRUE)]

classify <- function(r01, r2, r34, d_21, d_32) {
  n <- length(r01)
  out <- rep("invariant", n)
  a21 <- abs(d_21); a32 <- abs(d_32)

  mono_up   <- !is.na(r01) & !is.na(r2) & !is.na(r34) &
               (r01 < r2) & (r2 < r34)
  mono_down <- !is.na(r01) & !is.na(r2) & !is.na(r34) &
               (r01 > r2) & (r2 > r34)

  emerging   <- !is.na(a21) & !is.na(a32) & (a21 > a32) &
                !is.na(d_21) & (d_21 > 0) &
                !is.na(r01) & (r01 < 0.3) &
                !is.na(r2)  & (r2  >= 0.5)
  dissolving <- !is.na(a21) & !is.na(a32) & (a21 > a32) &
                !is.na(d_21) & (d_21 < 0) &
                !is.na(r01) & (r01 >= 0.5) &
                !is.na(r2)  & (r2  < 0.3)
  transient  <- !is.na(r01) & !is.na(r2) & !is.na(r34) &
                (((r2 - r01) > 0.2 & (r2 - r34) > 0.2) |
                 ((r01 - r2) > 0.2 & (r34 - r2) > 0.2))
  f34_spec   <- !is.na(a21) & !is.na(a32) & (a32 > a21)

  out[f34_spec]   <- "F34_specific"
  out[mono_up]    <- "progressive_up"
  out[mono_down]  <- "progressive_down"
  out[transient]  <- "transient_F2"
  out[emerging]   <- "F2_emerging"
  out[dissolving] <- "F2_dissolving"
  out
}

merged[, emergence_stage := classify(r_F01, r_F2, r_F34,
                                     delta_F2_vs_F01, delta_F34_vs_F2)]

message("[", Sys.time(), "] Emergence stage counts (pre-FDR):")
print(merged[, .N, by = emergence_stage])

message("[", Sys.time(), "] Building permutation null")
dge <- readRDS(DGE_PATH)
nodes <- fread(NODE_PATH)
strip_version <- function(x) sub("\\..*$", "", x)
sym_map <- setNames(nodes$human_symbol, strip_version(nodes$ensembl_id))

rn <- strip_version(rownames(dge))
keep <- rn %in% names(sym_map)
dge <- dge[keep, , keep.lib.sizes = FALSE]
row_sym <- sym_map[strip_version(rownames(dge))]
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
rownames(logcpm) <- row_sym
logcpm <- logcpm[!duplicated(rownames(logcpm)), , drop = FALSE]

genes_present <- intersect(rownames(logcpm),
                           unique(c(merged$gene_a, merged$gene_b)))
logcpm <- logcpm[genes_present, , drop = FALSE]

UNIFIED_META <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
unified_meta <- fread(UNIFIED_META)
meta <- data.table(sample_id = colnames(dge))
meta <- merge(meta, unified_meta[, .(sample_id, dataset, fibrosis_stage, sex)], by = "sample_id", all.x = TRUE, sort = FALSE)
setnames(meta, "dataset", "cohort")
setnames(meta, "sex", "inferred_sex")
stopifnot(all(meta$sample_id == colnames(dge)))
fs_raw <- as.character(meta$fibrosis_stage)
norm_stage <- function(x) {
  x <- toupper(trimws(x))
  x[x %in% c("0", "F0", "1", "F1")] <- "F01"
  x[x %in% c("2", "F2")]            <- "F2"
  x[x %in% c("3", "F3", "4", "F4")] <- "F34"
  x[!(x %in% c("F01", "F2", "F34"))] <- NA_character_
  x
}
stage <- norm_stage(fs_raw)
cohort <- as.character(meta$cohort)
sex    <- as.character(meta$inferred_sex)

keep_s <- !is.na(stage)
stage <- stage[keep_s]; cohort <- cohort[keep_s]; sex <- sex[keep_s]
logcpm <- logcpm[, keep_s, drop = FALSE]

strata <- paste(cohort, sex, sep = "|")
stratum_idx <- split(seq_along(strata), strata)

set.seed(42)
n_pairs_total <- nrow(merged)
pair_sample_idx <- sample.int(n_pairs_total,
                              size = min(N_PAIR_NULL, n_pairs_total))
pair_sample <- merged[pair_sample_idx, .(gene_a, gene_b)]

ga_in <- pair_sample$gene_a %in% rownames(logcpm)
gb_in <- pair_sample$gene_b %in% rownames(logcpm)
pair_sample <- pair_sample[ga_in & gb_in]
message("  Pair sample usable for null: ", nrow(pair_sample))

A <- logcpm[pair_sample$gene_a, , drop = FALSE]
B <- logcpm[pair_sample$gene_b, , drop = FALSE]

row_cor <- function(A, B, idx) {
  a <- A[, idx, drop = FALSE]
  b <- B[, idx, drop = FALSE]
  am <- rowMeans(a); bm <- rowMeans(b)
  a <- a - am; b <- b - bm
  num <- rowSums(a * b)
  den <- sqrt(rowSums(a * a) * rowSums(b * b))
  num / den
}

perm_once <- function(seed) {
  set.seed(seed)
  new_stage <- stage
  for (s in stratum_idx) {
    if (length(s) > 1) new_stage[s] <- sample(stage[s])
  }
  i01 <- which(new_stage == "F01")
  i2  <- which(new_stage == "F2")
  i34 <- which(new_stage == "F34")
  if (length(i01) < 3 || length(i2) < 3 || length(i34) < 3) return(NULL)
  r01 <- row_cor(A, B, i01)
  r2  <- row_cor(A, B, i2)
  r34 <- row_cor(A, B, i34)
  d21 <- r2 - r01
  d31 <- r34 - r01
  pmax(abs(d21), abs(d31), na.rm = TRUE)
}

message("[", Sys.time(), "] Running ", N_PERM,
        " permutations with ", N_CORES, " cores")
null_list <- mclapply(seq_len(N_PERM), perm_once, mc.cores = N_CORES)
null_list <- null_list[!vapply(null_list, is.null, logical(1))]
null_vals <- unlist(null_list, use.names = FALSE)
null_vals <- null_vals[is.finite(null_vals)]
message("  Null delta_max: n=", length(null_vals),
        " mean=", round(mean(null_vals), 4),
        " q95=", round(quantile(null_vals, 0.95), 4),
        " q99=", round(quantile(null_vals, 0.99), 4))

fwrite(data.table(delta_max_null = null_vals),
       file.path(STGDIR, "f2_delta_null_distribution.csv"))

obs <- merged$delta_max
null_sorted <- sort(null_vals)
n_null <- length(null_sorted)
ranks <- findInterval(obs, null_sorted)
p_perm <- (n_null - ranks + 1L) / (n_null + 1L)
p_perm[is.na(obs)] <- NA_real_
merged[, p_perm := p_perm]
merged[, fdr := p.adjust(p_perm, method = "BH")]

summary_tbl <- merged[, .(
  n_total = .N,
  n_fdr_sig = sum(!is.na(fdr) & fdr < FDR_THRESH, na.rm = TRUE),
  n_passing = sum(!is.na(fdr) & fdr < FDR_THRESH &
                  !is.na(delta_max) & delta_max >= DELTA_FLOOR,
                  na.rm = TRUE)
), by = emergence_stage]
fwrite(summary_tbl, file.path(STGDIR, "f2_delta_summary.csv"))
message("[", Sys.time(), "] Summary:")
print(summary_tbl)

keep_edge <- !is.na(merged$fdr) & merged$fdr < FDR_THRESH &
             !is.na(merged$delta_max) & merged$delta_max >= DELTA_FLOOR
out_edges <- merged[keep_edge,
  .(gene_a, gene_b, r_F01, r_F2, r_F34,
    delta_F2_vs_F01, delta_F34_vs_F2, delta_F34_vs_F01,
    delta_max, emergence_stage, p_perm, fdr)]

out_path <- file.path(NETDIR, "edges_d_f2.csv")
fwrite(out_edges, out_path)
message("[", Sys.time(), "] Wrote ", nrow(out_edges),
        " D-F2 edges -> ", out_path)

message("[", Sys.time(), "] Script 281 done in ",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
        " min")
