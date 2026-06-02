#!/usr/bin/env Rscript
#SBATCH --partition=cpu
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --job-name=net_281c_f2_continuous
#SBATCH --output=logs/net_281c_f2_continuous_%j.out
#SBATCH --error=logs/net_281c_f2_continuous_%j.err
# ===========================================================================
# Script 281c: F2 delta edges — CONTINUOUS variant.
# ---------------------------------------------------------------------------
# Same inputs as 281 / 281b, but WITHOUT the hard emergence_stage cutpoints
# (r_F01 < 0.3 & r_F2 >= 0.5, etc). The portal UI slider thresholds
# |delta_max| at display time; the edge table just needs to be FDR-clean.
#
# - Keeps 281's permutation FDR framework; reuses the null distribution
#   already cached at stage/f2_delta_null_distribution.csv (1M samples).
# - Drops the emergence_stage hard gates. Instead, assigns a soft
#   `direction_at_F2` label from the sign of delta_F2_vs_F01.
# - |delta_max| floor lowered from 0.2 -> 0.1 so the UI slider has room.
# - Output columns match edges_d_f2_loco.csv so existing atlas merge works.
#
# Output: edges_d_f2_continuous.csv
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
})

t0 <- Sys.time()

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

NETDIR    <- file.path(PROJ, "RNA-seq/results/network")
STGDIR    <- file.path(NETDIR, "stage")
NODE_PATH <- file.path(NETDIR, "network_nodes.csv")
NULL_PATH <- file.path(STGDIR, "f2_delta_null_distribution.csv")
OUT       <- file.path(NETDIR, "edges_d_f2_continuous.csv")

FDR_THRESH   <- 0.05
DELTA_FLOOR  <- 0.10

message("[281c] Continuous F2-delta edges starting at ", format(t0))

stopifnot(file.exists(NODE_PATH),
          file.exists(NULL_PATH),
          file.exists(file.path(STGDIR, "edges_coexpr_F01.csv")),
          file.exists(file.path(STGDIR, "edges_coexpr_F2.csv")),
          file.exists(file.path(STGDIR, "edges_coexpr_F34.csv")))

nodes <- fread(NODE_PATH)
V <- unique(nodes$human_symbol)
message("  Node set V: ", length(V), " genes")

message("[", Sys.time(), "] Loading stage tables")
e01 <- fread(file.path(STGDIR, "edges_coexpr_F01.csv"))
setnames(e01, "r", "r_F01"); e01[, n_samples := NULL]
e2 <- fread(file.path(STGDIR, "edges_coexpr_F2.csv"))
setnames(e2, "r", "r_F2"); e2[, n_samples := NULL]
e34 <- fread(file.path(STGDIR, "edges_coexpr_F34.csv"))
setnames(e34, "r", "r_F34"); e34[, n_samples := NULL]

message("[", Sys.time(), "] Pre-filter to V on both endpoints")
e01 <- e01[gene_a %in% V & gene_b %in% V]
e2  <- e2 [gene_a %in% V & gene_b %in% V]
e34 <- e34[gene_a %in% V & gene_b %in% V]

message("[", Sys.time(), "] Outer-join")
merged <- merge(e01, e2,  by = c("gene_a", "gene_b"), all = TRUE)
merged <- merge(merged, e34, by = c("gene_a", "gene_b"), all = TRUE)
rm(e01, e2, e34); invisible(gc(verbose = FALSE))
message("  Merged pairs: ", nrow(merged))

merged[, delta_F2_vs_F01  := r_F2  - r_F01]
merged[, delta_F34_vs_F2  := r_F34 - r_F2]
merged[, delta_F34_vs_F01 := r_F34 - r_F01]
merged[, delta_max := pmax(abs(delta_F2_vs_F01),
                           abs(delta_F34_vs_F01),
                           abs(delta_F34_vs_F2), na.rm = TRUE)]

message("[", Sys.time(), "] Loading cached null distribution")
null_vals <- fread(NULL_PATH)[[1]]
null_vals <- null_vals[is.finite(null_vals)]
message("  Null n=", length(null_vals),
        " mean=", round(mean(null_vals), 4),
        " q95=", round(quantile(null_vals, 0.95), 4))

message("[", Sys.time(), "] Computing permutation p-values and BH FDR")
null_sorted <- sort(null_vals)
n_null <- length(null_sorted)
ranks <- findInterval(merged$delta_max, null_sorted)
merged[, p_perm := (n_null - ranks + 1L) / (n_null + 1L)]
merged[is.na(delta_max), p_perm := NA_real_]
merged[, fdr := p.adjust(p_perm, method = "BH")]

message("[", Sys.time(), "] Applying FDR <= ", FDR_THRESH,
        " and |delta_max| >= ", DELTA_FLOOR)
keep <- !is.na(merged$fdr) & merged$fdr <= FDR_THRESH &
        !is.na(merged$delta_max) & merged$delta_max >= DELTA_FLOOR
out <- merged[keep]
message("  Passing edges: ", nrow(out),
        sprintf("  (%.3f%% of merged)", 100 * nrow(out) / nrow(merged)))

# Soft direction label at F2 (not a hard class; just helpful for UI)
out[, direction_at_F2 := fifelse(
  !is.na(delta_F2_vs_F01),
  fifelse(delta_F2_vs_F01 > 0, "F2_up", "F2_down"),
  NA_character_
)]

# Canonicalize order (undirected)
swap <- out$gene_a > out$gene_b
a2 <- out$gene_a; b2 <- out$gene_b
out[, gene_a := ifelse(swap, b2, a2)]
out[, gene_b := ifelse(swap, a2, b2)]
out <- out[gene_a != gene_b]
out <- unique(out, by = c("gene_a", "gene_b"))

emit <- out[, .(
  gene_a, gene_b,
  r_F01, r_F2, r_F34,
  delta_F2_vs_F01, delta_F34_vs_F2, delta_F34_vs_F01, delta_max,
  direction_at_F2,
  p_perm, fdr,
  raw_score = delta_max
)]

setorder(emit, -delta_max, gene_a, gene_b)

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
fwrite(emit, OUT)
message("[", Sys.time(), "] Wrote ", nrow(emit),
        " continuous D-F2 edges -> ", OUT)

# ---- Diagnostics ----
message("  |delta_max| quantiles:")
print(round(quantile(emit$delta_max, c(0.5, 0.75, 0.9, 0.95, 0.99)), 3))
message("  Direction breakdown:")
print(emit[, .N, by = direction_at_F2])
message("  Genes touched: ",
        length(unique(c(emit$gene_a, emit$gene_b))),
        " / ", length(V))
message("[281c] Elapsed: ",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
        " min")
