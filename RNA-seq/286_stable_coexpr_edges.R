#!/usr/bin/env Rscript
#SBATCH --partition=cpu
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --job-name=net_286_dstable
#SBATCH --output=logs/net_286_dstable_%j.out
#SBATCH --error=logs/net_286_dstable_%j.err
# ===========================================================================
# Script 286: D-STABLE edges -- stable co-expression across fibrosis stages.
# ---------------------------------------------------------------------------
# Complement of D-F2 (script 281). Reads the three per-stage correlation
# tables produced by 280 and selects pairs that are strongly AND stably
# co-expressed across F01, F2, F34. No permutation test: stability is
# defined geometrically (all three stage correlations close, all strong).
#
# Selection criteria (conjunction, all required):
#   r_mean       >= 0.6           (strong average co-expression)
#   delta_max    <  0.15          (stage-to-stage wobble bounded)
#   min_r        >= 0.40          (not rescued by a single outlier stage)
#   max_r - min_r < 0.25          (range constraint, extra guard)
#
# Inputs:
#   RNA-seq/results/network/stage/edges_coexpr_F{01,2,34}.csv
#   RNA-seq/results/network/network_nodes.csv          (node set V, symbols)
#
# Output:
#   RNA-seq/results/network/edges_d_stable.csv
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
OUT       <- file.path(NETDIR, "edges_d_stable.csv")

R_MEAN_MIN  <- 0.60
DELTA_CEIL  <- 0.15
MIN_R_FLOOR <- 0.40
RANGE_CEIL  <- 0.25

message("[286] D-STABLE edge construction starting at ", format(t0))

stopifnot(file.exists(NODE_PATH),
          file.exists(file.path(STGDIR, "edges_coexpr_F01.csv")),
          file.exists(file.path(STGDIR, "edges_coexpr_F2.csv")),
          file.exists(file.path(STGDIR, "edges_coexpr_F34.csv")))

nodes <- fread(NODE_PATH)
V <- unique(nodes$human_symbol)
message("  Node set V: ", length(V), " genes (symbols)")

message("[", Sys.time(), "] Loading per-stage tables")
e01 <- fread(file.path(STGDIR, "edges_coexpr_F01.csv"))
setnames(e01, "r", "r_F01"); e01[, n_samples := NULL]
message("  F01 rows: ", nrow(e01))

e2 <- fread(file.path(STGDIR, "edges_coexpr_F2.csv"))
setnames(e2, "r", "r_F2"); e2[, n_samples := NULL]
message("  F2 rows:  ", nrow(e2))

e34 <- fread(file.path(STGDIR, "edges_coexpr_F34.csv"))
setnames(e34, "r", "r_F34"); e34[, n_samples := NULL]
message("  F34 rows: ", nrow(e34))

message("[", Sys.time(), "] Restricting to V BEFORE merge (saves memory)")
e01 <- e01[gene_a %in% V & gene_b %in% V]
e2  <- e2 [gene_a %in% V & gene_b %in% V]
e34 <- e34[gene_a %in% V & gene_b %in% V]
message("  F01 in V: ", nrow(e01))
message("  F2 in V:  ", nrow(e2))
message("  F34 in V: ", nrow(e34))

message("[", Sys.time(), "] Inner-join across stages (a stable edge needs r in all three)")
merged <- merge(e01, e2,  by = c("gene_a", "gene_b"), all = FALSE)
merged <- merge(merged, e34, by = c("gene_a", "gene_b"), all = FALSE)
message("  Pairs with r in all three stages: ", nrow(merged))
rm(e01, e2, e34); invisible(gc(verbose = FALSE))

merged[, r_mean  := (r_F01 + r_F2 + r_F34) / 3]
merged[, min_r   := pmin(r_F01, r_F2, r_F34)]
merged[, max_r   := pmax(r_F01, r_F2, r_F34)]
merged[, range_r := max_r - min_r]
merged[, delta_F2_vs_F01  := r_F2  - r_F01]
merged[, delta_F34_vs_F2  := r_F34 - r_F2]
merged[, delta_F34_vs_F01 := r_F34 - r_F01]
merged[, delta_max := pmax(abs(delta_F2_vs_F01),
                           abs(delta_F34_vs_F2),
                           abs(delta_F34_vs_F01))]

message("[", Sys.time(), "] Applying stability filters")
keep_pos <- merged[r_mean >= R_MEAN_MIN &
                   delta_max < DELTA_CEIL &
                   min_r >= MIN_R_FLOOR &
                   range_r < RANGE_CEIL]

keep_neg <- merged[r_mean <= -R_MEAN_MIN &
                   delta_max < DELTA_CEIL &
                   max_r <= -MIN_R_FLOOR &
                   range_r < RANGE_CEIL]

stable <- rbind(keep_pos, keep_neg)
rm(merged, keep_pos, keep_neg); invisible(gc(verbose = FALSE))
message("  Stable pairs (positive |r| >= ", R_MEAN_MIN, "): ",
        nrow(stable))
message("  Direction breakdown (r_mean sign):")
print(stable[, .N, by = sign(r_mean)])

# Canonicalize (gene_a <= gene_b alphabetically, undirected)
swap <- stable$gene_a > stable$gene_b
a2 <- stable$gene_a; b2 <- stable$gene_b
stable[, gene_a := ifelse(swap, b2, a2)]
stable[, gene_b := ifelse(swap, a2, b2)]
stable <- stable[gene_a != gene_b]
stable <- unique(stable, by = c("gene_a", "gene_b"))

# Emit in edges_d_* schema
out <- stable[, .(
  gene_a, gene_b,
  r_F01, r_F2, r_F34,
  r_mean, min_r, max_r, range_r,
  delta_F2_vs_F01, delta_F34_vs_F2, delta_F34_vs_F01, delta_max,
  raw_score = abs(r_mean)
)]

setorder(out, -raw_score, gene_a, gene_b)

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
fwrite(out, OUT)
message("[", Sys.time(), "] Wrote ", nrow(out),
        " D-STABLE edges -> ", OUT)

# ---- Diagnostics ----
message("  Mean |r_mean|: ", round(mean(abs(out$r_mean)), 4))
message("  Mean delta_max: ", round(mean(out$delta_max), 4))
message("  Genes touched: ",
        length(unique(c(out$gene_a, out$gene_b))),
        " / ", length(V))
message("  Top 10 hub genes by degree:")
deg <- rbind(out[, .(g = gene_a)], out[, .(g = gene_b)])
print(deg[, .N, by = g][order(-N)][1:10])

message("[286] Elapsed: ",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
        " min")
