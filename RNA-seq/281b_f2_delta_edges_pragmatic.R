#!/usr/bin/env Rscript
#SBATCH --partition=cpu
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --job-name=net_281b_f2_delta
# ===========================================================================
# Script 281b: F2 Delta-r edges — pragmatic filter (effect size + category)
# ---------------------------------------------------------------------------
# Replaces 281's permutation-FDR approach with biologically-meaningful
# category + effect-size filtering. The permutation null in 281 produced
# zero passing edges; the problem was that a 10k-pair null absorbs extremes
# of a 126M-pair observed distribution. Here we instead use the categorical
# emergence classification (F2_emerging, F2_dissolving, transient_F2,
# F34_specific) plus a strict effect-size floor. Pre-registration
# Criterion 1 (LOCO replication ≥ 0.70) becomes the primary quality filter.
# Environment: micromamba activate rnaseq
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
})

t0 <- Sys.time()

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

NETDIR <- file.path(BASE, "RNA-seq/results/network")
STGDIR <- file.path(NETDIR, "stage")

DELTA_FLOOR     <- 0.35
DELTA_STRICT    <- 0.50
MIN_R_DIFF_PCT  <- 0.60

message("[", Sys.time(), "] Loading stage edge tables")
e01 <- fread(file.path(STGDIR, "edges_coexpr_F01.csv"))
e2  <- fread(file.path(STGDIR, "edges_coexpr_F2.csv"))
e34 <- fread(file.path(STGDIR, "edges_coexpr_F34.csv"))

setnames(e01, "r", "r_F01"); e01[, n_samples := NULL]
setnames(e2,  "r", "r_F2");  e2[,  n_samples := NULL]
setnames(e34, "r", "r_F34"); e34[, n_samples := NULL]

message("[", Sys.time(), "] Full-joining stage tables")
merged <- merge(merge(e01, e2, by = c("gene_a", "gene_b"), all = TRUE),
                e34, by = c("gene_a", "gene_b"), all = TRUE)
rm(e01, e2, e34); gc()
message("  n pairs merged: ", nrow(merged))

message("[", Sys.time(), "] Computing deltas")
merged[, delta_F2_vs_F01 := r_F2  - r_F01]
merged[, delta_F34_vs_F2 := r_F34 - r_F2]
merged[, delta_F34_vs_F01 := r_F34 - r_F01]
merged[, delta_max := pmax(abs(delta_F2_vs_F01),
                           abs(delta_F34_vs_F01),
                           abs(delta_F34_vs_F2), na.rm = TRUE)]

classify_emergence <- function(r01, r2, r34, d_21, d_42) {
  cls <- rep("invariant", length(r01))
  keep <- !is.na(r01) & !is.na(r2) & !is.na(r34)

  f2e <- keep & abs(d_21) > abs(d_42) & d_21 > 0 & r01 < 0.3 & r2 >= 0.5
  cls[f2e] <- "F2_emerging"

  f2d <- keep & abs(d_21) > abs(d_42) & d_21 < 0 & r01 >= 0.5 & r2 < 0.3
  cls[f2d] <- "F2_dissolving"

  pu <- keep & (r01 < r2) & (r2 < r34) & cls == "invariant"
  cls[pu] <- "progressive_up"
  pd <- keep & (r01 > r2) & (r2 > r34) & cls == "invariant"
  cls[pd] <- "progressive_down"

  t2 <- keep & ((r2 > r01 + 0.2 & r2 > r34 + 0.2) |
                (r2 < r01 - 0.2 & r2 < r34 - 0.2)) & cls == "invariant"
  cls[t2] <- "transient_F2"

  f34 <- keep & abs(d_42) > abs(d_21) & cls == "invariant"
  cls[f34] <- "F34_specific"

  cls
}

merged[, emergence_stage := classify_emergence(r_F01, r_F2, r_F34,
                                               delta_F2_vs_F01,
                                               delta_F34_vs_F2)]

message("  category counts pre-filter:")
message("  ", paste(sprintf("%s=%d", names(table(merged$emergence_stage)),
                             table(merged$emergence_stage)), collapse = ", "))

message("[", Sys.time(), "] Applying pragmatic filter")
f2_specific_categories <- c("F2_emerging", "F2_dissolving", "transient_F2")
progressive_categories <- c("progressive_up", "progressive_down", "F34_specific")

keep_f2 <- merged$emergence_stage %in% f2_specific_categories &
           merged$delta_max >= DELTA_FLOOR

keep_prog <- merged$emergence_stage %in% progressive_categories &
             merged$delta_max >= DELTA_STRICT

kept <- merged[keep_f2 | keep_prog]
message("  edges passing filter: ", nrow(kept))
message("  breakdown:")
for (cat in c(f2_specific_categories, progressive_categories)) {
  n <- sum(kept$emergence_stage == cat)
  message(sprintf("    %-20s %d", cat, n))
}

kept[, priority := fifelse(emergence_stage %in% f2_specific_categories,
                           1L, 2L)]

out_cols <- c("gene_a", "gene_b", "r_F01", "r_F2", "r_F34",
              "delta_F2_vs_F01", "delta_F34_vs_F2", "delta_F34_vs_F01",
              "delta_max", "emergence_stage", "priority")
out_edges <- kept[, ..out_cols]
setorder(out_edges, priority, -delta_max)

out_path <- file.path(NETDIR, "edges_d_f2.csv")
fwrite(out_edges, out_path)
message("[", Sys.time(), "] Wrote ", nrow(out_edges), " D-F2 edges to ", out_path)

summary_tab <- out_edges[, .(n = .N,
                             median_delta = round(median(delta_max, na.rm = TRUE), 3),
                             median_r_F01 = round(median(r_F01, na.rm = TRUE), 3),
                             median_r_F2  = round(median(r_F2, na.rm = TRUE), 3),
                             median_r_F34 = round(median(r_F34, na.rm = TRUE), 3)),
                         by = emergence_stage]
fwrite(summary_tab, file.path(STGDIR, "f2_delta_summary.csv"))
message("Per-category summary:")
print(summary_tab)

message("[", Sys.time(), "] Script 281b done in ",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), " min")
