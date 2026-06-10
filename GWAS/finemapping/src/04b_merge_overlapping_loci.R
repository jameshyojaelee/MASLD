#!/usr/bin/env Rscript
# 04b_merge_overlapping_loci.R
# Post-hoc locus merging: identifies overlapping lead SNPs within each study
# and selects a single representative locus per group.
#
# Problem: FinnGen (and some other studies) report many GWAS lead SNPs within
# 100-500kb of each other (e.g., 7 lead SNPs on chr19:19-20Mb). These produce
# redundant credible sets covering the same causal signal.
#
# Algorithm:
#   1. Parse locus IDs (chr.bp) from SuSiE result filenames
#   2. Within each study, cluster loci on the same chromosome where
#      lead SNP positions are within MERGE_DIST (500kb) of each other
#      (single-linkage: loci A,B merge if any pair is <500kb)
#   3. Post-hoc split: break apart any chained group where:
#      a. Total span exceeds 1Mb, OR
#      b. Converged loci have ZERO shared credible set variants
#      Split point: largest gap between consecutive loci (recursive)
#   4. Within each cluster, keep the locus with:
#      a. The most variants in the 95% credible set, OR
#      b. The highest max PIP if tied on CS size
#   5. Output mapping table + per-study summary (with split_reason column)
#
# Outputs:
#   results/merged_loci_map.csv      — per-locus decisions
#   results/independent_loci_summary.csv — de-duplicated counts per study
#
# Usage: Rscript src/04b_merge_overlapping_loci.R
# Environment: finemapping conda env with data.table

suppressPackageStartupMessages({
  library(data.table)
})

# ── Configuration ──────────────────────────────────────────────────────────
BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR  <- file.path(BASE_DIR, "GWAS/finemapping")
setwd(FM_DIR)

MERGE_DIST <- 500000  # 500kb merge window
RESULTS_DIR <- file.path(FM_DIR, "results")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("Post-hoc Locus Merging\n")
cat("Merge distance:", MERGE_DIST / 1000, "kb\n")
cat("============================================================\n\n")

# ── Step 1: Discover all SuSiE result files ────────────────────────────────
# Glob both converged (*_adjustedvar.tsv) and not-converged (*_notconverged.tsv)
susie_files <- Sys.glob(file.path(FM_DIR, "output", "*", "*Mb", "susie", "*_cov0.95_*.tsv"))

if (length(susie_files) == 0) {
  stop("No SuSiE result files found. Run finemapping first.")
}
cat("Found", length(susie_files), "SuSiE result files\n\n")

# ── Step 2: Parse metadata from file paths and read CS/PIP stats ───────────
parse_one_file <- function(fpath) {
  # Extract study name from path: output/{study}/{pop}_{window}Mb/susie/...
  parts <- strsplit(fpath, "/")[[1]]
  susie_idx <- which(parts == "susie")
  study <- parts[susie_idx - 2]

  # Extract locus from filename: e.g. UKBB_19.19393714_cov0.95_notconverged.tsv
  fname <- parts[susie_idx + 1]
  # Locus pattern: digits.digits (chr.bp) — find it in the filename

  locus_match <- regmatches(fname, regexpr("\\d+\\.\\d+", fname))
  if (length(locus_match) == 0) return(NULL)
  locus <- locus_match[1]

  # Determine convergence status from filename
  converged <- !grepl("notconverged", fname)

  # Read the file to get CS and PIP statistics
  dt <- tryCatch(fread(fpath, select = c("PIP", "CS")), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) {
    return(data.table(
      study = study, locus = locus,
      chr = as.integer(sub("\\..*", "", locus)),
      bp = as.integer(sub(".*\\.", "", locus)),
      n_cs_variants = 0L, max_pip = 0, n_variants = 0L,
      converged = converged, file = fpath
    ))
  }

  # Count variants in credible set (CS > 0) and max PIP
  n_cs <- sum(dt$CS > 0, na.rm = TRUE)
  max_pip <- max(dt$PIP, na.rm = TRUE)
  n_var <- nrow(dt)

  data.table(
    study = study, locus = locus,
    chr = as.integer(sub("\\..*", "", locus)),
    bp = as.integer(sub(".*\\.", "", locus)),
    n_cs_variants = n_cs, max_pip = max_pip, n_variants = n_var,
    converged = converged, file = fpath
  )
}

cat("Parsing SuSiE result files...\n")
locus_list <- rbindlist(lapply(susie_files, parse_one_file))
locus_list <- locus_list[!is.null(study)]

# De-duplicate: if a locus appears in both converged + notconverged (shouldn't happen,
# but be safe), keep the converged version
locus_list <- locus_list[order(study, locus, -converged)]
locus_list <- locus_list[!duplicated(paste(study, locus))]

# Restrict to the active GWAS portfolio in the registry (cirrhosis/HCC GWAS were
# dropped 2026-06 to mirror the eQTL-COLOC portfolio cut). output/ still contains
# stale per-GWAS dirs (cirrhosis/HCC + retired 2023_36653562_*/Anstee/Pazoki runs)
# that the glob in Step 1 picks up; filter them out so the summary tables reflect
# the canonical 23-study portfolio.
REGISTRY <- fread(file.path(FM_DIR, "config/gwas_registry.tsv"))
active_studies <- unique(REGISTRY$study_name)
n_before <- length(unique(locus_list$study))
locus_list <- locus_list[study %in% active_studies]
cat("Restricted to registry portfolio:", n_before, "->",
    length(unique(locus_list$study)), "studies\n")

cat("Parsed", nrow(locus_list), "unique study-locus entries across",
    length(unique(locus_list$study)), "studies\n\n")

# ── Step 3: Cluster overlapping loci per study ─────────────────────────────
# Single-linkage clustering: two loci on the same chromosome merge if their
# lead SNP positions are within MERGE_DIST bp of each other.

assign_clusters <- function(positions, merge_dist) {
  n <- length(positions)
  if (n <= 1) return(rep(1L, n))

  ord <- order(positions)
  sorted <- positions[ord]
  cluster <- integer(n)
  cl <- 1L
  cluster[ord[1]] <- cl

  for (i in 2:n) {
    if (sorted[i] - sorted[i - 1] <= merge_dist) {
      cluster[ord[i]] <- cl
    } else {
      cl <- cl + 1L
      cluster[ord[i]] <- cl
    }
  }
  cluster
}

# Process each study x chromosome
locus_list[, merge_group := {
  cl <- assign_clusters(bp, MERGE_DIST)
  paste0(chr, "_g", cl)
}, by = .(study, chr)]

# ── Step 3b: Post-hoc split of chained groups ─────────────────────────────
# Single-linkage can chain distant loci through intermediate non-converged loci.
# Split any merged group where:
#   (a) total span exceeds MAX_SPAN (1Mb), OR
#   (b) converged loci within the group share ZERO credible set variants
# Split point: largest gap between consecutive loci (sorted by bp).

MAX_SPAN <- 1000000  # 1Mb

# Helper: read CS variant positions from a SuSiE TSV file
read_cs_positions <- function(fpath) {
  dt <- tryCatch(fread(fpath, select = c("position", "CS")), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(integer(0))
  dt[CS > 0, position]
}

# Helper: check whether any pair of converged loci in a group shares >= 1 CS variant
check_cs_overlap <- function(files, converged_flags) {
  conv_idx <- which(converged_flags)
  if (length(conv_idx) < 2) return(TRUE)  # <2 converged → nothing to compare → keep merged

  # Read CS positions for each converged locus
  cs_list <- lapply(files[conv_idx], read_cs_positions)

  # Check all pairs — if ANY pair shares a variant, overlap exists

  for (i in seq_along(cs_list)[-length(cs_list)]) {
    for (j in (i + 1):length(cs_list)) {
      if (length(intersect(cs_list[[i]], cs_list[[j]])) > 0) return(TRUE)
    }
  }
  return(FALSE)
}

# Helper: split a group at the largest gap, returns new sub-group labels
split_at_largest_gap <- function(positions) {
  n <- length(positions)
  if (n <= 1) return(rep(1L, n))
  ord <- order(positions)
  sorted <- positions[ord]
  gaps <- diff(sorted)
  split_idx <- which.max(gaps)  # index in sorted order where gap is largest
  # Loci 1..split_idx go to subgroup 1; split_idx+1..n go to subgroup 2
  sub <- integer(n)
  sub[ord[1:split_idx]] <- 1L
  sub[ord[(split_idx + 1):n]] <- 2L
  sub
}

cat("Post-hoc split check for chained groups...\n")

# Identify groups with >1 locus that may need splitting
multi_groups <- locus_list[, .(
  n_loci = .N,
  span = max(bp) - min(bp),
  n_converged = sum(converged)
), by = .(study, merge_group)]
multi_groups <- multi_groups[n_loci > 1]

locus_list[, split_reason := NA_character_]
n_splits <- 0L

for (i in seq_len(nrow(multi_groups))) {
  mg <- multi_groups[i]
  rows_idx <- which(locus_list$study == mg$study & locus_list$merge_group == mg$merge_group)
  grp <- locus_list[rows_idx]

  needs_split <- FALSE
  reason <- NA_character_

  # Check (a): span > 1Mb
  if (mg$span > MAX_SPAN) {
    needs_split <- TRUE
    reason <- sprintf("span_%.0fkb_exceeds_1Mb", mg$span / 1000)
  }

  # Check (b): converged loci share zero CS variants
  if (!needs_split && mg$n_converged >= 2) {
    has_overlap <- check_cs_overlap(grp$file, grp$converged)
    if (!has_overlap) {
      needs_split <- TRUE
      reason <- "zero_cs_overlap_between_converged"
    }
  }

  if (needs_split) {
    n_splits <- n_splits + 1L
    sub <- split_at_largest_gap(grp$bp)
    old_grp <- mg$merge_group
    new_grps <- paste0(old_grp, LETTERS[sub])

    set(locus_list, rows_idx, "merge_group", new_grps)
    set(locus_list, rows_idx, "split_reason", reason)

    cat(sprintf("  SPLIT %s | %s: %d loci, span=%.0fkb -> subgroups %s (%s)\n",
                mg$study, old_grp, mg$n_loci, mg$span / 1000,
                paste(unique(new_grps), collapse = "+"), reason))

    # Recursively check the new sub-groups — they might still need splitting
    # (e.g., 3 clusters chained together). We handle this by re-checking after
    # the loop below.
  }
}

# Second pass: re-check any newly created sub-groups that still span >1Mb or lack CS overlap
# (handles 3+-way chains)
repeat {
  did_split <- FALSE
  multi_groups2 <- locus_list[, .(
    n_loci = .N,
    span = max(bp) - min(bp),
    n_converged = sum(converged)
  ), by = .(study, merge_group)]
  multi_groups2 <- multi_groups2[n_loci > 1]

  for (i in seq_len(nrow(multi_groups2))) {
    mg <- multi_groups2[i]
    rows_idx <- which(locus_list$study == mg$study & locus_list$merge_group == mg$merge_group)
    grp <- locus_list[rows_idx]

    needs_split <- FALSE
    reason <- NA_character_

    if (mg$span > MAX_SPAN) {
      needs_split <- TRUE
      reason <- sprintf("span_%.0fkb_exceeds_1Mb", mg$span / 1000)
    }

    if (!needs_split && mg$n_converged >= 2) {
      # Only re-check if this sub-group doesn't already have a split_reason
      # (i.e., it was created by a prior split and hasn't been checked yet)
      has_overlap <- check_cs_overlap(grp$file, grp$converged)
      if (!has_overlap) {
        needs_split <- TRUE
        reason <- "zero_cs_overlap_between_converged"
      }
    }

    if (needs_split) {
      did_split <- TRUE
      n_splits <- n_splits + 1L
      sub <- split_at_largest_gap(grp$bp)
      old_grp <- mg$merge_group
      # Append letters to create unique sub-group names
      new_grps <- paste0(old_grp, LETTERS[sub])

      set(locus_list, rows_idx, "merge_group", new_grps)
      set(locus_list, rows_idx, "split_reason", reason)

      cat(sprintf("  RE-SPLIT %s | %s: %d loci, span=%.0fkb -> %s (%s)\n",
                  mg$study, old_grp, mg$n_loci, mg$span / 1000,
                  paste(unique(new_grps), collapse = "+"), reason))
    }
  }
  if (!did_split) break
}

cat(sprintf("Post-hoc splitting: %d groups split\n\n", n_splits))

# ── Step 4: Within each group, pick the best locus ─────────────────────────
# Priority: (1) most CS variants, (2) highest max PIP
locus_list[, group_size := .N, by = .(study, merge_group)]

# Rank within group: descending n_cs_variants, then descending max_pip
locus_list[, rank_in_group := frank(
  list(-n_cs_variants, -max_pip),
  ties.method = "first"
), by = .(study, merge_group)]

locus_list[, kept := (rank_in_group == 1L)]

# Assign reason
locus_list[, reason := fifelse(
  group_size == 1L, "singleton",
  fifelse(kept, "best_in_group", "merged_out")
)]

# Create merged_locus: for kept loci, use their own locus; for merged-out, point to the kept one
locus_list[, merged_locus := {
  kept_locus <- locus[rank_in_group == 1L][1]
  kept_locus
}, by = .(study, merge_group)]

# ── Step 5: Report ─────────────────────────────────────────────────────────
cat("=== Merge Report ===\n\n")

# Per-study summary
study_stats <- locus_list[, .(
  original_loci = .N,
  independent_loci = sum(kept),
  merged_out = sum(!kept)
), by = study][order(-merged_out)]

cat("Per-study locus counts (before -> after merging):\n")
for (i in 1:nrow(study_stats)) {
  s <- study_stats[i]
  flag <- if (s$merged_out > 0) " ***" else ""
  cat(sprintf("  %-45s %3d -> %3d  (-%d)%s\n",
              s$study, s$original_loci, s$independent_loci, s$merged_out, flag))
}

total_before <- sum(study_stats$original_loci)
total_after  <- sum(study_stats$independent_loci)
total_merged <- sum(study_stats$merged_out)
cat(sprintf("\nTotal: %d loci -> %d independent (%d merged out, %.1f%% reduction)\n",
            total_before, total_after, total_merged,
            100 * total_merged / total_before))

# Which studies were most affected?
affected <- study_stats[merged_out > 0][order(-merged_out)]
if (nrow(affected) > 0) {
  cat("\nMost affected studies:\n")
  for (i in 1:min(10, nrow(affected))) {
    s <- affected[i]
    cat(sprintf("  %d. %s: %d loci merged out (%.0f%% reduction)\n",
                i, s$study, s$merged_out,
                100 * s$merged_out / s$original_loci))
  }
}

# Show merge groups with >1 locus
merged_groups <- locus_list[group_size > 1][order(study, merge_group, rank_in_group)]
if (nrow(merged_groups) > 0) {
  cat("\nMerge groups (overlapping loci within", MERGE_DIST / 1000, "kb):\n")
  for (sg in unique(merged_groups[, paste(study, merge_group)])) {
    rows <- merged_groups[paste(study, merge_group) == sg]
    cat(sprintf("\n  %s | group %s | chr%s:%s-%s\n",
                rows$study[1], rows$merge_group[1], rows$chr[1],
                format(min(rows$bp), big.mark = ","),
                format(max(rows$bp), big.mark = ",")))
    for (j in 1:nrow(rows)) {
      r <- rows[j]
      tag <- if (r$kept) "[KEEP]" else "[DROP]"
      cat(sprintf("    %s locus %s (bp=%s, CS=%d, maxPIP=%.3f)\n",
                  tag, r$locus, format(r$bp, big.mark = ","), r$n_cs_variants, r$max_pip))
    }
  }
}

# ── Step 6: Write outputs ──────────────────────────────────────────────────
# Mapping table
map_out <- locus_list[, .(study, original_locus = locus, merged_locus, kept, reason,
                          chr, bp, n_cs_variants, max_pip, n_variants, converged,
                          merge_group, group_size, split_reason)]
fwrite(map_out, file.path(RESULTS_DIR, "merged_loci_map.csv"))
cat("\nWrote:", file.path(RESULTS_DIR, "merged_loci_map.csv"), "\n")

# Independent loci summary
summary_out <- study_stats[order(study)]
fwrite(summary_out, file.path(RESULTS_DIR, "independent_loci_summary.csv"))
cat("Wrote:", file.path(RESULTS_DIR, "independent_loci_summary.csv"), "\n")

cat("\n============================================================\n")
cat("Done. Use merged_loci_map.csv to filter downstream analyses.\n")
cat("============================================================\n")
