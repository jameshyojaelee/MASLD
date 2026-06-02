#!/usr/bin/env Rscript
# 30_spatial_ccc_consensus.R
#
# Analysis I1 (v1) — Multi-method spatial CCC consensus.
#
# Strategy (v1, uses existing COMMOT + squidpy ligrec + LIANA scRNA):
#   1. Load COMMOT differential pathway communication (Healthy vs Steatotic).
#   2. Load squidpy differential LR pairs.
#   3. Load scRNA LIANA differential LR pairs (agg across cell-type pairs).
#   4. Cross-reference: LR pairs appearing in >=2 methods are "consensus".
#   5. Prioritize: scRNA MASLD-enriched + spatial disease-emergent + spatial
#      neighborhood-enriched = triple validated.
#
# Env: rnaseq
# Outputs: Analysis/Spatial/results/spatial_ccc_multimethod/

suppressPackageStartupMessages({library(data.table)})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUTDIR <- file.path(BASE, "Analysis/Spatial/results/spatial_ccc_multimethod")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

message("[1] Loading spatial squidpy differential LR pairs...")
sq <- fread(file.path(BASE, "Analysis/Spatial/results/communication/differential_lr_pairs.csv"))
# Parse LR tuple
sq[, ligand   := sub("^\\('?([^',]+)'?,.*", "\\1", lr_pair)]
sq[, receptor := sub(".*,\\s*'?([^')]+)'?\\)$", "\\1", lr_pair)]
sq[, ligand := trimws(ligand)]
sq[, receptor := trimws(receptor)]
message(sprintf("  squidpy LR pairs: %d", nrow(sq)))

message("[2] Loading COMMOT pathway communication shifts...")
commot <- fread(file.path(BASE, "Analysis/Spatial/results/commot/differential_communication.csv"))
message(sprintf("  COMMOT pathways: %d", nrow(commot)))

message("[3] Loading scRNA LIANA differential LR pairs (aggregated)...")
liana <- fread(file.path(BASE, "Analysis/SingleCell/results_gpu_v2/fig2_data/liana_differential_interactions.csv"))
# F114: LIANA multi-subunit complexes use underscore-joined subunits in
# ligand_complex/receptor_complex (e.g. receptor_complex='ITGAL_ITGB2'), whereas
# squidpy lr_pair tuples are single genes. A raw-string merge silently drops
# every heteromeric pair. Expand each complex into its constituent subunits so a
# single-gene squidpy ligand/receptor can match one subunit; retain the original
# complex name for provenance.
n_lig_cplx <- sum(grepl("_", liana$ligand_complex))
n_rec_cplx <- sum(grepl("_", liana$receptor_complex))
message(sprintf("  LIANA rows carrying a complex (ligand|receptor): %d | %d (of %d)",
                n_lig_cplx, n_rec_cplx, nrow(liana)))
liana[, abs_sd := abs(score_diff)]
liana[, row_id := .I]
# Explode underscore complexes on both sides into one row per subunit.
# Step 1: one row per ligand subunit (carry all metadata via the grouping key).
liana_exp <- liana[, .(ligand = unlist(strsplit(ligand_complex, "_", fixed = TRUE))),
                   by = .(row_id, receptor_complex, score_diff, abs_sd,
                          source, target, ligand_complex)]
# Step 2: one row per receptor subunit.
liana_exp <- liana_exp[, .(receptor = unlist(strsplit(receptor_complex, "_", fixed = TRUE))),
                       by = .(row_id, ligand, score_diff, abs_sd, source, target,
                              ligand_complex, receptor_complex)]
# Aggregate per (subunit ligand, subunit receptor) — keep max |score_diff|.
liana_exp <- liana_exp[order(-abs_sd)]
liana_agg <- liana_exp[!duplicated(paste(ligand, receptor))]
liana_agg <- liana_agg[, .(ligand, receptor,
                            liana_score_diff = score_diff, liana_source = source,
                            liana_target = target,
                            liana_ligand_complex = ligand_complex,
                            liana_receptor_complex = receptor_complex)]
message(sprintf("  LIANA unique LR pairs (subunit-expanded): %d", nrow(liana_agg)))

message("[4] Merging: LR pair consensus across methods...")
consensus <- merge(liana_agg, sq[, .(ligand, receptor,
                                      sq_category = category,
                                      sq_delta_expr = delta_mean_expr,
                                      sq_source = source, sq_target = target,
                                      sq_pair_type = pair_type)],
                   by = c("ligand","receptor"), all = FALSE)
message(sprintf("  LIANA ∩ squidpy spatial: %d LR pairs", nrow(consensus)))

# Annotate
# F113/F118: three-level direction. score_diff == 0 carries no direction (no
# differential), so it is 'no_change' and excluded from concordance — assigning
# it MASLD_down (the old fifelse) inflated one margin. The 'both_concordant'
# rate is NOT a standalone validation rate: LIANA score_diff is overwhelmingly
# positive and squidpy delta_mean_expr overwhelmingly negative, so the two
# marginals are skewed in OPPOSITE directions and concordance is suppressed by
# construction. We therefore benchmark the observed concordant fraction against
# a label-permutation null (shuffle scRNA_direction within the consensus set)
# and report both margins so the asymmetry is visible.
consensus[, scRNA_direction := fifelse(liana_score_diff > 0, "MASLD_up",
                                fifelse(liana_score_diff < 0, "MASLD_down", "no_change"))]
consensus[, sq_direction := fifelse(sq_delta_expr > 0, "spatial_up",
                            fifelse(sq_delta_expr < 0, "spatial_down", "no_change"))]
consensus[, both_concordant := (scRNA_direction == "MASLD_up" & sq_delta_expr > 0) |
                                (scRNA_direction == "MASLD_down" & sq_delta_expr < 0)]

fwrite(consensus, file.path(OUTDIR, "scRNA_spatial_LR_consensus.csv"))

# Summary stats
n_consensus <- nrow(consensus)
n_concordant <- sum(consensus$both_concordant, na.rm = TRUE)

# F113: permutation null for the concordant fraction. Restrict to directionally
# defined pairs (drop 'no_change' on either axis), then shuffle scRNA_direction
# labels and recompute the concordant fraction; the empirical p-value tells
# whether observed concordance exceeds chance GIVEN the marginal skews.
dir_ok <- consensus[scRNA_direction != "no_change" & sq_direction != "no_change"]
n_dir <- nrow(dir_ok)
obs_frac <- if (n_dir > 0) mean(dir_ok$both_concordant) else NA_real_
set.seed(42)
N_PERM <- 10000
perm_fracs <- if (n_dir > 0) {
  sq_pos <- dir_ok$sq_delta_expr > 0
  sc_dir <- dir_ok$scRNA_direction
  vapply(seq_len(N_PERM), function(i) {
    sh <- sample(sc_dir)
    mean((sh == "MASLD_up" & sq_pos) | (sh == "MASLD_down" & !sq_pos))
  }, numeric(1))
} else numeric(0)
perm_p <- if (n_dir > 0) (sum(perm_fracs >= obs_frac) + 1) / (N_PERM + 1) else NA_real_
perm_mean <- if (length(perm_fracs)) mean(perm_fracs) else NA_real_

# Marginal direction counts (the source of the structural asymmetry).
sc_up <- sum(consensus$scRNA_direction == "MASLD_up")
sc_dn <- sum(consensus$scRNA_direction == "MASLD_down")
sc_nc <- sum(consensus$scRNA_direction == "no_change")
sq_up <- sum(consensus$sq_direction == "spatial_up")
sq_dn <- sum(consensus$sq_direction == "spatial_down")

# COMMOT: map pathway-level MASLD-increased communication
commot_up <- commot[total_communication_delta > 0][order(-total_communication_delta)]
commot_dn <- commot[total_communication_delta < 0][order(total_communication_delta)]
fwrite(commot_up[, .(pathway, total_communication_delta, total_communication_fc)],
       file.path(OUTDIR, "commot_masld_up_pathways.csv"))
fwrite(commot_dn[, .(pathway, total_communication_delta, total_communication_fc)],
       file.path(OUTDIR, "commot_masld_down_pathways.csv"))

summary_lines <- c(
  sprintf("LIANA scRNA LR pairs: %d", nrow(liana_agg)),
  sprintf("squidpy spatial differential LR pairs: %d", nrow(sq)),
  sprintf("COMMOT pathways: %d", nrow(commot)),
  "",
  sprintf("LR pairs in BOTH scRNA (LIANA) + spatial (squidpy): %d", n_consensus),
  sprintf("  Both-concordant (same direction scRNA ↔ spatial): %d (%.1f%%)",
          n_concordant, 100 * n_concordant / max(1, n_consensus)),
  "",
  "MARGINAL DIRECTION COUNTS (F113 — the two axes are skewed in OPPOSITE",
  "directions, so the raw concordance fraction is structurally suppressed and",
  "is NOT a standalone validation rate):",
  sprintf("  scRNA LIANA direction: MASLD_up=%d, MASLD_down=%d, no_change=%d",
          sc_up, sc_dn, sc_nc),
  sprintf("  spatial squidpy direction: up=%d, down=%d", sq_up, sq_dn),
  sprintf("  Directionally-defined consensus pairs (no_change dropped): %d", n_dir),
  sprintf("  Concordant fraction among defined pairs: %.1f%% (observed)",
          100 * (if (is.na(obs_frac)) 0 else obs_frac)),
  sprintf("  Permutation null (shuffle scRNA labels, %d perms): mean=%.1f%%, empirical p=%.4g",
          N_PERM, 100 * (if (is.na(perm_mean)) 0 else perm_mean), perm_p),
  "",
  "CAVEAT (T1.14, 2026-04-22): COMMOT pathway-level LR detection operates on",
  "Visium 55 μm spots, which capture multiple cells per spot and cannot",
  "resolve autocrine vs paracrine signaling. Top spatial CCC pathways (e.g.",
  "FGF21) should be interpreted as aggregate spot-level communication; single-",
  "cell CCC validation (LIANA/CellChat on scRNA) is required to confirm the",
  "cell-type pairs driving the signal. TODO: random-spot permutation null.",
  "",
  "=== Top 20 COMMOT MASLD-upregulated pathways ===",
  capture.output(print(commot_up[, .(pathway, total_communication_delta, total_communication_fc)][1:20])),
  "",
  "=== Top 15 scRNA + spatial consensus LR pairs ===",
  capture.output(print(consensus[order(-abs(liana_score_diff))][1:15,
                       .(ligand, receptor, liana_source, liana_target,
                         liana_score_diff, sq_source, sq_target,
                         sq_category, sq_delta_expr, both_concordant)],
                       nrows = 15))
)
writeLines(summary_lines, file.path(OUTDIR, "spatial_ccc_consensus_summary.txt"))
writeLines(summary_lines)

message("Done. Outputs in: ", OUTDIR)
