# 246b_or_equality_lookelsewhere.R
# Phase 2.5b — OR-equality bootstrap test + look-elsewhere correction.
#
# Dropping the misleading "F2→F3 has 4 orders of magnitude stronger COLOC
# enrichment than F1→F2" framing. The original p-value gap is dominated by the
# fact that the F2→F3 set is ~30% larger than F1→F2 (3,346 vs 2,595). With
# unmatched set sizes, p-values are not comparable. We do two things:
#
# 1) OR-equality bootstrap. Subsample the F2→F3 set to N=2,595 (matching F1→F2)
#    and bootstrap 1,000× with replacement. Report 95% CI for both the
#    F2→F3-subsampled OR and the OR-ratio (F2F3/F1F2). Verdict: pre-registered
#    pass = OR-ratio CI excludes 1.0 at PP4≥0.5.
#
# 2) Look-elsewhere correction. Enumerate the post-hoc claims of the pivot:
#    pathway-uniqueness × transitions, cell-type × transitions, NMF k=6
#    program × stage cells, sex × transition cells, COLOC × transition × PP4
#    threshold, sc bimodality + axis annotation. Total family size is fixed
#    pre-correction. Apply Bonferroni (α=0.05) and BH-FDR.
#
# Inputs:
#   results/granular_staging/substate_genetic_enrichment.csv  (existing OR table)
#   results/granular_staging/substate_genetic_top_genes.csv   (existing gene-level)
#   results/multi_evidence/multi_evidence_atlas.csv           (33,943 × 218 cols)
#   results/integration/merged_dge.rds, meta_matched.rds      (DGE input for resampling)
#
# Outputs (results/granular_staging/):
#   or_equality_bootstrap.csv
#   look_elsewhere_correction.csv
#   or_equality_lookelsewhere_summary.md

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

set.seed(42)

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/results/granular_staging")
ATLAS <- file.path(PROJECT_ROOT, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")

cat("== 246b OR-equality bootstrap + look-elsewhere correction ==\n")
cat(sprintf("Started: %s\n", Sys.time()))

# ---------------------------------------------------------------------------
# Load atlas + COLOC hit sets
# ---------------------------------------------------------------------------
atlas <- read.csv(ATLAS, stringsAsFactors = FALSE)
cat(sprintf("Atlas: %d genes × %d cols\n", nrow(atlas), ncol(atlas)))

coloc_cols <- intersect(c("coloc_best_susie_pp4", "coloc_pp4", "coloc_abf_best_pp4",
                          "broadaway_coloc_pp4", "ast_coloc_pp4", "ggt_coloc_pp4",
                          "pdff_coloc_pp4", "ukbb_alt_coloc_pp4",
                          "finngen_nafld_coloc_pp4", "finngen_nash_coloc_pp4",
                          "finngen_hcc_coloc_pp4", "bbj_alt_coloc_pp4"),
                        colnames(atlas))
cat("COLOC columns:", paste(coloc_cols, collapse = ", "), "\n")

coloc_mat <- as.matrix(atlas[, coloc_cols, drop = FALSE])
coloc_mat[is.na(coloc_mat)] <- 0
atlas$any_pp4_05 <- rowSums(coloc_mat >= 0.5) > 0
atlas$any_pp4_08 <- rowSums(coloc_mat >= 0.8) > 0
atlas$any_pp4_09 <- rowSums(coloc_mat >= 0.9) > 0

hit_sets <- list(
  "PP4>=0.5" = atlas$human_symbol[atlas$any_pp4_05],
  "PP4>=0.8" = atlas$human_symbol[atlas$any_pp4_08],
  "PP4>=0.9" = atlas$human_symbol[atlas$any_pp4_09]
)
for (cls in names(hit_sets)) {
  cat(sprintf("  %s: %d genes\n", cls, length(hit_sets[[cls]])))
}

# ---------------------------------------------------------------------------
# Re-derive transition gene sets (mirror 246)
# ---------------------------------------------------------------------------
atlas_sym <- atlas[, c("ensembl_id", "human_symbol")]
atlas_sym$ensembl_clean <- sub("\\..*", "", atlas_sym$ensembl_id)
ensg_to_sym <- function(ensg) {
  ensg_clean <- sub("\\..*", "", ensg)
  m <- atlas_sym$human_symbol[match(ensg_clean, atlas_sym$ensembl_clean)]
  ifelse(is.na(m) | m == "", ensg_clean, m)
}

dge <- readRDS(file.path(INT_DIR, "merged_dge.rds"))
meta <- readRDS(file.path(INT_DIR, "meta_matched.rds"))
meta <- meta[match(colnames(dge), meta$sample_id), ]

transition_contrast <- function(idx_a, idx_b) {
  group <- factor(c(rep("A", length(idx_a)), rep("B", length(idx_b))), levels = c("A", "B"))
  cohort <- factor(meta$dataset[c(idx_a, idx_b)])
  dge_sub <- dge[, c(idx_a, idx_b)]
  keep <- filterByExpr(dge_sub, group = group, min.count = 5)
  dge_sub <- dge_sub[keep, , keep.lib.sizes = FALSE]
  dge_sub <- calcNormFactors(dge_sub)
  design <- if (length(unique(cohort)) > 1) model.matrix(~ group + cohort) else model.matrix(~ group)
  v <- voom(dge_sub, design)
  fit <- eBayes(lmFit(v, design))
  res <- topTable(fit, coef = "groupB", number = Inf, sort.by = "none")
  res$gene_id <- rownames(res)
  res$gene_symbol <- ensg_to_sym(res$gene_id)
  res
}

f1 <- which(meta$fibrosis_stage == 1)
f2 <- which(meta$fibrosis_stage == 2)
f3 <- which(meta$fibrosis_stage == 3)

cat(sprintf("\nF1: %d, F2: %d, F3: %d samples\n", length(f1), length(f2), length(f3)))
cat("Computing F1→F2 contrast...\n"); res_f1f2 <- transition_contrast(f1, f2)
cat("Computing F2→F3 contrast...\n"); res_f2f3 <- transition_contrast(f2, f3)

# Same gene-set definition as 246: padj<0.05 (or top-500 fallback if fewer than 100)
extract_set <- function(res) {
  sig <- res[!is.na(res$adj.P.Val) & res$adj.P.Val < 0.05, ]
  if (nrow(sig) < 100) sig <- res[order(res$P.Value), ][1:500, ]
  unique(sig$gene_symbol[!is.na(sig$gene_symbol) & sig$gene_symbol != ""])
}

# Track per-gene padj for ranked subsampling
rank_set_with_padj <- function(res) {
  res <- res[!is.na(res$adj.P.Val) & res$adj.P.Val < 0.05, ]
  res <- res[!is.na(res$gene_symbol) & res$gene_symbol != "", ]
  res <- res[!duplicated(res$gene_symbol), ]
  res[order(res$P.Value), c("gene_symbol", "P.Value", "adj.P.Val")]
}

set_f1f2 <- extract_set(res_f1f2)
set_f2f3 <- extract_set(res_f2f3)
ranked_f2f3 <- rank_set_with_padj(res_f2f3)

cat(sprintf("\nF1→F2 set size: %d (target = match)\n", length(set_f1f2)))
cat(sprintf("F2→F3 set size: %d (will be subsampled to %d)\n", length(set_f2f3), length(set_f1f2)))

N_TARGET <- length(set_f1f2)  # subsample target = F1F2 set size

# ---------------------------------------------------------------------------
# Fisher utilities
# ---------------------------------------------------------------------------
fisher_or <- function(set_genes, hit_set, atlas_genes) {
  set_in <- atlas_genes %in% set_genes
  hit_in <- atlas_genes %in% hit_set
  set_hits <- sum(set_in & hit_in)
  set_nonhits <- sum(set_in & !hit_in)
  nonset_hits <- sum(!set_in & hit_in)
  nonset_nonhits <- sum(!set_in & !hit_in)
  tab <- matrix(c(set_hits, set_nonhits, nonset_hits, nonset_nonhits), nrow = 2)
  ft <- tryCatch(fisher.test(tab, alternative = "greater"), error = function(e) NULL)
  list(
    OR = if (!is.null(ft)) as.numeric(ft$estimate) else NA_real_,
    p = if (!is.null(ft)) ft$p.value else NA_real_
  )
}

atlas_genes <- atlas$human_symbol

# ---------------------------------------------------------------------------
# Analysis 1: OR-equality bootstrap (1,000 reps, with replacement)
# ---------------------------------------------------------------------------
N_BOOT <- 1000
boot_rows <- list()

for (cls in names(hit_sets)) {
  hit_set <- hit_sets[[cls]]
  cat(sprintf("\nBootstrap PP4 %s (N_TARGET=%d, N_BOOT=%d)\n", cls, N_TARGET, N_BOOT))

  # baseline OR for F1F2 (no resampling — the smaller set is the anchor)
  or_f1f2 <- fisher_or(set_f1f2, hit_set, atlas_genes)$OR

  # bootstrap F2F3 subsample (with replacement) to N_TARGET, then unique
  # since duplicates collapse in set membership: we sample with replacement,
  # take unique genes, and accept variable effective set size — this matches
  # the "with replacement" specification directly. Atlas-membership testing
  # is set-based, so duplicate draws don't double-count.
  or_boot <- numeric(N_BOOT)
  for (b in seq_len(N_BOOT)) {
    samp <- unique(sample(set_f2f3, size = N_TARGET, replace = TRUE))
    or_boot[b] <- fisher_or(samp, hit_set, atlas_genes)$OR
  }

  finite_or <- or_boot[is.finite(or_boot)]
  med_or <- median(finite_or)
  ci <- quantile(finite_or, probs = c(0.025, 0.975), na.rm = TRUE)

  # OR-ratio CI: bootstrap distribution of (OR_f2f3_b / OR_f1f2_const)
  ratio_dist <- finite_or / or_f1f2
  ratio_ci <- quantile(ratio_dist, probs = c(0.025, 0.975), na.rm = TRUE)

  pass <- ratio_ci[1] > 1.0  # CI excludes 1.0 from below

  boot_rows[[cls]] <- data.frame(
    pp4_threshold = cls,
    n_target = N_TARGET,
    n_boot = N_BOOT,
    or_f1f2 = or_f1f2,
    or_f2f3_subsampled_median = med_or,
    or_f2f3_ci_low = unname(ci[1]),
    or_f2f3_ci_high = unname(ci[2]),
    or_ratio_median = median(ratio_dist),
    or_ratio_ci_low = unname(ratio_ci[1]),
    or_ratio_ci_high = unname(ratio_ci[2]),
    pass_real_difference = pass,
    stringsAsFactors = FALSE
  )
  cat(sprintf("  OR_f1f2=%.3f, OR_f2f3_sub median=%.3f [%.3f, %.3f]\n",
              or_f1f2, med_or, ci[1], ci[2]))
  cat(sprintf("  OR-ratio median=%.3f [%.3f, %.3f]  pass=%s\n",
              median(ratio_dist), ratio_ci[1], ratio_ci[2], pass))
}

bootstrap_df <- do.call(rbind, boot_rows)
write.csv(bootstrap_df, file.path(OUT_DIR, "or_equality_bootstrap.csv"), row.names = FALSE)
cat("\nWrote or_equality_bootstrap.csv\n")

# ---------------------------------------------------------------------------
# Analysis 2: Look-elsewhere correction
# ---------------------------------------------------------------------------
# Build a claims registry. Each claim has uncorrected_p; we then apply
# Bonferroni and BH-FDR over the full family.
#
# Family composition (per plan):
#   - Pathway-uniqueness × transitions:
#       268 pathways × 4 transitions = 1,072 cells. We test only the
#       "uniquely up at this transition" claims, i.e. 14 unique-pathway
#       claims (10 F1→F2 + 4 F2→F3) — but the family is the discovery
#       space (1,072 hypothesis cells).
#   - Cell-type × transition: 22 cell types × 4 transitions = 88 tests
#   - NMF k=6 program × stage: 6 programs × 5 stages = 30 cells
#   - Sex × transition: 4 transitions × 2 sexes = 8 cells
#   - COLOC × transition × PP4 threshold: 4 transitions × 3 thresholds = 12
#   - sc bimodality + axis annotation: ~6 tests
# Per-plan total ≈ 144 tests. To avoid double-counting, we use the
# explicitly enumerated post-hoc claims listed in the plan (144 family).

# Read existing OR enrichment table to source COLOC × transition p-values
enrich <- read.csv(file.path(OUT_DIR, "substate_genetic_enrichment.csv"), stringsAsFactors = FALSE)

# Read pathway / cell-type / NMF / sex / sc tables to source p-values
read_tbl <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE)
}
pathway_tbl <- read_tbl(file.path(OUT_DIR, "two_transition_pathway_decomposition.csv"))
celltype_tbl <- read_tbl(file.path(OUT_DIR, "two_transition_celltype_by_stage.csv"))
nmf_tbl <- read_tbl(file.path(OUT_DIR, "two_transition_nmf_program_by_stage.csv"))
sex_tbl <- read_tbl(file.path(OUT_DIR, "two_transition_sex_stratified.csv"))
vu_tbl <- read_tbl(file.path(OUT_DIR, "vu_spot_bimodality.csv"))
sc_tbl <- read_tbl(file.path(OUT_DIR, "f3_substate_sc_donor_diagnostics.csv"))

# Helper: pull the relevant p-value column generously
extract_p <- function(df, p_cols = c("p_value", "padj", "P.Value", "p", "BH.padj", "BH_padj", "bh_padj")) {
  if (is.null(df)) return(NULL)
  for (col in p_cols) if (col %in% colnames(df)) return(df[[col]])
  NULL
}

claims <- list()

# 1) COLOC × transition × PP4 threshold (12) — anchored on existing OR table
for (i in seq_len(nrow(enrich))) {
  claims[[length(claims) + 1L]] <- data.frame(
    claim_id = sprintf("coloc_%s", gsub("[^A-Za-z0-9]+", "_", enrich$transition[i])),
    family = "COLOC_transition_PP4",
    description = sprintf("COLOC enrichment %s (OR=%.2f, set_size=%d)",
                          enrich$transition[i], enrich$OR[i], enrich$set_size[i]),
    uncorrected_p = enrich$p_value[i],
    stringsAsFactors = FALSE
  )
}

# 2) Pathway × transition (use raw padj from pathway decomposition) — count all rows
if (!is.null(pathway_tbl)) {
  pcol <- intersect(c("p_value", "pval", "P.Value", "p", "padj", "BH"), colnames(pathway_tbl))[1]
  if (!is.na(pcol) && !is.null(pcol)) {
    for (i in seq_len(nrow(pathway_tbl))) {
      pv <- pathway_tbl[[pcol]][i]
      if (is.na(pv)) next
      tag <- if ("transition" %in% colnames(pathway_tbl)) pathway_tbl$transition[i] else "X"
      pw <- if ("pathway" %in% colnames(pathway_tbl)) pathway_tbl$pathway[i] else
            if ("Term" %in% colnames(pathway_tbl)) pathway_tbl$Term[i] else paste0("p", i)
      claims[[length(claims) + 1L]] <- data.frame(
        claim_id = sprintf("pathway_%s_%s", tag, gsub("[^A-Za-z0-9]+", "_", pw)),
        family = "Pathway_transition",
        description = sprintf("Pathway %s @ %s", pw, tag),
        uncorrected_p = pv,
        stringsAsFactors = FALSE
      )
    }
  }
}

# 3) Cell-type × transition
if (!is.null(celltype_tbl)) {
  pcol <- intersect(c("p_value", "p", "P.Value", "padj", "BH"), colnames(celltype_tbl))[1]
  if (!is.na(pcol) && !is.null(pcol)) {
    for (i in seq_len(nrow(celltype_tbl))) {
      pv <- celltype_tbl[[pcol]][i]
      if (is.na(pv)) next
      ct <- if ("celltype" %in% colnames(celltype_tbl)) celltype_tbl$celltype[i] else
            if ("cell_type" %in% colnames(celltype_tbl)) celltype_tbl$cell_type[i] else paste0("c", i)
      tag <- if ("transition" %in% colnames(celltype_tbl)) celltype_tbl$transition[i] else "X"
      claims[[length(claims) + 1L]] <- data.frame(
        claim_id = sprintf("celltype_%s_%s", tag, gsub("[^A-Za-z0-9]+", "_", ct)),
        family = "Celltype_transition",
        description = sprintf("Cell-type %s @ %s", ct, tag),
        uncorrected_p = pv,
        stringsAsFactors = FALSE
      )
    }
  }
}

# 4) NMF program × stage
if (!is.null(nmf_tbl)) {
  pcol <- intersect(c("p_value", "p", "P.Value", "padj", "BH"), colnames(nmf_tbl))[1]
  if (!is.na(pcol) && !is.null(pcol)) {
    for (i in seq_len(nrow(nmf_tbl))) {
      pv <- nmf_tbl[[pcol]][i]
      if (is.na(pv)) next
      pg <- if ("program" %in% colnames(nmf_tbl)) nmf_tbl$program[i] else paste0("P", i)
      tag <- if ("stage" %in% colnames(nmf_tbl)) nmf_tbl$stage[i] else "X"
      claims[[length(claims) + 1L]] <- data.frame(
        claim_id = sprintf("nmf_%s_%s", pg, tag),
        family = "NMF_program_stage",
        description = sprintf("NMF %s @ %s", pg, tag),
        uncorrected_p = pv,
        stringsAsFactors = FALSE
      )
    }
  }
}

# 5) Sex × transition
if (!is.null(sex_tbl)) {
  pcol <- intersect(c("p_value", "p", "P.Value", "padj", "BH"), colnames(sex_tbl))[1]
  if (!is.na(pcol) && !is.null(pcol)) {
    for (i in seq_len(nrow(sex_tbl))) {
      pv <- sex_tbl[[pcol]][i]
      if (is.na(pv)) next
      sx <- if ("sex" %in% colnames(sex_tbl)) sex_tbl$sex[i] else paste0("s", i)
      tag <- if ("transition" %in% colnames(sex_tbl)) sex_tbl$transition[i] else "X"
      claims[[length(claims) + 1L]] <- data.frame(
        claim_id = sprintf("sex_%s_%s", tag, sx),
        family = "Sex_transition",
        description = sprintf("Sex %s @ %s", sx, tag),
        uncorrected_p = pv,
        stringsAsFactors = FALSE
      )
    }
  }
}

# 6) sc bimodality + axis annotation (a small set, ~6 claims)
if (!is.null(sc_tbl)) {
  for (i in seq_len(nrow(sc_tbl))) {
    test_name <- if ("test" %in% colnames(sc_tbl)) sc_tbl$test[i]
                 else if ("metric" %in% colnames(sc_tbl)) sc_tbl$metric[i]
                 else paste0("sc_", i)
    pv <- NA_real_
    for (col in c("p_value", "p", "pval", "BH_padj", "padj")) {
      if (col %in% colnames(sc_tbl)) {
        pv <- sc_tbl[[col]][i]
        if (!is.na(pv)) break
      }
    }
    if (is.na(pv)) next
    claims[[length(claims) + 1L]] <- data.frame(
      claim_id = sprintf("sc_%s", gsub("[^A-Za-z0-9]+", "_", test_name)),
      family = "sc_bimodality",
      description = sprintf("sc test %s", test_name),
      uncorrected_p = pv,
      stringsAsFactors = FALSE
    )
  }
}

# Vu spatial bimodality if available
if (!is.null(vu_tbl)) {
  for (i in seq_len(nrow(vu_tbl))) {
    test_name <- if ("test" %in% colnames(vu_tbl)) vu_tbl$test[i]
                 else if ("axis" %in% colnames(vu_tbl)) vu_tbl$axis[i]
                 else paste0("vu_", i)
    pv <- NA_real_
    for (col in c("p_value", "p", "pval", "dip_mc_p", "BH_padj", "padj")) {
      if (col %in% colnames(vu_tbl)) {
        pv <- vu_tbl[[col]][i]
        if (!is.na(pv)) break
      }
    }
    if (is.na(pv)) next
    claims[[length(claims) + 1L]] <- data.frame(
      claim_id = sprintf("vu_%s", gsub("[^A-Za-z0-9]+", "_", test_name)),
      family = "vu_spatial_bimodality",
      description = sprintf("Vu spatial bimodality %s", test_name),
      uncorrected_p = pv,
      stringsAsFactors = FALSE
    )
  }
}

claims_df <- do.call(rbind, claims)

# Family-size-aware correction:
# Use the pre-registered total family size from the plan (≈144) as a
# CONSERVATIVE Bonferroni baseline (because the actual collected list may be
# smaller than the planned discovery family — applying the planned family to
# all collected claims is the right "look-elsewhere" answer).
PLANNED_FAMILY_SIZE <- 144L

n_observed <- nrow(claims_df)
m_for_bonferroni <- max(n_observed, PLANNED_FAMILY_SIZE)
cat(sprintf("\nLook-elsewhere family: observed = %d, planned = %d, used m = %d\n",
            n_observed, PLANNED_FAMILY_SIZE, m_for_bonferroni))

claims_df$bonferroni_p <- pmin(claims_df$uncorrected_p * m_for_bonferroni, 1.0)
claims_df$bh_fdr_p <- p.adjust(claims_df$uncorrected_p, method = "BH")
# pre-registered Bonferroni alpha = 0.05 → corrected p-threshold = 0.05/144 ≈ 3.47e-4
BONF_ALPHA_PER_TEST <- 0.05 / m_for_bonferroni
claims_df$survives_bonferroni <- claims_df$uncorrected_p < BONF_ALPHA_PER_TEST
claims_df$survives_bh_fdr <- claims_df$bh_fdr_p < 0.05

write.csv(claims_df, file.path(OUT_DIR, "look_elsewhere_correction.csv"), row.names = FALSE)
cat(sprintf("Wrote look_elsewhere_correction.csv (%d claims)\n", nrow(claims_df)))

# ---------------------------------------------------------------------------
# Markdown summary
# ---------------------------------------------------------------------------
sink(file.path(OUT_DIR, "or_equality_lookelsewhere_summary.md"))
cat("# 246b OR-equality bootstrap + look-elsewhere correction\n\n")
cat(sprintf("Date: %s\n\n", Sys.Date()))
cat("## Verdict\n\n")

# Pre-registered acceptance #1: PP4≥0.5 OR-ratio CI excludes 1.0
pp05_row <- bootstrap_df[bootstrap_df$pp4_threshold == "PP4>=0.5", ]
v1 <- pp05_row$pass_real_difference
cat(sprintf("**Pre-registered acceptance #1 (F2→F3 stronger COLOC than F1→F2 with set-size matching):**\n\n"))
if (v1) {
  cat(sprintf("PASS — at PP4≥0.5, OR-ratio 95%% CI = [%.3f, %.3f] excludes 1.0.\n",
              pp05_row$or_ratio_ci_low, pp05_row$or_ratio_ci_high))
} else {
  cat(sprintf("FAIL — at PP4≥0.5, OR-ratio 95%% CI = [%.3f, %.3f] does NOT exclude 1.0. ",
              pp05_row$or_ratio_ci_low, pp05_row$or_ratio_ci_high))
  cat("The original p-value gap of '4 orders of magnitude' is dominated by set-size asymmetry, not a real OR difference. The two transitions have indistinguishable COLOC enrichment when matched.\n")
}
cat("\n")

cat("## OR-equality bootstrap (1,000 reps, F2→F3 subsampled to F1→F2 size)\n\n")
cat("| PP4 threshold | OR(F1→F2) | OR(F2→F3 sub) median [95% CI] | OR-ratio median [95% CI] | Pass |\n")
cat("|---|---:|---|---|:---:|\n")
for (i in seq_len(nrow(bootstrap_df))) {
  r <- bootstrap_df[i, ]
  cat(sprintf("| %s | %.3f | %.3f [%.3f, %.3f] | %.3f [%.3f, %.3f] | %s |\n",
              r$pp4_threshold, r$or_f1f2,
              r$or_f2f3_subsampled_median, r$or_f2f3_ci_low, r$or_f2f3_ci_high,
              r$or_ratio_median, r$or_ratio_ci_low, r$or_ratio_ci_high,
              ifelse(r$pass_real_difference, "PASS", "FAIL")))
}
cat("\n")

cat("## Look-elsewhere correction\n\n")
cat(sprintf("Claims collected: %d (pre-registered family size %d → Bonferroni α/test = %.2e)\n\n",
            nrow(claims_df), PLANNED_FAMILY_SIZE, BONF_ALPHA_PER_TEST))

n_bonf <- sum(claims_df$survives_bonferroni, na.rm = TRUE)
n_bh <- sum(claims_df$survives_bh_fdr, na.rm = TRUE)
cat(sprintf("- Bonferroni-surviving claims: **%d / %d**\n", n_bonf, nrow(claims_df)))
cat(sprintf("- BH-FDR-surviving claims: **%d / %d**\n\n", n_bh, nrow(claims_df)))

cat("### Survives Bonferroni (uncorrected p < 0.05/m)\n\n")
surv <- claims_df[claims_df$survives_bonferroni & !is.na(claims_df$survives_bonferroni), ]
if (nrow(surv) > 0) {
  surv <- surv[order(surv$uncorrected_p), ]
  cat("| Claim | Family | Uncorrected p | Bonferroni p | BH FDR |\n")
  cat("|---|---|---:|---:|---:|\n")
  for (i in seq_len(min(50, nrow(surv)))) {
    r <- surv[i, ]
    cat(sprintf("| %s | %s | %.2e | %.2e | %.2e |\n",
                substr(r$description, 1, 70), r$family, r$uncorrected_p,
                r$bonferroni_p, r$bh_fdr_p))
  }
  if (nrow(surv) > 50) cat(sprintf("\n... %d more in CSV.\n", nrow(surv) - 50))
} else {
  cat("None.\n")
}
cat("\n")

cat("### Family-level breakdown\n\n")
fam_breakdown <- aggregate(survives_bonferroni ~ family, data = claims_df,
                           FUN = function(x) c(n = length(x), n_bonf = sum(x, na.rm = TRUE)))
fam_summary <- data.frame(
  family = fam_breakdown$family,
  n_total = fam_breakdown$survives_bonferroni[, "n"],
  n_bonf = fam_breakdown$survives_bonferroni[, "n_bonf"]
)
cat("| Family | Total | Survives Bonferroni |\n|---|---:|---:|\n")
for (i in seq_len(nrow(fam_summary))) {
  cat(sprintf("| %s | %d | %d |\n", fam_summary$family[i],
              fam_summary$n_total[i], fam_summary$n_bonf[i]))
}
cat("\n")

cat("## Pre-registered claim outcomes\n\n")
cat("**Headline COLOC claim**: 'F2→F3 has 4 orders of magnitude stronger COLOC enrichment than F1→F2'.\n")
if (v1) {
  cat(sprintf("- Set-size matched OR-equality bootstrap **PASSES** at PP4≥0.5 (CI [%.3f, %.3f]).\n",
              pp05_row$or_ratio_ci_low, pp05_row$or_ratio_ci_high))
} else {
  cat(sprintf("- Set-size matched OR-equality bootstrap **FAILS** at PP4≥0.5 (CI [%.3f, %.3f]). RETIRE this framing.\n",
              pp05_row$or_ratio_ci_low, pp05_row$or_ratio_ci_high))
}
# Bonferroni for COLOC claims under family ≈ 144
coloc_rows <- claims_df[claims_df$family == "COLOC_transition_PP4", ]
n_coloc_surv <- sum(coloc_rows$survives_bonferroni, na.rm = TRUE)
cat(sprintf("- COLOC × transition × PP4 cells surviving Bonferroni: %d / %d (corrected α = %.2e).\n",
            n_coloc_surv, nrow(coloc_rows), BONF_ALPHA_PER_TEST))
sink()

cat("\n== 246b complete. ==\n")
cat(sprintf("Finished: %s\n", Sys.time()))
