#!/usr/bin/env Rscript
# 242_resilient_de.R
#
# Healthy-control characterization pipeline, step 3 of 6.
#
# 3-way differential expression: resilient (obese-no-MASLD) vs lean-healthy
# vs obese-MASLD.
#
# Two analyses, written to disk in parallel:
#   (A) Primary, within-GSE126848 (Suppli 2019). 4-arm cohort by design:
#       Control (lean, n=14), Control_Obese (resilient, n=12),
#       NAFL (n=15), NASH (n=16). NO cohort confounding.
#   (B) Cross-cohort (resilient_combined ∪ lean_healthy_explicit ∪ obese_MASLD)
#       with ~ cohort + sex + age + group3 design. Lean_healthy is
#       GSE126848-only ⇒ cohort × group is partially aliased; reported as
#       sensitivity, not primary.
#
# Per-cohort leave-one-out is run on (B) only.
#
# Spec: docs/superpowers/specs/2026-04-27-healthy-control-audit-design.md
# Env: rnaseq

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})

BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
                     "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
HCDIR  <- file.path(BASE, "RNA-seq/results/audit_sensitivity/healthy_control_audit")
DEDIR  <- file.path(HCDIR, "resilience_de")
INT    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
DGE_RDS <- file.path(INT, "merged_dge.rds")

dir.create(DEDIR, recursive = TRUE, showWarnings = FALSE)

cat("Loading controls_definitions.csv...\n")
defs <- fread(file.path(HCDIR, "controls_definitions.csv"))

cat("Loading merged_dge.rds...\n")
dge_full <- readRDS(DGE_RDS)
cat(sprintf("  Counts shape: %d genes x %d samples\n",
            nrow(dge_full$counts), ncol(dge_full$counts)))
cat(sprintf("  Sample columns: %s\n", paste(head(colnames(dge_full$samples), 6), collapse = ", ")))

# Build the 3-way factor for analysis B (cross-cohort)
defs_B <- defs[is_resilient | is_lean_healthy | is_obese_MASLD,
                .(sample_id, dataset, sex, age, condition,
                  group3 = fcase(
                    is_resilient,    "resilient",
                    is_lean_healthy, "lean_healthy",
                    is_obese_MASLD,  "obese_MASLD"))]
defs_B <- defs_B[!is.na(sample_id)]
cat(sprintf("\n[Analysis B] Cross-cohort 3-way pool: %d subjects\n", nrow(defs_B)))
print(table(defs_B$group3, defs_B$dataset))

# --- Analysis (A): Within-GSE126848 ---
defs_A <- defs[dataset == "GSE126848" &
                 condition %in% c("Control", "Control_Obese", "NAFL", "NASH"),
                .(sample_id, sex, age, condition,
                  group3 = fcase(
                    condition == "Control",       "lean_healthy",
                    condition == "Control_Obese", "resilient",
                    condition %in% c("NAFL", "NASH"), "obese_MASLD"))]
defs_A <- defs_A[!is.na(sample_id)]
cat(sprintf("\n[Analysis A] Within-GSE126848 4-arm pool: %d subjects\n", nrow(defs_A)))
print(table(defs_A$group3, defs_A$condition))

# Helper: run limma-voom DE for a defined sample subset
run_de <- function(sample_meta, design_formula, contrast_spec, name_tag,
                    min_count = 10, min_total = 15) {
  # sample_meta has: sample_id, group3, sex, age, dataset (optional)
  ids <- intersect(sample_meta$sample_id, colnames(dge_full$counts))
  if (length(ids) < length(sample_meta$sample_id)) {
    cat(sprintf("  [warn] %d / %d sample_ids not in counts matrix; dropped\n",
                length(sample_meta$sample_id) - length(ids), length(sample_meta$sample_id)))
  }
  if (length(ids) < 10) {
    cat(sprintf("  [%s] Insufficient samples (%d); skipping.\n", name_tag, length(ids)))
    return(NULL)
  }
  sample_meta <- sample_meta[sample_id %in% ids]
  setkey(sample_meta, sample_id)
  sample_meta <- sample_meta[ids]  # reorder

  cnts <- dge_full$counts[, ids, drop = FALSE]
  dge <- DGEList(counts = cnts)

  # Filter low-count genes by group
  keep <- filterByExpr(dge, group = sample_meta$group3,
                        min.count = min_count, min.total.count = min_total)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)

  # Build design
  group <- factor(sample_meta$group3, levels = c("lean_healthy", "resilient", "obese_MASLD"))
  sex   <- factor(sample_meta$sex, levels = c("F", "M"))
  age   <- as.numeric(sample_meta$age)
  cohort <- if ("dataset" %in% names(sample_meta)) factor(sample_meta$dataset) else NULL

  # Robust age handling: detect whether age is usable
  age_usable <- sum(!is.na(age)) >= 0.7 * length(age) &&
                isTRUE(sd(age, na.rm = TRUE) > 0)
  if (age_usable) {
    age[is.na(age)] <- median(age, na.rm = TRUE)
  }
  # Drop subjects with missing sex/age covariates if used
  sex_usable <- nlevels(droplevels(sex)) > 1
  drop_idx <- rep(FALSE, length(group))
  if (sex_usable) drop_idx <- drop_idx | is.na(sex)
  if (age_usable) drop_idx <- drop_idx | is.na(age)
  if (any(drop_idx)) {
    cat(sprintf("  [%s] Dropping %d subjects with NA covariates\n",
                name_tag, sum(drop_idx)))
    sample_meta <- sample_meta[!drop_idx]
    group <- group[!drop_idx]
    sex   <- droplevels(sex[!drop_idx])
    age   <- age[!drop_idx]
    if (!is.null(cohort)) cohort <- droplevels(cohort[!drop_idx])
    dge   <- dge[, !drop_idx]
    sex_usable <- nlevels(sex) > 1
  }

  # Build no-intercept design ~ 0 + group [+ sex] [+ age] [+ cohort]
  # so that contrasts can be named directly as "groupresilient - grouplean_healthy"
  if (design_formula == "cohort_sex_age_group3") {
    cohort_usable <- !is.null(cohort) && nlevels(cohort) > 1
    rhs <- "0 + group"
    if (sex_usable) rhs <- paste(rhs, "+ sex")
    if (age_usable) rhs <- paste(rhs, "+ age")
    if (cohort_usable) rhs <- paste(rhs, "+ cohort")
    design <- model.matrix(as.formula(paste("~", rhs)))
  } else if (design_formula == "sex_age_group3") {
    rhs <- "0 + group"
    if (sex_usable) rhs <- paste(rhs, "+ sex")
    if (age_usable) rhs <- paste(rhs, "+ age")
    design <- model.matrix(as.formula(paste("~", rhs)))
  } else stop("Unknown design")

  cat(sprintf("  [%s] design rank: %d, samples: %d, genes after filter: %d\n",
              name_tag, qr(design)$rank, ncol(dge), nrow(dge)))

  v <- voom(dge, design)
  fit <- lmFit(v, design)
  # Contrasts
  if (length(contrast_spec) > 0) {
    contr <- makeContrasts(contrasts = contrast_spec, levels = design)
    fit2  <- contrasts.fit(fit, contr)
    fit2  <- eBayes(fit2)
    rs <- list()
    for (cn in colnames(contr)) {
      tt <- topTable(fit2, coef = cn, number = Inf, sort.by = "none")
      tt <- as.data.table(tt, keep.rownames = "gene_id")
      tt[, contrast := cn]
      tt[, analysis := name_tag]
      rs[[cn]] <- tt
    }
    pairs <- rbindlist(rs, fill = TRUE)
    # F-test (across all contrasts)
    ftest <- topTable(fit2, number = Inf, sort.by = "none")
    ftest <- as.data.table(ftest, keep.rownames = "gene_id")
    ftest[, analysis := name_tag]
    return(list(pairs = pairs, ftest = ftest, n_samples = ncol(dge),
                n_genes = nrow(dge)))
  }
  return(NULL)
}

# === Analysis A: within-GSE126848 ===
cat("\n=== Analysis (A): Within-GSE126848 4-arm DE ===\n")
contrast_A <- c(
  "groupresilient - grouplean_healthy",
  "groupresilient - groupobese_MASLD"
)
res_A <- run_de(defs_A, "sex_age_group3", contrast_A, "GSE126848_within_cohort")
if (!is.null(res_A)) {
  fwrite(res_A$pairs, file.path(DEDIR, "resilient_de_pairwise_GSE126848.csv"))
  fwrite(res_A$ftest,  file.path(DEDIR, "resilient_de_ftest_GSE126848.csv"))
  cat(sprintf("  Wrote: resilient_de_pairwise_GSE126848.csv (%d rows)\n", nrow(res_A$pairs)))
  cat("  Top genes by contrast (FDR<0.05):\n")
  for (cn in unique(res_A$pairs$contrast)) {
    n_sig <- sum(res_A$pairs$contrast == cn & res_A$pairs$adj.P.Val < 0.05, na.rm = TRUE)
    cat(sprintf("    %s: %d FDR<0.05 / %d genes\n", cn, n_sig, sum(res_A$pairs$contrast == cn)))
  }
}

# === Analysis B: Cross-cohort ===
cat("\n=== Analysis (B): Cross-cohort 3-way DE (cohort + sex + age + group3) ===\n")
contrast_B <- c(
  "groupresilient - grouplean_healthy",
  "groupresilient - groupobese_MASLD"
)
res_B <- run_de(defs_B, "cohort_sex_age_group3", contrast_B, "cross_cohort")
if (!is.null(res_B)) {
  fwrite(res_B$pairs, file.path(DEDIR, "resilient_de_pairwise_crosscohort.csv"))
  fwrite(res_B$ftest,  file.path(DEDIR, "resilient_de_ftest_crosscohort.csv"))
  cat(sprintf("  Wrote: resilient_de_pairwise_crosscohort.csv (%d rows)\n", nrow(res_B$pairs)))
  for (cn in unique(res_B$pairs$contrast)) {
    n_sig <- sum(res_B$pairs$contrast == cn & res_B$pairs$adj.P.Val < 0.05, na.rm = TRUE)
    cat(sprintf("    %s: %d FDR<0.05 / %d genes\n", cn, n_sig, sum(res_B$pairs$contrast == cn)))
  }
}

# === Per-cohort leave-one-out (Analysis B only) ===
cat("\n=== Per-cohort LOO on cross-cohort design ===\n")
cohorts_in_B <- unique(defs_B$dataset)
loo_list <- list()
for (cohort_drop in cohorts_in_B) {
  defs_loo <- defs_B[dataset != cohort_drop]
  ng <- table(defs_loo$group3)
  if (any(ng < 5)) {
    cat(sprintf("  [skip] LOO %s: under-represented group (min n=%d)\n",
                cohort_drop, min(ng)))
    next
  }
  # Adapt contrasts: skip lean_healthy contrasts if that group is empty
  contrast_loo <- if ("lean_healthy" %in% names(ng) && ng["lean_healthy"] >= 5)
    contrast_B
  else {
    cat(sprintf("  [%s] no lean_healthy in remainder; only resilient vs obese_MASLD\n",
                cohort_drop))
    "groupresilient - groupobese_MASLD"
  }
  res_loo <- tryCatch(
    run_de(defs_loo, "cohort_sex_age_group3", contrast_loo,
           sprintf("LOO_%s", cohort_drop)),
    error = function(e) {
      cat(sprintf("  [LOO_%s] error: %s\n", cohort_drop, conditionMessage(e)))
      NULL
    }
  )
  if (!is.null(res_loo)) {
    res_loo$pairs[, dropped_cohort := cohort_drop]
    loo_list[[cohort_drop]] <- res_loo$pairs
  }
}
if (length(loo_list) > 0) {
  loo_all <- rbindlist(loo_list, use.names = TRUE)
  fwrite(loo_all, file.path(DEDIR, "resilient_de_loo.csv"))
  cat(sprintf("  Wrote: resilient_de_loo.csv (%d rows across %d LOO iterations)\n",
              nrow(loo_all), length(loo_list)))
}

cat("\nDone.\n")
