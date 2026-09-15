#!/usr/bin/env Rscript
# KEY MESSAGE: The F0 reference in the adjacent-fibrosis contrasts is a mixture of
# source controls and F0 disease biopsies, and the mixture is nested in cohort, so
# the cohort term cannot absorb it. This runs the same contrasts three ways.
#
# The F0 group is 104 participants:
#
#   GSE130970   23 controls,  0 disease     (fibrosis_stage == 0 IS the control label)
#   GSE135251    8 controls, 38 disease
#   GSE162694    0 controls, 35 disease     (its controls carry fibrosis_stage = NA
#                                            and are excluded upstream)
#
# So F0 -> F1 is a control-versus-disease contrast inside GSE130970 and a
# disease-stage contrast inside GSE162694. Because control status is nested in
# cohort, `dataset` cannot separate the two. The 1,471 F0-to-F1 genes and the
# GNMT and MAT1A F0-to-F1 effects are therefore partly disease-versus-control.
#
# This is a perspective the existing analysis does not report, not a refutation of
# it. Arm A reproduces the published contrasts unchanged. Arm B removes the source
# controls. Arm C keeps them and adds a control-status term. All three are written
# side by side so a reader can see what the reference is doing.
#
# It also fixes a real defect in the adopted producer: at
# scripts/figures/build_pi_stage_extensions.R:112, `critical <- qt(0.975, df =
# fit$df.total)` returns an UNNAMED vector, and the next line indexes it by gene
# name (`critical[genes]`). Every confidence interval in the adopted stage
# extension output is therefore NA -- 257,070 of them. Point estimates,
# p-values and BH values were unaffected.
#
# Env: MASLD_PROJECT_ROOT, FIGURE_CANDIDATE_ROOT (must not exist), STAGE_RELEASE_ROOT.

suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma)})
options(digits = 17, scipen = 999)
set.seed(42)

fail <- function(...) stop(..., call. = FALSE)
write_tsv_once <- function(x, path) {
  if (file.exists(path)) fail("Refusing overwrite: ", path)
  fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
}

project_root <- normalizePath(Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"), mustWork = TRUE)
candidate_root <- Sys.getenv("FIGURE_CANDIDATE_ROOT", "")
stage_root <- normalizePath(Sys.getenv("STAGE_RELEASE_ROOT", ""), mustWork = TRUE)
if (!nzchar(candidate_root)) fail("FIGURE_CANDIDATE_ROOT is required")
if (file.exists(candidate_root)) fail("FIGURE_CANDIDATE_ROOT must not exist: ", candidate_root)

# same contract as the adopted producer
mm <- fread(file.path(stage_root, "model_input_manifest.tsv"))
if (nrow(mm) != 1L || mm$n_cohorts != 5L || mm$n_samples != 844L || mm$n_genes != 23370L)
  fail("Five-cohort model-input contract drift")
dge_all <- readRDS(normalizePath(mm$dge_path, mustWork = TRUE))
if (!inherits(dge_all, "DGEList") || !identical(dim(dge_all), c(23370L, 844L)))
  fail("Five-cohort DGE identity drift")
ann <- fread(normalizePath(mm$gene_metadata_path, mustWork = TRUE))[, .(
  gene_id_versioned = gene_id, gene_id_base = ensembl_base, gene_name)]
ann <- ann[match(rownames(dge_all), gene_id_versioned)]
if (anyNA(ann$gene_id_versioned)) fail("gene annotation join failed")

elig <- as.data.table(readRDS(normalizePath(mm$meta_path, mustWork = TRUE)))[
  sample_id %in% colnames(dge_all)]
elig <- elig[!is.na(inferred_sex) & !is.na(dataset)]
if (nrow(elig) != 844L) fail("Biological-unit identity drift")

out_root <- file.path(candidate_root, "analysis", "f0_arms")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------- the F0 audit
stage <- elig[dataset %in% c("GSE130970", "GSE135251", "GSE162694") & fibrosis_stage %in% 0:4,
              .(sample_id, dataset, fibrosis_stage, inferred_sex, group_binary, nas_score)]
stage[, fib_group := paste0("F", fibrosis_stage)]
audit <- stage[, .N, by = .(fib_group, dataset, group_binary)][order(fib_group, dataset)]
write_tsv_once(audit, file.path(out_root, "f0_composition_audit.tsv"))
cat("\n=== composition of every fibrosis group by cohort and control status ===\n")
print(dcast(audit, fib_group + dataset ~ group_binary, value.var = "N", fill = 0L))

# ---------------------------------------------------------------- fitting
fit_contrast <- function(tab, ref, cmp, contrast_id, arm, add_group_binary = FALSE) {
  sel <- tab[fib_group %in% c(ref, cmp)]
  sel[, group := factor(fib_group, levels = c(ref, cmp))]
  if (uniqueN(sel$group) != 2L || uniqueN(sel$dataset) < 2L) {
    cat(sprintf("  [skip] %s / %s: %d groups, %d cohorts\n",
                arm, contrast_id, uniqueN(sel$group), uniqueN(sel$dataset))); return(NULL)
  }
  setorder(sel, sample_id)
  dge <- dge_all[, match(sel$sample_id, colnames(dge_all))]
  info <- data.frame(dataset = droplevels(factor(sel$dataset)),
                     inferred_sex = droplevels(factor(sel$inferred_sex)),
                     group = sel$group, row.names = sel$sample_id)
  covs <- "dataset+inferred_sex"
  if (add_group_binary) {
    gb <- droplevels(factor(sel$group_binary))
    if (nlevels(gb) < 2L) {
      cat(sprintf("  [skip] %s / %s: control status is constant\n", arm, contrast_id)); return(NULL)
    }
    info$control_status <- gb
    design <- model.matrix(~ dataset + inferred_sex + control_status + group, info)
    covs <- "dataset+inferred_sex+control_status"
  } else {
    design <- model.matrix(~ dataset + inferred_sex + group, info)
  }
  coef_name <- paste0("group", cmp)
  if (!coef_name %in% colnames(design) || qr(design)$rank != ncol(design)) {
    cat(sprintf("  [skip] %s / %s: rank-deficient design (control status is nested in cohort)\n",
                arm, contrast_id)); return(NULL)
  }
  v <- voomWithQualityWeights(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  ci <- match(coef_name, colnames(design))
  tt <- topTable(fit, coef = ci, number = Inf, sort.by = "none")
  se <- fit$stdev.unscaled[, ci] * sqrt(fit$s2.post)
  # FIX: qt() returns an unnamed vector; the adopted producer indexes it by gene
  # name and gets NA for every interval. Name it so the intervals are real.
  critical <- qt(0.975, df = fit$df.total)
  names(critical) <- rownames(fit)
  g <- rownames(dge)
  res <- cbind(ann, data.table(
    logFC = tt[g, "logFC"], SE = se[g],
    CI_low  = tt[g, "logFC"] - critical[g] * se[g],
    CI_high = tt[g, "logFC"] + critical[g] * se[g],
    t = tt[g, "t"], P.Value = tt[g, "P.Value"],
    FDR = p.adjust(tt[g, "P.Value"], method = "BH"),
    AveExpr = tt[g, "AveExpr"],
    arm = arm, contrast = contrast_id, reference = ref, comparison = cmp,
    n_reference = sum(sel$group == ref), n_comparison = sum(sel$group == cmp),
    n_ref_control = sum(sel$group == ref & sel$group_binary == "Control"),
    n_ref_disease = sum(sel$group == ref & sel$group_binary == "Disease"),
    n_cohorts = uniqueN(sel$dataset), covariates = covs,
    bh_family_size = nrow(dge)))
  if (nrow(res) != nrow(dge_all) ||
      max(abs(res$FDR - p.adjust(res$P.Value, "BH"))) > 1e-12) fail("BH family check failed")
  if (anyNA(res$CI_low) || anyNA(res$CI_high)) fail("CI repair failed: intervals are still NA")
  res
}

adjacent <- list(c("F0","F1"), c("F1","F2"), c("F2","F3"), c("F3","F4"))
vs0 <- paste0("F", 1:4)

arms <- list(
  A_all_F0            = list(tab = copy(stage), gb = FALSE),
  B_disease_only_F0   = list(tab = stage[!(fib_group == "F0" & group_binary == "Control")], gb = FALSE),
  C_all_F0_plus_term  = list(tab = copy(stage), gb = TRUE))

all_res <- list()
for (arm in names(arms)) {
  cat(sprintf("\n=== arm %s ===\n", arm))
  a <- arms[[arm]]
  cat(sprintf("  F0 n = %d (%d control, %d disease)\n",
      sum(a$tab$fib_group == "F0"),
      sum(a$tab$fib_group == "F0" & a$tab$group_binary == "Control"),
      sum(a$tab$fib_group == "F0" & a$tab$group_binary == "Disease")))
  for (p in adjacent) {
    r <- fit_contrast(a$tab, p[1], p[2], paste0(p[1], "_to_", p[2]), arm, a$gb)
    if (!is.null(r)) all_res[[length(all_res) + 1L]] <- r
  }
  for (cmp in vs0) {
    r <- fit_contrast(a$tab, "F0", cmp, paste0(cmp, "_vs_F0"), arm, a$gb)
    if (!is.null(r)) all_res[[length(all_res) + 1L]] <- r
  }
}
res <- rbindlist(all_res)
fwrite(res, file.path(out_root, "f0_arm_results.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")

# ---------------------------------------------------------------- summaries
summ <- res[, .(n_reference = n_reference[1], n_comparison = n_comparison[1],
                n_ref_control = n_ref_control[1], n_ref_disease = n_ref_disease[1],
                n_cohorts = n_cohorts[1], covariates = covariates[1],
                n_BH_sig = sum(FDR < 0.05),
                n_up = sum(FDR < 0.05 & logFC > 0),
                n_down = sum(FDR < 0.05 & logFC < 0),
                n_CI_NA = sum(is.na(CI_low))),
            by = .(arm, contrast)][order(contrast, arm)]
write_tsv_once(summ, file.path(out_root, "f0_arm_summary.tsv"))
cat("\n=== DEG counts by arm ===\n"); print(summ)

heroes <- c("GNMT", "MAT1A", "CYP2C19")
hero <- res[gene_name %in% heroes,
            .(gene_name, arm, contrast, logFC = round(logFC, 4),
              CI_low = round(CI_low, 4), CI_high = round(CI_high, 4),
              FDR = signif(FDR, 3))][order(gene_name, contrast, arm)]
write_tsv_once(hero, file.path(out_root, "f0_arm_hero_genes.tsv"))
cat("\n=== worked-example genes, F0->F1 ===\n")
print(hero[contrast == "F0_to_F1"])

# effect agreement between arms, all genes
w <- dcast(res[contrast == "F0_to_F1"], gene_id_versioned ~ arm, value.var = "logFC")
cc <- data.table(contrast = "F0_to_F1",
                 rho_A_vs_B = if ("B_disease_only_F0" %in% names(w))
                   cor(w$A_all_F0, w$B_disease_only_F0, method = "spearman", use = "complete.obs") else NA_real_,
                 rho_A_vs_C = if ("C_all_F0_plus_term" %in% names(w))
                   cor(w$A_all_F0, w$C_all_F0_plus_term, method = "spearman", use = "complete.obs") else NA_real_)
write_tsv_once(cc, file.path(out_root, "f0_arm_effect_agreement.tsv"))
cat("\n=== all-gene effect agreement, F0->F1 ===\n"); print(cc)

writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))
cat("\nF0_ARMS_COMPLETE: ", out_root, "\n", sep = "")
