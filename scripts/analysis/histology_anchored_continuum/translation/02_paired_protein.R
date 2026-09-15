#!/usr/bin/env Rscript
# Matched liver/plasma association, not tissue origin or treatment response.
source("scripts/analysis/histology_anchored_continuum/translation/lib_translation.R")
suppressPackageStartupMessages(library(limma))
args <- commandArgs(TRUE); stopifnot(length(args) == 1L)
out <- file.path(args[[1]], "paired_protein")
stopifnot(!dir.exists(out)); dir.create(out, recursive = TRUE)
set.seed(20260906)
writeLines(c("unit=unique_participant_initial_sample", "join=exact_source_metadata_filenames",
  "assay=unambiguous_single_gene_DIA_MS_protein_groups",
  "normalization=log2_positive_intensities_then_assay_separate_quantile_normalization",
  "filter=at_least_50pct_observed_among_matched_participants_per_assay",
  "duplicate_gene_groups=median_log2_abundance_per_sample", "imputation=none",
  "primary=plasma_z~liver_z+factor(fibrosis)+age+bmi+sex+liver_batch+plasma_batch",
  "sensitivity=primary_without_fibrosis", "minimum_complete_participants=20",
  "minimum_residual_df=10", "BH_family=all_shared_genes_after_condition_blind_filter",
  "claim=within_participant_cross_compartment_association_not_liver_source_or_causality",
  "relation_to_bulk_RNA=unpaired", "source_owned_public_cohort=true"),
  file.path(out, "analysis_contract.txt"))
m <- fread("data/PXD051911/meta_data.txt")
l <- fread("data/PXD051911/liver_protein_quant.txt")
p <- fread("data/PXD051911/plasma_protein_quant.txt")
m[, `:=`(liver_available = liver_proteomics_filename %in% names(l),
         plasma_available = plasma_proteomics_filename %in% names(p))]
m[, included := liver_available & plasma_available & sample_group == "initial_sample"]
fwrite(m[, .(unique_identifier, patient_name, sample_group, saf_diagnosis,
  liver_available, plasma_available, included, liver_proteomics_filename,
  plasma_proteomics_filename, plasma_batch_effects)],
  file.path(out, "participant_linkage.tsv"), sep = "\t")
matched <- m[included == TRUE]
stopifnot(!anyDuplicated(matched$patient_name), !anyDuplicated(matched$unique_identifier),
  !anyDuplicated(matched$liver_proteomics_filename),
  !anyDuplicated(matched$plasma_proteomics_filename))
native <- function(raw, samples) {
  r <- raw[!is.na(Genes) & nzchar(Genes) & !grepl(";", Genes, fixed = TRUE)]
  x <- as.matrix(r[, ..samples]); storage.mode(x) <- "double"
  x[!is.finite(x) | x <= 0] <- NA_real_
  x <- log2(x); rownames(x) <- r$Genes
  x <- x[rowMeans(is.finite(x)) >= .5, , drop = FALSE]
  x <- normalizeBetweenArrays(x, method = "quantile")
  idx <- split(seq_len(nrow(x)), rownames(x))
  z <- t(vapply(idx, function(i) apply(x[i, , drop = FALSE], 2, median, na.rm = TRUE),
                numeric(ncol(x))))
  colnames(z) <- matched$unique_identifier
  z[!is.finite(z)] <- NA_real_
  z
}
liv <- native(l, matched$liver_proteomics_filename)
pla <- native(p, matched$plasma_proteomics_filename)
genes <- sort(intersect(rownames(liv), rownames(pla)))
stopifnot(length(genes) > 0L)
matched[, `:=`(fibrosis = suppressWarnings(as.numeric(sub("^F", "", kleiner_fibrosis_grade))),
  age = as.numeric(alder), bmi = as.numeric(bmi), sex = factor(gender),
  liver_batch = factor(ifelse(grepl("^2019", liver_proteomics_filename), "2019", "2020")),
  plasma_batch = factor(plasma_batch_effects))]
fit_rows <- list()
for (gene in genes) {
  for (arm in c("stage_adjusted", "without_stage")) {
    d <- copy(matched)
    d[, `:=`(liver_z = as.numeric(scale(liv[gene, ])), plasma_z = as.numeric(scale(pla[gene, ])))]
    cols <- c("liver_z", "age", "bmi", "sex", "liver_batch", "plasma_batch")
    if (arm == "stage_adjusted") cols <- c(cols, "fibrosis")
    d <- droplevels(d[complete.cases(d[, c("plasma_z", cols), with = FALSE])])
    rhs <- cols
    rhs <- rhs[vapply(rhs, function(k) !is.factor(d[[k]]) || nlevels(d[[k]]) > 1L, logical(1))]
    rhs[rhs == "fibrosis"] <- "factor(fibrosis)"
    r <- data.table(gene_name = gene, model = arm, n_participants = nrow(d),
      n_control = sum(d$saf_diagnosis == "No_MASLD"), n_masld = sum(d$saf_diagnosis != "No_MASLD"),
      beta = NA_real_, se = NA_real_, p_value = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      residual_df = NA_integer_, evidence_state = "untestable", reason = "insufficient_or_singular_design")
    if (nrow(d) >= 20) {
      f <- tryCatch(lm(reformulate(rhs, "plasma_z"), data = d), error = function(e) NULL)
      if (!is.null(f) && f$rank == ncol(model.matrix(f)) && df.residual(f) >= 10 &&
          "liver_z" %in% rownames(coef(summary(f)))) {
        b <- coef(summary(f))["liver_z", ]; t <- qt(.975, df.residual(f))
        r[, `:=`(beta = b[[1]], se = b[[2]], p_value = b[[4]],
          ci_low = b[[1]] - t * b[[2]], ci_high = b[[1]] + t * b[[2]],
          residual_df = df.residual(f), evidence_state = "tested", reason = "complete_participant_native_assay_model")]
      }
    }
    fit_rows[[length(fit_rows) + 1L]] <- r
  }
}
result <- rbindlist(fit_rows)
result[, BH_q := tr_bh(p_value, length(genes)), by = model]
result[evidence_state == "tested", evidence_state := fifelse(BH_q < .05,
  "matched_tissue_plasma_association", "tested_without_BH_support")]
result[, `:=`(family_size = length(genes), effect_unit = "SD_plasma_log2_protein_per_SD_liver_log2_protein",
  matched_to_bulk_RNA = FALSE, tissue_origin_established = FALSE)]
fwrite(result, file.path(out, "paired_tissue_plasma_effects.tsv"), sep = "\t")
fwrite(data.table(n_matched_participants = nrow(matched), n_shared_genes = length(genes),
  n_control = sum(matched$saf_diagnosis == "No_MASLD"),
  n_masld = sum(matched$saf_diagnosis != "No_MASLD")),
  file.path(out, "denominators.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
cat("Matched protein participants:", nrow(matched), "; complete shared-gene family:", length(genes), "\n")
print(result[, .N, by = .(model, evidence_state)])
