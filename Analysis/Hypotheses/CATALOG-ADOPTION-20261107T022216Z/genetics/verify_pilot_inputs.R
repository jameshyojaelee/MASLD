#!/usr/bin/env Rscript
# Lightweight inspection of saved fit inputs; no fitting/eigendecomposition.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
out <- normalizePath(args[[1]], mustWork = TRUE)
script_arg <- commandArgs(FALSE)[startsWith(commandArgs(FALSE), "--file=")]
here <- dirname(normalizePath(sub("--file=", "", script_arg)))
targets <- read.delim(file.path(here, "pilot_targets.tsv"), stringsAsFactors = FALSE)
rows <- list()
for (i in seq_len(nrow(targets))) {
  target <- targets[i, ]
  prefix <- paste(target$gwas_name, target$ensembl, sep = "__")
  object <- readRDS(file.path(out, "exports", paste0(prefix, "__inputs_and_fits.rds")))
  variants <- object$aligned_summary
  ids <- paste0(variants$chromosome, ":", variants$position)
  unique_ids <- !anyDuplicated(ids)
  stopifnot(identical(object$snp_ids, ids),
    identical(rownames(object$R_regularized), ids), identical(colnames(object$R_regularized), ids),
    identical(colnames(object$gwas_fit$lbf_variable), ids),
    isTRUE(all.equal(unname(object$z), variants$beta / variants$standard_error, tolerance = 1e-12)),
    isTRUE(object$gwas_fit$converged), isTRUE(object$eqtl_fit$converged),
    isTRUE(all.equal(object$R_regularized, t(object$R_regularized), tolerance = 1e-12)),
    all(is.finite(object$R_regularized)), all(is.finite(object$z)),
    object$L == 10L, object$ld_diagonal == 1e-3, object$residual_variance == 1,
    object$p1 == 1e-4, object$p2 == 1e-4, object$p12 == 5e-6, object$seed == 42,
    object$build == "GRCh37/hg19")
  # The export precedes coloc and must preserve the actual source eQTL object.
  stopifnot(identical(object$eqtl_fit, readRDS(object$eqtl_fit_path)))
  # SNP includes LD reference allele order, not just chromosome/position.
  allele_ids <- paste(variants$chromosome, variants$position,
                      variants$allele1, variants$allele2, sep = ":")
  stopifnot(identical(allele_ids, variants$SNP))
  credible <- read.delim(file.path(out, "exports", paste0(prefix, "__gwas_credible_sets.tsv")),
                        stringsAsFactors = FALSE)
  stopifnot(all(credible$variant %in% ids), all(credible$alpha >= 0), all(credible$alpha <= 1))
  rows[[i]] <- data.frame(gene = target$gene, study = target$gwas_name,
    chromosome = target$chr, n_variants = length(ids), n_eff = object$n_eff,
    gwas_credible_sets = length(object$gwas_fit$sets$cs),
    eqtl_credible_sets = length(object$eqtl_fit$sets$cs),
    lambda_s_locus = object$lambda_s_locus,
    variant_order_alleles_z_and_LD_dimensions_match = TRUE,
    eqtl_source_object_identical = TRUE, convergence_both = TRUE,
    method_parameters_match = TRUE,
    source_matched_LD_available = FALSE,
    variant_identity_unique = unique_ids,
    terminal_state = if (unique_ids) "original_reference_inputs_reproduced" else
      "unqualified_duplicate_variant_identity",
    interpretation = "original-reference reconstruction only; LD adequacy unresolved")
  if (!unique_ids) {
    collision <- variants[duplicated(ids) | duplicated(ids, fromLast = TRUE), ]
    write.table(collision, file.path(out, "confirmed_variant_identity_collisions.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    collision_rows <- which(duplicated(ids) | duplicated(ids, fromLast = TRUE))
    write.table(object$R_regularized[collision_rows, collision_rows],
      file.path(out, "confirmed_variant_identity_collision_LD.tsv"), sep = "\t", quote = FALSE)
  }
  rm(object); gc(FALSE)
}
result <- do.call(rbind, rows)
destination <- file.path(out, "independent_input_checks.tsv")
stopifnot(!file.exists(destination))
write.table(result, destination, sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "independent_sessionInfo.txt"))
print(result)
if (!all(result$variant_identity_unique)) quit(status = 2L)
