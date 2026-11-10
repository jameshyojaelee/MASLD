# Unique variant identities are required by position-keyed coloc SuSiE fits.
# Distinct source variants at one coordinate cannot safely share CHR:POS IDs.
# Return all collisions; never select an allele representation by its result.
COLOC_AMBIGUOUS_VARIANT_METHOD <- "susie_untestable_ambiguous_variant_identity"
COLOC_INVALID_REFERENCE_LD_METHOD <- "susie_untestable_invalid_reference_LD"

# Producers add a positive diagonal before this check. Such a regularized
# correlation matrix must admit Cholesky factorization; a failed check must not
# be silently repaired inside an LD diagnostic or accepted with check_R=FALSE.
coloc_reference_ld_valid <- function(R) {
  if (!is.matrix(R) || nrow(R)!=ncol(R) || !all(is.finite(R)) ||
      !isTRUE(all.equal(unname(R),unname(t(R)),tolerance=1e-10))) return(FALSE)
  tryCatch({chol(R);TRUE},error=function(e) FALSE)
}

coloc_variant_identity_collisions <- function(ids) {
  which(duplicated(ids) | duplicated(ids, fromLast = TRUE))
}

coloc_variant_identity_guard <- function() {
  c(
    '          identity_collisions <- coloc_variant_identity_collisions(snp_ids)',
    '          eqtl_identity_collisions <- coloc_variant_identity_collisions(eqtl_snp_ids)',
    '          if (length(identity_collisions) || length(eqtl_identity_collisions)) {',
    '            susie_method <- COLOC_AMBIGUOUS_VARIANT_METHOD',
    '            if (length(identity_collisions)) {',
    '              collision_rows <- as.data.table(ss_ld)[identity_collisions]',
    '              collision_rows[, coloc_variant_id := snp_ids[identity_collisions]]',
    '              fwrite(collision_rows, file.path(OUT_DIR, paste0(',
    '                "variant_identity_collisions_chr", chr_num, "_", gene_id, "_gwas.tsv")), sep = "\t")',
    '              fwrite(merged[eqtl_pos %in% ss_ld$position[identity_collisions]],',
    '                file.path(OUT_DIR, paste0("variant_identity_collisions_chr", chr_num,',
    '                  "_", gene_id, "_harmonized.tsv")), sep = "\t")',
    '            }',
    '            if (length(eqtl_identity_collisions)) {',
    '              fwrite(data.table(coloc_variant_id = eqtl_snp_ids[eqtl_identity_collisions]),',
    '                file.path(OUT_DIR, paste0("variant_identity_collisions_chr", chr_num,',
    '                  "_", gene_id, "_eqtl.tsv")), sep = "\t")',
    '            }',
    '            stop("Ambiguous variant identity: duplicate CHR:POS SuSiE columns; ABF retained")',
    '          }',
    '          if (!coloc_reference_ld_valid(R_reg)) {',
    '            susie_method <- COLOC_INVALID_REFERENCE_LD_METHOD',
    '            stop("Invalid regularized reference LD; SuSiE not fitted; ABF retained separately")',
    '          }'
  )
}

coloc_identity_guarded_runner <- function(lines) {
  anchor <- '          snp_ids <- paste0(chr_num, ":", ss_ld$position)'
  hit <- which(lines == anchor)
  stopifnot(length(hit) == 1L)
  append(lines, coloc_variant_identity_guard(), after = hit)
}

# The eQTL producer uses the same LD matching code. Refuse its collisions before
# fitting, rather than saving a posterior whose columns cannot identify variants.
coloc_identity_guarded_eqtl_runner <- function(lines) {
  anchor <- '  snp_ids <- paste0(CHR, ":", ss_filtered$position)'
  hit <- which(lines == anchor)
  stopifnot(length(hit) == 1L)
  guard <- c(
    '  identity_collisions <- coloc_variant_identity_collisions(snp_ids)',
    '  if (length(identity_collisions)) {',
    '    collision_rows <- as.data.table(ss_filtered)[identity_collisions]',
    '    collision_rows[, coloc_variant_id := snp_ids[identity_collisions]]',
    '    fwrite(collision_rows, file.path(chr_out_dir, paste0(ensg_id, "_identity_collisions.tsv")), sep="\t")',
    '    fwrite(data.table(gene=gene_name, ensg=ensg_id, chr=CHR, n_snps=n_snps,',
    '      converged=NA, n_cs=0L, max_pip=NA_real_, status="ambiguous_variant_identity"),',
    '      summary_file, sep="\t", append=TRUE, col.names=FALSE)',
    '    n_failed <- n_failed + 1L',
    '    next',
    '  }',
    '  if (!coloc_reference_ld_valid(R_mat + diag(LD_REGULARIZE,nrow(R_mat)))) {',
    '    fwrite(data.table(gene=gene_name, ensg=ensg_id, chr=CHR, n_snps=n_snps,',
    '      converged=NA, n_cs=0L, max_pip=NA_real_, status="invalid_reference_LD"),',
    '      summary_file, sep="\t", append=TRUE, col.names=FALSE)',
    '    n_failed <- n_failed + 1L',
    '    next',
    '  }',
    '  set.seed(42)'
  )
  append(lines, guard, after=hit)
}
