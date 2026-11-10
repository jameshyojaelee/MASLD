# Unique variant identities are required by position-keyed coloc SuSiE fits.
# Distinct source variants at one coordinate cannot safely share CHR:POS IDs.
# Return all collisions; never select an allele representation by its result.
COLOC_AMBIGUOUS_VARIANT_METHOD <- "susie_untestable_ambiguous_variant_identity"

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
    '          }'
  )
}

coloc_identity_guarded_runner <- function(lines) {
  anchor <- '          snp_ids <- paste0(chr_num, ":", ss_ld$position)'
  hit <- which(lines == anchor)
  stopifnot(length(hit) == 1L)
  append(lines, coloc_variant_identity_guard(), after = hit)
}
