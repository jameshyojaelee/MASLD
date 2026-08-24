# Helpers copied into the isolated replay directory by 07_prepare_genetic_replay.py.

export_coloc_susie_result <- function(
    susie_res, merged, gwas_name, gene_symbol, gene_id, chr_num, export_root) {
  if (is.null(susie_res$summary) || is.null(susie_res$results)) {
    stop("coloc.susie result lacks summary or variant posterior tables")
  }
  signal_summary <- data.table::as.data.table(data.table::copy(susie_res$summary))
  posterior_wide <- data.table::as.data.table(data.table::copy(susie_res$results))
  if (!"snp" %in% names(posterior_wide)) {
    stop("coloc.susie variant posterior table lacks snp")
  }
  posterior_cols <- setdiff(names(posterior_wide), "snp")
  if (length(posterior_cols) != nrow(signal_summary)) {
    stop("coloc.susie signal-pair and posterior-column counts differ")
  }
  expected_posterior_cols <- if (nrow(signal_summary) == 1L) {
    "SNP.PP.H4.abf"
  } else {
    paste0("SNP.PP.H4.row", seq_len(nrow(signal_summary)))
  }
  if (!identical(posterior_cols, expected_posterior_cols)) {
    stop("coloc.susie posterior columns do not align to signal-summary rows")
  }
  if (!all(c("hit1", "hit2", "PP.H4.abf") %in% names(signal_summary))) {
    stop("coloc.susie signal summary schema is incomplete")
  }
  if (anyDuplicated(posterior_wide$snp) || anyNA(posterior_wide$snp)) {
    stop("coloc.susie variant posterior identifiers are not unique and complete")
  }

  signal_summary[, signal_pair_index := seq_len(.N)]
  signal_summary[, posterior_column := posterior_cols]
  signal_summary[, `:=`(
    gwas_name = gwas_name,
    gene = gene_symbol,
    ensembl = gene_id,
    chr = chr_num,
    gwas_signal = hit1,
    eqtl_signal = hit2
  )]
  posterior_long <- data.table::melt(
    posterior_wide,
    id.vars = "snp",
    measure.vars = posterior_cols,
    variable.name = "posterior_column",
    value.name = "SNP.PP.H4"
  )
  posterior_long[, signal_pair_index := match(posterior_column, posterior_cols)]
  allele_lookup <- unique(merged[, .(
    snp = merge_key,
    hg19_position = eqtl_pos,
    allele1 = gwas_a1,
    allele2 = gwas_a2
  )])
  allele_counts <- allele_lookup[, .N, by = snp]
  if (any(allele_counts$N != 1L)) {
    stop("variant posterior has ambiguous hg19 allele metadata")
  }
  posterior_long <- merge(
    posterior_long, allele_lookup, by = "snp", all.x = TRUE, sort = FALSE
  )
  if (posterior_long[, anyNA(hg19_position) || anyNA(allele1) || anyNA(allele2)]) {
    stop("variant posterior lacks a complete hg19 allele lookup")
  }
  posterior_long[, `:=`(
    gwas_name = gwas_name,
    gene = gene_symbol,
    ensembl = gene_id,
    chr = chr_num
  )]
  posterior_sums <- posterior_long[, .(posterior_sum = sum(SNP.PP.H4)),
    by = signal_pair_index]
  if (posterior_long[, any(!is.finite(SNP.PP.H4)) || any(SNP.PP.H4 < 0) ||
      any(SNP.PP.H4 > 1)]) {
    stop("coloc.susie SNP.PP.H4 contains an invalid probability")
  }
  if (nrow(posterior_sums) != nrow(signal_summary) ||
      any(!is.finite(posterior_sums$posterior_sum)) ||
      any(abs(posterior_sums$posterior_sum - 1) > 1e-8)) {
    stop("coloc.susie SNP.PP.H4 does not sum to one for every signal pair")
  }

  gene_export <- file.path(
    export_root, gwas_name, paste0("chr", chr_num), gene_id
  )
  if (dir.exists(gene_export) || file.exists(gene_export)) {
    stop("Refusing to overwrite replay export: ", gene_export)
  }
  dir.create(gene_export, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(
    signal_summary, file.path(gene_export, "signal_pairs.tsv"), sep = "\t"
  )
  data.table::fwrite(
    posterior_long,
    file.path(gene_export, "variant_posteriors.tsv.gz"),
    sep = "\t"
  )
  invisible(list(
    n_signal_pairs = nrow(signal_summary),
    n_variant_posterior_rows = nrow(posterior_long)
  ))
}
