# coloc_signal_eligibility.R — per-signal-pair shared-posterior check for coloc.susie()
# Sourced by src/06_susie_coloc.R and by the synthetic test in
# Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/tests/.
#
# coloc.bf_bf() (coloc 5.2.3) drops a signal pair when either signal keeps less
# than overlap.min (default 0.5) of its posterior on the SNPs the two fits share.
# It measures that share on the FULL lbf_variable rows, via
# logbf_to_pp(..., last_is_null = TRUE), which treats the last column as the
# null (SuSiE fits have no null column, so the last SNP gets prior 1 - n*p).
# The check therefore only works when coloc.susie() receives the full fits;
# subsetting both fits to the shared SNPs first makes every share exactly 1.
# The functions below replicate the 5.2.3 calculation so that every pair,
# dropped or kept, is recorded with its share.

COLOC_OVERLAP_MIN <- 0.5

# Share of each signal's (row's) posterior on `shared`, coloc 5.2.3 convention,
# plus the top variant of the full posterior (coloc's all-dropped `hit` rule).
signal_shared_prop <- function(bf, shared, prior) {
  if ("null" %in% colnames(bf)) bf <- bf - matrix(bf[, "null"], nrow(bf), ncol(bf))
  pp <- coloc:::logbf_to_pp(bf, prior, last_is_null = TRUE)
  list(prop = rowSums(pp[, shared, drop = FALSE]) /
              rowSums(pp[, setdiff(colnames(pp), "null"), drop = FALSE]),
       hit  = colnames(pp)[apply(pp, 1, which.max)])
}

# Same share under a uniform prior with no null column (the normalised SuSiE
# single-effect posterior). Diagnostic only; eligibility uses the coloc share.
signal_shared_prop_uniform <- function(bf, shared) {
  pp <- exp(bf - apply(bf, 1, max))
  rowSums(pp[, shared, drop = FALSE]) / rowSums(pp)
}

# One row per credible-set pair of two SuSiE fits, joined to the coloc.susie()
# result computed on the SAME full fits. Returns the pair table, the kept coloc
# rows, and whether our eligibility call equals coloc's (kept == eligible).
coloc_susie_pair_table <- function(s1, s2, susie_res) {
  bf_formals <- formals(coloc:::coloc.bf_bf)
  cs1 <- s1$sets$cs_index
  cs2 <- s2$sets$cs_index
  bf1 <- s1$lbf_variable[cs1, , drop = FALSE]
  bf2 <- s2$lbf_variable[cs2, , drop = FALSE]
  shared <- setdiff(intersect(colnames(bf1), colnames(bf2)), "null")
  sp1 <- signal_shared_prop(bf1, shared, eval(bf_formals$p1))
  sp2 <- signal_shared_prop(bf2, shared, eval(bf_formals$p2))
  u1 <- signal_shared_prop_uniform(bf1, shared)
  u2 <- signal_shared_prop_uniform(bf2, shared)

  pt <- data.table::CJ(k1 = seq_along(cs1), k2 = seq_along(cs2))
  pt[, `:=`(idx1 = cs1[k1], idx2 = cs2[k2],
            prop_gwas = sp1$prop[k1], prop_eqtl = sp2$prop[k2],
            prop_gwas_uniform = u1[k1], prop_eqtl_uniform = u2[k2])]
  pt[, eligible := !(prop_gwas < COLOC_OVERLAP_MIN | prop_eqtl < COLOC_OVERLAP_MIN)]

  # coloc returns only the kept pairs, or every pair with NA PP.H4 when all
  # were dropped.
  kept <- data.table::as.data.table(susie_res$summary)[!is.na(PP.H4.abf)]
  hit <- match(paste(pt$idx1, pt$idx2), paste(kept$idx1, kept$idx2))
  pt[, `:=`(hit1  = ifelse(is.na(hit), sp1$hit[k1], kept$hit1[hit]),
            hit2  = ifelse(is.na(hit), sp2$hit[k2], kept$hit2[hit]),
            PP.H4 = kept$PP.H4.abf[hit],
            PP.H3 = kept$PP.H3.abf[hit])]
  matches <- identical(pt$eligible %in% TRUE, !is.na(hit)) &&
    isTRUE(all.equal(eval(bf_formals$overlap.min), COLOC_OVERLAP_MIN))
  list(pairs = pt[, .(hit1, hit2, idx1, idx2, PP.H4, PP.H3, prop_gwas, prop_eqtl,
                      prop_gwas_uniform, prop_eqtl_uniform, eligible)],
       kept = kept,
       matches_coloc = matches)
}

# Gene-level SuSiE status from the kept coloc rows. When every pair is dropped
# the gene is untestable on these SNPs, not PP4-negative.
coloc_susie_method <- function(kept) {
  if (nrow(kept) > 0) "susie" else "susie_untestable_insufficient_shared_posterior"
}
