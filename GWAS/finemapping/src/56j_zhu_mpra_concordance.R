#!/usr/bin/env Rscript
# 56j_zhu_mpra_concordance.R  --  external experimental corroboration of the GWAS-ATAC arm
#
# Cross-validate our computational regulatory-variant predictions (motif disruption,
# ABC/cCRE/caQTL variant->gene) against the INDEPENDENT experimental MPRA + sc-CRISPRi
# of Zhu / Hu et al. 2026 (Nat Genet; DOI 10.1038/s41588-026-02617-8; the repo's
# seqfunc "Hu 2025" ingestion). Analog of 56i_currin_caqtl_concordance.R, one modality
# over (experimental reporter/CRISPRi vs bulk caQTL). APPLY-ONLY: writes only to
# results/gwas_atac/, never the multi-evidence atlas or the convergence score.
#
# ---------------------------------------------------------------------------------
# IMPORTANT (method-modality / orientation):
#   A SIGNED MPRA-direction concordance (sign(mpra_log2fc) vs our sign(alleleDiff)) is
#   NOT computed. The on-disk Zhu DAV tables (parsed by 92_hu2025_mpra_benchmark.R from
#   Supp Tables S2-S5) carry only the tested oligo INTERVAL + rsid + mpra_log2fc -- NO
#   effect allele. mpra_log2fc's sign is therefore relative to Zhu's (unstated) allele,
#   not our ALT, so a signed comparison is not well-defined and would return a
#   meaningless ~50% agreement (cf. the below-chance-caQTL trap in 56i's F117 fix).
#   The oligo is centered on the variant, so the POSITION mapping is precise; only the
#   sign is unoriented. We therefore report two orientation-free tests:
#     (A) UNSIGNED MPRA-activity enrichment: are our motif-disrupted regulatory variants
#         enriched for experimentally-significant MPRA differential activity (DAV) vs the
#         MPRA-tested finemapping-substrate background? (Fisher)
#     (B) CRISPRi target-gene SET-MEMBERSHIP: for Zhu's 20 sc-CRISPRi loci, does our
#         variant->gene assignment fall within their candidate target-gene set?
#   A signed arm can be added later IF Zhu's per-variant effect-allele table is fetched.
# ---------------------------------------------------------------------------------
#
# Inputs (read-only):
#   GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv   (SNP_id hg19, alleleDiff, max_pip)
#   GWAS/finemapping/results/gwas_atac/allele_concordance.csv        (SNP_id, linked_gene_resolved, coloc_best_pp4)
#   GWAS/finemapping/results/gwas_atac/gene_level_gwas_atac.csv      (assigned_gene)
#   GWAS/finemapping/results/gwas_atac/abc_variant_to_gene.csv       (variant_id, abc_target_gene)
#   GWAS/finemapping/results/seqfunc/mpra_benchmark/mpra_substrate_truth.tsv  (variant_id_hg19, dav, mpra_log2fc, ...)
#   GWAS/finemapping/results/seqfunc/mpra_benchmark/cropseq_candidate_truth.tsv (rsid, candidate_target_genes)
#
# Outputs (write-only):
#   GWAS/finemapping/results/gwas_atac/zhu_mpra_overlap.csv           (per-variant: our motif-disrupt x MPRA-DAV)
#   GWAS/finemapping/results/gwas_atac/zhu_mpra_overlap_summary.csv   (Fisher enrichment, overall + per context)
#   GWAS/finemapping/results/gwas_atac/zhu_crispri_v2g_agreement.csv  (per-locus set-membership of our V2G)
# Env: rnaseq

suppressPackageStartupMessages(library(data.table))

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT <- file.path(BASE, "GWAS/finemapping/results/gwas_atac")
SF  <- file.path(BASE, "GWAS/finemapping/results/seqfunc/mpra_benchmark")
num <- function(v) suppressWarnings(as.numeric(v))

md_f   <- file.path(OUT, "motif_disruption_scores.csv")
ac_f   <- file.path(OUT, "allele_concordance.csv")
gl_f   <- file.path(OUT, "gene_level_gwas_atac.csv")
abc_f  <- file.path(OUT, "abc_variant_to_gene.csv")
sub_f  <- file.path(SF, "mpra_substrate_truth.tsv")
crop_f <- file.path(SF, "cropseq_candidate_truth.tsv")
stopifnot(file.exists(md_f), file.exists(ac_f), file.exists(sub_f), file.exists(crop_f))

cat("============================================================\n")
cat("56j_zhu_mpra_concordance.R -- GWAS-ATAC x Zhu/Hu 2026 MPRA/CRISPRi\n")
cat("============================================================\n\n")

md  <- fread(md_f);  ac <- fread(ac_f);  sub <- fread(sub_f)
gl  <- if (file.exists(gl_f))  fread(gl_f)  else data.table(assigned_gene = character())
abc <- if (file.exists(abc_f)) fread(abc_f) else data.table(abc_target_gene = character())

# ─────────────────────────────────────────────────────────────────────────────
# (A) UNSIGNED MPRA-activity enrichment
# ─────────────────────────────────────────────────────────────────────────────
# Collapse the MPRA substrate to one row per tested variant (DAV in ANY context).
sub[, dav := num(dav)]
sub_v <- sub[, .(
  is_dav_any     = as.integer(any(dav == 1, na.rm = TRUE)),
  best_abs_log2fc= suppressWarnings(max(abs(num(mpra_log2fc)), na.rm = TRUE)),
  best_coloc_pp4 = suppressWarnings(max(num(coloc_pp4_best), na.rm = TRUE)),
  max_pip        = suppressWarnings(max(num(max_pip), na.rm = TRUE)),
  n_contexts_dav = sum(dav == 1, na.rm = TRUE)
), by = .(variant_id_hg19)]
sub_v[!is.finite(best_abs_log2fc), best_abs_log2fc := NA_real_]
sub_v[!is.finite(best_coloc_pp4),  best_coloc_pp4  := NA_real_]

motif_v <- unique(md$SNP_id)                 # our motif-disrupted regulatory variants (hg19 chr:pos:ref:alt)
sub_v[, our_motif_disrupted := variant_id_hg19 %in% motif_v]

n_tested        <- nrow(sub_v)
n_motif         <- sum(sub_v$our_motif_disrupted)
n_dav           <- sum(sub_v$is_dav_any == 1)
n_motif_and_dav <- sum(sub_v$our_motif_disrupted & sub_v$is_dav_any == 1)

cat("(A) MPRA-activity enrichment among MPRA-tested finemapping variants\n")
cat("  MPRA-tested substrate variants:            ", n_tested, "\n", sep = "")
cat("  ... that are our motif-disrupted variants: ", n_motif, "\n", sep = "")
cat("  ... that are MPRA-DAVs (any context):      ", n_dav, "\n", sep = "")
cat("  ... motif-disrupted AND DAV:               ", n_motif_and_dav, "\n", sep = "")

# Fisher 2x2: motif-disrupted (Y/N) x DAV (Y/N) over MPRA-tested variants
tab <- table(motif = sub_v$our_motif_disrupted, dav = sub_v$is_dav_any == 1)
ft  <- tryCatch(fisher.test(tab), error = function(e) NULL)
dav_rate_motif    <- if (n_motif > 0) n_motif_and_dav / n_motif else NA_real_
dav_rate_nonmotif <- {
  nm <- n_tested - n_motif; if (nm > 0) (n_dav - n_motif_and_dav) / nm else NA_real_
}
cat(sprintf("  DAV rate | motif-disrupted:   %.1f%% (%d/%d)\n",
            100*dav_rate_motif, n_motif_and_dav, n_motif))
cat(sprintf("  DAV rate | not motif-disrupt: %.1f%% (%d/%d)\n",
            100*dav_rate_nonmotif, n_dav - n_motif_and_dav, n_tested - n_motif))
if (!is.null(ft)) cat(sprintf("  Fisher OR = %.2f  [%.2f, %.2f]  p = %.3g\n",
            ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))

fwrite(sub_v[order(-is_dav_any, -best_abs_log2fc)], file.path(OUT, "zhu_mpra_overlap.csv"))

# per-context DAV rate among our motif-disrupted variants
ctx_summ <- sub[variant_id_hg19 %in% motif_v, .(
  n_tested_motif = uniqueN(variant_id_hg19),
  n_dav_motif    = uniqueN(variant_id_hg19[dav == 1])
), by = context][order(context)]
ovl_summ <- rbind(
  data.table(context = "ANY", n_tested_motif = n_motif, n_dav_motif = n_motif_and_dav),
  ctx_summ, fill = TRUE)
ovl_summ[, dav_rate_pct := round(100 * n_dav_motif / n_tested_motif, 1)]
ovl_summ[, fisher_or := NA_real_][context == "ANY",
         `:=`(fisher_or = if (!is.null(ft)) round(unname(ft$estimate),2) else NA_real_,
              fisher_p  = if (!is.null(ft)) signif(ft$p.value,3) else NA_real_,
              n_tested_all = n_tested, n_dav_all = n_dav)]
fwrite(ovl_summ, file.path(OUT, "zhu_mpra_overlap_summary.csv"))
cat("  Wrote zhu_mpra_overlap.csv + zhu_mpra_overlap_summary.csv\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# (B) CRISPRi target-gene SET-MEMBERSHIP
# ─────────────────────────────────────────────────────────────────────────────
crop <- fread(crop_f)
# our variant->gene universe (regulatory-nominated genes), symbol-level
our_v2g <- unique(c(
  if ("assigned_gene" %in% names(gl))          gl$assigned_gene            else character(),
  if ("abc_target_gene" %in% names(abc))       abc$abc_target_gene         else character(),
  if ("linked_gene_resolved" %in% names(ac))   ac$linked_gene_resolved     else character()
))
our_v2g <- our_v2g[!is.na(our_v2g) & our_v2g != "" & our_v2g != "-"]
cat("(B) CRISPRi candidate-set membership vs our variant->gene universe (", length(our_v2g), " genes)\n", sep = "")

crop[, cand := candidate_target_genes]
res <- crop[, {
  cand_genes <- trimws(unlist(strsplit(as.character(cand), ";", fixed = TRUE)))
  cand_genes <- cand_genes[cand_genes != ""]
  hit <- intersect(cand_genes, our_v2g)
  .(n_candidate_genes = length(cand_genes),
    candidate_target_genes = paste(cand_genes, collapse = ";"),
    our_v2g_hit = length(hit) > 0,
    our_v2g_match_gene = paste(hit, collapse = ";"))
}, by = .(rsid)]

n_loci    <- nrow(res)
n_agree   <- sum(res$our_v2g_hit)
# crude expected-by-chance: mean per-locus P(>=1 of k candidates in our set) under a
# background of all protein-coding genes (~20000) -- documents the candidate-set-size caveat.
N_PC <- 20000
res[, p_chance := 1 - choose(N_PC - length(our_v2g), n_candidate_genes) /
                      choose(N_PC, n_candidate_genes)]
exp_by_chance <- sum(res$p_chance, na.rm = TRUE)
cat(sprintf("  Loci with our V2G in Zhu's candidate set: %d / %d (%.0f%%)\n",
            n_agree, n_loci, 100*n_agree/n_loci))
cat(sprintf("  Expected by chance (sum per-locus P, |our set|=%d, ~%d PC genes): %.1f loci\n",
            length(our_v2g), N_PC, exp_by_chance))
cat("  CAVEAT: candidate sets are Zhu's CRISPRi candidate gene SETS (not narrowed single\n")
cat("          targets); this is set-membership agreement, not exact-target agreement.\n")
setorder(res, -our_v2g_hit, -n_candidate_genes)
fwrite(res, file.path(OUT, "zhu_crispri_v2g_agreement.csv"))
print(res[our_v2g_hit == TRUE, .(rsid, n_candidate_genes, our_v2g_match_gene)])
cat("  Wrote zhu_crispri_v2g_agreement.csv\n")

cat("\n[", as.character(Sys.time()), "] Done.\n", sep = "")
