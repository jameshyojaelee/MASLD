#!/usr/bin/env Rscript
# Apply the 90-gene eQTL parity refit to the CONF1 candidate Catalog and rerun the
# genetics-by-RNA overlap test. Writes a new candidate directory; reads everything else.
suppressPackageStartupMessages(library(data.table))
set.seed(42)
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
run <- file.path(root, "Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z")
src <- file.path(root, "RNA-seq/results/evidence_passports/candidates/gene-catalog-conf1-20261110T023225Z")
out <- Sys.getenv("OUT_DIR"); stopifnot(nzchar(out), !dir.exists(out))
dir.create(out, recursive = TRUE)

t <- fread(file.path(src, "catalog_annotated_gene_table_conf1.tsv.gz"), na.strings = c("", "NA"))
cmp <- fread(file.path(run, "results/eqtl_dedup_gene_comparison.tsv"))
stopifnot(nrow(cmp) == 90, all(cmp$ensembl %in% t$ensembl_id))

# ---- overlap test (adopted definition, reconstructed) ----
overlap <- function(tab, label) {
  u <- tab[gen_susie_state_t12 %in% c("positive", "tested_le_0.5") & r8_rna_testability_state == "testable"]
  u[, gen := gen_susie_state_t12 == "positive"]
  u[, rna := rna_canonical_call == TRUE]
  br <- unique(quantile(u$rna_ave_log2cpm, probs = seq(0, 1, 0.1), na.rm = TRUE))
  u[, bin := cut(rna_ave_log2cpm, br, include.lowest = TRUE)]
  stopifnot(!anyNA(u$bin))
  arr <- table(factor(u$gen, c(TRUE, FALSE)), factor(u$rna, c(TRUE, FALSE)), u$bin)
  mh <- mantelhaen.test(arr, exact = TRUE)
  expct <- u[, .(e = sum(gen) * sum(rna) / .N), by = bin][, sum(e)]
  data.table(set = label, n = nrow(u), genetic = sum(u$gen), rna = sum(u$rna),
             shared = sum(u$gen & u$rna), expected = expct,
             or = unname(mh$estimate), ci_low = mh$conf.int[1], ci_high = mh$conf.int[2],
             p_exact = mh$p.value)
}
base <- overlap(t, "adopted_r2")

# ---- apply refit ----
new <- copy(t)
new[, eqtl_refit_state := NA_character_]
new[cmp, on = .(ensembl_id = ensembl), eqtl_refit_state := i.new_state]
pos <- cmp[new_state == "positive"]
new[pos, on = .(ensembl_id = ensembl), `:=`(gen_coloc_best_susie_pp4 = i.new_best_susie,
                                           gen_coloc_best_susie_study = i.new_driver)]
lost <- cmp[new_state == "lost", ensembl]
unt <- cmp[grepl("^untestable", new_state), ensembl]
stopifnot(length(lost) == 1, length(unt) == 6)
stopifnot(all(new[ensembl_id %in% c(lost, unt), primary_evidence_class == "genetically_anchored" &
                  candidate_class == "genetically_anchored" & coding_fine_mapped == FALSE &
                  rna_canonical_call == FALSE]))
new[ensembl_id %in% lost, `:=`(
  genetics_call_state = "indeterminate", genetics_testability_state = "testable",
  gen_susie_state_t12 = "tested_le_0.5",
  gen_coloc_best_susie_pp4 = cmp[new_state == "lost", new_best_susie],
  gen_coloc_best_susie_study = cmp[new_state == "lost", new_driver],
  primary_evidence_class = "unresolved", candidate_class = "unresolved",
  class_display_r9 = "unresolved")]
new[ensembl_id %in% unt, `:=`(
  genetics_call_state = "untestable", genetics_testability_state = "eqtl_refit_untestable",
  gen_susie_state_t12 = "untestable", gen_coloc_best_susie_pp4 = NA_real_,
  gen_coloc_best_susie_study = NA_character_,
  primary_evidence_class = "untested", candidate_class = "untested",
  class_display_r9 = "untested")]
new[ensembl_id %in% c(lost, unt), `:=`(
  regulatory_colocalization = FALSE, genetic_mechanism_candidate = "none",
  disease_liability_candidate = "not_assessed", liability_sidecar_class = NA_character_,
  class_changed = TRUE)]

refit <- overlap(new, "eqtl_refit")
ov <- rbind(base, refit)
fwrite(ov, file.path(out, "overlap_test.tsv"), sep = "\t")

counts <- merge(t[, .(conf1 = .N), by = .(class = candidate_class)],
                new[, .(refit = .N), by = .(class = candidate_class)], by = "class", all = TRUE)
adopted <- merge(t[, .(adopted = .N), by = .(class = primary_evidence_class)],
                 new[, .(adopted_refit = .N), by = .(class = primary_evidence_class)], by = "class", all = TRUE)
fwrite(merge(adopted, counts, by = "class", all = TRUE), file.path(out, "class_counts.tsv"), sep = "\t")
liab <- merge(t[candidate_class %in% c("genetically_anchored", "both_supported"), .(conf1 = .N), by = disease_liability_candidate],
              new[candidate_class %in% c("genetically_anchored", "both_supported"), .(refit = .N), by = disease_liability_candidate],
              by = "disease_liability_candidate", all = TRUE)
fwrite(liab, file.path(out, "liability_counts.tsv"), sep = "\t")
fwrite(new[ensembl_id %in% c(lost, unt), .(ensembl_id, symbol, primary_evidence_class, candidate_class,
          genetics_call_state, eqtl_refit_state)], file.path(out, "reclassified_genes.tsv"), sep = "\t")
fwrite(new, file.path(out, "catalog_annotated_gene_table_conf1_eqtlfix.tsv.gz"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(ov); print(merge(adopted, counts, by = "class", all = TRUE)); print(liab)
