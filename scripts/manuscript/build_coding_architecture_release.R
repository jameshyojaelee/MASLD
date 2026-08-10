#!/usr/bin/env Rscript

# Quantify coding versus noncoding posterior mass within credible sets linked to
# primary Tier-1/2 SuSiE-colocalized genes. This bounds the result to analyzed
# colocalized loci and replaces the one-top-variant architecture headline.

suppressPackageStartupMessages(library(data.table))

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
RELEASE_ID <- Sys.getenv("MANUSCRIPT_RELEASE_ID", "2026-07-15-r2")
OUT <- file.path(BASE, "RNA-seq/results/manuscript_release", RELEASE_ID)
DOC <- file.path(BASE, "docs/manuscript/release")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

COLOC_FILE <- file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv")
TIER_FILE <- file.path(BASE, "GWAS/finemapping/config/gwas_trait_tier.tsv")
CS_FILE <- file.path(BASE, "GWAS/finemapping/results/credible_sets.csv")
ANN_FILE <- file.path(BASE,
  "GWAS/finemapping/results/credible_set_variant_consequences.csv")

coloc <- fread(COLOC_FILE)
tiers <- fread(TIER_FILE)
cs <- fread(CS_FILE)
ann <- fread(ANN_FILE)

coloc <- merge(coloc, tiers, by.x = "gwas_name", by.y = "study_name", all = FALSE)
coloc <- coloc[
  placement == "main" & tier %in% c(1L, 2L) &
  !is.na(gene) & gene != "" &
  !is.na(`PP.H4.susie`) & `PP.H4.susie` > 0.5 &
  !is.na(top_snp) & top_snp != ""
]
setorder(coloc, gwas_name, gene, -`PP.H4.susie`)
coloc <- coloc[!duplicated(coloc, by = c("gwas_name", "gene"))]
coloc[, c("top_chr", "top_pos") := tstrsplit(top_snp, ":", fixed = TRUE)]
coloc[, `:=`(
  top_chr = as.integer(sub("^chr", "", top_chr)),
  top_pos = as.integer(top_pos)
)]

cs <- cs[
  susie_converged %in% TRUE & susie_reliable %in% TRUE &
  !is.na(susie_cs) & susie_cs > 0 & either_in_cs %in% TRUE
]
cs[, variant_chrpos := paste0(chromosome, ":", position)]

# The COLOC top SNP identifies the credible set supporting each gene-study pair.
hits <- merge(
  coloc[, .(
    gwas_name, gene, ensembl, tier, tier_label, ancestry,
    pp4_susie = `PP.H4.susie`, top_chr, top_pos
  )],
  cs[, .(study, chromosome, position, locus, susie_cs)],
  by.x = c("gwas_name", "top_chr", "top_pos"),
  by.y = c("study", "chromosome", "position"),
  all = FALSE,
  allow.cartesian = TRUE
)
keys <- unique(hits[, .(
  study = gwas_name, locus, susie_cs, tier, tier_label, ancestry
)])
members <- merge(
  cs,
  keys,
  by = c("study", "locus", "susie_cs", "ancestry"),
  all = FALSE,
  allow.cartesian = TRUE
)

gene_map <- hits[, .(
  coloc_genes = paste(sort(unique(gene)), collapse = ";"),
  max_coloc_pp4 = max(pp4_susie, na.rm = TRUE)
), by = .(study = gwas_name, locus, susie_cs)]
members <- merge(members, gene_map, by = c("study", "locus", "susie_cs"), all.x = TRUE)

ann[, chrpos := sub("^([^:]+:[^:]+):.*$", "\\1", variant_id)]
class_rank <- c(
  coding_protein_altering = 1L,
  coding_synonymous = 2L,
  splice_region = 3L,
  noncoding = 4L
)
ann[, class_order := class_rank[class]]
setorder(ann, chrpos, class_order)
ann_chrpos <- ann[!duplicated(chrpos), .(chrpos, consequence_class = class, Consequence)]
members <- merge(
  members,
  ann_chrpos,
  by.x = "variant_chrpos",
  by.y = "chrpos",
  all.x = TRUE
)
members[, is_coding := consequence_class %in% c(
  "coding_protein_altering", "coding_synonymous"
)]
members[, is_protein_altering := consequence_class == "coding_protein_altering"]
members[, pip := as.numeric(recommended_pip)]
members[!is.finite(pip) | pip < 0, pip := 0]
members[, analysis_release_id := RELEASE_ID]

locus <- members[, {
  total_pip <- sum(pip, na.rm = TRUE)
  coding_pip <- sum(pip[is_coding %in% TRUE], na.rm = TRUE)
  altering_pip <- sum(pip[is_protein_altering %in% TRUE], na.rm = TRUE)
  .(
    n_members = .N,
    n_annotated = sum(!is.na(consequence_class)),
    total_recommended_pip = total_pip,
    coding_posterior_mass = if (total_pip > 0) coding_pip / total_pip else NA_real_,
    protein_altering_posterior_mass = if (total_pip > 0) altering_pip / total_pip else NA_real_,
    any_coding_variant = any(is_coding %in% TRUE),
    any_coding_pip_gt_0_1 = any(is_coding %in% TRUE & pip > 0.1),
    any_coding_pip_gt_0_5 = any(is_coding %in% TRUE & pip > 0.5),
    coloc_genes = first(coloc_genes),
    max_coloc_pp4 = first(max_coloc_pp4),
    tier = first(tier),
    tier_label = first(tier_label),
    ancestry = first(ancestry)
  )
}, by = .(analysis_release_id, study, locus, susie_cs)]
fwrite(locus, file.path(OUT, "credible_set_coding_architecture.tsv"), sep = "\t", na = "")

summary <- locus[, .(
  n_credible_sets = .N,
  n_coloc_genes = uniqueN(unlist(strsplit(coloc_genes, ";", fixed = TRUE))),
  median_coding_posterior_mass = median(coding_posterior_mass, na.rm = TRUE),
  mean_coding_posterior_mass = mean(coding_posterior_mass, na.rm = TRUE),
  n_any_coding_variant = sum(any_coding_variant),
  n_any_coding_pip_gt_0_1 = sum(any_coding_pip_gt_0_1),
  n_any_coding_pip_gt_0_5 = sum(any_coding_pip_gt_0_5),
  n_majority_coding_mass = sum(coding_posterior_mass > 0.5, na.rm = TRUE)
), by = .(analysis_release_id, trait_scope = tier_label)]
summary <- rbind(
  locus[, .(
    n_credible_sets = .N,
    n_coloc_genes = uniqueN(unlist(strsplit(coloc_genes, ";", fixed = TRUE))),
    median_coding_posterior_mass = median(coding_posterior_mass, na.rm = TRUE),
    mean_coding_posterior_mass = mean(coding_posterior_mass, na.rm = TRUE),
    n_any_coding_variant = sum(any_coding_variant),
    n_any_coding_pip_gt_0_1 = sum(any_coding_pip_gt_0_1),
    n_any_coding_pip_gt_0_5 = sum(any_coding_pip_gt_0_5),
    n_majority_coding_mass = sum(coding_posterior_mass > 0.5, na.rm = TRUE)
  ), by = .(analysis_release_id)][, trait_scope := "all"],
  summary,
  fill = TRUE
)
summary[, `:=`(
  n_primary_susie_genes = uniqueN(coloc$gene),
  pct_primary_genes_linked = 100 * n_coloc_genes / uniqueN(coloc$gene)
)]
setcolorder(summary, c("analysis_release_id", "trait_scope"))
fwrite(summary, file.path(OUT, "credible_set_coding_summary.tsv"), sep = "\t", na = "")

# Sensitivity: classify the single top COLOC variant per gene with coding
# precedence. This is retained only to explain the historical "coding-led" count.
top_gene <- coloc[order(-`PP.H4.susie`)][!duplicated(gene)]
top_gene[, chrpos := paste0(top_chr, ":", top_pos)]
top_gene <- merge(top_gene, ann_chrpos, by = "chrpos", all.x = TRUE)
top_gene[, coding_led := consequence_class %in% c(
  "coding_protein_altering", "coding_synonymous"
)]
fwrite(top_gene[, .(
  analysis_release_id = RELEASE_ID,
  gene, ensembl, gwas_name, tier_label, ancestry,
  pp4_susie = `PP.H4.susie`, top_snp, consequence_class, Consequence, coding_led
)], file.path(OUT, "top_coloc_variant_sensitivity.tsv"), sep = "\t", na = "")

top_summary <- top_gene[, .(
  n_genes = .N,
  n_annotated = sum(!is.na(consequence_class)),
  n_coding_led = sum(coding_led, na.rm = TRUE),
  pct_coding_led = 100 * mean(coding_led, na.rm = TRUE)
)]
fwrite(top_summary, file.path(OUT, "top_coloc_variant_summary.tsv"), sep = "\t", na = "")

ledger_path <- file.path(DOC, "claim_ledger.tsv")
if (file.exists(ledger_path)) {
  ledger <- fread(ledger_path)
  all_summary <- summary[trait_scope == "all"][1]
  coverage_ok <- all_summary$pct_primary_genes_linked >= 80
  wording <- sprintf(
    "Credible-set coding architecture was estimable for %d of %d primary SuSiE-colocalized genes (%.1f%% coverage), spanning %d reliable credible sets. This coverage is insufficient for a genome-wide architecture claim.",
    all_summary$n_coloc_genes,
    all_summary$n_primary_susie_genes,
    all_summary$pct_primary_genes_linked,
    all_summary$n_credible_sets
  )
  ledger[claim_id == "C_CODING_SCOPE", `:=`(
    status = if (coverage_ok) "supported" else "insufficient_credible_set_mapping",
    allowed_wording = wording,
    source_output = "credible_set_coding_architecture.tsv; credible_set_coding_summary.tsv"
  )]
  fwrite(ledger, ledger_path, sep = "\t", na = "")

  gates_path <- file.path(OUT, "acceptance_gates.tsv")
  if (file.exists(gates_path)) {
    gates <- fread(gates_path)
    gates <- gates[gate != "coding_architecture_coverage"]
    gates <- rbind(gates, data.table(
      analysis_release_id = RELEASE_ID,
      gate = "coding_architecture_coverage",
      passed = coverage_ok,
      criterion = ">=80% of primary SuSiE-colocalized genes link unambiguously to reliable credible sets",
      detail = sprintf("linked=%d/%d (%.1f%%)", all_summary$n_coloc_genes,
                       all_summary$n_primary_susie_genes,
                       all_summary$pct_primary_genes_linked)
    ), fill = TRUE)
    fwrite(gates, gates_path, sep = "\t", na = "")
  }
}

print(summary)
print(top_summary)
