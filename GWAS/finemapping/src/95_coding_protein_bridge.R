#!/usr/bin/env Rscript
# Apply-only bridge from predicted coding consequence to measured liver protein
# abundance and external liver-pQTL evidence. Absence of protein evidence is
# reported as untested/undetected, never as a benign functional call.
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SF <- file.path(ROOT, "GWAS/finemapping/results/seqfunc")
OUT <- file.path(SF, "coding_protein_bridge")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
coding <- fread(file.path(SF, "coding_hardening_v2.tsv"))
prot <- fread(file.path(ROOT, "Analysis/Proteomics/results/protein_transcript_concordance_v3.csv"))
pqtl <- fread(file.path(ROOT, "RNA-seq/results/multi_evidence/external_pqtl/pqtl_external_annotation.tsv"))
prot[, protein_padj := suppressWarnings(as.numeric(protein_padj))]
prot[, protein_logFC := suppressWarnings(as.numeric(protein_logFC))]
prot_summary <- prot[, .(
  protein_tested = TRUE,
  n_protein_datasets = uniqueN(dataset),
  min_protein_padj = suppressWarnings(min(protein_padj, na.rm = TRUE)),
  strongest_protein_logFC = protein_logFC[which.max(abs(protein_logFC))],
  any_protein_fdr05 = any(protein_padj < 0.05, na.rm = TRUE),
  all_measured_directions_agree = uniqueN(sign(protein_logFC[is.finite(protein_logFC)])) <= 1
), by = gene]
prot_summary[!is.finite(min_protein_padj), min_protein_padj := NA_real_]
pqkeep <- intersect(c("gene", "gobeil_pqtl_coloc_pph4", "layer_relationship",
                      "gobeil_pqtl_hit", "concordance_call"), names(pqtl))
out <- merge(coding, prot_summary, by = "gene", all.x = TRUE)
out <- merge(out, pqtl[, ..pqkeep], by = "gene", all.x = TRUE)
out[, protein_tested := !is.na(protein_tested) & protein_tested]
out[, protein_abundance_context := fifelse(
  protein_tested & any_protein_fdr05, "gene_has_disease_associated_protein_abundance",
  fifelse(protein_tested, "measured_no_FDR05_abundance_change", "protein_not_measured"))]
out[, interpretation := fifelse(
  protein_abundance_context == "protein_not_measured", "unknown_not_benign",
  fifelse(!is.na(gobeil_pqtl_hit) & gobeil_pqtl_hit,
          "coding_prediction_with_gene_level_pQTL_coloc_context_not_variant_linked",
          "coding_prediction_with_gene_level_abundance_context_only"))]
out[, same_variant_protein_mechanism_validated := FALSE]
out[, pqtl_supported := !is.na(gobeil_pqtl_hit) & gobeil_pqtl_hit]
out[, abundance_supported := !is.na(any_protein_fdr05) & any_protein_fdr05]
setorder(out, -pqtl_supported, -abundance_supported, gene)
fwrite(out, file.path(OUT, "coding_variant_protein_evidence.tsv"), sep = "\t", na = "")
summary <- list(
  status = "complete_apply_only_v1", apply_only_firewall = TRUE,
  n_coding_variants = nrow(out), n_genes = uniqueN(out$gene),
  n_protein_measured = out[protein_tested == TRUE, .N],
  n_protein_fdr05 = out[any_protein_fdr05 == TRUE, .N],
  n_external_liver_pqtl = out[gobeil_pqtl_hit == TRUE, .N],
  caveat = paste("Bulk disease proteomics and germline pQTL answer different questions and are not summed.",
    "Both joins are gene-level context; neither demonstrates that the nominated coding allele changes protein abundance.",
    "MR direction fields are deliberately excluded under the project-wide no-MR rule.")
)
write_json(summary, file.path(OUT, "summary.json"), pretty = TRUE, auto_unbox = TRUE)
message("[95] wrote coding-to-protein bridge for ", nrow(out), " variants")
