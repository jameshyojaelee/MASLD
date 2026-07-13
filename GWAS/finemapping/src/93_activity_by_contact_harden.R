#!/usr/bin/env Rscript

# Predicted-link variant-to-gene hypotheses for the seqfunc layer.
#
# Authority order is explicit:
#   1. Nasser et al. liver/HepG2 ABC: published ABC predictions using activity
#      and AvgHiC/contact estimates; these are not direct measurements of each
#      variant-gene relationship.
#   2. MASLD multiome SCENIC+ peak-gene correlation: disease-context support,
#      but not a direct contact measurement and therefore never called "ABC".
# AlphaGenome O/E contact-map argmax is retained only as a diagnostic in the
# historical output and is not used here. Outputs are apply-only.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SF <- file.path(ROOT, "GWAS/finemapping/results/seqfunc")
ABC <- file.path(ROOT, "GWAS/finemapping/results/gwas_atac/abc_variant_to_gene.csv")
SUB <- file.path(SF, "variant_substrate_hg38.tsv")
SCENIC <- file.path(ROOT, "Analysis/ATAC/Human_Multiome/scenic_plus")
OUT <- file.path(SF, "activity_by_contact")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(ABC), file.exists(SUB))
sub <- fread(SUB)
sub[, variant_key := sub("^chr", "", variant_id_hg19)]
sub[, max_pip := suppressWarnings(as.numeric(max_pip))]
sub[, pos_hg38 := as.integer(pos_hg38)]
sub1 <- sub[order(-max_pip), .SD[1], by = variant_key]

# Source 1: published ABC links in the source genome build (hg19).
abc <- fread(ABC)
abc[, variant_key := sub("^chr", "", variant_id)]
abc[, `:=`(abc_score = as.numeric(abc_score), max_rec_pip = as.numeric(max_rec_pip))]
abc <- merge(abc, sub1[, .(variant_key, pos_hg38, var_class, consequence,
                            substrate_gene = gene, substrate_pip = max_pip)],
             by = "variant_key", all = FALSE)
abc[, `:=`(
  source = "Nasser2021_ABC_AvgHiC",
  target_gene = abc_target_gene,
  cell_context = abc_celltype,
  link_strength = abc_score,
  pip_weighted_link = pmax(substrate_pip, max_rec_pip, na.rm = TRUE) * abc_score,
  direct_link_measured = FALSE,
  evidence_class = "published_ABC_prediction",
  disease_context = FALSE
)]

# Source 2: adult MASLD multiome peak-gene links. Coordinates are hg38.
link_files <- list.files(SCENIC, pattern = "enhancer_gene_links\\.csv$", full.names = TRUE)
ctx_from_file <- function(p) {
  x <- sub("_enhancer_gene_links\\.csv$", "", basename(p))
  if (x == "enhancer_gene_links.csv") "global" else x
}
multi <- rbindlist(lapply(link_files, function(p) {
  x <- fread(p)
  x[, cell_context := ctx_from_file(p)]
  x
}), fill = TRUE)
multi[, `:=`(enhancer_start = as.integer(enhancer_start),
             enhancer_end = as.integer(enhancer_end),
             correlation_rna_atac = as.numeric(correlation_rna_atac),
             tf_enhancer_corr = as.numeric(tf_enhancer_corr))]
setkey(multi, enhancer_chr, enhancer_start, enhancer_end)
q <- sub1[, .(variant_key, chr, start = pos_hg38, end = pos_hg38,
              pos_hg38, var_class, consequence, substrate_gene = gene,
              substrate_pip = max_pip)]
setkey(q, chr, start, end)
ov <- foverlaps(q, multi, by.x = c("chr", "start", "end"),
                by.y = c("enhancer_chr", "enhancer_start", "enhancer_end"),
                type = "within", nomatch = 0L)
sc <- ov[, .(
  variant_key, pos_hg38, var_class, consequence, substrate_gene, substrate_pip,
  target_gene, cell_context, tf_name,
  correlation_rna_atac, tf_enhancer_corr, distance_to_tss
)]
sc[, `:=`(
  source = "MASLD_multiome_SCENICplus",
  link_strength = abs(correlation_rna_atac),
  pip_weighted_link = substrate_pip * abs(correlation_rna_atac),
  direct_link_measured = FALSE,
  evidence_class = "disease_context_eGRN_prediction",
  disease_context = TRUE
)]

common <- c("variant_key", "pos_hg38", "var_class", "consequence",
            "substrate_gene", "substrate_pip", "target_gene", "cell_context",
            "source", "evidence_class", "link_strength", "pip_weighted_link", "direct_link_measured",
            "disease_context")
links <- rbindlist(list(abc[, ..common], sc[, ..common]), fill = TRUE)
links[, coding_ineligible := var_class == "coding"]
links <- links[coding_ineligible == FALSE]
links[, source_rank_within_variant := frank(-pip_weighted_link, ties.method = "min"),
      by = .(variant_key, source)]
links[, cross_method_support := uniqueN(source) > 1L,
      by = .(variant_key, target_gene)]
links[, effector_gene_validated := FALSE]
setorder(links, variant_key, -cross_method_support, source,
         source_rank_within_variant, -pip_weighted_link)

fwrite(links, file.path(OUT, "measured_variant_gene_links.tsv"), sep = "\t", na = "")
summary <- links[, .(
  n_links = .N,
  n_variants = uniqueN(variant_key),
  n_genes = uniqueN(target_gene),
  n_cross_method_variant_gene = uniqueN(paste(variant_key, target_gene)[cross_method_support])
), by = source]
fwrite(summary, file.path(OUT, "link_summary.tsv"), sep = "\t")

contract <- list(
  status = "complete_predicted_link_v2",
  apply_only_firewall = TRUE,
  primary_link_prediction = "Nasser2021 liver/HepG2 ABC AvgHiC prediction",
  disease_context_prediction = "MASLD adult-liver multiome SCENIC+ region-to-gene prediction",
  excluded_from_effector_call = "AlphaGenome O/E contact-map argmax",
  score_interpretation = paste(
    "Link-strength values rank predictions within a source. Legacy PIP x link",
    "values are descriptive only; ABC and SCENIC+ scores are not on a common",
    "scale, are not independent validation, and are never summed."
  ),
  limitations = c(
    "Published ABC liver is not disease-state or cell-type resolved.",
    "HepG2 is a carcinoma line and is reported separately from primary liver.",
    "SCENIC+ correlation is not a physical contact measurement.",
    "Neither source experimentally validates a specific variant-gene link."
  )
)
write_json(contract, file.path(OUT, "contract.json"), pretty = TRUE, auto_unbox = TRUE)
message("[93] wrote ", nrow(links), " noncoding predicted/contextual links to ", OUT)
