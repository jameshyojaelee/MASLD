# SS S3d reproduction mocks (not SteatoSITE). Writes S2-layout inputs (S2_CONTRACT.md) from existing cohorts so that
# s3d_programs.R can be checked against stored results:
#  fig4e_<cohort>/  GSE130970, GSE135251, GSE162694 (all records of the cohort in the 774-unit release DGE, release
#                   RLE factors kept) and GSE174478 (all its records in the BULK-PROGRAM-MAP-v9 nonholdout DGE, the
#                   substrate the Fig 4E four-cohort arm scored it on). Stage file = recorded fibrosis stage; units.tsv
#                   adds nas_score. batch = one level.
#  h54_A91/         GSE281797 set A91 of PRESPEC_H54 (91 sex-concordant participants, sex_expressed), counts as
#                   deposited (fractional), order as in h54_lib.R build_set_meta(); nas_score added.
# Usage: Rscript s3d_mock_inputs.R <new out dir>
suppressPackageStartupMessages({library(edgeR); library(data.table)})
out <- commandArgs(trailingOnly = TRUE)[1]; stopifnot(!dir.exists(out)); dir.create(out, recursive = TRUE)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
Sys.setenv(MASLD_PROJECT_ROOT = ROOT)
source(file.path(ROOT, "scripts/manuscript/program_observability_map/config.R"))
REL <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")

write_mock <- function(dir, dge, meta) {
  # meta: sample_id, sex_final, fibrosis_stage, nas_score
  dir.create(dir); stopifnot(identical(colnames(dge), meta$sample_id))
  lib <- colSums(dge$counts)
  u <- data.table(sample_alias = meta$sample_id, analysis_unit_id = meta$sample_id, batch = "mock",
                  inferred_sex = meta$sex_final, assigned_fragments = lib, log10_assigned = log10(lib),
                  unique_pct = NA_real_, intron_exon_ratio = NA_real_, q1_pass = TRUE, q2_pass = TRUE, q3_pass = TRUE,
                  include = TRUE, exclusion_reason = NA_character_, nas_score = meta$nas_score)
  fwrite(u, file.path(dir, "units.tsv"), sep = "\t", na = "NA")
  dge$samples <- cbind(dge$samples[, c("group", "lib.size", "norm.factors")], as.data.frame(u))
  saveRDS(dge, file.path(dir, "dge_included.rds"))
  saveRDS(dge$counts, file.path(dir, "counts_all_genes.rds"))
  st <- paste0(c("sample_accession\tsample_alias\tnash_crn_kleiner_fibrosis_stage",
                 sprintf("%s\t%s\t%s", meta$sample_id, meta$sample_id, ifelse(is.na(meta$fibrosis_stage), "NA", meta$fibrosis_stage))), "\r")
  writeLines(st, file.path(dir, "stage.tsv"))
  cat(dir, ncol(dge), "units,", nrow(dge), "rows\n")
}

d <- readRDS(file.path(REL, "merged_dge.rds")); m <- as.data.table(readRDS(file.path(REL, "meta_matched.rds")))
m <- m[match(colnames(d), sample_id)]
for (cid in c("GSE130970", "GSE135251", "GSE162694")) {
  j <- which(m$dataset == cid)
  write_mock(file.path(out, paste0("fig4e_", cid)), d[, j], m[j, .(sample_id, sex_final, fibrosis_stage, nas_score)])
}
rm(d); invisible(gc())
od <- readRDS(nonholdout_dge); om <- as.data.table(readRDS(nonholdout_meta)); om <- om[match(colnames(od), sample_id)]
j <- which(om$dataset == "GSE174478")
write_mock(file.path(out, "fig4e_GSE174478"), od[, j], om[j, .(sample_id, sex_final, fibrosis_stage, nas_score)])
rm(od); invisible(gc())

a <- fread(file.path(ROOT, "data/GSE281797/metadata/sample_phenotypes_resolved.tsv"), colClasses = "character", na.strings = "")
a <- a[sex_concordant == "True"]; stopifnot(nrow(a) == 91L)
cts <- fread(file.path(ROOT, "data/GSE281797/GSE281797_counts_ensembl_stable.tsv.gz"), header = TRUE)
M <- as.matrix(cts[, -1, with = FALSE]); rownames(M) <- sub("\\..*$", "", cts[[1]]); storage.mode(M) <- "double"
M <- M[, a$matrix_column]; colnames(M) <- a$sample_id
write_mock(file.path(out, "h54_A91"), DGEList(M),
           a[, .(sample_id, sex_final = sex_expressed, fibrosis_stage = as.integer(fibrosis), nas_score = as.numeric(nas))])
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
