#!/usr/bin/env Rscript
# C3c: Model C external-test targets from TCGA-LIHC non-tumor liver (PRESPEC_C Amendment 6).
# Same scorer as c3_targets.R (collapse_symbols -> score_programs, coverage 0.80) and the
# same two negative-control gene lists, with the C1b-matched TCGA non-tumor samples as the
# cohort for within-cohort z-scores. Input: UCSC Xena GDC-hub STAR TPM, already
# log2(TPM + 1). No image feature is read here.
# Usage: Rscript c3c_tcga_targets.R <c1b_exec_dir> <star_tpm.tsv.gz> <out_dir>
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
c1b <- args[[1]]; tpm_file <- args[[2]]; out <- args[[3]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
repo <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
here <- file.path(repo, "scripts/analysis/he_molecular")
owd <- setwd(here); stopifnot(system2("sha256sum", c("-c", "--quiet", "PRESPEC.sha256")) == 0); setwd(owd)
source(file.path(repo, "scripts/manuscript/program_observability_map/analysis_lib.R"))
rel <- file.path(repo, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10")
annotation <- fread(file.path(rel, "inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"), select = c("gene_id", "gene_name"))
hot <- file.path(repo, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot")
registry <- fread(file.path(hot, "program_registry_v2.tsv")); membership <- fread(file.path(hot, "program_membership_v2.tsv"))

cov <- fread(file.path(c1b, "liver_slides_matched.tsv"))[, .(SUBJID = PatientID, sample, tss, gender, age)]
hdr <- names(fread(tpm_file, nrows = 0))
stopifnot(all(cov$sample %in% hdr))
tpm <- fread(tpm_file, select = c(hdr[1], cov$sample))
lib <- as.matrix(tpm[, -1]); rownames(lib) <- tpm[[1]]
lib <- lib[rowMeans(lib >= 1) >= 0.10, , drop = FALSE]  # TPM >= 1 in at least 10% of samples
sm <- collapse_symbols(lib, rownames(lib), annotation)
colnames(sm) <- cov$SUBJID[match(colnames(sm), cov$sample)]
meta <- data.table(sample_id = colnames(sm), dataset = "TCGA_LIHC_nontumor")
s <- score_programs(sm, meta, membership, registry, 0.80)
scores <- t(s$primary)

Z <- zscore_rows(sm); Z <- Z[rowSums(is.finite(Z)) == ncol(Z), , drop = FALSE]
neg <- sapply(c(hypoxia = "negative_control_hypoxia_hallmark.txt", immediate_early = "negative_control_immediate_early.txt"),
              function(f) {g <- intersect(toupper(readLines(file.path(here, f))), rownames(Z)); colMeans(Z[g, , drop = FALSE])})

tab <- merge(cov, data.table(SUBJID = rownames(scores), scores), by = "SUBJID")
tab <- merge(tab, data.table(SUBJID = rownames(neg), neg_hypoxia = neg[, 1], neg_immediate_early = neg[, 2]), by = "SUBJID")
fwrite(tab, file.path(out, "targets_tcga.tsv"), sep = "\t")
fwrite(s$coverage, file.path(out, "program_coverage.tsv"), sep = "\t")
cat("donors", nrow(tab), "| programs scored", sum(colSums(is.finite(scores)) > 0), "\n")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_c3c_tcga_targets.txt"))
