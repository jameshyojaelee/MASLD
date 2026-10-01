#!/usr/bin/env Rscript
# B2: Model B development inputs (PRESPEC_B sections 2-3 and Amendment 1).
# One row per participant (first biopsy) in the six development cohorts:
# histology labels, sex, age, and the 117 pinned program scores from
# score_programs() (coverage 0.80), computed per cohort from
#   corr: ALT-corrected counts (fold_alt_copies; ALT/patch rows dropped)   [primary]
#   orig: the counts the Resource used, all rows                           [sensitivity]
# Genes kept per cohort: CPM >= 1 in >= 10% of its participants. GSE267145 uses its
# deposited counts in both versions. GSE281797 (sealed) is not read here.
# Usage: Rscript b2_build_inputs.R <multi_counts_dir> <out_dir>
suppressPackageStartupMessages({library(data.table); library(edgeR)})
args <- commandArgs(trailingOnly = TRUE)
multi_dir <- args[[1]]; out <- args[[2]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
repo <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(repo, "scripts/manuscript/program_observability_map/analysis_lib.R"))
source(file.path(repo, "scripts/analysis/alt_contig_recount/fold_lib.R"))
rel <- file.path(repo, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10")
prov <- fread(file.path(repo, "results/remediation/bg001/bg001-fragment-v211-gencode49-20260807T195243Z/arms/F_five/provenance/effective_count_sources.tsv"))
annotation <- fread(file.path(rel, "inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"), select = c("gene_id", "gene_name"))
hot <- file.path(repo, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot")
registry <- fread(file.path(hot, "program_registry_v2.tsv")); membership <- fread(file.path(hot, "program_membership_v2.tsv"))
stopifnot(sha256_file(file.path(hot, "program_registry_v2.tsv")) == "b134b5f8e46a716e857af22271aa54e1294432b9f9ab735559d7f603cdd29fe7",
          sha256_file(file.path(hot, "program_membership_v2.tsv")) == "feec3fc9ccaaffaa9abc1ac5d593cc06845ed8140460c8873a51bc8652d7336b",
          sha256_file(file.path(repo, "scripts/manuscript/program_observability_map/analysis_lib.R")) == "a6439390506398be9283eb335b3866c45e9a510fa87aec07f1bf7fc21accb71b")
genes <- gene_table("/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz")
meta_all <- fread(file.path(repo, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"), na.strings = "")

# ---------------------------------------------------------------- labels
lab <- list()
sra <- fread(file.path(repo, "RNA-seq/Human/Patient_Cohorts/metadata/source/GSE130970/GSE130970_SraRunTable.csv"))
lab$GSE130970 <- sra[, .(sample = Run, participant = Run, steatosis = steatosis_grade,
                         ballooning = cytological_ballooning_grade, inflammation = lobular_inflammation_grade,
                         nas = NA_integer_, fib_lo = fibrosis_stage, fib_hi = fibrosis_stage,
                         sex = toupper(substr(sex, 1, 1)), age = as.numeric(age_at_biopsy))]
pj <- fread(file.path(repo, "Analysis/MASLD_Model_Benchmark/executions/gse267145-authoritative-join-21064930/participant_join.tsv"))
lab$GSE267145 <- pj[, .(sample = participant_id, participant = participant_id, steatosis, ballooning,
                        inflammation = lobular_inflammation, nas = NA_integer_, fib_lo = fibrosis, fib_hi = fibrosis,
                        sex, age = NA_real_)]
for (c in c("GSE135251", "GSE162694")) {
  lab[[c]] <- meta_all[dataset == c, .(sample = sample_id, participant = sample_id, steatosis = NA_integer_,
                                       ballooning = NA_integer_, inflammation = NA_integer_, nas = as.integer(nas_score),
                                       fib_lo = as.integer(fibrosis_stage), fib_hi = as.integer(fibrosis_stage), sex, age)]
}
cw <- fread(file.path(rel, "workstreams/BULK-STAGE/PREFLIGHT-v1/audits/gse193066_participant_crosswalk.tsv"))
m193 <- merge(meta_all[dataset == "GSE193066"], cw[is_first_biopsy == TRUE, .(sample_id = run_id, participant_token)], by = "sample_id")
lab$GSE193066 <- m193[, .(sample = sample_id, participant = participant_token, steatosis = NA_integer_, ballooning = NA_integer_,
                          inflammation = NA_integer_, nas = as.integer(nas_score), fib_lo = as.integer(fibrosis_stage),
                          fib_hi = as.integer(fibrosis_stage), sex, age)]
m213 <- meta_all[dataset == "GSE213621" & condition != "Control"]
iv <- list(Fibrosis_F0F1 = c(0L, 1L), Fibrosis_F2 = c(2L, 2L), Fibrosis_F3F4 = c(3L, 4L))
stopifnot(all(m213$condition %in% names(iv)))
lab$GSE213621 <- m213[, .(sample = sample_id, participant = sample_id, steatosis = NA_integer_, ballooning = NA_integer_,
                          inflammation = NA_integer_, nas = NA_integer_,
                          fib_lo = vapply(condition, function(x) iv[[x]][1], 1L),
                          fib_hi = vapply(condition, function(x) iv[[x]][2], 1L), sex, age)]
lab$GSE240729 <- meta_all[dataset == "GSE240729", .(sample = sample_id, participant = sample_id, steatosis = NA_integer_,
                          ballooning = NA_integer_, inflammation = NA_integer_, nas = NA_integer_,
                          fib_lo = as.integer(fibrosis_stage), fib_hi = as.integer(fibrosis_stage), sex, age)]

# ---------------------------------------------------------------- scores
score_cohort <- function(counts, cohort) {
  keep <- rowMeans(cpm(counts) >= 1) >= 0.10
  lib <- cpm(DGEList(counts[keep, , drop = FALSE]), log = TRUE, prior.count = 1)
  sm <- collapse_symbols(lib, rownames(lib), annotation)
  s <- score_programs(sm, data.table(sample_id = colnames(sm), dataset = cohort), membership, registry, 0.80)
  list(scores = t(s$primary), coverage = s$cohort_coverage)
}

res <- list(corr = list(), orig = list()); cover <- list(); replaced <- list()
for (c in names(lab)) {
  L <- lab[[c]][!is.na(participant)]
  if (c == "GSE267145") {
    x <- fread(file.path(repo, "Analysis/MASLD_Model_Benchmark/executions/gse267145-matrix-audit-21064333/raw/GSE269412_RNA.txt.gz"))
    m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]
    m <- m[, intersect(colnames(m), L$sample), drop = FALSE]
    versions <- list(corr = m, orig = m)
  } else {
    orig <- read_featurecounts(prov[dataset == c, count_path])
    multi <- read_featurecounts(file.path(multi_dir, paste0(c, ".multi.txt")))
    runs <- intersect(colnames(orig), L$sample)
    f <- fold_alt_copies(orig[, runs, drop = FALSE], multi[, runs, drop = FALSE], genes)
    replaced[[c]] <- length(f$replaced)
    versions <- list(corr = f$counts, orig = orig[, runs, drop = FALSE])
  }
  for (v in names(versions)) {
    s <- score_cohort(versions[[v]], c)
    d <- data.table(sample = rownames(s$scores), s$scores)
    res[[v]][[c]] <- merge(L, d, by = "sample")[, cohort := c]
    if (v == "corr") cover[[c]] <- s$coverage
  }
  cat(c, ": participants", nrow(res$corr[[c]]), "\n")
}
for (v in names(res)) fwrite(rbindlist(res[[v]], use.names = TRUE, fill = TRUE), file.path(out, paste0("b_inputs_", v, ".tsv.gz")), sep = "\t")
fwrite(rbindlist(cover), file.path(out, "program_cohort_coverage_corr.tsv"), sep = "\t")
writeLines(jsonlite::toJSON(list(replaced_genes = replaced, participants = lapply(res$corr, nrow)), auto_unbox = TRUE, pretty = TRUE),
           file.path(out, "b2_summary.json"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_b2.txt"))
